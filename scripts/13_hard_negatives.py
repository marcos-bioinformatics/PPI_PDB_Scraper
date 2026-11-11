#!/usr/bin/env python3
"""Step 13 - Hard negatives: chains of the same assembly that are not in contact.

Input : paths.chain_metadata (06), paths.entry_metadata (04), PISA JSON (01),
        paths.clustered (09), paths.test and paths.test_reserved (11),
        paths.trainval_splits (12)
Output: paths.hard_negatives                 one row per negative pair, with its split
        reports_dir/13_hard_negatives.json   candidates left after each filter

A hard negative is a pair of protein chains of the same biological assembly
with no PISA interface between them (of any size). Candidates are filtered in
this order:
  1. hetero pairs only, unless hard_negatives.include_homo
  2. both chain sequences present in the clustered positives (step 09), so the
     negative shares the protein universe, cluster IDs and splits of the positives
  3. entry and chain quality: the same method, R-free, resolution, coverage and
     length thresholds and the same scope (release date, membrane proteins,
     antibodies) as step 07 (interface thresholds do not apply)
  4. no conflict: the pair of chain clusters (hard_negatives.conflict_level) is
     not in contact in any downloaded assembly
  5. structure: no forbidden special residues (same rules as step 08) and a
     heavy-atom distance between the chains of at least hard_negatives.min_distance
     (and at most hard_negatives.max_distance, if set, to keep only near-miss pairs)
  6. one negative per (entry, sequence pair)
  7. at most hard_negatives.max_per_pair_cluster per pair cluster
     (hard_negatives.pair_cluster_key) across entries, then at most
     hard_negatives.max_per_entry per entry (random, seeded)
  8. IntAct (hard_negatives.intact): every UniProt pair of a negative is looked
     up in IntAct (PSICQUIC, MITAB 2.7). intact_evidences counts its evidences;
     intact_direct_evidences counts the binary ones (no spoke or matrix
     expansion) of type direct interaction (MI:0407) or any of its descendants
     in PSI-MI (read from OLS). Two chains of one complex without contact are
     expected to have association or physical association evidence; direct
     evidence contradicts the negative, which is dropped if
     hard_negatives.intact.drop_direct. Answers are cached, with their date, in
     paths.intact_cache.
Split, for every split scheme of step 12 (splits.schemes, e.g. split_c60, split_c40):
"test" if a chain appears only in test positives or belongs to a chain cluster
held out of training by step 11 (the unseen proteins of C2 and C3); else "val" if a chain appears only in
validation positives of that scheme; else "train". Negatives of the same pair
cluster at hard_negatives.family_level (c40 by default) all get the most
restrictive split of the group (test > val > train), so a family of negatives
is never split between sets. Proteins shared by train and test (as
in C1) do not pull a negative into the test set. The
column `split` repeats the scheme chosen in splits.validation_split.
Missing or unreadable PISA files stop the step (they would hide contacts).
"""

from __future__ import annotations

import re
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple
from urllib.parse import quote

import pandas as pd
import requests
from tqdm import tqdm

from common import ensure_parent, get_logger, make_session, parse_args, pisa_subdir, read_json, seq_uid, write_json
from structure_checks import contact_summary, load_chains, special_residues

log = get_logger("13_hard_negatives")


def pisa_contacts(task: Tuple[str, int, str]) -> Tuple[str, int, List[Tuple[str, str]], str]:
    """Protein chain pairs with any PISA interface in one assembly (and an error, if any)."""
    pdb_id, asm, path = task
    try:
        data = read_json(path)
        block = (data.get(pdb_id, data) or {}).get("assembly") or {}
        out = []
        for itf in block.get("interfaces") or []:
            mols = [m for m in itf.get("molecules") or [] if str(m.get("molecule_class", "")).lower() == "protein"]
            if len(mols) == 2:
                out.append(tuple(sorted((str(mols[0]["chain_id"]), str(mols[1]["chain_id"])))))
        return pdb_id, asm, out, ""
    except Exception as exc:  # reported by main(): a missing file would hide contacts
        return pdb_id, asm, [], f"{path}: {type(exc).__name__}: {exc}"


def check_assembly(task: Dict[str, Any]) -> List[Dict[str, Any]]:
    needed = {c for pair in task["pairs"] for c in pair}
    try:
        chains = load_chains(task["path"], needed)
    except Exception:  # unreadable file: its candidates are dropped
        return []
    special_ok = {c: not special_residues(chains[c], task["allowed"], task["ignore"])[1] for c in chains}
    out = []
    for a, b in task["pairs"]:
        if special_ok.get(a) and special_ok.get(b):
            out.append({"pdb_id": task["pdb_id"], "assembly_id": task["assembly_id"],
                        "molA_chain_id": a, "molB_chain_id": b,
                        **contact_summary(chains[a], chains[b], task["cutoff"])})
    return out


MI_0407_DESCENDANTS = ("https://www.ebi.ac.uk/ols4/api/ontologies/mi/terms/"
                       "http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FMI_0407/hierarchicalDescendants?size=1000")


def intact_crosscheck(cand: pd.DataFrame, icfg: Dict, cache_path: Path) -> Tuple[pd.Series, pd.Series, Dict[str, int]]:
    """IntAct evidences of each negative (all evidences, binary direct ones); NaN if a query failed."""
    session = make_session(pool_size=4)
    today = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    cache = read_json(cache_path) if cache_path.is_file() else {}
    cache.setdefault("pairs", {})
    if not cache.get("direct_types"):
        r = session.get(MI_0407_DESCENDANTS, timeout=(10, 120))
        r.raise_for_status()
        ids = {t["obo_id"] for t in r.json()["_embedded"]["terms"]} | {"MI:0407"}
        cache["direct_types"] = {"date": today, "source": "OLS, PSI-MI descendants of MI:0407", "ids": sorted(ids)}
    direct_types = set(cache["direct_types"]["ids"])

    def acc(x) -> List[str]:
        return [a for a in str(x).split(";") if a and a != "nan"]
    row_pairs = [sorted({"|".join(sorted((a, b))) for a in acc(x) for b in acc(y) if a != b})
                 for x, y in zip(cand["molA_uniprot_ids"], cand["molB_uniprot_ids"])]
    todo = sorted({p for ps in row_pairs for p in ps} - set(cache["pairs"]))

    def query(pair: str):
        a, b = pair.split("|")
        url = icfg["url"].rstrip("/") + "/" + quote(f"id:{a} AND id:{b}") + "?format=tab27"
        for attempt in range(4):
            try:
                r = session.get(url, timeout=(10, 120))
                if r.status_code == 200:
                    rows = [line.split("\t") for line in r.text.splitlines() if line.strip()]
                    return pair, [[(re.findall(r"MI:\d{4}", x[11]) or [""])[0],
                                   (re.findall(r"MI:\d{4}", x[6]) or [""])[0],
                                   x[15] if len(x) > 15 else "-"] for x in rows]
            except requests.RequestException:
                pass
            time.sleep(2 ** attempt)
        return pair, None

    failures = 0
    with ThreadPoolExecutor(max_workers=4) as ex:
        for pair, rows in ex.map(query, todo):
            if rows is None:
                failures += 1
            else:
                cache["pairs"][pair] = {"date": today, "rows": rows}
    write_json(ensure_parent(cache_path), cache)

    def count(pairs, direct_only):
        if any(p not in cache["pairs"] for p in pairs):
            return float("nan")
        rows = [r for p in pairs for r in cache["pairs"][p]["rows"]]
        if direct_only:
            rows = [r for r in rows if r[0] in direct_types and r[2] in ("-", "")]
        return len(rows)
    evid = pd.Series([count(ps, False) if ps else float("nan") for ps in row_pairs], index=cand.index)
    direct = pd.Series([count(ps, True) if ps else float("nan") for ps in row_pairs], index=cand.index)
    stats = {"intact_pairs_queried": len(todo), "intact_query_failures": failures,
             "negatives_without_uniprot_pair": int(sum(not ps for ps in row_pairs)),
             "negatives_with_intact_evidence": int((evid > 0).sum()),
             "negatives_with_direct_evidence": int((direct > 0).sum())}
    return evid, direct, stats


def main() -> None:
    cfg = parse_args(__doc__)
    paths, hn, qf = cfg["paths"], cfg["hard_negatives"], cfg["quality_filter"]
    workers = int(cfg["workers"])
    counts: Dict[str, int] = {}

    chains = pd.read_csv(paths["chain_metadata"], keep_default_na=False, na_values=[""],
                         dtype={"pdb_id": str, "chain_id": str, "entity_id": str})
    entries = pd.read_csv(paths["entry_metadata"], keep_default_na=False, na_values=[""], dtype={"pdb_id": str})
    chains["uid"] = chains["seq"].map(seq_uid)
    seq_of = {(p, a, c): u for p, a, c, u in zip(chains["pdb_id"], chains["assembly_id"], chains["chain_id"], chains["uid"])}

    # --- cluster maps from the positives (step 09): uid -> representative per threshold
    pos = pd.read_csv(paths["clustered"], keep_default_na=False, na_values=[""], low_memory=False)
    tags = sorted({c.split("_")[-1] for c in pos.columns if c.startswith("molA_cluster_c")}, reverse=True)
    cluster_of: Dict[str, Dict[str, str]] = {t: {} for t in tags}
    for side in ("molA", "molB"):
        uids = pos[f"{side}_seq"].map(seq_uid)
        for t in tags:
            cluster_of[t].update(zip(uids, pos[f"{side}_cluster_{t}"]))
    known = set(cluster_of[tags[0]])
    level = str(hn["conflict_level"])

    def conflict_key(u: str) -> str:
        return u if level == "exact" else cluster_of[level].get(u, u)

    # --- 1) all PISA contacts of the downloaded assemblies
    pisa_dir = Path(paths["pisa_dir"])
    assemblies = sorted(set(zip(chains["pdb_id"], chains["assembly_id"])))
    tasks = [(p, a, str(pisa_dir / pisa_subdir(p) / f"{p}_assembly{a}_interfaces.json.gz")) for p, a in assemblies]
    contacts: Dict[Tuple[str, int], Set[Tuple[str, str]]] = {}
    errors = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for p, a, pairs, err in tqdm(ex.map(pisa_contacts, tasks, chunksize=200), total=len(tasks),
                                     desc="PISA contacts"):
            contacts[(p, a)] = set(pairs)
            if err:
                errors.append(err)
    if errors:
        raise SystemExit(f"{len(errors)} PISA files are missing or unreadable (e.g. {errors[0]}); "
                         "rerun step 01 (snakemake --forcerun download_pisa).")
    touching = set()
    for (p, a), pairs in contacts.items():
        for x, y in pairs:
            if (p, a, x) in seq_of and (p, a, y) in seq_of:
                touching.add(tuple(sorted((conflict_key(seq_of[(p, a, x)]), conflict_key(seq_of[(p, a, y)])))))

    # --- 2) candidates: chain pairs of the same assembly without PISA interface
    rows = []
    for (p, a), g in chains.groupby(["pdb_id", "assembly_id"]):
        recs = list(g[["chain_id", "uid"]].itertuples(index=False))
        for (ca, ua), (cb, ub) in combinations(sorted(recs), 2):
            if (ca, cb) not in contacts.get((p, a), set()):
                rows.append((p, a, ca, cb, ua, ub, len(recs)))
    cand = pd.DataFrame(rows, columns=["pdb_id", "assembly_id", "molA_chain_id", "molB_chain_id",
                                       "uidA", "uidB", "n_protein_chains_assembly"])
    counts["candidates_no_pisa_interface"] = len(cand)

    if not hn["include_homo"]:
        cand = cand[cand["uidA"] != cand["uidB"]]
    counts["after_homo_filter"] = len(cand)
    cand = cand[cand["uidA"].isin(known) & cand["uidB"].isin(known)]
    counts["after_known_proteins"] = len(cand)

    # --- 3) entry and chain quality (same thresholds as step 07)
    cols = ["pdb_id", "assembly_id", "chain_id", "entity_id", "seq", "seq_len", "coverage", "coverage_pct",
            "internal_coverage", "missing_nterm", "missing_cterm", "uniprot_ids", "membrane_sources", "is_antibody"]
    for side in ("molA", "molB"):
        cand = cand.merge(chains[cols].rename(columns={c: f"{side}_{c}" for c in cols[2:]}),
                          on=["pdb_id", "assembly_id", f"{side}_chain_id"], how="left", validate="many_to_one")
    cand = cand.merge(entries, on="pdb_id", how="left", validate="many_to_one")
    ok = cand["r_free"] < qf["max_r_free"]
    if qf.get("exptl_method") is not None:
        ok &= cand["exptl_method"] == qf["exptl_method"]
    if qf.get("max_resolution") is not None:
        ok &= cand["resolution"] <= qf["max_resolution"]
    for side in ("molA", "molB"):
        ok &= (cand[f"{side}_{qf['coverage_metric']}"] > qf["min_coverage"])
        ok &= (cand[f"{side}_seq_len"] > qf["min_chain_length"]) & (cand[f"{side}_seq_len"] <= qf["max_chain_length"])
    # scope, as in step 07: release-date cut-off, membrane proteins, antibodies (on request)
    sc = cfg.get("scope") or {}
    if sc.get("max_release_date") is not None:
        ok &= pd.to_datetime(cand["initial_release_date"]) <= pd.Timestamp(str(sc["max_release_date"]))
    sources = set(sc.get("exclude_membrane") or [])
    for side in ("molA", "molB"):
        cand[f"{side}_is_antibody"] = cand[f"{side}_is_antibody"].fillna(False).astype(bool)
        if sources:
            ok &= ~cand[f"{side}_membrane_sources"].fillna("").str.split(";").map(lambda v: bool(sources & set(v)))
        if sc.get("exclude_antibodies"):
            ok &= ~cand[f"{side}_is_antibody"]
    cand = cand[ok.fillna(False)]
    counts["after_quality"] = len(cand)

    # --- 4) conflicts
    keys = [tuple(sorted((conflict_key(a), conflict_key(b)))) for a, b in zip(cand["uidA"], cand["uidB"])]
    cand = cand[[k not in touching for k in keys]]
    counts["after_conflicts"] = len(cand)

    # --- 5) structure checks: special residues and minimum distance
    asm_dir = Path(paths["assembly_dir"])
    allowed = {x.upper() for x in cfg["special_residues"]["allowed"]}
    ignore = {x.upper() for x in cfg["special_residues"].get("ignore", [])}
    by_asm = defaultdict(list)
    for p, a, ca, cb in cand[["pdb_id", "assembly_id", "molA_chain_id", "molB_chain_id"]].itertuples(index=False):
        by_asm[(p, a)].append((ca, cb))
    tasks = [{"pdb_id": p, "assembly_id": a, "pairs": prs, "allowed": allowed, "ignore": ignore,
              "cutoff": float(cfg["contacts"]["cutoff"]),
              "path": str(asm_dir / pisa_subdir(p) / f"{p}-assembly{a}.cif.gz")}
             for (p, a), prs in sorted(by_asm.items())]
    checked = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for res in tqdm(ex.map(check_assembly, tasks, chunksize=5), total=len(tasks), desc="Structure checks"):
            checked.extend(res)
    geo = pd.DataFrame(checked, columns=["pdb_id", "assembly_id", "molA_chain_id", "molB_chain_id",
                                         "n_atom_contacts", "n_residue_contacts", "min_interchain_distance"])
    cand = cand.merge(geo, on=["pdb_id", "assembly_id", "molA_chain_id", "molB_chain_id"], how="inner")
    counts["after_special_residues"] = len(cand)
    cand = cand[cand["min_interchain_distance"] >= float(hn["min_distance"])]
    if hn.get("max_distance") is not None:
        cand = cand[cand["min_interchain_distance"] <= float(hn["max_distance"])]
    counts["after_min_distance"] = len(cand)

    # --- 6) one negative per (entry, sequence pair): the best-covered copy
    cand = cand.assign(_pair=[tuple(sorted(x)) for x in zip(cand["uidA"], cand["uidB"])],
                       _cov=cand[[f"molA_{qf['coverage_metric']}", f"molB_{qf['coverage_metric']}"]].min(axis=1))
    cand = (cand.sort_values(["pdb_id", "_cov", "assembly_id", "molA_chain_id", "molB_chain_id"],
                             ascending=[True, False, True, True, True])
                .drop_duplicates(["pdb_id", "_pair"]))
    counts["after_dedup_per_entry"] = len(cand)

    # --- clusters of both chains and of the pair
    for side, u in (("molA", "uidA"), ("molB", "uidB")):
        for t in tags:
            cand[f"{side}_cluster_{t}"] = cand[u].map(cluster_of[t])
    for t in tags:
        a, b = cand[f"molA_cluster_{t}"], cand[f"molB_cluster_{t}"]
        cand[f"pair_cluster_{t}"] = a.where(a <= b, b) + "||" + b.where(a <= b, a)

    # --- 7) caps across entries (same complex deposited many times) and per entry
    seed = int(hn["seed"])
    cand = (cand.sample(frac=1.0, random_state=seed)
                .groupby(f"pair_cluster_{hn['pair_cluster_key']}", sort=False).head(int(hn["max_per_pair_cluster"])))
    counts["after_cap_per_pair_cluster"] = len(cand)
    cand = (cand.sample(frac=1.0, random_state=seed)
                .groupby("pdb_id", sort=False).head(int(hn["max_per_entry"]))
                .sort_values(["pdb_id", "assembly_id", "molA_chain_id", "molB_chain_id"]))
    counts["final"] = len(cand)

    # --- 8) IntAct cross-check
    icfg = hn.get("intact") or {}
    if icfg.get("enabled"):
        evid, direct, stats = intact_crosscheck(cand, icfg, Path(paths["intact_cache"]))
        cand["intact_evidences"], cand["intact_direct_evidences"] = evid, direct
        counts.update(stats)
        if stats["intact_query_failures"]:
            log.warning("%d IntAct queries failed; their negatives have no IntAct counts (rerun to retry)",
                        stats["intact_query_failures"])
        if icfg.get("drop_direct"):
            cand = cand[~(cand["intact_direct_evidences"] > 0)]
        counts["after_intact"] = len(cand)

    test = pd.read_csv(paths["test"], keep_default_na=False, na_values=[""], low_memory=False)
    tv = pd.read_csv(paths["trainval_splits"], keep_default_na=False, na_values=[""], low_memory=False)
    split_col = cfg["splits"]["validation_split"]

    def uids(df: pd.DataFrame) -> Set[str]:
        return set(df["molA_seq"].map(seq_uid)) | set(df["molB_seq"].map(seq_uid))

    test_uids = uids(test)
    level_tag = f"c{int(round(float(cfg['test_set']['level']) * 100)):02d}"
    reserved = set(pd.read_csv(paths["test_reserved"], dtype={"cluster": str})["cluster"])
    held = [cluster_of[level_tag].get(a) in reserved or cluster_of[level_tag].get(b) in reserved
            for a, b in zip(cand["uidA"], cand["uidB"])]
    rank = {"train": 0, "val": 1, "test": 2}
    group_col = f"pair_cluster_{hn['family_level']}"
    if group_col not in cand.columns:
        raise SystemExit(f"hard_negatives.family_level must be one of {', '.join(tags)}")
    schemes = [c for c in tv.columns if c.startswith("split_")]
    for scheme in schemes:
        val_uids, train_uids = uids(tv[tv[scheme] == "val"]), uids(tv[tv[scheme] == "train"])
        test_only = test_uids - train_uids - val_uids
        val_only = val_uids - train_uids - test_uids
        row_split = pd.Series(["test" if h or a in test_only or b in test_only
                               else "val" if a in val_only or b in val_only else "train"
                               for a, b, h in zip(cand["uidA"], cand["uidB"], held)], index=cand.index)
        # whole pair-cluster families go to the most restrictive set
        worst = row_split.map(rank).groupby(cand[group_col]).transform("max")
        row_split = worst.map({v: k for k, v in rank.items()})
        cand[scheme] = row_split
        counts.update({f"{scheme}_{k}": int(v) for k, v in cand[scheme].value_counts().items()})
    cand["split"] = cand[split_col]
    cand.insert(4, "label", 0)
    # coverage_pct repeats coverage x 100; the positives do not export it either
    cand = cand.drop(columns=["uidA", "uidB", "_pair", "_cov", "molA_coverage_pct", "molB_coverage_pct"])

    dest = ensure_parent(paths["hard_negatives"])
    tmp = dest.with_name(dest.name + ".part")
    cand.to_csv(tmp, index=False)
    tmp.replace(dest)
    write_json(Path(paths["reports_dir"]) / "13_hard_negatives.json",
               {"config": hn, "candidates_left_after_each_filter": counts})
    log.info("Hard negatives: %s", counts)


if __name__ == "__main__":
    main()
