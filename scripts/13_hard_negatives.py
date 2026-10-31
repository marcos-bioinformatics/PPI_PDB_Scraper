#!/usr/bin/env python3
"""Step 13 - Hard negatives: chains of the same assembly that are not in contact.

Input : paths.chain_metadata, paths.entry_metadata (06), PISA JSON (01),
        paths.clustered (09), paths.test (11), paths.trainval_splits (12)
Output: paths.hard_negatives                 one row per negative pair, with its split
        reports_dir/13_hard_negatives.json   candidates left after each filter

A hard negative is a pair of protein chains of the same biological assembly
with no PISA interface between them (of any size). Candidates are filtered in
this order:
  1. hetero pairs only, unless hard_negatives.include_homo
  2. both chain sequences present in the clustered positives (step 09), so the
     negative shares the protein universe, cluster IDs and splits of the positives
  3. entry and chain quality: the same method, R-free, coverage and length
     thresholds as step 07 (interface thresholds do not apply)
  4. no conflict: the pair of chain clusters (hard_negatives.conflict_level) is
     not in contact in any downloaded assembly
  5. one negative per (entry, sequence pair)
  6. structure: no forbidden special residues and a minimum heavy-atom distance
     between the chains of hard_negatives.min_distance
  7. at most hard_negatives.max_per_entry per entry (random, seeded)
Split: "test" if a chain is in the test set or in a cluster held out for the
unseen-chain (C3) buckets; else "val" if a chain appears only in validation
positives; else "train".
"""

from __future__ import annotations

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple

import pandas as pd
from tqdm import tqdm

from common import ensure_parent, get_logger, parse_args, pisa_subdir, read_json, seq_uid, write_json
from structure_checks import contact_summary, heavy_atoms, load_model, special_residues

log = get_logger("13_hard_negatives")


def pisa_contacts(task: Tuple[str, int, str]) -> Tuple[str, int, List[Tuple[str, str]]]:
    """Protein chain pairs with any PISA interface in one assembly."""
    pdb_id, asm, path = task
    if not Path(path).is_file():
        return pdb_id, asm, []
    data = read_json(path)
    block = (data.get(pdb_id, data) or {}).get("assembly") or {}
    out = []
    for itf in block.get("interfaces") or []:
        mols = [m for m in itf.get("molecules") or [] if str(m.get("molecule_class", "")).lower() == "protein"]
        if len(mols) == 2:
            out.append(tuple(sorted((str(mols[0]["chain_id"]), str(mols[1]["chain_id"])))))
    return pdb_id, asm, out


def check_assembly(task: Dict[str, Any]) -> List[Dict[str, Any]]:
    try:
        model = load_model(task["path"])
    except Exception:  # unreadable file: its candidates are dropped
        return []
    atoms, special_ok = {}, {}
    for chain_id in {c for pair in task["pairs"] for c in pair}:
        if chain_id in model:
            special_ok[chain_id] = not special_residues(model[chain_id], task["allowed"])[1]
            atoms[chain_id] = heavy_atoms(model[chain_id])
    out = []
    for a, b in task["pairs"]:
        if special_ok.get(a) and special_ok.get(b):
            out.append({"pdb_id": task["pdb_id"], "assembly_id": task["assembly_id"],
                        "molA_chain_id": a, "molB_chain_id": b,
                        **contact_summary(atoms[a], atoms[b], task["cutoff"])})
    return out


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
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for p, a, pairs in tqdm(ex.map(pisa_contacts, tasks, chunksize=200), total=len(tasks), desc="PISA contacts"):
            contacts[(p, a)] = set(pairs)
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
    cols = ["pdb_id", "assembly_id", "chain_id", "entity_id", "seq", "seq_len", "coverage", "coverage_pct"]
    for side in ("molA", "molB"):
        cand = cand.merge(chains[cols].rename(columns={c: f"{side}_{c}" for c in cols[2:]}),
                          on=["pdb_id", "assembly_id", f"{side}_chain_id"], how="left", validate="many_to_one")
    cand = cand.merge(entries, on="pdb_id", how="left", validate="many_to_one")
    ok = (cand["exptl_method"] == qf["exptl_method"]) & (cand["r_free"] < qf["max_r_free"])
    for side in ("molA", "molB"):
        ok &= (cand[f"{side}_coverage"] > qf["min_coverage"])
        ok &= (cand[f"{side}_seq_len"] > qf["min_chain_length"]) & (cand[f"{side}_seq_len"] <= qf["max_chain_length"])
    cand = cand[ok.fillna(False)]
    counts["after_quality"] = len(cand)

    # --- 4) conflicts and 5) one per (entry, sequence pair)
    keys = [tuple(sorted((conflict_key(a), conflict_key(b)))) for a, b in zip(cand["uidA"], cand["uidB"])]
    cand = cand[[k not in touching for k in keys]]
    counts["after_conflicts"] = len(cand)
    cand = cand.assign(_pair=[tuple(sorted(x)) for x in zip(cand["uidA"], cand["uidB"])],
                       _cov=cand[["molA_coverage", "molB_coverage"]].min(axis=1))
    cand = (cand.sort_values(["pdb_id", "_cov", "assembly_id", "molA_chain_id", "molB_chain_id"],
                             ascending=[True, False, True, True, True])
                .drop_duplicates(["pdb_id", "_pair"]))
    counts["after_dedup_per_entry"] = len(cand)

    # --- 6) structure checks: special residues and minimum distance
    asm_dir = Path(paths["assembly_dir"])
    allowed = {x.upper() for x in cfg["special_residues"]["allowed"]}
    by_asm = defaultdict(list)
    for p, a, ca, cb in cand[["pdb_id", "assembly_id", "molA_chain_id", "molB_chain_id"]].itertuples(index=False):
        by_asm[(p, a)].append((ca, cb))
    tasks = [{"pdb_id": p, "assembly_id": a, "pairs": prs, "allowed": allowed,
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
    counts["after_min_distance"] = len(cand)

    # --- 7) cap per entry
    cand = (cand.sample(frac=1.0, random_state=int(hn["seed"]))
                .groupby("pdb_id", sort=False).head(int(hn["max_per_entry"]))
                .sort_values(["pdb_id", "assembly_id", "molA_chain_id", "molB_chain_id"]))
    counts["final"] = len(cand)

    # --- clusters and split
    for side, u in (("molA", "uidA"), ("molB", "uidB")):
        for t in tags:
            cand[f"{side}_cluster_{t}"] = cand[u].map(cluster_of[t])
    for t in tags:
        a, b = cand[f"molA_cluster_{t}"], cand[f"molB_cluster_{t}"]
        cand[f"pair_cluster_{t}"] = a.where(a <= b, b) + "||" + b.where(a <= b, a)

    test = pd.read_csv(paths["test"], keep_default_na=False, na_values=[""], low_memory=False)
    tv = pd.read_csv(paths["trainval_splits"], keep_default_na=False, na_values=[""], low_memory=False)
    split_col = cfg["splits"]["validation_split"]

    def uids(df: pd.DataFrame) -> Set[str]:
        return set(df["molA_seq"].map(seq_uid)) | set(df["molB_seq"].map(seq_uid))

    test_uids, val_uids, train_uids = uids(test), uids(tv[tv[split_col] == "val"]), uids(tv[tv[split_col] == "train"])
    held_out = {t: set() for t in ("c40", "c50")}
    for t in held_out:
        c3 = test[test["test_subbucket"] == f"unseen_chain_{t}"]
        held_out[t] = {cluster_of[t].get(u) for u in uids(c3)} - {None}

    def assign(ua: str, ub: str) -> str:
        if ua in test_uids or ub in test_uids or any(
                cluster_of[t].get(u) in held_out[t] for t in held_out for u in (ua, ub)):
            return "test"
        if any(u in val_uids and u not in train_uids for u in (ua, ub)):
            return "val"
        return "train"

    cand["split"] = [assign(a, b) for a, b in zip(cand["uidA"], cand["uidB"])]
    cand.insert(4, "label", 0)
    cand = cand.drop(columns=["uidA", "uidB", "_pair", "_cov"])
    counts.update({f"split_{k}": int(v) for k, v in cand["split"].value_counts().items()})

    dest = ensure_parent(paths["hard_negatives"])
    tmp = dest.with_name(dest.name + ".part")
    cand.to_csv(tmp, index=False)
    tmp.replace(dest)
    write_json(Path(paths["reports_dir"]) / "13_hard_negatives.json",
               {"config": hn, "candidates_left_after_each_filter": counts})
    log.info("Hard negatives: %s", counts)


if __name__ == "__main__":
    main()
