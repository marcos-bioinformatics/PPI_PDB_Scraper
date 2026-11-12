#!/usr/bin/env python3
"""Step 14 - Audit homology leakage between the evaluation sets and the data used for training.

Input : paths.test (11), paths.trainval_splits (12), paths.hard_negatives (13, if enabled)
Output: paths.leakage_audit                     summary per set at its design threshold
        reports_dir/14_leakage_by_threshold.tsv  every set at every clustering threshold,
                                                 with and without the coverage rule

Chains are searched with MMseqs2 in both directions (as in step 09) and linked
with the identity, coverage and E-value criteria of step 09. Sets and the data
they are compared with:
  test buckets (C1, C2, C3)       all train + validation positives
  test buckets, "_vs_negatives"   the train + validation hard negatives
  validation of each split_cXX     train positives of that scheme
  hard negatives in the test set   train + validation positives and negatives
For each set and threshold it reports:
  chains_with_homologue  chains with a reference chain at >= threshold identity
  rows_any_chain         pairs with at least one such chain
  rows_pair_homologue    pairs for which a reference pair links homologues of both chains
The design threshold of each set is where leakage should be zero (L =
test_set.level, 0.40 by default):
  C3                       chain level at L, both chains (proteins unseen in training)
  C2_unseen_chain          chain level at L, the unseen chain of each C2 pair
  C2, C1                   pair level at L (no training pair joins homologues of both chains)
  validation split_cXX     pair level at 0.XX
  hard negatives (test)    pair level at the conflict level of step 13
The design-threshold summary uses the linking rule of step 09 (identity and
coverage). The audit is repeated without the coverage rule (any alignment
that reaches the identity threshold links two chains, e.g. a single shared
domain); that second summary is reported, not enforced.
The audit only reports; it does not change any dataset.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

from common import ensure_parent, get_logger, homology_linked, parse_args, seq_uid, write_json

log = get_logger("14_leakage_audit")



def search(query: Dict[str, str], target: Dict[str, str], work: Path, ccfg: Dict, threads: int) -> pd.DataFrame:
    """All hits query -> target above the lowest threshold, with a column
    `linked` that says whether they satisfy the link criterion of step 09."""
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    for name, seqs in (("q.fa", query), ("t.fa", target)):
        (work / name).write_text("".join(f">{u}\n{s}\n" for u, s in seqs.items()), encoding="utf-8")
    cmd = ["mmseqs", "easy-search", str(work / "q.fa"), str(work / "t.fa"), str(work / "hits.m8"), str(work / "tmp"),
           "--min-seq-id", str(min(ccfg["thresholds"])), "-c", "0", "-e", str(ccfg["evalue"]),
           "-s", str(ccfg["search_sensitivity"]), "--max-seqs", str(ccfg["max_seqs"]), "--threads", str(threads),
           "--format-output", "query,target,fident,qcov,tcov"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"MMseqs2 search failed:\n{p.stderr}")
    hits = pd.read_csv(work / "hits.m8", sep="\t", header=None, dtype={0: str, 1: str},
                       names=["query", "target", "fident", "qcov", "tcov"])
    shutil.rmtree(work)
    hits["linked"] = homology_linked(hits, ccfg)
    return hits[["query", "target", "fident", "linked"]]


def uid_pairs(df: pd.DataFrame) -> List[Tuple[str, str]]:
    return list(zip(df["molA_seq"].map(seq_uid), df["molB_seq"].map(seq_uid)))


def evaluation_sets(cfg: Dict) -> Tuple[List[tuple], Dict[str, str]]:
    """Every audited set as (name, pairs, reference pairs, design level, design
    threshold, chains to check or None), and the sequence of every chain uid.
    Also used by analysis/leakage_control.py, so that it audits the same sets."""
    paths = cfg["paths"]
    test = pd.read_csv(paths["test"], keep_default_na=False, na_values=[""], low_memory=False)
    tv = pd.read_csv(paths["trainval_splits"], keep_default_na=False, na_values=[""], low_memory=False)
    schemes = [c for c in tv.columns if c.startswith("split_")]
    negatives: Optional[pd.DataFrame] = None
    if cfg["hard_negatives"]["enabled"] and Path(paths["hard_negatives"]).is_file():
        negatives = pd.read_csv(paths["hard_negatives"], keep_default_na=False, na_values=[""], low_memory=False)

    # (set name, its pairs, reference pairs, design level, design threshold, chains to check or None)
    sets = []
    neg_ref = (uid_pairs(negatives[negatives[cfg["splits"]["validation_split"]] != "test"]) if negatives is not None else [])
    L = float(cfg["test_set"]["level"])
    for bucket, g in test.groupby("test_bucket"):
        level = "chain" if bucket == "C3" else "pair"
        variants = [(bucket, level, None)]
        if bucket == "C2":
            unseen = {seq_uid(s) for side, col in (("A", "molA_seq"), ("B", "molB_seq"))
                      for s in g.loc[g["test_unseen_chains"].fillna("").str.contains(side), col]}
            variants.append(("C2_unseen_chain", "chain", unseen))
        for name, lv, chains in variants:
            sets.append((name, uid_pairs(g), uid_pairs(tv), lv, L, chains))
            if neg_ref:
                sets.append((f"{name}_vs_negatives", uid_pairs(g), neg_ref, lv, L, chains))
    for scheme in schemes:
        m = re.search(r"(\d+)$", scheme)
        sets.append((f"validation_{scheme}", uid_pairs(tv[tv[scheme] == "val"]),
                     uid_pairs(tv[tv[scheme] == "train"]), "pair", int(m.group(1)) / 100 if m else None, None))
    if negatives is not None:
        level = str(cfg["hard_negatives"]["conflict_level"])
        thr = int(level[1:]) / 100 if re.fullmatch(r"c\d+", level) else None
        ref = uid_pairs(tv) + uid_pairs(negatives[negatives[cfg["splits"]["validation_split"]] != "test"])
        sets.append(("hard_negatives_test", uid_pairs(negatives[negatives[cfg["splits"]["validation_split"]] == "test"]), ref, "pair", thr, None))

    seq_of: Dict[str, str] = {}
    for df in [test, tv] + ([negatives] if negatives is not None else []):
        for col in ("molA_seq", "molB_seq"):
            seq_of.update({seq_uid(s): s for s in df[col]})
    return sets, seq_of


def main() -> None:
    cfg = parse_args(__doc__)
    paths, ccfg = cfg["paths"], cfg["clustering"]
    work = Path(paths["reports_dir"]) / "14_search"
    sets, seq_of = evaluation_sets(cfg)
    query_uids = {u for _, pairs, _, _, _, _ in sets for p in pairs for u in p}
    ref_uids = {u for _, _, ref, _, _, _ in sets for p in ref for u in p}
    log.info("Searching %d evaluation chains against %d reference chains (both directions)",
             len(query_uids), len(ref_uids))
    q = {u: seq_of[u] for u in query_uids}
    t = {u: seq_of[u] for u in ref_uids}
    fwd = search(q, t, work, ccfg, int(cfg["threads"]))
    rev = search(t, q, work, ccfg, int(cfg["threads"])).rename(columns={"query": "target", "target": "query"})
    both = pd.concat([fwd, rev], ignore_index=True)
    hit_sets = {
        "step09": both[both["linked"]].groupby(["query", "target"], as_index=False)["fident"].max(),
        "no_coverage": both.groupby(["query", "target"], as_index=False)["fident"].max(),
    }

    rows, summary, summary_nocov = [], {}, {}
    for (name, pairs, ref, level, design_thr, only), (criterion, hits) in (
            (s, h) for h in hit_sets.items() for s in sets):
        ref_uids_set = {u for p in ref for u in p}
        pairs_of: Dict[str, Set[int]] = defaultdict(set)
        for i, (a, b) in enumerate(ref):
            pairs_of[a].add(i)
            pairs_of[b].add(i)
        chains = {u for p in pairs for u in p} if only is None else set(only)
        sub = hits[hits["query"].isin(chains) & hits["target"].isin(ref_uids_set)]
        for thr in sorted(ccfg["thresholds"], reverse=True):
            homologues: Dict[str, Set[str]] = defaultdict(set)
            for qq, tt in zip(*[sub.loc[sub["fident"] >= thr, c] for c in ("query", "target")]):
                homologues[qq].add(tt)

            def pair_leaks(a: str, b: str) -> bool:
                ha, hb = homologues.get(a, set()), homologues.get(b, set())
                if not ha or not hb:
                    return False
                candidates = set().union(*(pairs_of[x] for x in ha)) & set().union(*(pairs_of[y] for y in hb))
                return any((x in ha and y in hb) or (y in ha and x in hb) for x, y in (ref[i] for i in candidates))

            leak_chain = {u for u in chains if homologues.get(u)}
            rec = {"set": name, "criterion": criterion, "threshold": thr, "rows": len(pairs), "chains": len(chains),
                   "chains_with_homologue": len(leak_chain),
                   "rows_any_chain": int(sum(a in leak_chain or b in leak_chain for a, b in pairs)),
                   "rows_pair_homologue": int(sum(pair_leaks(a, b) for a, b in pairs))}
            rows.append(rec)
            if design_thr is not None and abs(thr - design_thr) < 1e-9:
                key = "chains_with_homologue" if level == "chain" else "rows_pair_homologue"
                (summary if criterion == "step09" else summary_nocov)[name] = {
                    "design_level": level, "design_threshold": thr,
                    "leaks": rec[key], "of": rec["chains"] if level == "chain" else rec["rows"]}

    pd.DataFrame(rows).to_csv(ensure_parent(Path(paths["reports_dir"]) / "14_leakage_by_threshold.tsv"),
                              sep="\t", index=False)
    write_json(paths["leakage_audit"], {
        "criteria": {k: ccfg[k] for k in ("coverage", "cov_mode", "evalue", "search_sensitivity", "max_seqs")},
        "design_threshold_leaks": summary,
        "design_threshold_leaks_without_coverage_rule": summary_nocov})
    for name, s in summary.items():
        log.info("%-36s %s level at %.2f: %d of %d leak (without coverage rule: %d)", name, s["design_level"],
                 s["design_threshold"], s["leaks"], s["of"], summary_nocov.get(name, {}).get("leaks", -1))


if __name__ == "__main__":
    main()
