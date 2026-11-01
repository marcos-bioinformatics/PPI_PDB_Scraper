#!/usr/bin/env python3
"""Step 14 - Audit homology leakage between the training set and the test / validation sets.

Input : paths.test (11), paths.trainval_splits (12)
Output: paths.leakage_audit                     summary per set at its design threshold
        reports_dir/14_leakage_by_threshold.tsv  every set at every clustering threshold

Every chain of the test buckets and of the validation set is searched with
MMseqs2 against the chains of the training set, with the identity and coverage
criteria of step 09. For each set and threshold it reports:
  chains_with_homologue  chains with a training chain at >= threshold identity
  rows_any_chain         pairs with at least one such chain
  rows_pair_homologue    pairs for which a training pair links homologues of both chains
The design threshold of each set is where leakage should be zero:
  unseen_chain_c40 / _c50  chain level at 0.40 / 0.50 (C3: proteins unseen in training)
  c70 / c80                pair level at 0.70 / 0.80 (C2: new pairs of known proteins)
  validation               pair level at the identity of the split column (split_c90 -> 0.90)
The audit only reports; it does not change any dataset.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pandas as pd

from common import ensure_parent, get_logger, parse_args, seq_uid, write_json

log = get_logger("14_leakage_audit")

DESIGN = {"unseen_chain_c40": ("chain", 0.40), "unseen_chain_c50": ("chain", 0.50),
          "c70": ("pair", 0.70), "c80": ("pair", 0.80)}


def search(query: Dict[str, str], target: Dict[str, str], work: Path, ccfg: Dict, threads: int) -> pd.DataFrame:
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    for name, seqs in (("q.fa", query), ("t.fa", target)):
        (work / name).write_text("".join(f">{u}\n{s}\n" for u, s in seqs.items()), encoding="utf-8")
    cmd = ["mmseqs", "easy-search", str(work / "q.fa"), str(work / "t.fa"), str(work / "hits.m8"), str(work / "tmp"),
           "--min-seq-id", str(min(ccfg["thresholds"])), "-c", "0", "-s", str(ccfg["search_sensitivity"]),
           "--max-seqs", str(ccfg["max_seqs"]), "--threads", str(threads),
           "--format-output", "query,target,fident,qcov,tcov"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"MMseqs2 search failed:\n{p.stderr}")
    hits = pd.read_csv(work / "hits.m8", sep="\t", header=None, dtype={0: str, 1: str},
                       names=["query", "target", "fident", "qcov", "tcov"])
    shutil.rmtree(work)
    c = float(ccfg["coverage"])
    # same link criterion as step 09: both coverages for cov_mode 0, either one otherwise
    covered = (hits["qcov"] >= c) & (hits["tcov"] >= c) if int(ccfg["cov_mode"]) == 0 else \
              (hits["qcov"] >= c) | (hits["tcov"] >= c)
    return hits[covered]


def main() -> None:
    cfg = parse_args(__doc__)
    paths, ccfg = cfg["paths"], cfg["clustering"]
    split_col = cfg["splits"]["validation_split"]

    test = pd.read_csv(paths["test"], keep_default_na=False, na_values=[""], low_memory=False)
    tv = pd.read_csv(paths["trainval_splits"], keep_default_na=False, na_values=[""], low_memory=False)
    for df in (test, tv):
        df["uidA"], df["uidB"] = df["molA_seq"].map(seq_uid), df["molB_seq"].map(seq_uid)
    train, val = tv[tv[split_col] == "train"], tv[tv[split_col] == "val"]

    sets = {b: g for b, g in test.groupby("test_subbucket")}
    sets["validation"] = val
    m = re.search(r"(\d+)$", split_col)
    design = dict(DESIGN, validation=("pair", int(m.group(1)) / 100 if m else 0.90))

    seq_of = dict(zip(pd.concat([tv["uidA"], tv["uidB"], test["uidA"], test["uidB"]]),
                      pd.concat([tv["molA_seq"], tv["molB_seq"], test["molA_seq"], test["molB_seq"]])))
    train_uids = set(train["uidA"]) | set(train["uidB"])
    query_uids = set().union(*({*g["uidA"], *g["uidB"]} for g in sets.values()))
    log.info("Searching %d test/validation chains against %d training chains", len(query_uids), len(train_uids))
    hits = search({u: seq_of[u] for u in query_uids}, {u: seq_of[u] for u in train_uids},
                  Path(paths["reports_dir"]) / "14_search", ccfg, int(cfg["threads"]))

    train_pairs: List[Tuple[str, str]] = list(zip(train["uidA"], train["uidB"]))
    pairs_of: Dict[str, Set[int]] = defaultdict(set)
    for i, (a, b) in enumerate(train_pairs):
        pairs_of[a].add(i)
        pairs_of[b].add(i)

    rows, summary = [], {}
    for thr in sorted(ccfg["thresholds"], reverse=True):
        homologues: Dict[str, Set[str]] = defaultdict(set)
        for q, t in zip(*[hits.loc[hits["fident"] >= thr, c] for c in ("query", "target")]):
            homologues[q].add(t)

        def pair_leaks(a: str, b: str) -> bool:
            ha, hb = homologues.get(a, set()), homologues.get(b, set())
            candidates = set().union(*(pairs_of[x] for x in ha)) & set().union(*(pairs_of[y] for y in hb)) \
                if ha and hb else set()
            return any((x in ha and y in hb) or (y in ha and x in hb)
                       for x, y in (train_pairs[i] for i in candidates))

        for name, g in sets.items():
            chains = set(g["uidA"]) | set(g["uidB"])
            leak_chain = {u for u in chains if homologues.get(u)}
            rec = {
                "set": name, "threshold": thr, "rows": len(g), "chains": len(chains),
                "chains_with_homologue": len(leak_chain),
                "rows_any_chain": int(sum(a in leak_chain or b in leak_chain for a, b in zip(g["uidA"], g["uidB"]))),
                "rows_pair_homologue": int(sum(pair_leaks(a, b) for a, b in zip(g["uidA"], g["uidB"]))),
            }
            rows.append(rec)
            level, dthr = design.get(name, (None, None))
            if dthr is not None and abs(thr - dthr) < 1e-9:
                key = "chains_with_homologue" if level == "chain" else "rows_pair_homologue"
                summary[name] = {"design_level": level, "design_threshold": thr,
                                 "leaks": rec[key], "of": rec["chains"] if level == "chain" else rec["rows"]}

    table = pd.DataFrame(rows)
    table.to_csv(ensure_parent(Path(paths["reports_dir"]) / "14_leakage_by_threshold.tsv"), sep="\t", index=False)
    write_json(paths["leakage_audit"], {"criteria": {k: ccfg[k] for k in ("coverage", "cov_mode", "search_sensitivity")},
                                        "design_threshold_leaks": summary})
    for name, s in summary.items():
        log.info("%-18s %s level at %.2f: %d of %d leak", name, s["design_level"], s["design_threshold"], s["leaks"], s["of"])


if __name__ == "__main__":
    main()
