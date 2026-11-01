#!/usr/bin/env python3
"""Sanity checks on the outputs of a finished run.

    python tests/check_outputs.py --config config/test.yaml

Checks that the final tables exist and are consistent with the config: no row
in both test and training, no duplicated rows, chains within the length limits,
no hard negative in contact, and zero homology leaks at every design threshold
of the leakage audit. Exits with status 1 if any check fails.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from common import parse_args  # noqa: E402

FAILED = []


def check(ok: bool, what: str) -> None:
    print(("PASS  " if ok else "FAIL  ") + what)
    if not ok:
        FAILED.append(what)


def main() -> None:
    cfg = parse_args(__doc__)
    paths, qf = cfg["paths"], cfg["quality_filter"]

    test = pd.read_csv(paths["test"], low_memory=False)
    tv = pd.read_csv(paths["trainval_splits"], low_memory=False)
    val = pd.read_csv(paths["validation"], low_memory=False)
    check(len(test) > 0 and len(tv) > 0 and len(val) > 0, "test, trainval and validation are not empty")
    check(test["row_id"].is_unique and tv["row_id"].is_unique, "no duplicated row_id")
    check(not set(test["row_id"]) & set(tv["row_id"]), "no row in both test and trainval")
    split = cfg["splits"]["validation_split"]
    check(set(val["row_id"]) == set(tv.loc[tv[split] == "val", "row_id"]), f"validation.csv == {split} 'val' rows")

    for name, df in (("test", test), ("trainval", tv)):
        lengths = pd.concat([df["molA_seq_len"], df["molB_seq_len"]])
        check(lengths.min() > qf["min_chain_length"] and lengths.max() <= qf["max_chain_length"],
              f"{name}: chain lengths within ({qf['min_chain_length']}, {qf['max_chain_length']}]")
        check(df[["molA_seq", "molB_seq", "pair_cluster_90"]].notna().all().all(),
              f"{name}: sequences and clusters present")
        check((df["n_residue_contacts"] >= cfg["contacts"]["min_residue_contacts"]).all(),
              f"{name}: every pair has atomic contacts")

    if cfg["hard_negatives"]["enabled"]:
        neg = pd.read_csv(paths["hard_negatives"], low_memory=False)
        check((neg["label"] == 0).all(), "hard negatives: label 0")
        check((neg["min_interchain_distance"] >= cfg["hard_negatives"]["min_distance"]).all(),
              "hard negatives: chains farther apart than min_distance")
        schemes = [c for c in neg.columns if c.startswith("split_")]
        check(all(neg[c].isin(["train", "val", "test"]).all() for c in schemes), "hard negatives: valid splits")

    audit = json.loads(Path(paths["leakage_audit"]).read_text())["design_threshold_leaks"]
    for name, s in audit.items():
        check(s["leaks"] == 0, f"leakage {name} at {s['design_threshold']:.2f} ({s['design_level']}): "
                               f"{s['leaks']} of {s['of']}")

    print(f"\n{'All checks passed' if not FAILED else f'{len(FAILED)} checks failed'}")
    sys.exit(1 if FAILED else 0)


if __name__ == "__main__":
    main()
