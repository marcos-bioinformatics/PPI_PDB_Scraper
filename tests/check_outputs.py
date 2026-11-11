#!/usr/bin/env python3
"""Sanity checks on the outputs of a finished run.

    python tests/check_outputs.py --config config/test.yaml

Checks that the final tables exist and are consistent with the config: no row
in both test and training, no duplicated rows, chains within the length limits,
no hard negative in contact or equal to a positive pair, no file dropped by the
parsers beyond max_failed_fraction, no failed IntAct query, every audited set
present, and zero homology leaks at every design threshold of the leakage
audit. Exits with status 1 if any check fails.
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
    reports = Path(paths["reports_dir"])
    level = f"c{round(100 * float(cfg['test_set']['level']))}"

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
        check(df[["molA_seq", "molB_seq", f"pair_cluster_{level}"]].notna().all().all(),
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
        key = lambda d: {(str(p).lower(), str(a), *sorted((str(x), str(y))))  # noqa: E731
                         for p, a, x, y in zip(d["pdb_id"], d["assembly_id"], d["molA_chain_id"], d["molB_chain_id"])}
        check(len(key(neg)) == len(neg), "hard negatives: no duplicated chain pair")
        check(not key(neg) & (key(test) | key(tv)), "hard negatives: no negative is also a positive pair")
        if (cfg["hard_negatives"].get("intact") or {}).get("enabled"):
            n13 = json.loads((reports / "13_hard_negatives.json").read_text())["candidates_left_after_each_filter"]
            check(n13.get("intact_query_failures", 0) == 0,
                  f"hard negatives: no failed IntAct query ({n13.get('intact_query_failures')} failed; rerun step 13 if any)")

    # files the parsers could not read are dropped with a report: they must stay below max_failed_fraction
    maxf = float(cfg.get("max_failed_fraction", 0.001))
    p02 = json.loads((reports / "02_parse_summary.json").read_text())
    check(p02["files_with_errors"] <= maxf * max(p02["assemblies"], 1),
          f"step 02: PISA files with errors {p02['files_with_errors']} of {p02['assemblies']} (max {maxf:g})")
    e06 = reports / "06_mmcif_parse_errors.tsv"
    n06 = sum(1 for line in e06.read_text().splitlines() if line.strip()) if e06.is_file() else 0
    n_files = max(sum(1 for _ in open(paths["structures_manifest"])) - 1, 1)
    check(n06 <= maxf * n_files, f"step 06: mmCIF files with errors {n06} of {n_files} (max {maxf:g})")

    if cfg.get("prodigy", {}).get("filter"):
        band = pd.read_csv(paths["prodigy_band"], dtype=str)
        final = set(test["row_id"]) | set(tv["row_id"])
        xtal = set(band.loc[band["prodigy_class"] == "XTAL", "row_id"])
        bio = set(band.loc[band["prodigy_class"] == "BIO", "row_id"])
        check(not xtal & final and len(bio & final) > 0,
              f"PRODIGY-cryst filter: no XTAL pair below {cfg['prodigy']['max_area']} A2 in test or trainval")

    audit = json.loads(Path(paths["leakage_audit"]).read_text())["design_threshold_leaks"]
    expected = set(test["test_bucket"].unique()) | {f"validation_{c}" for c in tv.columns if c.startswith("split_")}
    if "C2" in expected:
        expected.add("C2_unseen_chain")
    if cfg["hard_negatives"]["enabled"]:
        expected |= {f"{b}_vs_negatives" for b in set(expected) if b.startswith("C")}
        if str(cfg["hard_negatives"]["conflict_level"]).startswith("c"):
            expected.add("hard_negatives_test")
    check(expected <= set(audit), f"leakage audit: every set present (missing: {sorted(expected - set(audit)) or 'none'})")
    for name, s in audit.items():
        check(s["leaks"] == 0, f"leakage {name} at {s['design_threshold']:.2f} ({s['design_level']}): "
                               f"{s['leaks']} of {s['of']}")

    print(f"\n{'All checks passed' if not FAILED else f'{len(FAILED)} checks failed'}")
    sys.exit(1 if FAILED else 0)


if __name__ == "__main__":
    main()
