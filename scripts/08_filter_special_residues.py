#!/usr/bin/env python3
"""Step 08 - Drop pairs whose chains contain non-standard amino acids.

Input : paths.quality_filtered (step 07), assembly mmCIF files (step 04)
Output: paths.special_residues                         rows where both chains pass
        reports_dir/08_special_chain_summary.csv       one row per scanned chain
        reports_dir/08_special_rows_excluded.csv       excluded rows and reason

Every observed amino-acid residue of each chain (Biopython MMCIFParser, first
model, author chain IDs as in PISA) is checked: the 20 standard residues are
always allowed, special_residues.allowed (MSE, SEM, SEC by default) are
tolerated and counted; any other amino-acid residue excludes the chain.
"""

from __future__ import annotations

import warnings
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
from Bio import BiopythonWarning
from Bio.PDB import MMCIFParser, is_aa
from tqdm import tqdm

from common import ensure_parent, get_logger, open_mmcif, parse_args, pisa_subdir, write_json

warnings.simplefilter("ignore", BiopythonWarning)
log = get_logger("08_filter_special_residues")

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


def fmt(counter: Counter) -> str:
    return ";".join(f"{k}:{counter[k]}" for k in sorted(counter))


def scan_file(task: Dict[str, Any]) -> List[Dict[str, Any]]:
    base = {"pdb_id": task["pdb_id"], "assembly_id": task["assembly_id"], "source_file": task["path"]}
    try:
        with open_mmcif(task["path"]) as fh:
            model = MMCIFParser(QUIET=True).get_structure("X", fh)[0]
    except Exception as exc:  # unreadable file: every requested chain is excluded
        return [{**base, "chain_id": c, "keep_chain": False, "error": f"parse_error:{type(exc).__name__}",
                 "allowed_special_count": 0, "allowed_special_types": "", "forbidden_special_types": ""}
                for c in task["chains"]]

    rows = []
    for chain_id in task["chains"]:
        if chain_id not in model:
            rows.append({**base, "chain_id": chain_id, "keep_chain": False, "error": "chain_not_found",
                         "allowed_special_count": 0, "allowed_special_types": "", "forbidden_special_types": ""})
            continue
        allowed, forbidden = Counter(), Counter()
        for res in model[chain_id].get_residues():
            if not is_aa(res, standard=False):
                continue
            name = res.get_resname().upper().strip()
            if name in STANDARD_AA:
                continue
            (allowed if name in task["allowed"] else forbidden)[name] += 1
        rows.append({**base, "chain_id": chain_id, "keep_chain": not forbidden, "error": "",
                     "allowed_special_count": sum(allowed.values()),
                     "allowed_special_types": fmt(allowed), "forbidden_special_types": fmt(forbidden)})
    return rows


def main() -> None:
    cfg = parse_args(__doc__)
    paths = cfg["paths"]
    allowed = {x.upper() for x in cfg["special_residues"]["allowed"]}
    reports = Path(paths["reports_dir"])

    df = pd.read_csv(paths["quality_filtered"], keep_default_na=False, na_values=[""],
                     dtype={"pdb_id": str, "molA_chain_id": str, "molB_chain_id": str})
    df = df.loc[:, ~df.columns.str.startswith(("Unnamed", "X."))]
    df["assembly_id"] = df["assembly_id"].astype(int)

    chains = defaultdict(set)
    for pid, asm, a, b in df[["pdb_id", "assembly_id", "molA_chain_id", "molB_chain_id"]].itertuples(index=False):
        chains[(pid, asm)].update((a, b))

    asm_dir, asu_dir = Path(paths["assembly_dir"]), Path(paths["asu_dir"])
    tasks = []
    for (pid, asm), chain_ids in sorted(chains.items()):
        path = asm_dir / pisa_subdir(pid) / f"{pid}-assembly{asm}.cif.gz"
        if not path.is_file() and asm == 1:  # assembly 1 can fall back to the ASU
            path = asu_dir / pisa_subdir(pid) / f"{pid}.cif.gz"
        tasks.append({"pdb_id": pid, "assembly_id": asm, "path": str(path),
                      "chains": sorted(chain_ids), "allowed": allowed})
    log.info("Scanning %d structures for %d rows", len(tasks), len(df))

    rows = []
    with ProcessPoolExecutor(max_workers=int(cfg["workers"])) as ex:
        for result in tqdm(ex.map(scan_file, tasks, chunksize=20), total=len(tasks), desc="Special residues"):
            rows.extend(result)
    summary = pd.DataFrame(rows)
    summary.to_csv(ensure_parent(reports / "08_special_chain_summary.csv"), index=False)

    # Per-side columns (same names as the original pipeline)
    keys = ["pdb_id", "assembly_id"]
    out = df
    for side in ("molA", "molB"):
        s = summary.rename(columns={"chain_id": f"{side}_chain_id"})
        s = s.assign(**{
            f"{side}_has_special": s["allowed_special_count"] > 0,
            f"{side}_number_special": s["allowed_special_count"],
            f"{side}_has_mse_sec": s["allowed_special_count"] > 0,
            f"{side}_special_residues": s["allowed_special_types"],
            f"{side}_special_status": s["keep_chain"].map({True: "ok", False: "forbidden_or_error"}),
            f"{side}_special_error": s["error"].where(s["error"] != "", s["forbidden_special_types"]),
            f"_{side}_keep": s["keep_chain"],
            f"_{side}_types": s["allowed_special_types"],
        })
        cols = keys + [f"{side}_chain_id"] + [c for c in s.columns if c.startswith((f"{side}_", f"_{side}_"))
                                              and c != f"{side}_chain_id"]
        out = out.merge(s[cols], on=keys + [f"{side}_chain_id"], how="left", validate="many_to_one")

    def combine(a: str, b: str) -> Counter:
        total = Counter()
        for part in f"{a};{b}".split(";"):
            if part:
                k, v = part.split(":")
                total[k] += int(v)
        return total

    combined = [combine(a, b) for a, b in zip(out["_molA_types"].fillna(""), out["_molB_types"].fillna(""))]
    out["has_special"] = [sum(c.values()) > 0 for c in combined]
    out["number_special"] = [sum(c.values()) for c in combined]
    out["has_mse_sec"] = out["has_special"]
    out["special_residues"] = [fmt(c) for c in combined]

    keep = out["_molA_keep"].fillna(False).astype(bool) & out["_molB_keep"].fillna(False).astype(bool)
    excluded = out.loc[~keep, keys + ["molA_chain_id", "molB_chain_id", "molA_special_error", "molB_special_error"]]
    excluded.to_csv(ensure_parent(reports / "08_special_rows_excluded.csv"), index=False)

    out = out.loc[keep].drop(columns=[c for c in out.columns if c.startswith("_")])
    dest = ensure_parent(paths["special_residues"])
    tmp = dest.with_name(dest.name + ".part")
    out.to_csv(tmp, index=False)
    tmp.replace(dest)

    write_json(reports / "08_special_residues_summary.json", {
        "rows_in": int(len(df)), "rows_kept": int(len(out)), "rows_excluded": int((~keep).sum()),
        "chains_scanned": int(len(summary)), "chains_excluded": int((~summary["keep_chain"]).sum()),
        "allowed_special": sorted(allowed),
    })
    log.info("Kept %d of %d rows (%d excluded)", len(out), len(df), (~keep).sum())


if __name__ == "__main__":
    main()
