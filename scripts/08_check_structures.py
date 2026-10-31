#!/usr/bin/env python3
"""Step 08 - Structure checks on the quality-filtered pairs.

Input : paths.quality_filtered (step 07), assembly mmCIF files (step 04)
Output: paths.checked                              rows that pass both checks
        reports_dir/08_chain_summary.csv           one row per scanned chain
        reports_dir/08_rows_excluded.csv           excluded rows and reason
        reports_dir/08_structure_checks.json       counts

For every pair, in its assembly file (first model, author chain IDs as in PISA):
  1. Special residues: the 20 standard amino acids are always allowed;
     special_residues.allowed (MSE, SEM, SEC by default) are tolerated and
     counted; any other amino-acid residue excludes the chain.
  2. Contacts: heavy-atom contacts between the two chains within
     contacts.cutoff A. A pair needs at least contacts.min_residue_contacts
     residue pairs in contact; a PISA interface without real contacts points
     to a chain-mapping problem.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
from tqdm import tqdm

from common import ensure_parent, get_logger, parse_args, pisa_subdir, write_json
from structure_checks import contact_summary, fmt_counter, heavy_atoms, load_model, special_residues

log = get_logger("08_check_structures")


def check_assembly(task: Dict[str, Any]) -> Dict[str, List[Dict[str, Any]]]:
    base = {"pdb_id": task["pdb_id"], "assembly_id": task["assembly_id"]}
    try:
        model = load_model(task["path"])
    except Exception as exc:  # unreadable file: every row of this assembly is excluded
        err = f"parse_error:{type(exc).__name__}"
        return {"chains": [{**base, "chain_id": c, "keep_chain": False, "error": err, "allowed_special_count": 0,
                            "allowed_special_types": "", "forbidden_special_types": ""} for c in task["chains"]],
                "pairs": []}

    chains, atoms = [], {}
    for chain_id in task["chains"]:
        if chain_id not in model:
            chains.append({**base, "chain_id": chain_id, "keep_chain": False, "error": "chain_not_found",
                           "allowed_special_count": 0, "allowed_special_types": "", "forbidden_special_types": ""})
            continue
        ok, bad = special_residues(model[chain_id], task["allowed"])
        chains.append({**base, "chain_id": chain_id, "keep_chain": not bad, "error": "",
                       "allowed_special_count": sum(ok.values()),
                       "allowed_special_types": fmt_counter(ok), "forbidden_special_types": fmt_counter(bad)})
        atoms[chain_id] = heavy_atoms(model[chain_id])

    pairs = []
    for a, b in task["pairs"]:
        if a in atoms and b in atoms:
            pairs.append({**base, "molA_chain_id": a, "molB_chain_id": b,
                          **contact_summary(atoms[a], atoms[b], task["cutoff"])})
    return {"chains": chains, "pairs": pairs}


def combine(a: str, b: str) -> Counter:
    total = Counter()
    for part in f"{a};{b}".split(";"):
        if part:
            k, v = part.split(":")
            total[k] += int(v)
    return total


def main() -> None:
    cfg = parse_args(__doc__)
    paths, contacts = cfg["paths"], cfg["contacts"]
    allowed = {x.upper() for x in cfg["special_residues"]["allowed"]}
    reports = Path(paths["reports_dir"])

    df = pd.read_csv(paths["quality_filtered"], keep_default_na=False, na_values=[""],
                     dtype={"pdb_id": str, "molA_chain_id": str, "molB_chain_id": str})
    df = df.loc[:, ~df.columns.str.startswith(("Unnamed", "X."))]
    df["assembly_id"] = df["assembly_id"].astype(int)

    chains, pairs = defaultdict(set), defaultdict(set)
    for pid, asm, a, b in df[["pdb_id", "assembly_id", "molA_chain_id", "molB_chain_id"]].itertuples(index=False):
        chains[(pid, asm)].update((a, b))
        pairs[(pid, asm)].add((a, b))

    asm_dir = Path(paths["assembly_dir"])
    tasks = []
    for (pid, asm), chain_ids in sorted(chains.items()):
        path = asm_dir / pisa_subdir(pid) / f"{pid}-assembly{asm}.cif.gz"
        tasks.append({"pdb_id": pid, "assembly_id": asm, "path": str(path), "chains": sorted(chain_ids),
                      "pairs": sorted(pairs[(pid, asm)]), "allowed": allowed, "cutoff": float(contacts["cutoff"])})
    log.info("Checking %d structures for %d rows", len(tasks), len(df))

    chain_rows, pair_rows = [], []
    with ProcessPoolExecutor(max_workers=int(cfg["workers"])) as ex:
        for result in tqdm(ex.map(check_assembly, tasks, chunksize=20), total=len(tasks), desc="Structure checks"):
            chain_rows.extend(result["chains"])
            pair_rows.extend(result["pairs"])
    summary = pd.DataFrame(chain_rows)
    summary.to_csv(ensure_parent(reports / "08_chain_summary.csv"), index=False)

    # Special-residue columns per side (same names as the original pipeline)
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
        cols = [c for c in s.columns if c.startswith((f"{side}_", f"_{side}_"))]
        out = out.merge(s[keys + cols], on=keys + [f"{side}_chain_id"], how="left", validate="many_to_one")

    combined = [combine(a, b) for a, b in zip(out["_molA_types"].fillna(""), out["_molB_types"].fillna(""))]
    out["has_special"] = [sum(c.values()) > 0 for c in combined]
    out["number_special"] = [sum(c.values()) for c in combined]
    out["has_mse_sec"] = out["has_special"]
    out["special_residues"] = [fmt_counter(c) for c in combined]

    # Contact columns
    out = out.merge(pd.DataFrame(pair_rows), on=keys + ["molA_chain_id", "molB_chain_id"], how="left",
                    validate="many_to_one")

    special_ok = out["_molA_keep"].fillna(False).astype(bool) & out["_molB_keep"].fillna(False).astype(bool)
    contacts_ok = out["n_residue_contacts"].fillna(0) >= int(contacts["min_residue_contacts"])
    keep = special_ok & contacts_ok
    excluded = out.loc[~keep, keys + ["molA_chain_id", "molB_chain_id", "molA_special_error",
                                      "molB_special_error", "n_residue_contacts"]].copy()
    excluded["reason"] = ["special_residues" if not s else "no_contacts"
                          for s in special_ok[~keep]]
    excluded.to_csv(ensure_parent(reports / "08_rows_excluded.csv"), index=False)

    out = out.loc[keep].drop(columns=[c for c in out.columns if c.startswith("_")])
    dest = ensure_parent(paths["checked"])
    tmp = dest.with_name(dest.name + ".part")
    out.to_csv(tmp, index=False)
    tmp.replace(dest)

    write_json(reports / "08_structure_checks.json", {
        "rows_in": int(len(df)), "rows_kept": int(len(out)),
        "excluded_special_residues": int((~special_ok).sum()),
        "excluded_no_contacts": int((special_ok & ~contacts_ok).sum()),
        "chains_scanned": int(len(summary)), "chains_with_forbidden_residues": int((~summary["keep_chain"]).sum()),
        "allowed_special": sorted(allowed), "contact_cutoff": contacts["cutoff"],
        "min_residue_contacts": contacts["min_residue_contacts"],
    })
    log.info("Kept %d of %d rows (%d special residues, %d without contacts)", len(out), len(df),
             (~special_ok).sum(), (special_ok & ~contacts_ok).sum())


if __name__ == "__main__":
    main()
