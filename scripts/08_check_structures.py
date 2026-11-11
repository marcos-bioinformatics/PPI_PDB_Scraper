#!/usr/bin/env python3
"""Step 08 - Structure checks on the quality-filtered pairs.

Input : paths.quality_filtered (step 07), assembly mmCIF files (step 04)
Output: paths.checked                              rows that pass both checks
        reports_dir/08_chain_summary.csv           one row per scanned chain
        reports_dir/08_rows_excluded.csv           excluded rows and reason
        reports_dir/08_structure_checks.json       counts

For every pair, in its assembly file (first model, author chain IDs as in PISA,
polymer residues only):
  1. Special residues: the 20 standard amino acids are always allowed;
     special_residues.allowed (MSE, SEC by default) are tolerated and
     counted, special_residues.ignore (terminal caps ACE, NH2) are skipped; any
     other polymer residue (UNK, D-amino acids, chromophores, PTMs, ...) excludes
     the chain. Ligands that share the chain ID are not considered.
  2. Contacts: heavy-atom contacts between the two chains within
     contacts.cutoff A. A pair needs at least contacts.min_residue_contacts
     residue pairs in contact; a PISA interface without real contacts points
     to a chain-mapping problem.
  3. Clashes: n_clash_atoms (atoms of either chain closer than
     contacts.clash_distance A to the other chain) and max_clash_fraction (the
     larger fraction of the two chains). Pairs above contacts.max_clash_fraction
     are excluded (AlphaFold 3 removes chains with more than 30 % of their atoms
     within 1.7 A of another chain); null = report only.
  4. Interface gaps (molA_/molB_interface_gaps): unresolved SEQRES residues
     within contacts.interface_gap_window positions of an interface residue.
     Reported for every pair; pairs with more than contacts.max_interface_gaps
     on either side are excluded only if that value is set (null = report only).
"""

from __future__ import annotations

from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
from tqdm import tqdm

from common import ensure_parent, get_logger, parse_args, pisa_subdir, write_json
from structure_checks import (clash_summary, contact_summary, fmt_counter, interface_gaps, load_chains,
                              special_residues)

log = get_logger("08_check_structures")


def check_assembly(task: Dict[str, Any]) -> Dict[str, List[Dict[str, Any]]]:
    base = {"pdb_id": task["pdb_id"], "assembly_id": task["assembly_id"]}
    try:
        model = load_chains(task["path"], task["chains"])
    except Exception as exc:  # unreadable file: every row of this assembly is excluded
        err = f"parse_error:{type(exc).__name__}"
        return {"chains": [{**base, "chain_id": c, "keep_chain": False, "error": err, "allowed_special_count": 0,
                            "allowed_special_types": "", "forbidden_special_types": ""} for c in task["chains"]],
                "pairs": []}

    chains = []
    for chain_id in task["chains"]:
        if chain_id not in model:
            chains.append({**base, "chain_id": chain_id, "keep_chain": False, "error": "chain_not_found",
                           "allowed_special_count": 0, "allowed_special_types": "", "forbidden_special_types": ""})
            continue
        ok, bad = special_residues(model[chain_id], task["allowed"], task["ignore"])
        chains.append({**base, "chain_id": chain_id, "keep_chain": not bad, "error": "",
                       "allowed_special_count": sum(ok.values()),
                       "allowed_special_types": fmt_counter(ok), "forbidden_special_types": fmt_counter(bad)})

    pairs = []
    for a, b in task["pairs"]:
        if a in model and b in model:
            pairs.append({**base, "molA_chain_id": a, "molB_chain_id": b,
                          **contact_summary(model[a], model[b], task["cutoff"]),
                          **interface_gaps(model[a], model[b], task["cutoff"], task["gap_window"]),
                          **clash_summary(model[a], model[b], task["clash_distance"])})
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
    ignore = {x.upper() for x in cfg["special_residues"].get("ignore", [])}
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
                      "pairs": sorted(pairs[(pid, asm)]), "allowed": allowed, "ignore": ignore,
                      "cutoff": float(contacts["cutoff"]), "gap_window": int(contacts["interface_gap_window"]),
                      "clash_distance": float(contacts["clash_distance"])})
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
            f"_{side}_err": s["error"],
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
    structure_error = (out["_molA_err"].fillna("missing") != "") | (out["_molB_err"].fillna("missing") != "")
    contacts_ok = out["n_residue_contacts"].fillna(0) >= int(contacts["min_residue_contacts"])
    max_gaps = contacts.get("max_interface_gaps")
    gaps_ok = (pd.Series(True, index=out.index) if max_gaps is None else
               (out[["molA_interface_gaps", "molB_interface_gaps"]].fillna(0).max(axis=1) <= int(max_gaps)))
    max_clash = contacts.get("max_clash_fraction")
    clash_ok = (pd.Series(True, index=out.index) if max_clash is None else
                out["max_clash_fraction"].fillna(0) <= float(max_clash))
    keep = special_ok & contacts_ok & clash_ok & gaps_ok
    excluded = out.loc[~keep, keys + ["molA_chain_id", "molB_chain_id", "molA_special_error",
                                      "molB_special_error", "n_residue_contacts"]].copy()
    excluded["reason"] = ["structure_error" if e else "special_residues" if not s else
                          "no_contacts" if not c else "clashes" if not k else "interface_gaps"
                          for e, s, c, k in zip(structure_error[~keep], special_ok[~keep], contacts_ok[~keep],
                                                clash_ok[~keep])]
    excluded.to_csv(ensure_parent(reports / "08_rows_excluded.csv"), index=False)

    out = out.loc[keep].drop(columns=[c for c in out.columns if c.startswith("_")])
    dest = ensure_parent(paths["checked"])
    tmp = dest.with_name(dest.name + ".part")
    out.to_csv(tmp, index=False)
    tmp.replace(dest)

    write_json(reports / "08_structure_checks.json", {
        "rows_in": int(len(df)), "rows_kept": int(len(out)),
        "excluded_structure_errors": int(structure_error.sum()),
        "excluded_special_residues": int((~special_ok & ~structure_error).sum()),
        "excluded_no_contacts": int((special_ok & ~contacts_ok).sum()),
        "excluded_clashes": int((special_ok & contacts_ok & ~clash_ok).sum()),
        "pairs_with_clashes": int((out["n_clash_atoms"].fillna(0) > 0).sum()),
        "clash_distance": contacts["clash_distance"], "max_clash_fraction": max_clash,
        "excluded_interface_gaps": int((special_ok & contacts_ok & clash_ok & ~gaps_ok).sum()),
        "pairs_with_interface_gaps": int((out[["molA_interface_gaps", "molB_interface_gaps"]].fillna(0)
                                          .max(axis=1) > 0).sum()),
        "interface_gap_window": contacts["interface_gap_window"], "max_interface_gaps": max_gaps,
        "chains_scanned": int(len(summary)), "chains_with_forbidden_residues": int((~summary["keep_chain"]).sum()),
        "allowed_special": sorted(allowed), "ignored": sorted(ignore), "contact_cutoff": contacts["cutoff"],
        "min_residue_contacts": contacts["min_residue_contacts"],
    })
    log.info("Kept %d of %d rows (%d special residues, %d without contacts)", len(out), len(df),
             (~special_ok).sum(), (special_ok & ~contacts_ok).sum())


if __name__ == "__main__":
    main()
