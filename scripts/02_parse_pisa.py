#!/usr/bin/env python3
"""Step 02 - Parse PISA JSON into one row per protein-protein interface.

Input : paths.pisa_manifest, paths.pisa_dir
Output: paths.master_pisa                    one row per candidate chain pair
        reports_dir/02_parse_summary.json     interfaces kept / dropped by reason
        reports_dir/02_parse_errors.tsv       unreadable files

A PISA interface is always formed by two molecules. It becomes a candidate
pair when both are proteins and
    interface_area >= pisa.min_interface_area   (A^2 buried per side)
    p_value        <= pisa.max_p_value          (PISA Delta-iG P-value, 0-1)
These are permissive pre-filters (they decide which structures are downloaded);
the final thresholds are applied in step 07 and can only be stricter.
Assembly-level metrics come from {id}_assembly{n}.json.gz.
"""

from __future__ import annotations

import csv
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
from tqdm import tqdm

from common import ensure_parent, get_logger, parse_args, pisa_subdir, read_json, write_json

log = get_logger("02_parse_pisa")

INTERFACE_INT_KEYS = [
    "number_interface_residues", "number_hydrogen_bonds", "number_salt_bridges",
    "number_disulfide_bonds", "number_covalent_bonds", "number_other_bonds",
]
ASSEMBLY_KEYS = [
    "dissociation_energy", "solvation_energy_gain", "entropy", "dissociation_area",
    "accessible_surface_area", "buried_surface_area", "formula", "composition", "interface_count",
]
COLUMNS = (
    ["pdb_id", "assembly_id", "interface_id", "pair_index_in_interface",
     "interface_area", "p_value", "solvation_energy", "stabilization_energy"]
    + INTERFACE_INT_KEYS
    + ["protein_molecule_count", "other_molecule_count",
       "molA_molecule_id", "molA_molecule_class", "molA_chain_id",
       "molB_molecule_id", "molB_molecule_class", "molB_chain_id"]
    + [f"assembly_{k}" for k in ASSEMBLY_KEYS]
    + ["source_file"]
)


def to_float(x: Any) -> Optional[float]:
    try:
        return None if x is None else float(x)
    except (TypeError, ValueError):
        return None


def to_int(x: Any) -> Optional[int]:
    try:
        return None if x is None else int(x)
    except (TypeError, ValueError):
        return None


def assembly_record(path: Path, pdb_id: str) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
    """Return (assembly block, pisa_version) of a PISA JSON file."""
    if not path.is_file():
        return None, None
    data = read_json(path)
    data = data.get(pdb_id, data) if isinstance(data, dict) else data
    if not isinstance(data, dict) or not isinstance(data.get("assembly"), dict):
        return None, None
    return data["assembly"], data.get("pisa_version")


def parse_assembly(task) -> Tuple[List[Dict[str, Any]], Counter, Optional[str]]:
    pdb_id, asm, itf_path, asm_path, min_area, max_p = task
    stats: Counter = Counter()
    try:
        itf_block, version = assembly_record(Path(itf_path), pdb_id)
        if itf_block is None:
            return [], stats, f"{itf_path}\tno assembly block"
        stats[f"pisa_version={version}"] += 1
        metrics = assembly_record(Path(asm_path), pdb_id)[0] or itf_block
        asm_cols = {f"assembly_{k}": metrics.get(k) for k in ASSEMBLY_KEYS}

        rows: List[Dict[str, Any]] = []
        for itf in itf_block.get("interfaces") or []:
            stats["interfaces"] += 1
            molecules = [m for m in itf.get("molecules") or [] if isinstance(m, dict)]
            proteins = [m for m in molecules if str(m.get("molecule_class", "")).strip().lower() == "protein"]
            area = to_float(itf.get("interface_area"))
            p_value = to_float(itf.get("p_value"))
            if len(molecules) != 2:
                stats["dropped_not_two_molecules"] += 1
            elif len(proteins) != 2:
                stats["dropped_not_protein_protein"] += 1
            elif area is None or p_value is None or not 0.0 <= p_value <= 1.0:
                stats["dropped_missing_or_invalid_values"] += 1
            elif area < min_area:
                stats["dropped_interface_area"] += 1
            elif p_value > max_p:
                stats["dropped_p_value"] += 1
            else:
                stats["kept"] += 1
                a, b = proteins
                rows.append({
                    "pdb_id": pdb_id,
                    "assembly_id": asm,
                    "interface_id": itf.get("interface_id"),
                    "pair_index_in_interface": 1,
                    "interface_area": area,
                    "p_value": p_value,
                    "solvation_energy": to_float(itf.get("solvation_energy")),
                    "stabilization_energy": to_float(itf.get("stabilization_energy")),
                    **{k: to_int(itf.get(k)) for k in INTERFACE_INT_KEYS},
                    "protein_molecule_count": len(proteins),
                    "other_molecule_count": len(molecules) - len(proteins),
                    "molA_molecule_id": a.get("molecule_id"),
                    "molA_molecule_class": a.get("molecule_class"),
                    "molA_chain_id": a.get("chain_id"),
                    "molB_molecule_id": b.get("molecule_id"),
                    "molB_molecule_class": b.get("molecule_class"),
                    "molB_chain_id": b.get("chain_id"),
                    **asm_cols,
                    "source_file": itf_path,
                })
        return rows, stats, None
    except (OSError, EOFError, ValueError, AttributeError) as exc:
        return [], stats, f"{itf_path}\t{type(exc).__name__}: {exc}"


def main() -> None:
    cfg = parse_args(__doc__)
    paths, pisa = cfg["paths"], cfg["pisa"]
    pisa_dir = Path(paths["pisa_dir"])
    reports = Path(paths["reports_dir"])

    manifest = pd.read_csv(paths["pisa_manifest"], sep="\t", dtype={"pdb_id": str})
    manifest = manifest[manifest["status"] == "ok"]
    tasks = []
    for pid, asm in zip(manifest["pdb_id"], manifest["assembly_id"].astype(int)):
        folder = pisa_dir / pisa_subdir(pid)
        tasks.append((pid, asm,
                      str(folder / f"{pid}_assembly{asm}_interfaces.json.gz"),
                      str(folder / f"{pid}_assembly{asm}.json.gz"),
                      float(pisa["min_interface_area"]), float(pisa["max_p_value"])))
    log.info("Parsing %d assemblies", len(tasks))

    out = ensure_parent(paths["master_pisa"])
    tmp = out.with_name(out.name + ".part")
    totals: Counter = Counter()
    errors: List[str] = []
    with tmp.open("w", newline="", encoding="utf-8") as fh, \
            ProcessPoolExecutor(max_workers=int(cfg["workers"])) as ex:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, lineterminator="\n")
        writer.writeheader()
        # map() keeps the input order, so the output is deterministic
        for rows, stats, err in tqdm(ex.map(parse_assembly, tasks, chunksize=200), total=len(tasks),
                                     desc="Parse PISA", unit="asm"):
            totals.update(stats)
            if err:
                errors.append(err)
            writer.writerows(rows)
    tmp.replace(out)

    ensure_parent(reports / "02_parse_errors.tsv").write_text("".join(e + "\n" for e in errors), encoding="utf-8")
    summary = {"assemblies": len(tasks), "files_with_errors": len(errors),
               "min_interface_area": pisa["min_interface_area"], "max_p_value": pisa["max_p_value"],
               **dict(sorted(totals.items()))}
    write_json(reports / "02_parse_summary.json", summary)
    log.info("%d candidate pairs from %d interfaces written to %s (%d files with errors)",
             totals["kept"], totals["interfaces"], out, len(errors))


if __name__ == "__main__":
    main()
