#!/usr/bin/env python3
"""Step 02 - Parse PISA JSON into one row per protein-protein chain pair.

Input : paths.pisa_manifest, paths.pisa_dir
Output: paths.master_pisa   (one row per pair of protein molecules in an interface)

An interface is kept as a candidate when interface_area >= pisa.min_interface_area,
p_value <= pisa.max_p_value and it involves at least two protein molecules. With
pisa.expand_multiprotein, interfaces with more than two proteins give one row
per protein pair. Assembly-level metrics come from {id}_assembly{n}.json.gz.
"""

from __future__ import annotations

import csv
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
from tqdm import tqdm

from common import ensure_parent, get_logger, parse_args, pisa_subdir, read_json

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
    ["pdb_id", "assembly_id", "interface_id", "pair_index_in_interface", "pisa_interacts",
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


def normalize_p_value(p: Optional[float]) -> Optional[float]:
    """Some PISA records store the p-value scaled by 1000 (565 -> 0.565)."""
    while p is not None and p > 1.0:
        p /= 1000.0
    return p


def assembly_block(path: Path, pdb_id: str) -> Optional[Dict[str, Any]]:
    if not path.is_file():
        return None
    data = read_json(path)
    data = data.get(pdb_id, data) if isinstance(data, dict) else data
    block = data.get("assembly") if isinstance(data, dict) else None
    return block if isinstance(block, dict) else None


def parse_assembly(task: Tuple[str, int, str, str, float, float, bool]) -> Tuple[List[Dict[str, Any]], Optional[str]]:
    pdb_id, asm, itf_path, asm_path, min_area, max_p, expand = task
    try:
        itf_block = assembly_block(Path(itf_path), pdb_id)
        if itf_block is None:
            return [], f"{itf_path}\tno assembly block"
        metrics = assembly_block(Path(asm_path), pdb_id) or itf_block
        asm_cols = {f"assembly_{k}": metrics.get(k) for k in ASSEMBLY_KEYS}

        rows: List[Dict[str, Any]] = []
        for itf in itf_block.get("interfaces") or []:
            area = to_float(itf.get("interface_area"))
            p_value = normalize_p_value(to_float(itf.get("p_value")))
            if area is None or area < min_area or p_value is None or p_value > max_p:
                continue
            molecules = [m for m in itf.get("molecules") or [] if isinstance(m, dict)]
            proteins = [m for m in molecules if str(m.get("molecule_class", "")).strip().lower() == "protein"]
            if len(proteins) < 2:
                continue
            pairs = list(combinations(proteins, 2)) if expand else [(proteins[0], proteins[1])]

            base = {
                "pdb_id": pdb_id,
                "assembly_id": asm,
                "interface_id": itf.get("interface_id"),
                "pisa_interacts": 1,
                "interface_area": area,
                "p_value": p_value,
                "solvation_energy": to_float(itf.get("solvation_energy")),
                "stabilization_energy": to_float(itf.get("stabilization_energy")),
                **{k: to_int(itf.get(k)) for k in INTERFACE_INT_KEYS},
                "protein_molecule_count": len(proteins),
                "other_molecule_count": len(molecules) - len(proteins),
                **asm_cols,
                "source_file": itf_path,
            }
            for idx, (a, b) in enumerate(pairs, start=1):
                rows.append({
                    **base,
                    "pair_index_in_interface": idx,
                    "molA_molecule_id": a.get("molecule_id"),
                    "molA_molecule_class": a.get("molecule_class"),
                    "molA_chain_id": a.get("chain_id"),
                    "molB_molecule_id": b.get("molecule_id"),
                    "molB_molecule_class": b.get("molecule_class"),
                    "molB_chain_id": b.get("chain_id"),
                })
        return rows, None
    except (OSError, EOFError, ValueError, AttributeError) as exc:
        return [], f"{itf_path}\t{type(exc).__name__}: {exc}"


def main() -> None:
    cfg = parse_args(__doc__)
    paths, pisa = cfg["paths"], cfg["pisa"]
    pisa_dir = Path(paths["pisa_dir"])

    manifest = pd.read_csv(paths["pisa_manifest"], sep="\t", dtype={"pdb_id": str})
    manifest = manifest[manifest["status"] == "ok"]
    tasks = []
    for pid, asm in zip(manifest["pdb_id"], manifest["assembly_id"].astype(int)):
        folder = pisa_dir / pisa_subdir(pid)
        tasks.append((pid, asm,
                      str(folder / f"{pid}_assembly{asm}_interfaces.json.gz"),
                      str(folder / f"{pid}_assembly{asm}.json.gz"),
                      float(pisa["min_interface_area"]), float(pisa["max_p_value"]),
                      bool(pisa["expand_multiprotein"])))
    log.info("Parsing %d assemblies", len(tasks))

    out = ensure_parent(paths["master_pisa"])
    tmp = out.with_name(out.name + ".part")
    errors_path = Path(paths["reports_dir"]) / "02_parse_errors.tsv"
    n_rows, errors = 0, []
    with tmp.open("w", newline="", encoding="utf-8") as fh, \
            ProcessPoolExecutor(max_workers=int(cfg["workers"])) as ex:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        # map() keeps the input order, so the output is deterministic
        for rows, err in tqdm(ex.map(parse_assembly, tasks, chunksize=200), total=len(tasks),
                              desc="Parse PISA", unit="asm"):
            if err:
                errors.append(err)
            writer.writerows(rows)
            n_rows += len(rows)
    tmp.replace(out)

    ensure_parent(errors_path).write_text("".join(e + "\n" for e in errors), encoding="utf-8")
    log.info("%d candidate pairs written to %s (%d files with errors)", n_rows, out, len(errors))


if __name__ == "__main__":
    main()
