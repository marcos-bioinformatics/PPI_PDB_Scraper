#!/usr/bin/env python3
"""Step 06 - Enrich every candidate pair with mmCIF data (entry + chain level).

Input : paths.master_selected (step 03), mmCIF files from step 04
Output: paths.enriched         master rows + sequence, coverage, homo/hetero and entry metadata
        paths.entry_metadata   one row per entry (ASU file: method, resolution, R-free, title)
        paths.chain_metadata   one row per protein chain of every needed assembly

Chains are matched by auth_asym_id: PISA's chain_id is the author chain ID of
the assembly file, including symmetry copies such as "A-2".

Per chain (assembly file, first model):
  seq              canonical sequence of its entity (_entity_poly)
  seqres_count     residues in _pdbx_poly_seq_scheme
  resolved_count   polymer residues with at least one non-hydrogen atom
                   (ATOM or HETATM, so MSE counts as resolved)
  coverage         resolved_count / seqres_count

is_homo = same entity OR positional sequence identity >= enrichment.homo_similarity.
Rows without entry metadata or without both chain sequences are dropped.
Parsed files are cached per chunk in a parts/ folder so an interrupted run resumes.
"""

from __future__ import annotations

import hashlib
import io
import re
import shutil
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from tqdm import tqdm

from common import ensure_parent, get_logger, open_mmcif, parse_args, pisa_subdir

log = get_logger("06_enrich_mmcif")

CHUNK = 500
CHAIN_FIELDS = ["entity_id", "seq", "seq_len", "seqres_count_chain", "resolved_count",
                "missing_residues", "coverage", "coverage_pct"]


# -----------------------------------------------------------------------------
# mmCIF parsing
# -----------------------------------------------------------------------------
def as_list(value: Any) -> List[str]:
    if value is None:
        return []
    return [str(v) for v in value] if isinstance(value, list) else [str(value)]


def first_str(mm: Dict[str, Any], key: str) -> Optional[str]:
    values = as_list(mm.get(key))
    return values[0] if values and values[0] not in (".", "?") else None


def first_float(mm: Dict[str, Any], key: str) -> Optional[float]:
    value = first_str(mm, key)
    try:
        return float(value) if value is not None else None
    except ValueError:
        return None


def read_header(path: str) -> Dict[str, Any]:
    """Parse only the categories before _atom_site (entry-level metadata)."""
    lines: List[str] = []
    with open_mmcif(path) as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                if lines and lines[-1].strip() == "loop_":
                    lines.pop()
                break
            lines.append(line)
    return MMCIF2Dict(io.StringIO("".join(lines)))


def parse_entry(task: Tuple[str, str]) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
    pdb_id, path = task
    try:
        mm = read_header(path)
        return {
            "pdb_id": pdb_id,
            "exptl_method": first_str(mm, "_exptl.method"),
            "resolution": first_float(mm, "_refine.ls_d_res_high"),
            "r_free": first_float(mm, "_refine.ls_R_factor_R_free"),
            "title": first_str(mm, "_struct.title"),
            "keywords": first_str(mm, "_struct_keywords.pdbx_keywords"),
        }, None
    except Exception as exc:  # malformed files are reported, never fatal
        return None, f"{pdb_id}\tasu\t{path}\t{type(exc).__name__}: {exc}"


def parse_assembly(task: Tuple[str, int, str, bool]) -> Tuple[List[Dict[str, Any]], Optional[str]]:
    pdb_id, assembly_id, path, ignore_h = task
    try:
        with open_mmcif(path) as fh:
            mm = MMCIF2Dict(fh)

        seq_by_entity = {}
        types = as_list(mm.get("_entity_poly.type"))
        for i, (ent, seq) in enumerate(zip(as_list(mm.get("_entity_poly.entity_id")),
                                           as_list(mm.get("_entity_poly.pdbx_seq_one_letter_code_can")))):
            if i < len(types) and "polypeptide" not in types[i].lower():
                continue
            seq_by_entity[ent] = re.sub(r"[^A-Za-z]", "", seq).upper()

        seqres = defaultdict(set)
        entity_of: Dict[str, str] = {}
        for chain, ent, seq_id in zip(as_list(mm.get("_pdbx_poly_seq_scheme.pdb_strand_id")),
                                      as_list(mm.get("_pdbx_poly_seq_scheme.entity_id")),
                                      as_list(mm.get("_pdbx_poly_seq_scheme.seq_id"))):
            seqres[chain].add(seq_id)
            entity_of.setdefault(chain, ent)

        models = as_list(mm.get("_atom_site.pdbx_PDB_model_num"))
        first_model = models[0] if models else None
        resolved = defaultdict(set)
        for i, (chain, seq_id, elem) in enumerate(zip(as_list(mm.get("_atom_site.auth_asym_id")),
                                                      as_list(mm.get("_atom_site.label_seq_id")),
                                                      as_list(mm.get("_atom_site.type_symbol")))):
            if seq_id in (".", "?") or (ignore_h and elem.upper() == "H"):
                continue
            if models and models[i] != first_model:
                continue
            resolved[chain].add(seq_id)

        rows = []
        for chain, ent in entity_of.items():
            seq = seq_by_entity.get(ent)
            if seq is None:
                continue
            n_seqres = len(seqres[chain]) or len(seq)
            n_res = len(resolved[chain]) if chain in resolved else None
            coverage = n_res / n_seqres if (n_res is not None and n_seqres) else None
            rows.append({
                "pdb_id": pdb_id, "assembly_id": assembly_id, "chain_id": chain,
                "entity_id": ent, "seq": seq, "seq_len": len(seq),
                "seqres_count_chain": n_seqres, "resolved_count": n_res,
                "missing_residues": max(0, n_seqres - n_res) if n_res is not None else None,
                "coverage": coverage,
                "coverage_pct": 100.0 * coverage if coverage is not None else None,
            })
        return rows, None
    except Exception as exc:  # malformed files are reported, never fatal
        return [], f"{pdb_id}\tassembly{assembly_id}\t{path}\t{type(exc).__name__}: {exc}"


# -----------------------------------------------------------------------------
# Chunked, resumable execution
# -----------------------------------------------------------------------------
def run_cached(name: str, func, tasks: List[tuple], parts_dir: Path, workers: int,
               flatten: bool) -> Tuple[pd.DataFrame, List[str]]:
    """Run ``func`` over ``tasks`` in chunks; each finished chunk is saved to
    parts_dir and reused if the same chunk is requested again."""
    frames, errors = [], []
    with ProcessPoolExecutor(max_workers=workers) as ex, \
            tqdm(total=len(tasks), desc=name, unit="file") as pbar:
        for start in range(0, len(tasks), CHUNK):
            chunk = tasks[start:start + CHUNK]
            key = hashlib.sha1(repr(chunk).encode()).hexdigest()[:16]
            part = parts_dir / f"{name}_{key}.pkl"
            err_part = parts_dir / f"{name}_{key}.err"
            if part.is_file():
                frames.append(pd.read_pickle(part))
                errors += err_part.read_text(encoding="utf-8").splitlines()
            else:
                records, chunk_errors = [], []
                for result, err in ex.map(func, chunk, chunksize=10):
                    if err:
                        chunk_errors.append(err)
                    elif flatten:
                        records.extend(result)
                    elif result is not None:
                        records.append(result)
                frame = pd.DataFrame.from_records(records)
                tmp = part.with_suffix(".tmp")
                frame.to_pickle(tmp)
                err_part.write_text("".join(e + "\n" for e in chunk_errors), encoding="utf-8")
                tmp.replace(part)
                frames.append(frame)
                errors += chunk_errors
            pbar.update(len(chunk))
    frames = [f for f in frames if not f.empty]
    return (pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()), errors


def positional_identity(a: Any, b: Any) -> Optional[float]:
    """Matches over the shared prefix divided by the longer length (no alignment)."""
    if not isinstance(a, str) or not isinstance(b, str) or not a or not b:
        return None
    matches = sum(x == y for x, y in zip(a, b))
    return matches / max(len(a), len(b))


def main() -> None:
    cfg = parse_args(__doc__)
    paths, enr = cfg["paths"], cfg["enrichment"]
    workers = int(cfg["workers"])
    parts_dir = Path(paths["enriched"]).parent / "06_parts"
    parts_dir.mkdir(parents=True, exist_ok=True)

    master = pd.read_csv(paths["master_selected"], keep_default_na=False, na_values=[""],
                         dtype={"pdb_id": str, "molA_chain_id": str, "molB_chain_id": str})
    master["pdb_id"] = master["pdb_id"].str.lower()
    master["assembly_id"] = master["assembly_id"].astype(int)
    log.info("%d candidate pairs from %d entries", len(master), master["pdb_id"].nunique())

    asu_dir, asm_dir = Path(paths["asu_dir"]), Path(paths["assembly_dir"])
    entries = sorted(master["pdb_id"].unique())
    pairs = sorted(set(zip(master["pdb_id"], master["assembly_id"])))
    entry_tasks = [(p, str(asu_dir / pisa_subdir(p) / f"{p}.cif.gz")) for p in entries]
    asm_tasks = [(p, a, str(asm_dir / pisa_subdir(p) / f"{p}-assembly{a}.cif.gz"), bool(enr["ignore_hydrogens"]))
                 for p, a in pairs]
    entry_tasks = [t for t in entry_tasks if Path(t[1]).is_file()]
    asm_tasks = [t for t in asm_tasks if Path(t[2]).is_file()]
    log.info("Files on disk: %d/%d ASU, %d/%d assemblies",
             len(entry_tasks), len(entries), len(asm_tasks), len(pairs))

    entry_df, err1 = run_cached("asu", parse_entry, entry_tasks, parts_dir, workers, flatten=False)
    chain_df, err2 = run_cached("assembly", parse_assembly, asm_tasks, parts_dir, workers, flatten=True)
    errors_path = Path(paths["reports_dir"]) / "06_mmcif_parse_errors.tsv"
    ensure_parent(errors_path).write_text("".join(e + "\n" for e in err1 + err2), encoding="utf-8")

    entry_df.to_csv(ensure_parent(paths["entry_metadata"]), index=False)
    chain_df.to_csv(ensure_parent(paths["chain_metadata"]), index=False)

    # --- merge: chain data for A and B, then entry metadata ------------------
    keys = ["pdb_id", "assembly_id"]
    out = master
    for side in ("molA", "molB"):
        side_df = chain_df.rename(columns={"chain_id": f"{side}_chain_id",
                                           **{c: f"{side}_{c}" for c in CHAIN_FIELDS}})
        side_df = side_df.rename(columns={f"{side}_entity_id": f"{side}_entity_id_mmcif"})
        out = out.merge(side_df, on=keys + [f"{side}_chain_id"], how="left", validate="many_to_one")
    out = out.merge(entry_df, on="pdb_id", how="left", validate="many_to_one")

    sim = [positional_identity(a, b) for a, b in zip(out["molA_seq"], out["molB_seq"])]
    out["seq_similarity_simple"] = sim
    out["is_homo_entity"] = (out["molA_entity_id_mmcif"].notna()
                             & (out["molA_entity_id_mmcif"] == out["molB_entity_id_mmcif"]))
    out["is_homo_80"] = [s is not None and s >= float(enr["homo_similarity"]) for s in sim]
    out["is_homo"] = out["is_homo_entity"] | out["is_homo_80"]
    out["has_entry_meta"] = out["exptl_method"].notna()
    out["has_chainA_asm"] = out["molA_seq"].notna()
    out["has_chainB_asm"] = out["molB_seq"].notna()

    keep = out["has_entry_meta"] & out["has_chainA_asm"] & out["has_chainB_asm"]
    log.info("Dropped rows: %d without entry metadata, %d without chain A, %d without chain B",
             (~out["has_entry_meta"]).sum(), (~out["has_chainA_asm"]).sum(), (~out["has_chainB_asm"]).sum())
    out = out[keep]

    dest = ensure_parent(paths["enriched"])
    tmp = dest.with_name(dest.name + ".part")
    out.to_csv(tmp, index=False)
    tmp.replace(dest)
    shutil.rmtree(parts_dir)
    log.info("%d enriched pairs from %d entries written to %s", len(out), out["pdb_id"].nunique(), dest)


if __name__ == "__main__":
    main()
