#!/usr/bin/env python3
"""Step 06 - Enrich every candidate pair with mmCIF data (entry + chain level).

Input : paths.master_selected (step 03), assembly mmCIF files,
        paths.entry_metadata and paths.entity_metadata (step 04)
Output: paths.enriched         master rows + sequence, coverage, homo/hetero and entry metadata
        paths.chain_metadata   one row per protein chain of every needed assembly

Chains are matched by auth_asym_id: PISA's chain_id is the author chain ID of
the assembly file, including symmetry copies such as "A-2".

Per chain (assembly file, first model):
  seq                canonical sequence of its entity (_entity_poly)
  seqres_count       residues in _pdbx_poly_seq_scheme
  resolved_count     polymer residues with at least one non-hydrogen atom of
                     occupancy > 0 (ATOM or HETATM, so MSE counts as resolved);
                     atoms with zero occupancy are not part of the
                     experimental model (PDBx groups them with unobserved residues)
  coverage           resolved_count / seqres_count
  internal_coverage  resolved_count / residues between the first and last
                     resolved ones (ignores disordered tags and tails)
  missing_nterm/_cterm  unresolved residues before the first / after the last
                     resolved residue
  uniprot_ids, membrane_sources, is_antibody   RCSB annotations of the chain's entity (step 04)

seq_identity: identical residues of a global alignment with free end gaps
(BLOSUM62) divided by the length of the shorter chain, so a truncated construct
of the same protein scores 1. is_homo = same entity OR seq_identity >=
enrichment.homo_similarity. Rows without entry metadata or without both chain
sequences are dropped.
Parsed files are cached per chunk in a parts/ folder so an interrupted run resumes.
"""

from __future__ import annotations

import hashlib
import re
import shutil
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from tqdm import tqdm

from common import ensure_parent, get_logger, open_mmcif, parse_args, pisa_subdir

log = get_logger("06_enrich_mmcif")

CHUNK = 500
CACHE_VERSION = 3  # bump when parse_assembly() output changes
CHAIN_FIELDS = ["entity_id", "seq", "seq_len", "seqres_count_chain", "resolved_count",
                "missing_residues", "coverage", "coverage_pct",
                "internal_coverage", "missing_nterm", "missing_cterm",
                "uniprot_ids", "membrane_sources", "is_antibody"]


# -----------------------------------------------------------------------------
# mmCIF parsing
# -----------------------------------------------------------------------------
def as_list(value: Any) -> List[str]:
    if value is None:
        return []
    return [str(v) for v in value] if isinstance(value, list) else [str(value)]


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
            # terminal X are caps such as ACE / NH2: dropping them gives capped and
            # uncapped constructs of a protein the same sequence
            clean = re.sub(r"[^A-Za-z]", "", seq).upper().strip("X")
            if clean:
                seq_by_entity[ent] = clean

        seqres = defaultdict(set)
        entity_of: Dict[str, str] = {}
        for chain, ent, seq_id in zip(as_list(mm.get("_pdbx_poly_seq_scheme.pdb_strand_id")),
                                      as_list(mm.get("_pdbx_poly_seq_scheme.entity_id")),
                                      as_list(mm.get("_pdbx_poly_seq_scheme.seq_id"))):
            seqres[chain].add(int(seq_id))
            entity_of.setdefault(chain, ent)

        models = as_list(mm.get("_atom_site.pdbx_PDB_model_num"))
        first_model = models[0] if models else None
        resolved = defaultdict(set)
        occupancy = as_list(mm.get("_atom_site.occupancy"))
        for i, (chain, seq_id, elem) in enumerate(zip(as_list(mm.get("_atom_site.auth_asym_id")),
                                                      as_list(mm.get("_atom_site.label_seq_id")),
                                                      as_list(mm.get("_atom_site.type_symbol")))):
            if seq_id in (".", "?") or (ignore_h and elem.upper() in ("H", "D")):
                continue
            if occupancy and occupancy[i] not in (".", "?") and float(occupancy[i]) <= 0:
                continue
            if models and models[i] != first_model:
                continue
            resolved[chain].add(int(seq_id))

        rows = []
        for chain, ent in entity_of.items():
            seq = seq_by_entity.get(ent)
            if seq is None:
                continue
            n_seqres = len(seqres[chain]) or len(seq)
            row = {"pdb_id": pdb_id, "assembly_id": assembly_id, "chain_id": chain,
                   "entity_id": ent, "seq": seq, "seq_len": len(seq), "seqres_count_chain": n_seqres}
            if resolved.get(chain):
                obs = resolved[chain]
                lo, hi = min(obs), max(obs)
                coverage = len(obs) / n_seqres
                row.update({
                    "resolved_count": len(obs),
                    "missing_residues": max(0, n_seqres - len(obs)),
                    "coverage": coverage,
                    "coverage_pct": 100.0 * coverage,
                    "internal_coverage": len(obs) / sum(1 for s in seqres[chain] if lo <= s <= hi),
                    "missing_nterm": sum(1 for s in seqres[chain] if s < lo),
                    "missing_cterm": sum(1 for s in seqres[chain] if s > hi),
                })
            rows.append(row)
        return rows, None
    except Exception as exc:  # malformed files are reported, never fatal
        return [], f"{pdb_id}\tassembly{assembly_id}\t{path}\t{type(exc).__name__}: {exc}"


# -----------------------------------------------------------------------------
# Chunked, resumable execution
# -----------------------------------------------------------------------------
def run_cached(name: str, func, tasks: List[tuple], parts_dir: Path, workers: int) -> Tuple[pd.DataFrame, List[str]]:
    """Run ``func`` over ``tasks`` in chunks; each finished chunk is saved to
    parts_dir and reused if the same chunk is requested again."""
    frames, errors = [], []
    with ProcessPoolExecutor(max_workers=workers) as ex, \
            tqdm(total=len(tasks), desc=name, unit="file") as pbar:
        for start in range(0, len(tasks), CHUNK):
            chunk = tasks[start:start + CHUNK]
            key = hashlib.sha1(repr((CACHE_VERSION, chunk)).encode()).hexdigest()[:16]
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
                    else:
                        records.extend(result)
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


_ALIGNER = None


def _to_matrix_alphabet(seq: str, alphabet: str) -> str:
    """Selenocysteine (U) and pyrrolysine (O) are not in BLOSUM62: score them as
    C and K; any other unknown letter as X."""
    seq = seq.replace("U", "C").replace("O", "K")
    return "".join(c if c in alphabet else "X" for c in seq)


def aligned_identity(pair: Tuple[str, str]) -> Optional[float]:
    """Identical residues of a global alignment with free end gaps (BLOSUM62,
    gap open -10 / extend -0.5) divided by the length of the shorter sequence."""
    global _ALIGNER
    if _ALIGNER is None:
        _ALIGNER = Align.PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5,
                                         end_gap_score=0, substitution_matrix=substitution_matrices.load("BLOSUM62"))
    alphabet = _ALIGNER.substitution_matrix.alphabet
    a, b = (_to_matrix_alphabet(x, alphabet) for x in pair)
    if not a or not b:
        return None
    try:
        aln = _ALIGNER.align(a, b)[0]
    except (ValueError, OverflowError):
        return None
    matches = sum(1 for x, y in zip(aln[0], aln[1]) if x == y and x != "-")
    return matches / min(len(a), len(b))


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

    asm_dir = Path(paths["assembly_dir"])
    pairs = sorted(set(zip(master["pdb_id"], master["assembly_id"])))
    asm_tasks = [(p, a, str(asm_dir / pisa_subdir(p) / f"{p}-assembly{a}.cif.gz"), bool(enr["ignore_hydrogens"]))
                 for p, a in pairs]
    asm_tasks = [t for t in asm_tasks if Path(t[2]).is_file()]
    log.info("Assembly files on disk: %d/%d", len(asm_tasks), len(pairs))

    entry_df = pd.read_csv(paths["entry_metadata"], keep_default_na=False, na_values=[""], dtype={"pdb_id": str})
    entity_df = pd.read_csv(paths["entity_metadata"], keep_default_na=False, na_values=[""],
                            dtype={"pdb_id": str, "entity_id": str, "uniprot_ids": str, "membrane_sources": str})
    chain_df, errors = run_cached("assembly", parse_assembly, asm_tasks, parts_dir, workers)
    errors_path = Path(paths["reports_dir"]) / "06_mmcif_parse_errors.tsv"
    ensure_parent(errors_path).write_text("".join(e + "\n" for e in errors), encoding="utf-8")
    chain_df = chain_df.merge(entity_df, on=["pdb_id", "entity_id"], how="left", validate="many_to_one")
    chain_df[["uniprot_ids", "membrane_sources"]] = chain_df[["uniprot_ids", "membrane_sources"]].fillna("")
    chain_df["is_antibody"] = chain_df["is_antibody"].fillna(False).astype(bool)
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

    out["is_homo_entity"] = (out["molA_entity_id_mmcif"].notna()
                             & (out["molA_entity_id_mmcif"] == out["molB_entity_id_mmcif"]))
    # same entity -> identical sequences; align only the pairs of different entities
    to_align = sorted({(a, b) for a, b, same in zip(out["molA_seq"], out["molB_seq"], out["is_homo_entity"])
                       if not same and isinstance(a, str) and isinstance(b, str)})
    with ProcessPoolExecutor(max_workers=workers) as ex:
        identity = dict(zip(to_align, tqdm(ex.map(aligned_identity, to_align, chunksize=50),
                                           total=len(to_align), desc="Align", unit="pair")))
    out["seq_identity"] = [1.0 if same else identity.get((a, b))
                           for a, b, same in zip(out["molA_seq"], out["molB_seq"], out["is_homo_entity"])]
    out["is_homo_80"] = out["seq_identity"].fillna(0) >= float(enr["homo_similarity"])
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
