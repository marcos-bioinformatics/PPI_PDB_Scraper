"""Structure-level checks shared by steps 08 and 13.

- Non-standard amino acids per chain (Biopython residues of the first model).
- Inter-chain contacts between heavy atoms with a KD-tree. The neighbour search
  is adapted from NeighborSearch.py (pairs_within_cutoff) and the heavy-atom
  criterion from Test_distance.py of the BeastAffinity ML preprocessing code.
"""

from __future__ import annotations

import warnings
from collections import Counter
from typing import Dict, Iterable, Tuple

import numpy as np
from Bio import BiopythonWarning
from Bio.PDB import MMCIFParser, is_aa
from scipy.spatial import cKDTree

from common import open_mmcif

warnings.simplefilter("ignore", BiopythonWarning)

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


def load_model(path: str):
    """First model of an mmCIF file (author chain IDs, as used by PISA)."""
    with open_mmcif(path) as fh:
        return MMCIFParser(QUIET=True).get_structure("X", fh)[0]


def fmt_counter(counter: Counter) -> str:
    return ";".join(f"{k}:{counter[k]}" for k in sorted(counter))


def special_residues(chain, allowed: Iterable[str]) -> Tuple[Counter, Counter]:
    """(allowed, forbidden) counts of non-standard amino-acid residues in a chain."""
    allowed = set(allowed)
    ok, bad = Counter(), Counter()
    for res in chain.get_residues():
        if not is_aa(res, standard=False):
            continue
        name = res.get_resname().upper().strip()
        if name not in STANDARD_AA:
            (ok if name in allowed else bad)[name] += 1
    return ok, bad


def heavy_atoms(chain) -> Tuple[np.ndarray, np.ndarray]:
    """Coordinates (n, 3) of the non-hydrogen atoms of a chain and the index of
    the residue each atom belongs to."""
    coords, residue_idx = [], []
    for i, res in enumerate(chain.get_residues()):
        for atom in res.get_atoms():
            if atom.element not in ("H", "D"):
                coords.append(atom.coord)
                residue_idx.append(i)
    return np.asarray(coords, dtype=np.float32).reshape(-1, 3), np.asarray(residue_idx, dtype=np.int32)


def pairs_within_cutoff(coords_a: np.ndarray, coords_b: np.ndarray, cutoff: float) -> Tuple[np.ndarray, np.ndarray]:
    """Index pairs (i, j) with ||a_i - b_j|| <= cutoff (KD-tree on b)."""
    if len(coords_a) == 0 or len(coords_b) == 0:
        return np.zeros(0, dtype=np.int32), np.zeros(0, dtype=np.int32)
    neighbours = cKDTree(coords_b).query_ball_point(coords_a, r=float(cutoff))
    i = np.repeat(np.arange(len(coords_a), dtype=np.int32), [len(n) for n in neighbours])
    j = np.fromiter((x for n in neighbours for x in n), dtype=np.int32, count=len(i))
    return i, j


def contact_summary(a: Tuple[np.ndarray, np.ndarray], b: Tuple[np.ndarray, np.ndarray],
                    cutoff: float) -> Dict[str, float]:
    """Heavy-atom contacts between two chains given by heavy_atoms()."""
    (ca, ra), (cb, rb) = a, b
    if len(ca) == 0 or len(cb) == 0:
        return {"n_atom_contacts": 0, "n_residue_contacts": 0, "min_interchain_distance": float("nan")}
    i, j = pairs_within_cutoff(ca, cb, cutoff)
    residue_pairs = set(zip(ra[i].tolist(), rb[j].tolist()))
    min_dist = float(cKDTree(cb).query(ca, k=1)[0].min())
    return {"n_atom_contacts": int(len(i)), "n_residue_contacts": len(residue_pairs),
            "min_interchain_distance": round(min_dist, 3)}
