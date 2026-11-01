"""Structure-level checks shared by steps 08 and 13.

Chains are read from the mmCIF _atom_site table (first model, first alternate
location). Only polymer residues (those with a label_seq_id) are used, so
ligands, waters and free amino acids that share the author chain ID never
count as special residues nor create contacts.

- Special residues: every polymer residue that is not one of the 20 standard
  amino acids, is not in `allowed` and is not in `ignore` (terminal caps such as
  ACE / NH2) makes the chain forbidden. This includes UNK, D-amino acids,
  chromophores (CRO, ...) and post-translational modifications.
- Contacts: heavy-atom contacts between two chains with a KD-tree. The
  neighbour search is adapted from NeighborSearch.py (pairs_within_cutoff) and
  the heavy-atom criterion from Test_distance.py of the BeastAffinity ML code.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from scipy.spatial import cKDTree

from common import open_mmcif

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


@dataclass
class Chain:
    residues: List[str]       # residue name of every observed polymer residue
    coords: np.ndarray        # (n, 3) heavy-atom coordinates of those residues
    residue_idx: np.ndarray   # (n,) index into `residues` of each atom


def _column(mm: Dict, key: str) -> np.ndarray:
    value = mm.get(key, [])
    return np.asarray(value if isinstance(value, list) else [value], dtype=object)


def load_chains(path: str, chain_ids: Optional[Iterable[str]] = None) -> Dict[str, Chain]:
    """Polymer residues and heavy atoms of the requested author chains."""
    with open_mmcif(path) as fh:
        mm = MMCIF2Dict(fh)
    chain = _column(mm, "_atom_site.auth_asym_id")
    seq_id = _column(mm, "_atom_site.label_seq_id")
    comp = _column(mm, "_atom_site.label_comp_id")
    alt = _column(mm, "_atom_site.label_alt_id")
    elem = _column(mm, "_atom_site.type_symbol")
    model = _column(mm, "_atom_site.pdbx_PDB_model_num")
    xyz = np.column_stack([_column(mm, f"_atom_site.Cartn_{a}").astype(float) for a in "xyz"])

    keep = (seq_id != ".") & (seq_id != "?") & ~np.isin(elem, ["H", "D"])
    keep &= np.isin(alt, [".", "?", "A", "1"])
    if len(model):
        keep &= model == model[0]
    if chain_ids is not None:
        keep &= np.isin(chain, list(chain_ids))

    out: Dict[str, Chain] = {}
    for c in dict.fromkeys(chain[keep]):  # first-seen order
        sel = keep & (chain == c)
        sids = seq_id[sel]
        order = list(dict.fromkeys(sids))  # residues in file order
        index = {s: i for i, s in enumerate(order)}
        first_comp: Dict[str, str] = {}
        for s, name in zip(sids, comp[sel]):
            first_comp.setdefault(s, str(name).upper())  # microheterogeneity: first residue name
        out[str(c)] = Chain(residues=[first_comp[s] for s in order],
                            coords=xyz[sel].astype(np.float32),
                            residue_idx=np.fromiter((index[s] for s in sids), dtype=np.int32, count=len(sids)))
    return out


def fmt_counter(counter: Counter) -> str:
    return ";".join(f"{k}:{counter[k]}" for k in sorted(counter))


def special_residues(chain: Chain, allowed: Iterable[str], ignore: Iterable[str] = ()) -> Tuple[Counter, Counter]:
    """(allowed, forbidden) counts of non-standard polymer residues in a chain."""
    allowed, ignore = set(allowed), set(ignore)
    ok, bad = Counter(), Counter()
    for name in chain.residues:
        if name in STANDARD_AA or name in ignore:
            continue
        (ok if name in allowed else bad)[name] += 1
    return ok, bad


def pairs_within_cutoff(coords_a: np.ndarray, coords_b: np.ndarray, cutoff: float) -> Tuple[np.ndarray, np.ndarray]:
    """Index pairs (i, j) with ||a_i - b_j|| <= cutoff (KD-tree on b)."""
    if len(coords_a) == 0 or len(coords_b) == 0:
        return np.zeros(0, dtype=np.int32), np.zeros(0, dtype=np.int32)
    neighbours = cKDTree(coords_b).query_ball_point(coords_a, r=float(cutoff))
    i = np.repeat(np.arange(len(coords_a), dtype=np.int32), [len(n) for n in neighbours])
    j = np.fromiter((x for n in neighbours for x in n), dtype=np.int32, count=len(i))
    return i, j


def contact_summary(a: Chain, b: Chain, cutoff: float) -> Dict[str, float]:
    """Heavy-atom contacts between two chains."""
    if len(a.coords) == 0 or len(b.coords) == 0:
        return {"n_atom_contacts": 0, "n_residue_contacts": 0, "min_interchain_distance": float("nan")}
    i, j = pairs_within_cutoff(a.coords, b.coords, cutoff)
    residue_pairs = set(zip(a.residue_idx[i].tolist(), b.residue_idx[j].tolist()))
    min_dist = float(cKDTree(b.coords).query(a.coords, k=1)[0].min())
    return {"n_atom_contacts": int(len(i)), "n_residue_contacts": len(residue_pairs),
            "min_interchain_distance": round(min_dist, 3)}
