"""Structure-level checks shared by steps 08 and 13.

Chains are read from the mmCIF _atom_site table: first model, heavy atoms with
occupancy > 0 and, where a residue has alternate locations, the conformer with
the highest mean occupancy (ties: the first in the file), so steps 06, 08 and 13
see the same observed residues. Only polymer residues (those with a
label_seq_id) are used, so ligands, waters and free amino acids that share the
author chain ID never count as special residues nor create contacts.

- Special residues: every polymer residue that is not one of the 20 standard
  amino acids, is not in `allowed` and is not in `ignore` (terminal caps such as
  ACE / NH2) makes the chain forbidden. This includes UNK, D-amino acids,
  chromophores (CRO, ...) and post-translational modifications.
- Contacts: heavy-atom contacts between two chains with a KD-tree. The
  neighbour search is adapted from NeighborSearch.py (pairs_within_cutoff) and
  the heavy-atom criterion from Test_distance.py of the PPI_PDB_Scraper ML code.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple

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
    seq_ids: List[int]        # label_seq_id of every observed residue (same order as `residues`)
    seqres: Set[int]          # label_seq_id of every SEQRES residue of the chain


def _column(mm: Dict, key: str) -> np.ndarray:
    value = mm.get(key, [])
    return np.asarray(value if isinstance(value, list) else [value], dtype=object)


def _other_conformers(keep: np.ndarray, chain: np.ndarray, seq_id: np.ndarray, alt: np.ndarray,
                      occ: np.ndarray) -> np.ndarray:
    """Atoms of alternate locations other than the one kept for their residue:
    the conformer with the highest mean occupancy (ties: first in the file)."""
    drop = np.zeros(len(keep), dtype=bool)
    idx = np.flatnonzero(keep & ~np.isin(alt, [".", "?"]))
    if not len(idx):
        return drop
    totals: Dict[Tuple[str, str, str], List[float]] = {}
    for i in idx:
        t = totals.setdefault((chain[i], seq_id[i], alt[i]), [0.0, 0])
        t[0] += occ[i] if len(occ) else 1.0
        t[1] += 1
    best: Dict[Tuple[str, str], Tuple[str, float]] = {}
    for (c, s, a), (total, n) in totals.items():  # first-seen order, so ties keep the first
        if (c, s) not in best or total / n > best[(c, s)][1]:
            best[(c, s)] = (a, total / n)
    for i in idx:
        drop[i] = best[(chain[i], seq_id[i])][0] != alt[i]
    return drop


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
    occ = np.array([1.0 if o in (".", "?") else float(o) for o in _column(mm, "_atom_site.occupancy")])
    xyz = np.column_stack([_column(mm, f"_atom_site.Cartn_{a}").astype(float) for a in "xyz"])

    keep = (seq_id != ".") & (seq_id != "?") & ~np.isin(elem, ["H", "D"])
    if len(occ):
        keep &= occ > 0
    if len(model):
        keep &= model == model[0]
    if chain_ids is not None:
        keep &= np.isin(chain, list(chain_ids))
    keep &= ~_other_conformers(keep, chain, seq_id, alt, occ)
    seqres: Dict[str, Set[int]] = {}
    for c, s in zip(_column(mm, "_pdbx_poly_seq_scheme.pdb_strand_id"), _column(mm, "_pdbx_poly_seq_scheme.seq_id")):
        seqres.setdefault(str(c), set()).add(int(s))

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
                            residue_idx=np.fromiter((index[s] for s in sids), dtype=np.int32, count=len(sids)),
                            seq_ids=[int(s) for s in order],
                            seqres=seqres.get(str(c), set()))
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


def interface_gaps(a: Chain, b: Chain, cutoff: float, window: int) -> Dict[str, int]:
    """Unresolved SEQRES residues of each chain within `window` positions of one
    of its interface residues (residues with a heavy atom within `cutoff` of the
    other chain): missing density at or next to the interface."""
    if len(a.coords) == 0 or len(b.coords) == 0:
        return {"molA_interface_gaps": 0, "molB_interface_gaps": 0}
    i, j = pairs_within_cutoff(a.coords, b.coords, cutoff)
    out = {}
    for side, chain, idx in (("molA", a, i), ("molB", b, j)):
        iface = {chain.seq_ids[k] for k in set(chain.residue_idx[idx].tolist())}
        near = {p + d for p in iface for d in range(-window, window + 1)}
        missing = chain.seqres - set(chain.seq_ids)
        out[f"{side}_interface_gaps"] = len(missing & near)
    return out


def clash_summary(a: Chain, b: Chain, distance: float) -> Dict[str, float]:
    """Atoms of each chain closer than `distance` to the other chain: their number
    (both chains) and the larger fraction of the two chains (AlphaFold 3 removes
    chains with more than 30 % of their atoms within 1.7 A of another chain)."""
    if len(a.coords) == 0 or len(b.coords) == 0:
        return {"n_clash_atoms": 0, "max_clash_fraction": 0.0}
    da = cKDTree(b.coords).query(a.coords, k=1)[0]
    db = cKDTree(a.coords).query(b.coords, k=1)[0]
    return {"n_clash_atoms": int((da < distance).sum() + (db < distance).sum()),
            "max_clash_fraction": round(float(max((da < distance).mean(), (db < distance).mean())), 4)}
