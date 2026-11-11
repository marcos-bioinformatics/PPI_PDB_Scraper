#!/usr/bin/env python3
"""Optional annotation - PRODIGY-cryst call (biological or crystal) for every positive pair.

PRODIGY-cryst needs Python 3.8 and scikit-learn 0.22, so it runs in its own
environment (not part of `all`):

    conda env create -f envs/prodigy_cryst.yml
    conda run -n prodigy_cryst pip install --no-deps prodigy-cryst==1.0.1
    snakemake --cores 8 prodigy_cryst

Input : paths.test, paths.trainval_splits, assembly mmCIF files (step 04)
Output: paths.prodigy_cryst, one row per positive (row_id):
  prodigy_class         BIO or XTAL
  prodigy_p_bio         probability of BIO given by the classifier
  prodigy_ics           intermolecular residue contacts (PRODIGY definition)
  prodigy_link_density  contacts / (interface residues of A x interface residues of B)
  prodigy_error         why a pair could not be scored
It is an annotation: no step filters on it (step 08b uses the same scoring on
small interfaces). Structures are prepared as PRODIGY-cryst does (first model,
one conformer per atom, no HETATM or water), with two differences: only the two
chains of the pair are kept, so other molecules of the assembly (e.g. nucleic
acids) do not stop the classification, and selenomethionine (MSE) and
selenocysteine (SEC), which steps 08 and 13 allow, are scored as their sulphur
parents (MET, CYS) instead of being removed with the HETATM records (a SeMet
chain would otherwise lose its methionines and their contacts). Pairs whose own
chains have a non-standard residue in ATOM records are not scored.
"""

from __future__ import annotations

import csv
import gzip
import logging
import shutil
import tempfile
import warnings
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from prodigy_cryst.interface_classifier import ProdigyCrystal

from common import ensure_parent, get_logger, parse_args, pisa_subdir

log = get_logger("prodigy_cryst")
# selenium residue -> (sulphur parent, selenium atom, sulphur atom)
SELENIUM_PARENT = {"MSE": ("MET", "SE", "SD"), "SEC": ("CYS", "SE", "SG")}
COLUMNS = ["row_id", "pdb_id", "assembly_id", "molA_chain_id", "molB_chain_id", "prodigy_class", "prodigy_p_bio",
           "prodigy_ics", "prodigy_link_density", "prodigy_error"]


def to_sulphur_parent(res) -> None:
    """Turn a HETATM MSE or SEC residue into the ATOM residue of its sulphur parent, in place."""
    parent_name, se_name, s_name = SELENIUM_PARENT[res.resname]
    chain = res.parent
    chain.detach_child(res.id)
    res.id = (" ", res.id[1], res.id[2])
    res.resname = parent_name
    if se_name in res:
        atom = res[se_name]
        res.detach_child(se_name)
        atom.id = atom.name = s_name
        atom.fullname = f" {s_name} "
        atom.element = "S"
        res.add(atom)
    chain.add(res)


def prepare(path: str, chains):
    """Parse an assembly file keeping the requested chains, as PRODIGY-cryst prepares a structure."""
    tmpd = tempfile.mkdtemp()
    try:
        tmp = Path(tmpd) / "s.cif"
        with gzip.open(path, "rb") as fi, tmp.open("wb") as fo:
            shutil.copyfileobj(fi, fo)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            s = MMCIFParser(QUIET=True).get_structure("s", str(tmp))
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)
    for m in list(s)[1:]:
        s.detach_child(m.id)
    model = s[0]
    for c in [c.id for c in model if c.id not in chains]:
        model.detach_child(c)
    for atom in list(s.get_atoms()):  # one conformer: the one Biopython selects (highest occupancy)
        if atom.is_disordered():
            residue, sel = atom.parent, atom.selected_child
            sel.altloc, sel.disordered_flag = " ", 0
            residue.detach_child(atom.id)
            residue.add(sel)
    bad = {}
    for res in list(s.get_residues()):
        if res.id[0].startswith("H_") and res.resname in SELENIUM_PARENT:
            to_sulphur_parent(res)
        elif res.id[0][0] in ("W", "H"):
            res.parent.detach_child(res.id)
        elif not is_aa(res, standard=True):
            bad[res.parent.id] = res.resname
    return s, bad


def score(task):
    path, recs = task
    out = []
    chains = {c for r in recs for c in (r["molA_chain_id"], r["molB_chain_id"])}
    try:
        s, bad = prepare(path, chains)
    except Exception as exc:  # unreadable file: every pair of this assembly is reported
        return [{**r, "prodigy_error": f"parse {type(exc).__name__}: {exc}"[:120]} for r in recs]
    for r in recs:
        a, b = r["molA_chain_id"], r["molB_chain_id"]
        if a in bad or b in bad:
            out.append({**r, "prodigy_error": f"non-standard residue {bad.get(a) or bad.get(b)}"})
            continue
        try:
            pc = ProdigyCrystal(s, [a, b])
            pc.predict()
            cls, p_bio, _ = pc.predicted_class
            out.append({**r, "prodigy_class": cls, "prodigy_p_bio": round(float(p_bio), 4),
                        "prodigy_ics": len(pc.ic_network), "prodigy_link_density": round(pc.link_density, 4)})
        except Exception as exc:
            out.append({**r, "prodigy_error": f"{type(exc).__name__}: {exc}"[:120]})
    return out


def main() -> None:
    cfg = parse_args(__doc__)
    paths = cfg["paths"]
    logging.getLogger("Prodigy").setLevel(logging.ERROR)
    csv.field_size_limit(1 << 30)
    by_asm = defaultdict(list)
    for name in ("test", "trainval_splits"):
        with open(paths[name], newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                rec = {k: r[k] for k in COLUMNS[:5]}
                pid, asm = r["pdb_id"].lower(), int(float(r["assembly_id"]))
                path = Path(paths["assembly_dir"]) / pisa_subdir(pid) / f"{pid}-assembly{asm}.cif.gz"
                by_asm[str(path)].append(rec)
    tasks = sorted(by_asm.items())
    log.info("Scoring %d pairs in %d assemblies", sum(len(v) for v in by_asm.values()), len(tasks))
    rows = []
    with Pool(int(cfg["workers"])) as pool:
        for i, res in enumerate(pool.imap_unordered(score, tasks, chunksize=8), 1):
            rows.extend(res)
            if i % 2000 == 0:
                log.info("%d of %d assemblies", i, len(tasks))
    rows.sort(key=lambda r: r["row_id"])
    dest = ensure_parent(paths["prodigy_cryst"])
    tmp = dest.with_name(dest.name + ".part")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(dest)
    scored = [r for r in rows if r.get("prodigy_class")]
    log.info("%d of %d pairs scored; %d BIO (%.1f %%)", len(scored), len(rows),
             sum(r["prodigy_class"] == "BIO" for r in scored),
             100 * sum(r["prodigy_class"] == "BIO" for r in scored) / max(1, len(scored)))


if __name__ == "__main__":
    main()
