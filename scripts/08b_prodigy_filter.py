#!/usr/bin/env python3
"""Step 08b - Drop small interfaces that PRODIGY-cryst calls crystal contacts.

Biological interfaces and crystal contacts overlap at small interface areas
(I-03). Among the pairs with interface_area < prodigy.max_area, those that
PRODIGY-cryst, a contact-based classifier independent of PISA, calls crystal
(XTAL) are removed. Pairs it cannot score (a non-standard residue in their own
chains, an unreadable file) are kept and counted.

Runs only when prodigy.filter is true, in the prodigy_cryst environment
(PRODIGY-cryst needs Python 3.8 and scikit-learn 0.22):

    conda env create -f envs/prodigy_cryst.yml
    conda run -n prodigy_cryst pip install --no-deps prodigy-cryst==1.0.1

Input : paths.checked (step 08), assembly mmCIF files (step 04)
Output: paths.prodigy_filtered            rows of step 08 without the XTAL pairs (read by step 09)
        paths.prodigy_band                PRODIGY-cryst call of every pair below max_area
        reports/08b_prodigy_filter.json   counts
"""

from __future__ import annotations

import csv
import logging
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

from annotate_prodigy_cryst import COLUMNS, score
from common import ensure_parent, get_logger, parse_args, pisa_subdir, write_json

log = get_logger("08b_prodigy_filter")


def row_key(r) -> str:
    """Same fields as the row_id of step 10."""
    return "|".join(r[k] for k in ("pdb_id", "assembly_id", "interface_id", "pair_index_in_interface",
                                   "molA_chain_id", "molB_chain_id"))


def main() -> None:
    cfg = parse_args(__doc__)
    paths, pc = cfg["paths"], cfg["prodigy"]
    if not pc.get("filter"):
        raise SystemExit("prodigy.filter is false: step 09 reads paths.checked and this step is not needed")
    max_area = float(pc["max_area"])
    logging.getLogger("Prodigy").setLevel(logging.ERROR)
    csv.field_size_limit(1 << 30)

    by_asm, info, n_in = defaultdict(list), {}, 0
    with open(paths["checked"], newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            n_in += 1
            if float(r["interface_area"]) >= max_area:
                continue
            key = row_key(r)
            info[key] = {"interface_area": r["interface_area"], "is_homo": r["is_homo"], "tier": r["tier"]}
            pid, asm = r["pdb_id"].lower(), int(float(r["assembly_id"]))
            path = Path(paths["assembly_dir"]) / pisa_subdir(pid) / f"{pid}-assembly{asm}.cif.gz"
            by_asm[str(path)].append({"row_id": key, "pdb_id": r["pdb_id"], "assembly_id": r["assembly_id"],
                                      "molA_chain_id": r["molA_chain_id"], "molB_chain_id": r["molB_chain_id"]})
    tasks = sorted(by_asm.items())
    log.info("%d rows; scoring the %d with interface_area < %g in %d assemblies",
             n_in, len(info), max_area, len(tasks))

    rows = []
    with Pool(int(cfg["workers"])) as pool:
        for i, res in enumerate(pool.imap_unordered(score, tasks, chunksize=8), 1):
            rows.extend(res)
            if i % 2000 == 0:
                log.info("%d of %d assemblies", i, len(tasks))
    rows.sort(key=lambda r: r["row_id"])
    for r in rows:
        r.update(info[r["row_id"]])
    xtal = {r["row_id"] for r in rows if r.get("prodigy_class") == "XTAL"}

    band = ensure_parent(paths["prodigy_band"])
    tmp = band.with_name(band.name + ".part")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS + ["interface_area", "is_homo", "tier"], extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(band)

    dest = ensure_parent(paths["prodigy_filtered"])
    tmp = dest.with_name(dest.name + ".part")
    n_out = 0
    with open(paths["checked"], newline="", encoding="utf-8") as fi, tmp.open("w", newline="", encoding="utf-8") as fo:
        reader = csv.DictReader(fi)
        w = csv.DictWriter(fo, fieldnames=reader.fieldnames)
        w.writeheader()
        for r in reader:
            if row_key(r) not in xtal:
                w.writerow(r)
                n_out += 1
    tmp.replace(dest)

    scored = [r for r in rows if r.get("prodigy_class")]
    removed = [r for r in rows if r["row_id"] in xtal]
    kind = lambda r: "homo" if r["is_homo"].strip().lower() in ("true", "1") else "hetero"  # noqa: E731
    report = {
        "max_area": max_area,
        "rows_in": n_in,
        "rows_below_max_area": len(rows),
        "scored": len(scored),
        "not_scored_kept": len(rows) - len(scored),
        "not_scored_reasons": dict(Counter(r["prodigy_error"].split(":")[0] for r in rows if not r.get("prodigy_class"))),
        "removed_xtal": len(removed),
        "removed_by_kind": dict(Counter(kind(r) for r in removed)),
        "removed_by_tier": dict(Counter(r["tier"] for r in removed)),
        "rows_out": n_out,
    }
    write_json(Path(paths["reports_dir"]) / "08b_prodigy_filter.json", report)
    log.info("Removed %d XTAL pairs of %d scored below %g A2 (%d not scored, kept); %d rows left",
             len(removed), len(scored), max_area, len(rows) - len(scored), n_out)


if __name__ == "__main__":
    main()
