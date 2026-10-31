#!/usr/bin/env python3
"""Step 04 - Download the assembly mmCIF files and entry metadata of the selected pairs.

Input : paths.needed_pairs   (pdb_id, assembly_id) pairs from step 03
Output: paths.assembly_dir/{subdir}/{id}-assembly{n}.cif.gz   every needed assembly
        paths.structures_manifest                             status of every file
        paths.entry_metadata                                  method, resolution, R-free, title, keywords

Only the biological assemblies that PISA analysed are downloaded. Files are
stored gzipped, fully checked for gzip integrity when downloaded and written
atomically; files already on disk are kept, so the step can be rerun to resume
or to retry failures.

Entry-level metadata comes from the RCSB Data API (GraphQL, batches of 1000
IDs). The values are those of the first _exptl and _refine records of the
entry's mmCIF, so the asymmetric-unit files do not need to be downloaded.
"""

from __future__ import annotations

import csv
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from common import (download_file, ensure_parent, get_logger, gzip_is_valid, gzip_looks_valid, make_session,
                    parse_args, pisa_subdir, rcsb_entries)

log = get_logger("04_download_structures")

ENTRY_FIELDS = ("exptl { method } refine { ls_d_res_high ls_R_factor_R_free } "
                "struct { title } struct_keywords { pdbx_keywords }")


def first(records, key):
    return (records or [{}])[0].get(key)


def main() -> None:
    cfg = parse_args(__doc__)
    paths, rcsb = cfg["paths"], cfg["rcsb"]
    base_url = rcsb["download_url"].rstrip("/")
    asm_dir = Path(paths["assembly_dir"])
    workers = int(cfg["download_workers"])

    pairs = pd.read_csv(paths["needed_pairs"], dtype={"pdb_id": str})
    pairs["pdb_id"] = pairs["pdb_id"].str.lower()
    pairs["assembly_id"] = pairs["assembly_id"].astype(int)
    pairs = pairs.drop_duplicates().sort_values(["pdb_id", "assembly_id"])
    entries = sorted(pairs["pdb_id"].unique())
    session = make_session(pool_size=workers)

    # 1) Entry metadata from the RCSB Data API
    rows = []
    for e in rcsb_entries(session, rcsb["graphql_url"], entries, ENTRY_FIELDS):
        rows.append({
            "pdb_id": e["rcsb_id"].lower(),
            "exptl_method": first(e.get("exptl"), "method"),
            "resolution": first(e.get("refine"), "ls_d_res_high"),
            "r_free": first(e.get("refine"), "ls_R_factor_R_free"),
            "title": (e.get("struct") or {}).get("title"),
            "keywords": (e.get("struct_keywords") or {}).get("pdbx_keywords"),
        })
    meta = pd.DataFrame(rows, columns=["pdb_id", "exptl_method", "resolution", "r_free", "title", "keywords"])
    meta.sort_values("pdb_id").to_csv(ensure_parent(paths["entry_metadata"]), index=False)
    log.info("Entry metadata for %d of %d entries", len(meta), len(entries))

    # 2) Assembly files
    tasks = []
    for pid, asm in zip(pairs["pdb_id"], pairs["assembly_id"]):
        name = f"{pid}-assembly{asm}.cif.gz"
        tasks.append((pid, asm, f"{base_url}/{name}", asm_dir / pisa_subdir(pid) / name))
    log.info("%d assembly files needed", len(tasks))

    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = [(t, ex.submit(download_file, session, t[2], t[3], gzip_is_valid, exists_check=gzip_looks_valid))
                   for t in tasks]
        results = [(t, fut.result()) for t, fut in tqdm(futures, desc="mmCIF", unit="file")]

    manifest = ensure_parent(paths["structures_manifest"])
    tmp = manifest.with_name(manifest.name + ".part")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["pdb_id", "assembly_id", "path", "status"])
        for (pid, asm, _url, dest), status in results:
            w.writerow([pid, asm, dest, status])
    tmp.replace(manifest)

    counts = pd.Series([s.split(":")[0] for _, s in results]).value_counts().to_dict()
    log.info("Download summary: %s", counts)
    if counts.get("error") or counts.get("missing"):
        log.warning("Some files are missing or failed (see %s); rerun this step to retry errors.", manifest)


if __name__ == "__main__":
    main()
