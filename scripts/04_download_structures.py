#!/usr/bin/env python3
"""Step 04 - Download the mmCIF files needed by the selected PISA pairs.

Input : paths.needed_pairs   (pdb_id, assembly_id) pairs from step 03
Output: paths.asu_dir/{subdir}/{id}.cif.gz                    one ASU per entry
        paths.assembly_dir/{subdir}/{id}-assembly{n}.cif.gz   every needed assembly
        paths.structures_manifest                             status of every file

Only what the pipeline will use is downloaded: the ASU file (entry metadata:
method, resolution, R-free) of every entry with a candidate pair, and the
biological assembly files that PISA analysed. Files are stored gzipped, checked
for gzip integrity and written atomically; valid files already on disk are
skipped, so the step can be rerun to resume or to retry failures.
"""

from __future__ import annotations

import csv
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from common import download_file, ensure_parent, get_logger, gzip_is_valid, make_session, parse_args, pisa_subdir

log = get_logger("04_download_structures")


def main() -> None:
    cfg = parse_args(__doc__)
    paths = cfg["paths"]
    base_url = cfg["structures"]["base_url"].rstrip("/")
    asu_dir, asm_dir = Path(paths["asu_dir"]), Path(paths["assembly_dir"])
    workers = int(cfg["download_workers"])

    pairs = pd.read_csv(paths["needed_pairs"], dtype={"pdb_id": str})
    pairs["pdb_id"] = pairs["pdb_id"].str.lower()
    pairs["assembly_id"] = pairs["assembly_id"].astype(int)
    pairs = pairs.drop_duplicates().sort_values(["pdb_id", "assembly_id"])

    tasks = []
    for pid in pairs["pdb_id"].unique():
        name = f"{pid}.cif.gz"
        tasks.append(("asu", pid, "", f"{base_url}/{name}", asu_dir / pisa_subdir(pid) / name))
    for pid, asm in zip(pairs["pdb_id"], pairs["assembly_id"]):
        name = f"{pid}-assembly{asm}.cif.gz"
        tasks.append(("assembly", pid, asm, f"{base_url}/{name}", asm_dir / pisa_subdir(pid) / name))
    log.info("%d ASU files and %d assembly files needed", pairs["pdb_id"].nunique(), len(pairs))

    session = make_session(pool_size=workers)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = [(t, ex.submit(download_file, session, t[3], t[4], gzip_is_valid)) for t in tasks]
        results = [(t, fut.result()) for t, fut in tqdm(futures, desc="mmCIF", unit="file")]

    manifest = ensure_parent(paths["structures_manifest"])
    tmp = manifest.with_name(manifest.name + ".part")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["kind", "pdb_id", "assembly_id", "path", "status"])
        for (kind, pid, asm, _url, dest), status in results:
            w.writerow([kind, pid, asm, dest, status])
    tmp.replace(manifest)

    counts = pd.Series([s.split(":")[0] for _, s in results]).value_counts().to_dict()
    log.info("Download summary: %s", counts)
    if counts.get("error") or counts.get("missing"):
        log.warning("Some files are missing or failed (see %s); rerun this step to retry errors.", manifest)


if __name__ == "__main__":
    main()
