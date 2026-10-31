#!/usr/bin/env python3
"""Step 01 - Read the PDB ID lists and download their PISA data from PDBe.

Input : every .txt / .txt.gz in `input_dir` (RCSB Advanced Search exports).
Output: paths.ids             canonical, de-duplicated ID list
        paths.ids_rejected    tokens that are not valid PDB IDs
        paths.pisa_dir/       {subdir}/{id}_assembly{n}.json.gz and
                              {id}_assembly{n}_interfaces.json.gz (PISA-lite,
                              gzipped on download: ~10x smaller)
        paths.pisa_manifest   one row per (pdb_id, assembly_id) + IDs without PISA data

PDBe stores the PISA bulk files by the two middle characters of the ID
(1abc -> ab/). Each needed folder listing is fetched once, so every assembly of
every listed entry is found without guessing assembly numbers. Downloads are
atomic and resumable: valid files already on disk are skipped.
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from tqdm import tqdm

from common import (atomic_write_text, download_file, ensure_parent, get_logger, json_is_valid,
                    make_session, parse_args, pisa_subdir, read_id_lists)

LISTING_RE = re.compile(r'href="([0-9a-z]{4})_assembly(\d+)(_interfaces)?\.json"')

log = get_logger("01_download_pisa")


def fetch_listing(session, base_url: str, subdir: str):
    """Return {(pdb_id, assembly_id): set(kinds)} for one PDBe sub-folder."""
    r = session.get(f"{base_url}/{subdir}/", timeout=(10, 120))
    if r.status_code == 404:
        return subdir, {}, None
    if r.status_code != 200:
        return subdir, {}, f"HTTP {r.status_code}"
    found = defaultdict(set)
    for pdb_id, asm, itf in LISTING_RE.findall(r.text):
        found[(pdb_id, int(asm))].add("interfaces" if itf else "assembly")
    return subdir, found, None


def main() -> None:
    cfg = parse_args(__doc__)
    paths = cfg["paths"]
    base_url = cfg["pisa"]["base_url"].rstrip("/")
    pisa_dir = Path(paths["pisa_dir"])
    workers = int(cfg["download_workers"])

    ids, rejected = read_id_lists(cfg["input_dir"])
    ids = sorted(ids)
    atomic_write_text(paths["ids"], "\n".join(ids) + "\n")
    atomic_write_text(paths["ids_rejected"],
                      "file\ttoken\n" + "".join(f"{f}\t{t}\n" for f, t in rejected))
    log.info("%d unique PDB IDs (%d tokens rejected)", len(ids), len(rejected))

    by_subdir = defaultdict(set)
    for pid in ids:
        by_subdir[pisa_subdir(pid)].add(pid)

    session = make_session(pool_size=workers)

    # 1) Folder listings -> which (pdb, assembly) files exist for our IDs
    available = {}
    listing_errors = {}
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = [ex.submit(fetch_listing, session, base_url, sd) for sd in sorted(by_subdir)]
        for fut in tqdm(futures, desc="PISA listings", unit="dir"):
            subdir, found, err = fut.result()
            if err:
                listing_errors[subdir] = err
                continue
            wanted = by_subdir[subdir]
            available.update({k: v for k, v in found.items() if k[0] in wanted})
    if listing_errors:
        raise SystemExit(f"Could not list {len(listing_errors)} PDBe folders "
                         f"(e.g. {next(iter(listing_errors.items()))}); rerun to retry.")

    # 2) Download both JSON files of every available assembly
    tasks = []
    for (pid, asm), kinds in sorted(available.items()):
        for kind in sorted(kinds):
            name = f"{pid}_assembly{asm}{'_interfaces' if kind == 'interfaces' else ''}.json"
            tasks.append((pid, asm, kind, f"{base_url}/{pisa_subdir(pid)}/{name}",
                          pisa_dir / pisa_subdir(pid) / f"{name}.gz"))

    status = {}
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(download_file, session, url, dest, json_is_valid, gzip_output=True): (pid, asm, kind)
                   for pid, asm, kind, url, dest in tasks}
        for fut in tqdm(futures, desc="PISA JSON", unit="file"):
            status[futures[fut]] = fut.result()

    # 3) Manifest (also the Snakemake "done" marker, so it is written last)
    rows = []
    for (pid, asm), kinds in sorted(available.items()):
        st_itf = status.get((pid, asm, "interfaces"), "not_in_pdbe")
        st_asm = status.get((pid, asm, "assembly"), "not_in_pdbe")
        ok = st_itf in ("ok", "skip")
        rows.append((pid, asm, "ok" if ok else "error", st_itf, st_asm))
    with_data = {pid for pid, _ in available}
    rows += [(pid, "", "no_pisa_data", "", "") for pid in ids if pid not in with_data]

    manifest = ensure_parent(paths["pisa_manifest"])
    tmp = manifest.with_name(manifest.name + ".part")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pdb_id", "assembly_id", "status", "interfaces_json", "assembly_json"])
        w.writerows(rows)
    tmp.replace(manifest)

    n_err = sum(1 for r in rows if r[2] == "error")
    log.info("%d entries with PISA data, %d assemblies, %d IDs without PISA data, %d download errors",
             len(with_data), len(available), len(ids) - len(with_data), n_err)
    if n_err:
        log.warning("Some downloads failed (status=error in %s); rerun this step to retry them.", manifest)


if __name__ == "__main__":
    main()
