#!/usr/bin/env python3
"""Step 01 - Read the PDB ID lists and download their PISA data from PDBe.

Input : every .txt / .txt.gz in `input_dir` (RCSB Advanced Search exports).
Output: paths.ids             canonical, de-duplicated ID list
        paths.ids_rejected    tokens that are not valid PDB IDs
        paths.pisa_dir/       {subdir}/{id}_assembly{n}_interfaces.json.gz (PISA-lite,
                              gzipped on download: ~10x smaller)
        paths.pisa_manifest   one row per (pdb_id, assembly_id) found in PDBe, with
                              status and the PDBe file date, plus IDs without PISA data

PDBe stores the PISA bulk files by the two middle characters of the ID
(1abc -> ab/). Each needed folder listing is fetched once, so every assembly of
every listed entry is found without guessing assembly numbers. The interfaces
file also carries the assembly-level metrics, so the separate
{id}_assembly{n}.json is not needed.

The PDBe dump can keep files of assemblies that no longer exist after an entry
was revised. With pisa.only_current_assemblies, the current assembly IDs of
each entry are read from the RCSB Data API (GraphQL, batches of 1000 IDs):
outdated assemblies are listed as "stale" and entries RCSB no longer serves as
"obsolete", and neither is downloaded.

Downloads are atomic and resumable. A file on disk is kept unless PDBe has a
newer version (PISA is recomputed when an entry is revised). The step stops
with an error if the ID lists give no IDs, if too many entries have no PISA
data (pisa.max_missing_fraction, e.g. a changed listing format) or if more than
max_failed_fraction of the downloads fail after retries.
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import requests

from common import (atomic_write_text, check_failures, download_file, ensure_parent, get_logger, json_is_valid,
                    make_session, parse_args, pisa_subdir, rcsb_entries, read_id_lists, run_threaded)

# Apache index rows: <a href="1abc_assembly1_interfaces.json">...</a> 2024-09-23 12:09
LISTING_RE = re.compile(r'href="([0-9a-z]{4})_assembly(\d+)_interfaces\.json".*?(\d{4}-\d{2}-\d{2} \d{2}:\d{2})')

log = get_logger("01_download_pisa")


def fetch_listing(session, base_url: str, subdir: str):
    """Return (subdir, {(pdb_id, assembly_id): file date}, error) for one PDBe sub-folder."""
    try:
        r = session.get(f"{base_url}/{subdir}/", timeout=(10, 120))
    except requests.RequestException as exc:
        return subdir, {}, f"{type(exc).__name__}: {exc}"
    if r.status_code == 404:
        return subdir, {}, None
    if r.status_code != 200:
        return subdir, {}, f"HTTP {r.status_code}"
    found = {(pdb_id, int(asm)): date for pdb_id, asm, date in LISTING_RE.findall(r.text)}
    if not found and "_interfaces.json" in r.text:
        return subdir, {}, "listing format not recognised"
    return subdir, found, None


def current_assemblies(session, url: str, ids: List[str]) -> Dict[str, Set[int]]:
    """Current assembly IDs of each entry according to the RCSB Data API."""
    out: Dict[str, Set[int]] = {}
    for entry in rcsb_entries(session, url, ids, "rcsb_entry_container_identifiers { assembly_ids }"):
        asm_ids = (entry.get("rcsb_entry_container_identifiers") or {}).get("assembly_ids") or []
        out[entry["rcsb_id"].lower()] = {int(a) for a in asm_ids if str(a).isdigit()}
    return out


def epoch(date: str) -> float:
    return datetime.strptime(date, "%Y-%m-%d %H:%M").replace(tzinfo=timezone.utc).timestamp()


def main() -> None:
    cfg = parse_args(__doc__)
    paths, pisa = cfg["paths"], cfg["pisa"]
    base_url = pisa["base_url"].rstrip("/")
    pisa_dir = Path(paths["pisa_dir"])
    workers = int(cfg["download_workers"])

    ids, rejected = read_id_lists(cfg["input_dir"])
    ids = sorted(ids)
    atomic_write_text(paths["ids"], "\n".join(ids) + "\n")
    atomic_write_text(paths["ids_rejected"],
                      "file\ttoken\n" + "".join(f"{f}\t{t}\n" for f, t in rejected))
    log.info("%d unique PDB IDs (%d tokens rejected)", len(ids), len(rejected))
    if not ids:
        raise SystemExit(f"No valid PDB IDs in {cfg['input_dir']}/ (see {paths['ids_rejected']}).")

    by_subdir = defaultdict(set)
    for pid in ids:
        by_subdir[pisa_subdir(pid)].add(pid)
    session = make_session(pool_size=workers)

    # 1) Folder listings -> which assemblies exist for our IDs, and their dates
    available: Dict[Tuple[str, int], str] = {}
    listing_errors = {}
    for subdir, found, err in run_threaded(lambda sd: fetch_listing(session, base_url, sd),
                                           sorted(by_subdir), workers, "PISA listings"):
        if err:
            listing_errors[subdir] = err
            continue
        available.update({k: v for k, v in found.items() if k[0] in by_subdir[subdir]})
    if listing_errors:
        raise SystemExit(f"Could not list {len(listing_errors)} PDBe folders "
                         f"(e.g. {next(iter(listing_errors.items()))}); rerun this step to retry.")

    # 2) Optionally skip assemblies (and entries) that RCSB no longer has
    current: Optional[Dict[str, Set[int]]] = None
    if pisa.get("only_current_assemblies", True):
        current = current_assemblies(session, cfg["rcsb"]["graphql_url"], sorted({p for p, _ in available}))

    def state_of(pid: str, asm: int) -> str:
        if current is None:
            return "current"
        if pid not in current:
            return "obsolete"
        return "current" if asm in current[pid] else "stale"

    # 3) Download the interfaces JSON of every current assembly
    todo = sorted(k for k in available if state_of(*k) == "current")

    def fetch(key: Tuple[str, int]) -> str:
        pid, asm = key
        name = f"{pid}_assembly{asm}_interfaces.json"
        return download_file(session, f"{base_url}/{pisa_subdir(pid)}/{name}",
                             pisa_dir / pisa_subdir(pid) / f"{name}.gz", json_is_valid,
                             gzip_output=True, remote_mtime=epoch(available[key]),
                             attempts=int(cfg["download_retries"]))

    status = dict(zip(todo, run_threaded(fetch, todo, workers, "PISA JSON")))

    # 4) Manifest (also the Snakemake "done" marker, so it is written last)
    rows = []
    for key in sorted(available):
        pid, asm = key
        state = state_of(pid, asm)
        if state == "current":
            state = "ok" if status[key] in ("ok", "skip") else "error"
        rows.append((pid, asm, state, status.get(key, ""), available[key]))
    with_data = {pid for pid, _ in available}
    rows += [(pid, "", "no_pisa_data", "", "") for pid in ids if pid not in with_data]

    manifest = ensure_parent(paths["pisa_manifest"])
    tmp = manifest.with_name(manifest.name + ".part")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["pdb_id", "assembly_id", "status", "download", "pdbe_date"])
        w.writerows(rows)
    tmp.replace(manifest)

    counts: Dict[str, int] = defaultdict(int)
    for r in rows:
        counts[r[2]] += 1
    log.info("%d IDs | PISA assemblies: %s", len(ids), dict(counts))
    missing_fraction = counts["no_pisa_data"] / len(ids)
    if missing_fraction > float(pisa["max_missing_fraction"]):
        raise SystemExit(f"{counts['no_pisa_data']} of {len(ids)} IDs have no PISA data (more than "
                         f"{float(pisa['max_missing_fraction']):.0%}); check the ID lists or the PDBe listing format.")
    check_failures(counts["error"], len(todo), float(cfg["max_failed_fraction"]), "PISA downloads", log)


if __name__ == "__main__":
    main()
