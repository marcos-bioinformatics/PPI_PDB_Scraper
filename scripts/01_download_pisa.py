#!/usr/bin/env python3
"""Step 01 - Read the PDB ID lists and download their PISA data from PDBe.

Input : every .txt / .txt.gz in `input_dir` (RCSB Advanced Search exports).
Output: paths.ids             canonical, de-duplicated ID list
        paths.ids_rejected    tokens that are not valid PDB IDs
        paths.pisa_dir/       {subdir}/{id}_assembly{n}.json.gz and
                              {id}_assembly{n}_interfaces.json.gz (PISA-lite,
                              gzipped on download: ~10x smaller)
        paths.pisa_manifest   one row per (pdb_id, assembly_id) found in PDBe, with
                              status and the PDBe file date, plus IDs without PISA data

PDBe stores the PISA bulk files by the two middle characters of the ID
(1abc -> ab/). Each needed folder listing is fetched once, so every assembly of
every listed entry is found without guessing assembly numbers.

The PDBe dump can keep files of assemblies that no longer exist after an entry
was revised. With pisa.only_current_assemblies, the current assembly IDs of
each entry are read from the RCSB Data API (GraphQL, batches of 1000 IDs) and
outdated assemblies are listed as "stale" instead of downloaded.

Downloads are atomic and resumable: valid files already on disk are skipped.
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import requests
from tqdm import tqdm

from common import (atomic_write_text, download_file, ensure_parent, get_logger, json_is_valid,
                    make_session, parse_args, pisa_subdir, read_id_lists)

# Apache index rows: <a href="1abc_assembly1_interfaces.json">...</a> 2024-09-23 12:09
LISTING_RE = re.compile(
    r'href="([0-9a-z]{4})_assembly(\d+)(_interfaces)?\.json".*?(\d{4}-\d{2}-\d{2} \d{2}:\d{2})')

GRAPHQL_QUERY = """query($ids: [String!]!) {
  entries(entry_ids: $ids) { rcsb_id rcsb_entry_container_identifiers { assembly_ids } }
}"""

log = get_logger("01_download_pisa")


def fetch_listing(session, base_url: str, subdir: str):
    """Return {(pdb_id, assembly_id): {kind: file_date}} for one PDBe sub-folder."""
    try:
        r = session.get(f"{base_url}/{subdir}/", timeout=(10, 120))
    except requests.RequestException as exc:
        return subdir, {}, f"{type(exc).__name__}: {exc}"
    if r.status_code == 404:
        return subdir, {}, None
    if r.status_code != 200:
        return subdir, {}, f"HTTP {r.status_code}"
    found = defaultdict(dict)
    for pdb_id, asm, itf, date in LISTING_RE.findall(r.text):
        found[(pdb_id, int(asm))]["interfaces" if itf else "assembly"] = date
    return subdir, found, None


def current_assemblies(session, url: str, ids: List[str], batch: int = 1000) -> Dict[str, Set[int]]:
    """Current assembly IDs of each entry according to the RCSB Data API."""
    out: Dict[str, Set[int]] = {}
    for start in tqdm(range(0, len(ids), batch), desc="RCSB assemblies", unit="batch"):
        chunk = ids[start:start + batch]
        r = session.post(url, json={"query": GRAPHQL_QUERY, "variables": {"ids": chunk}}, timeout=(10, 300))
        if r.status_code != 200:
            raise SystemExit(f"RCSB GraphQL returned HTTP {r.status_code}; rerun this step to retry.")
        payload = r.json()
        if payload.get("errors"):
            raise SystemExit(f"RCSB GraphQL error: {payload['errors'][:1]}")
        for entry in payload["data"]["entries"] or []:
            asm_ids = (entry.get("rcsb_entry_container_identifiers") or {}).get("assembly_ids") or []
            out[entry["rcsb_id"].lower()] = {int(a) for a in asm_ids if str(a).isdigit()}
    return out


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

    by_subdir = defaultdict(set)
    for pid in ids:
        by_subdir[pisa_subdir(pid)].add(pid)

    session = make_session(pool_size=workers)

    # 1) Folder listings -> which (pdb, assembly) files exist for our IDs, and their dates
    available: Dict[Tuple[str, int], Dict[str, str]] = {}
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
                         f"(e.g. {next(iter(listing_errors.items()))}); rerun this step to retry.")

    # 2) Optionally skip assemblies that no longer exist in the current entry
    current: Optional[Dict[str, Set[int]]] = None
    if pisa.get("only_current_assemblies", True):
        current = current_assemblies(session, pisa["rcsb_graphql_url"], sorted({p for p, _ in available}))

    def is_current(pid: str, asm: int) -> bool:
        return current is None or asm in current.get(pid, set())

    # 3) Download both JSON files of every (current) assembly
    tasks = []
    for (pid, asm), kinds in sorted(available.items()):
        if not is_current(pid, asm):
            continue
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

    # 4) Manifest (also the Snakemake "done" marker, so it is written last)
    rows = []
    for (pid, asm), kinds in sorted(available.items()):
        if not is_current(pid, asm):
            state = "stale"
        else:
            state = "ok" if status.get((pid, asm, "interfaces"), "").startswith(("ok", "skip")) else "error"
        rows.append((pid, asm, state,
                     status.get((pid, asm, "interfaces"), ""), status.get((pid, asm, "assembly"), ""),
                     kinds.get("interfaces", ""), kinds.get("assembly", "")))
    with_data = {pid for pid, _ in available}
    rows += [(pid, "", "no_pisa_data", "", "", "", "") for pid in ids if pid not in with_data]

    manifest = ensure_parent(paths["pisa_manifest"])
    tmp = manifest.with_name(manifest.name + ".part")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["pdb_id", "assembly_id", "status", "interfaces_download", "assembly_download",
                    "interfaces_pdbe_date", "assembly_pdbe_date"])
        w.writerows(rows)
    tmp.replace(manifest)

    counts = defaultdict(int)
    for r in rows:
        counts[r[2]] += 1
    log.info("%d IDs | PISA assemblies: %s", len(ids), dict(counts))
    if counts["error"]:
        log.warning("Some downloads failed (status=error in %s); rerun this step to retry them.", manifest)


if __name__ == "__main__":
    main()
