#!/usr/bin/env python3
"""Step 04 - Download the assembly mmCIF files and entry metadata of the selected pairs.

Input : paths.needed_pairs   (pdb_id, assembly_id) pairs from step 03
Output: paths.assembly_dir/{subdir}/{id}-assembly{n}.cif.gz   every needed assembly
        paths.structures_manifest                             status of every file
        paths.entry_metadata                                  method, resolution, R-free, dates, revision
        paths.entity_metadata                                 membrane and antibody annotations per entity

Only the biological assemblies that PISA analysed are downloaded. Files are
stored gzipped, fully checked for gzip integrity when downloaded and written
atomically; files already on disk are kept, so the step can be rerun to resume
or to retry failures. Failed transfers are retried (download_retries) and the
step stops if more than max_failed_fraction of the files fail.

Entry-level metadata comes from the RCSB Data API (GraphQL, batches of 1000
IDs), so the asymmetric-unit files do not need to be downloaded. An entry can
be solved by several methods (e.g. joint X-ray and neutron refinement): if one
of them is quality_filter.exptl_method, that is the entry's method and the
resolution and R-free come from the _refine record of that method
(_refine.pdbx_refine_id), as RCSB does; otherwise the first method and the
first _refine record are used. exptl_methods keeps every method. The release
dates and the current revision of the entry are also recorded.

Entity-level annotations (RCSB, one row per polymer entity):
  uniprot_ids        UniProt accessions of the entity (SIFTS mapping, ";"-separated)
  membrane_sources   membrane-protein annotation sources (PDBTM, OPM, MemProtMD, mpstruc)
  is_antibody        the entity has a SAbDab antibody annotation
"""

from __future__ import annotations

import csv
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from common import (check_failures, download_file, ensure_parent, get_logger, gzip_is_valid, gzip_looks_valid,
                    make_session, parse_args, pisa_subdir, rcsb_entries, run_threaded)

log = get_logger("04_download_structures")

ENTRY_FIELDS = ("exptl { method } refine { pdbx_refine_id ls_d_res_high ls_R_factor_R_free } "
                "struct { title } struct_keywords { pdbx_keywords } "
                "rcsb_accession_info { deposit_date initial_release_date revision_date major_revision minor_revision } "
                "polymer_entities { rcsb_polymer_entity_container_identifiers { entity_id uniprot_ids } "
                "rcsb_polymer_entity_annotation { type } rcsb_polymer_entity_feature { type } }")
ENTRY_COLUMNS = ["pdb_id", "exptl_method", "exptl_methods", "resolution", "r_free", "title", "keywords",
                 "deposit_date", "initial_release_date", "revision_date", "revision"]
MEMBRANE_SOURCES = {"PDBTM", "OPM", "MemProtMD", "mpstruc"}


def pick_method(entry, wanted):
    """(method, all methods, refine record) of an entry: the wanted method and
    its _refine record if the entry has it, else the first of each."""
    methods = [m.get("method") for m in entry.get("exptl") or [] if m.get("method")]
    refines = entry.get("refine") or []
    method = wanted if wanted in methods else (methods[0] if methods else None)
    match = [r for r in refines if r.get("pdbx_refine_id") == method]
    if match:
        refine = match[0]
    elif len(refines) == 1 or (refines and len(methods) <= 1):
        refine = refines[0]  # single method: the refinement is of that method
    else:
        refine = {}
    return method, ";".join(methods), refine


def day(value):
    return value[:10] if value else None


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

    # 1) Entry and entity metadata from the RCSB Data API
    wanted = cfg["quality_filter"].get("exptl_method")
    rows, ent_rows = [], []
    for e in rcsb_entries(session, rcsb["graphql_url"], entries, ENTRY_FIELDS):
        pid = e["rcsb_id"].lower()
        method, methods, refine = pick_method(e, wanted)
        acc = e.get("rcsb_accession_info") or {}
        rows.append({
            "pdb_id": pid,
            "exptl_method": method,
            "exptl_methods": methods,
            "resolution": refine.get("ls_d_res_high"),
            "r_free": refine.get("ls_R_factor_R_free"),
            "title": (e.get("struct") or {}).get("title"),
            "keywords": (e.get("struct_keywords") or {}).get("pdbx_keywords"),
            "deposit_date": day(acc.get("deposit_date")),
            "initial_release_date": day(acc.get("initial_release_date")),
            "revision_date": day(acc.get("revision_date")),
            "revision": (f"{acc['major_revision']}.{acc['minor_revision']}"
                         if acc.get("major_revision") is not None else None),
        })
        for pe in e.get("polymer_entities") or []:
            annotations = {a.get("type") for a in pe.get("rcsb_polymer_entity_annotation") or []}
            features = {f.get("type") or "" for f in pe.get("rcsb_polymer_entity_feature") or []}
            ids = pe["rcsb_polymer_entity_container_identifiers"]
            ent_rows.append({
                "pdb_id": pid,
                "entity_id": str(ids["entity_id"]),
                "uniprot_ids": ";".join(ids.get("uniprot_ids") or []),
                "membrane_sources": ";".join(sorted(annotations & MEMBRANE_SOURCES)),
                "is_antibody": any(f.startswith("SABDAB_ANTIBODY") for f in features),
            })
    meta = pd.DataFrame(rows, columns=ENTRY_COLUMNS)
    meta.sort_values("pdb_id").to_csv(ensure_parent(paths["entry_metadata"]), index=False)
    ent = pd.DataFrame(ent_rows, columns=["pdb_id", "entity_id", "uniprot_ids", "membrane_sources", "is_antibody"])
    ent.sort_values(["pdb_id", "entity_id"]).to_csv(ensure_parent(paths["entity_metadata"]), index=False)
    multi = meta["exptl_methods"].str.contains(";", na=False)
    log.info("Entry metadata for %d of %d entries (%d solved by several methods); %d polymer entities "
             "(%d with UniProt, %d membrane, %d antibody)", len(meta), len(entries), multi.sum(), len(ent),
             (ent["uniprot_ids"] != "").sum(), (ent["membrane_sources"] != "").sum(), ent["is_antibody"].sum())

    # 2) Assembly files
    tasks = []
    for pid, asm in zip(pairs["pdb_id"], pairs["assembly_id"]):
        name = f"{pid}-assembly{asm}.cif.gz"
        tasks.append((pid, asm, f"{base_url}/{name}", asm_dir / pisa_subdir(pid) / name))
    log.info("%d assembly files needed", len(tasks))

    attempts = int(cfg["download_retries"])
    statuses = run_threaded(lambda t: download_file(session, t[2], t[3], gzip_is_valid,
                                                    exists_check=gzip_looks_valid, attempts=attempts),
                            tasks, workers, "mmCIF")
    results = list(zip(tasks, statuses))

    manifest = ensure_parent(paths["structures_manifest"])
    tmp = manifest.with_name(manifest.name + ".part")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["pdb_id", "assembly_id", "path", "status", "file_date"])
        for (pid, asm, _url, dest), status in results:
            # file_date: when the file on disk was downloaded (UTC)
            date = (datetime.fromtimestamp(dest.stat().st_mtime, timezone.utc).strftime("%Y-%m-%d %H:%M")
                    if dest.is_file() else "")
            w.writerow([pid, asm, dest, status, date])
    tmp.replace(manifest)

    counts = pd.Series([s.split(":")[0] for _, s in results]).value_counts().to_dict()
    log.info("Download summary: %s", counts)
    if counts.get("missing"):
        log.warning("%d assembly files do not exist at RCSB (see %s).", counts["missing"], manifest)
    check_failures(counts.get("error", 0), len(tasks), float(cfg["max_failed_fraction"]), "mmCIF downloads", log)


if __name__ == "__main__":
    main()
