#!/usr/bin/env python3
"""Step 05 - Audit that every file the next steps need is on disk.

Input : paths.ids, paths.pisa_manifest, paths.master_pisa, paths.needed_pairs,
        paths.structures_manifest, paths.entry_metadata
Output: paths.audit (JSON summary) and reports_dir/05_missing_structures.tsv

The audit follows the IDs through the acquisition steps (list -> PISA ->
candidate pairs -> mmCIF files and entry metadata) and checks on disk that every
needed assembly file exists. It never stops the pipeline: rows whose structure
or metadata is missing are dropped and counted by step 06.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from common import get_logger, parse_args, read_ids_file, write_json

log = get_logger("05_audit_downloads")


def main() -> None:
    cfg = parse_args(__doc__)
    paths = cfg["paths"]

    ids = read_ids_file(paths["ids"])
    pisa = pd.read_csv(paths["pisa_manifest"], sep="\t", dtype={"pdb_id": str})
    master = pd.read_csv(paths["master_pisa"], usecols=["pdb_id", "assembly_id"], dtype={"pdb_id": str})
    needed = pd.read_csv(paths["needed_pairs"], dtype={"pdb_id": str})
    files = pd.read_csv(paths["structures_manifest"], sep="\t", dtype={"pdb_id": str, "assembly_id": str})
    meta = pd.read_csv(paths["entry_metadata"], dtype={"pdb_id": str})
    selection = cfg["assemblies"]["selection"]
    if selection == "first":
        first = pisa.dropna(subset=["first_assembly"])[["pdb_id", "first_assembly"]].drop_duplicates()
        first = first.rename(columns={"first_assembly": "assembly_id"}).astype({"assembly_id": int})
        expected = master.merge(first, on=["pdb_id", "assembly_id"])
    else:
        expected = master

    files["on_disk"] = [Path(p).is_file() for p in files["path"]]
    missing = files[~files["on_disk"]]
    missing_path = Path(paths["reports_dir"]) / "05_missing_structures.tsv"
    missing.to_csv(missing_path, sep="\t", index=False)

    report = {
        "ids_in_list": len(ids),
        "pisa": {
            "entries_with_data": int(pisa.loc[pisa["status"] == "ok", "pdb_id"].nunique()),
            "entries_without_data": int((pisa["status"] == "no_pisa_data").sum()),
            "assemblies_ok": int((pisa["status"] == "ok").sum()),
            "assemblies_stale": int((pisa["status"] == "stale").sum()),
            "entries_obsolete": int(pisa.loc[pisa["status"] == "obsolete", "pdb_id"].nunique()),
            "assemblies_error": int((pisa["status"] == "error").sum()),
        },
        "candidate_pairs": {
            "rows": int(len(master)),
            "entries": int(master["pdb_id"].nunique()),
            "assemblies": int(len(master.drop_duplicates())),
        },
        "needed_pairs": int(len(needed)),
        "assembly_selection": selection,
        # entries whose candidate pairs are all in assemblies other than the first
        "entries_without_candidates_in_first_assembly": int(len(set(master["pdb_id"]) - set(expected["pdb_id"]))),
        # entries with candidate pairs that step 03 should have kept but did not (should be 0)
        "entries_lost_in_step_03": int(len(set(expected["pdb_id"]) - set(needed["pdb_id"]))),
        "assembly_files": {
            "needed": int(len(files)),
            "on_disk": int(files["on_disk"].sum()),
            "status": files["status"].str.split(":").str[0].value_counts().to_dict(),
        },
        "entry_metadata": {
            "needed": int(needed["pdb_id"].nunique()),
            "found": int(meta["pdb_id"].isin(set(needed["pdb_id"])).sum()),
        },
        "missing_files_tsv": str(missing_path),
    }
    write_json(paths["audit"], report)
    if report["entries_lost_in_step_03"]:
        log.warning("%d entries with candidate pairs were dropped by step 03", report["entries_lost_in_step_03"])
    log.info("Audit written to %s (%d needed files missing on disk)", paths["audit"], len(missing))


if __name__ == "__main__":
    main()
