#!/usr/bin/env python3
"""Provenance record of a run (written after the final outputs).

Output: paths.provenance (JSON) with
  run_date      when this record was written (UTC)
  code          git commit of the pipeline and whether it had uncommitted changes
  input_query   scope.input_query (the JSON is embedded when it is a file) and its date
  input_lists   every ID list of input_dir with its SHA-256 and number of valid IDs
  config        the resolved configuration of the run
  software      versions of Python, R, MMseqs2 and the main Python packages

Per-file provenance stays in the manifests: PISA file dates (paths.pisa_manifest),
download dates of the assembly files (paths.structures_manifest), and release
dates and current revision of every entry (paths.entry_metadata).
"""

from __future__ import annotations

import hashlib
import json
import platform
import re
import subprocess
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
from typing import Any, Dict, List, Optional

from common import _read_text, get_logger, normalize_pdb_id, parse_args, write_json

log = get_logger("provenance")

PACKAGES = ["pandas", "numpy", "scipy", "biopython", "requests", "PyYAML", "tqdm", "snakemake"]


def run(cmd: List[str]) -> Optional[str]:
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        return (out.stdout or out.stderr).strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def query_record(scope: Dict[str, Any]) -> Dict[str, Any]:
    query = scope.get("input_query")
    rec: Dict[str, Any] = {"query": query, "date": None if scope.get("input_query_date") is None
                           else str(scope["input_query_date"])}
    if query and Path(str(query)).is_file():
        text = Path(str(query)).read_text(encoding="utf-8")
        try:
            rec["query"] = {"file": str(query), "content": json.loads(text)}
        except json.JSONDecodeError:
            rec["query"] = {"file": str(query), "content": text}
    return rec


def main() -> None:
    cfg = parse_args(__doc__)
    input_dir = Path(cfg["input_dir"])
    lists = []
    for path in sorted(p for p in input_dir.iterdir() if p.name.endswith((".txt", ".txt.gz"))):
        ids = {normalize_pdb_id(t) for t in re.split(r"[,;\s]+", _read_text(path)) if t}
        lists.append({"file": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                      "ids": len(ids - {None})})

    commit = run(["git", "rev-parse", "HEAD"])
    dirty = run(["git", "status", "--porcelain", "--untracked-files=no"])
    mmseqs = run(["mmseqs", "version"])
    r_version = run(["Rscript", "-e", "cat(R.version.string)"])
    packages = {}
    for name in PACKAGES:
        try:
            packages[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            packages[name] = None

    record = {
        "run_date": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "code": {"git_commit": commit, "uncommitted_changes": bool(dirty) if commit else None},
        "input_query": query_record(cfg.get("scope") or {}),
        "input_lists": lists,
        "config": cfg,
        "software": {"python": platform.python_version(), "R": r_version, "mmseqs2": mmseqs,
                     "packages": packages},
    }
    write_json(cfg["paths"]["provenance"], json.loads(json.dumps(record, default=str)))
    log.info("Provenance written to %s", cfg["paths"]["provenance"])


if __name__ == "__main__":
    main()
