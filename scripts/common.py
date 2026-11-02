"""Shared helpers for the BeastAffinity pipeline scripts.

Every numbered script imports this module (it lives next to them in
``scripts/``), so configuration loading, PDB ID handling, HTTP downloads and
atomic file writes behave the same way in every step.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import logging
import os
import re
import time
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import requests
import yaml
from requests.adapters import HTTPAdapter
from tqdm import tqdm
from urllib3.util.retry import Retry

DEFAULT_CONFIG = "config/config.yaml"

# Classic 4-character PDB IDs (first character is a digit).
# Extended IDs (pdb_XXXXXXXX) will be handled in v2.
PDB_ID_RE = re.compile(r"^[0-9][a-z0-9]{3}$")
# Entity / assembly / chain forms exported by RCSB searches: 4HHB_1, 4HHB-1, 4HHB.A
PDB_ID_SUFFIX_RE = re.compile(r"^([0-9][a-z0-9]{3})[_\-.][a-z0-9]+$")
NUMBER_RE = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
def _deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    out = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = value
    return out


def _coerce_numbers(obj: Any) -> Any:
    """YAML 1.1 reads values such as 1e-3 as strings: turn numeric strings into numbers."""
    if isinstance(obj, dict):
        return {k: _coerce_numbers(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_coerce_numbers(v) for v in obj]
    if isinstance(obj, str) and NUMBER_RE.match(obj.strip()):
        return float(obj) if re.search(r"[.eE]", obj) else int(obj)
    return obj


def _load_raw(path: str | Path, seen: Tuple[str, ...] = ()) -> Dict[str, Any]:
    path = Path(path)
    if str(path.resolve()) in seen:
        raise SystemExit(f"Config inheritance cycle: {' -> '.join(seen)} -> {path}")
    with path.open("r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh) or {}
    base = cfg.pop("base", None)
    if not base:
        return cfg
    base_path = Path(base)
    if not base_path.is_file() and (path.parent / base).is_file():
        base_path = path.parent / base  # relative to the config file instead of the working directory
    return _deep_merge(_load_raw(base_path, seen + (str(path.resolve()),)), cfg)


def load_config(path: str | Path = DEFAULT_CONFIG) -> Dict[str, Any]:
    """Load a YAML config. A ``base:`` key points to a parent config whose
    values are inherited and can be overridden (used by config/test.yaml).
    ``{data_dir}`` in paths is replaced by the configured data_dir."""
    cfg = _coerce_numbers(_load_raw(path))
    data_dir = str(cfg.get("data_dir", "data"))
    cfg["paths"] = {k: str(v).replace("{data_dir}", data_dir) for k, v in cfg.get("paths", {}).items()}
    return cfg


def parse_args(description: str) -> Dict[str, Any]:
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--config", default=DEFAULT_CONFIG, help="YAML config file")
    args = parser.parse_args()
    return load_config(args.config)


def get_logger(name: str) -> logging.Logger:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )
    return logging.getLogger(name)


# -----------------------------------------------------------------------------
# PDB IDs
# -----------------------------------------------------------------------------
def normalize_pdb_id(token: str) -> Optional[str]:
    s = token.strip().strip("\"'").lower()
    if PDB_ID_RE.match(s):
        return s
    m = PDB_ID_SUFFIX_RE.match(s)
    return m.group(1) if m else None


def _read_text(path: Path) -> str:
    """Text of an ID list (.txt or .txt.gz), handling UTF-8 BOMs and UTF-16 files
    (e.g. saved from PowerShell)."""
    raw = gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes()
    if raw.startswith((b"\xff\xfe", b"\xfe\xff")):
        return raw.decode("utf-16", errors="replace")
    return raw.decode("utf-8-sig", errors="replace")


def read_id_lists(input_dir: str | Path) -> Tuple[List[str], List[Tuple[str, str]]]:
    """Read every .txt / .txt.gz file in ``input_dir`` (RCSB Advanced Search
    exports: IDs separated by commas, semicolons, whitespace or one per line;
    entity or assembly forms such as 4HHB_1 are reduced to the entry ID).

    Returns (unique IDs in first-seen order, rejected (file, token) pairs).
    """
    input_dir = Path(input_dir)
    files = sorted(p for p in input_dir.iterdir() if p.name.endswith((".txt", ".txt.gz")))
    if not files:
        raise SystemExit(f"No .txt or .txt.gz ID lists found in {input_dir}/")

    seen: Dict[str, None] = {}
    rejected: List[Tuple[str, str]] = []
    for path in files:
        for token in re.split(r"[,;\s]+", _read_text(path)):
            if not token:
                continue
            pid = normalize_pdb_id(token)
            if pid is None:
                rejected.append((path.name, token))
            else:
                seen.setdefault(pid, None)
    return list(seen), rejected


def read_ids_file(path: str | Path) -> List[str]:
    with Path(path).open("r", encoding="utf-8") as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def seq_uid(seq: str) -> str:
    """Stable identifier of a sequence: first 16 hex characters of its SHA-1."""
    return hashlib.sha1(seq.encode("utf-8")).hexdigest()[:16]


def homology_linked(hits, ccfg: Dict[str, Any]):
    """Boolean mask of MMseqs2 hits (columns qcov, tcov; searched with -c 0) that
    meet clustering.coverage. Shared by steps 09 and 14 so that clustering and
    audit link the same pairs. cov_mode 0 needs both coverages; cov_mode 1
    (target) and 2 (query) need either, because every pair is searched both ways
    and MMseqs2 may keep only a shorter local alignment in one direction."""
    c = float(ccfg["coverage"])
    if int(ccfg["cov_mode"]) == 0:
        return (hits["qcov"] >= c) & (hits["tcov"] >= c)
    return (hits["qcov"] >= c) | (hits["tcov"] >= c)


def pisa_subdir(pdb_id: str) -> str:
    """PDBe splits its archives by the two middle characters (1abc -> 'ab')."""
    return pdb_id[1:3]


# -----------------------------------------------------------------------------
# Files
# -----------------------------------------------------------------------------
def ensure_parent(path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def atomic_write_text(path: str | Path, text: str) -> None:
    path = ensure_parent(path)
    tmp = path.with_name(path.name + ".part")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def write_json(path: str | Path, obj: Any) -> None:
    atomic_write_text(path, json.dumps(obj, indent=2, ensure_ascii=False, default=str) + "\n")


def open_mmcif(path: str | Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("r", encoding="utf-8", errors="ignore")


def gzip_looks_valid(path: Path) -> bool:
    """Cheap check for files already on disk: non-trivial size and gzip magic bytes."""
    try:
        with path.open("rb") as fh:
            return path.stat().st_size > 20 and fh.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def gzip_is_valid(path: Path) -> bool:
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1 << 20):
                pass
        return True
    except (OSError, EOFError, zlib.error):
        return False


def read_json(path: str | Path) -> Any:
    path = Path(path)
    data = gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes()
    return json.loads(data)


def json_is_valid(path: Path) -> bool:
    try:
        read_json(path)
        return True
    except (OSError, EOFError, ValueError, zlib.error):
        return False


# -----------------------------------------------------------------------------
# HTTP
# -----------------------------------------------------------------------------
def make_session(pool_size: int, user_agent: str = "BeastAffinity-pipeline/1.0") -> requests.Session:
    retries = Retry(
        total=6,
        backoff_factor=1.0,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods={"GET", "HEAD", "POST"},  # the GraphQL POST queries are read-only
        respect_retry_after_header=True,
        raise_on_status=False,
    )
    session = requests.Session()
    adapter = HTTPAdapter(max_retries=retries, pool_connections=pool_size, pool_maxsize=pool_size)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    # Ask for the raw bytes: with gzip transfer-encoding the Content-Length
    # header refers to the compressed body and size checks become meaningless.
    session.headers.update({"User-Agent": user_agent, "Accept-Encoding": "identity"})
    return session


def _download_once(session: requests.Session, url: str, dest: Path, validator, timeout,
                   gzip_output: bool) -> str:
    tmp = dest.with_name(dest.name + ".part")
    try:
        with session.get(url, stream=True, timeout=timeout) as r:
            if r.status_code == 404:
                return "missing"
            if r.status_code != 200:
                return f"error:HTTP {r.status_code}"
            expected = r.headers.get("Content-Length")
            written = 0
            with tmp.open("wb") as fh:
                for chunk in r.iter_content(chunk_size=1 << 17):
                    fh.write(chunk)
                    written += len(chunk)
        if expected is not None and int(expected) != written:
            tmp.unlink(missing_ok=True)
            return f"error:size {written} != {expected}"
        if not validator(tmp):
            tmp.unlink(missing_ok=True)
            return "error:invalid content"
        if gzip_output:
            gz_tmp = dest.with_name(dest.name + ".gzpart")
            gz_tmp.write_bytes(gzip.compress(tmp.read_bytes()))
            tmp.unlink()
            tmp = gz_tmp
        os.replace(tmp, dest)
        return "ok"
    except (requests.RequestException, OSError) as exc:
        tmp.unlink(missing_ok=True)
        return f"error:{type(exc).__name__}: {exc}"


def download_file(session: requests.Session, url: str, dest: Path, validator, timeout=(10, 120),
                  gzip_output: bool = False, exists_check=None, remote_mtime: Optional[float] = None,
                  attempts: int = 4) -> str:
    """Download ``url`` to ``dest`` atomically (gzip-compressing it if ``gzip_output``).

    Returns "skip" (dest already valid and up to date), "ok", "missing" (HTTP
    404) or an error string starting with "error:". ``validator(path) -> bool``
    checks a fresh download; ``exists_check`` (default: validator) decides
    whether a file already on disk can be kept. With ``remote_mtime`` (epoch
    seconds of the server copy) an older local file is downloaded again, and a
    new file gets that modification time. Failed transfers (including bodies
    cut mid-way, which urllib3 does not retry) are retried ``attempts`` times.
    """
    if dest.is_file() and (exists_check or validator)(dest):
        if remote_mtime is None or dest.stat().st_mtime >= remote_mtime - 60:
            return "skip"
    ensure_parent(dest)
    status = "error:not attempted"
    for attempt in range(attempts):
        status = _download_once(session, url, dest, validator, timeout, gzip_output)
        if not status.startswith("error"):
            break
        time.sleep(2 ** attempt)
    if status == "ok" and remote_mtime is not None:
        os.utime(dest, (remote_mtime, remote_mtime))
    return status


def run_threaded(func: Callable, items: Sequence, workers: int, desc: str, chunk: int = 5000) -> List[Any]:
    """``func`` over ``items`` with a thread pool, in input order, submitting at
    most ``chunk`` tasks at a time (keeps memory flat for hundreds of thousands
    of downloads)."""
    out: List[Any] = []
    with ThreadPoolExecutor(max_workers=workers) as ex, tqdm(total=len(items), desc=desc, unit="file") as bar:
        for start in range(0, len(items), chunk):
            for result in ex.map(func, items[start:start + chunk]):
                out.append(result)
                bar.update(1)
    return out


def check_failures(n_failed: int, n_total: int, max_fraction: float, what: str, log: logging.Logger) -> None:
    """Stop the step (non-zero exit) when too many downloads failed; a rerun resumes."""
    if n_failed and n_total and n_failed / n_total > max_fraction:
        raise SystemExit(f"{n_failed} of {n_total} {what} failed (more than {max_fraction:.2%}); "
                         "rerun this step to retry them.")
    if n_failed:
        log.warning("%d of %d %s failed (within the tolerated %.2f%%).", n_failed, n_total, what, 100 * max_fraction)


def rcsb_entries(session: requests.Session, url: str, ids: List[str], fields: str,
                 batch: int = 1000, attempts: int = 4) -> List[Dict[str, Any]]:
    """Query the RCSB Data API (GraphQL) for ``fields`` of many entries, in batches.
    Results are not returned in input order: match them by ``rcsb_id``. Entries
    that RCSB no longer serves (obsolete IDs) are simply absent."""
    query = "query($ids: [String!]!) { entries(entry_ids: $ids) { rcsb_id %s } }" % fields
    out: List[Dict[str, Any]] = []
    for start in range(0, len(ids), batch):
        body = {"query": query, "variables": {"ids": ids[start:start + batch]}}
        problem = ""
        for attempt in range(attempts):
            try:
                r = session.post(url, json=body, timeout=(10, 300))
                payload = r.json() if r.status_code == 200 else None
                if payload is not None and not payload.get("errors"):
                    out.extend(payload["data"]["entries"] or [])
                    break
                problem = f"HTTP {r.status_code}" if payload is None else f"GraphQL error {payload['errors'][:1]}"
            except (requests.RequestException, ValueError) as exc:
                problem = f"{type(exc).__name__}: {exc}"
            time.sleep(2 ** attempt)
        else:
            raise SystemExit(f"RCSB Data API failed ({problem}); rerun this step to retry.")
    return out
