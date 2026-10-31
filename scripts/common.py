"""Shared helpers for the BeastAffinity pipeline scripts.

Every numbered script imports this module (it lives next to them in
``scripts/``), so configuration loading, PDB ID handling, HTTP downloads and
atomic file writes behave the same way in every step.
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import requests
import yaml
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

DEFAULT_CONFIG = "config/config.yaml"

# Classic 4-character PDB IDs (first character is a digit).
# Extended IDs (pdb_XXXXXXXX) will be handled in v2.
PDB_ID_RE = re.compile(r"^[0-9][a-z0-9]{3}$")


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


def _load_raw(path: str | Path) -> Dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh) or {}
    base = cfg.pop("base", None)
    return _deep_merge(_load_raw(base), cfg) if base else cfg


def load_config(path: str | Path = DEFAULT_CONFIG) -> Dict[str, Any]:
    """Load a YAML config. A ``base:`` key points to a parent config whose
    values are inherited and can be overridden (used by config/test.yaml).
    ``{data_dir}`` in paths is replaced by the configured data_dir."""
    cfg = _load_raw(path)
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
    s = token.strip().lower()
    return s if PDB_ID_RE.match(s) else None


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def read_id_lists(input_dir: str | Path) -> Tuple[List[str], List[Tuple[str, str]]]:
    """Read every .txt / .txt.gz file in ``input_dir`` (RCSB Advanced Search
    exports: IDs separated by commas, whitespace or one per line).

    Returns (unique IDs in first-seen order, rejected (file, token) pairs).
    """
    input_dir = Path(input_dir)
    files = sorted(p for p in input_dir.iterdir() if p.name.endswith((".txt", ".txt.gz")))
    if not files:
        raise SystemExit(f"No .txt or .txt.gz ID lists found in {input_dir}/")

    seen: Dict[str, None] = {}
    rejected: List[Tuple[str, str]] = []
    for path in files:
        with _open_text(path) as fh:
            for line in fh:
                for token in line.replace(",", " ").split():
                    pid = normalize_pdb_id(token)
                    if pid is None:
                        rejected.append((path.name, token))
                    else:
                        seen.setdefault(pid, None)
    return list(seen), rejected


def read_ids_file(path: str | Path) -> List[str]:
    with Path(path).open("r", encoding="utf-8") as fh:
        return [ln.strip() for ln in fh if ln.strip()]


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


def gzip_is_valid(path: Path) -> bool:
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1 << 20):
                pass
        return True
    except (OSError, EOFError):
        return False


def read_json(path: str | Path) -> Any:
    path = Path(path)
    data = gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes()
    return json.loads(data)


def json_is_valid(path: Path) -> bool:
    try:
        read_json(path)
        return True
    except (OSError, EOFError, ValueError):
        return False


# -----------------------------------------------------------------------------
# HTTP
# -----------------------------------------------------------------------------
def make_session(pool_size: int, user_agent: str = "BeastAffinity-pipeline/1.0") -> requests.Session:
    retries = Retry(
        total=6,
        backoff_factor=1.0,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods={"GET", "HEAD"},
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


def download_file(session: requests.Session, url: str, dest: Path, validator, timeout=(10, 120),
                  gzip_output: bool = False) -> str:
    """Download ``url`` to ``dest`` atomically (gzip-compressing it if
    ``gzip_output``).

    Returns "skip" (dest already valid), "ok", "missing" (HTTP 404) or an error
    string starting with "error:". ``validator(path) -> bool`` checks content.
    """
    if dest.is_file() and validator(dest):
        return "skip"
    ensure_parent(dest)
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
    except requests.RequestException as exc:
        tmp.unlink(missing_ok=True)
        return f"error:{type(exc).__name__}: {exc}"
