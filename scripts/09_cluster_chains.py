#!/usr/bin/env python3
"""Step 09 - Cluster all chain sequences with MMseqs2 at several identity thresholds.

Input : paths.checked (step 08)
Output: paths.clustered        input rows + molA/molB/pair cluster IDs per threshold
        paths.clustering_dir/  unique_sequences.fasta, cXX/clusters.tsv, summary.txt/json

Unique sequences are identified by sha1(sequence)[:16]. For each threshold XX:
  molA_cluster_cXX / molB_cluster_cXX  representative seq_uid of the chain's cluster
  pair_cluster_cXX                     order-invariant "min||max" of both chain clusters
"""

from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

from common import ensure_parent, get_logger, parse_args, seq_uid

log = get_logger("09_cluster_chains")


def thr_tag(thr: float) -> str:
    return f"c{int(round(thr * 100)):02d}"


def run(cmd: List[str]) -> None:
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\n{p.stdout}\n{p.stderr}")


def mmseqs_cluster(fasta: Path, out_dir: Path, thr: float, ccfg: Dict, threads: int) -> Path:
    """Cluster ``fasta`` at identity ``thr``; returns clusters.tsv (rep<TAB>member)."""
    if out_dir.exists():
        shutil.rmtree(out_dir)
    tmp = out_dir / "tmp"
    tmp.mkdir(parents=True)
    db, clu, tsv = out_dir / "seqdb", out_dir / "clu", out_dir / "clusters.tsv"
    run(["mmseqs", "createdb", str(fasta), str(db)])
    run(["mmseqs", ccfg["engine"], str(db), str(clu), str(tmp),
         "--min-seq-id", str(thr), "-c", str(ccfg["coverage"]),
         "--cov-mode", str(ccfg["cov_mode"]), "--threads", str(threads)])
    run(["mmseqs", "createtsv", str(db), str(db), str(clu), str(tsv)])
    # keep only the TSV: the MMseqs2 databases and tmp/ are large and disposable
    for p in out_dir.iterdir():
        if p != tsv:
            shutil.rmtree(p) if p.is_dir() else p.unlink()
    return tsv


@dataclass
class ThrSummary:
    thr: float
    n_clusters: int
    singleton_clusters: int
    singleton_clusters_pct: float
    max_unique: int
    p50: int
    p90: int
    p99: int
    top10_share_pct: float
    max_instances: int


def summarize(clusters: pd.DataFrame, instances: Dict[str, int], thr: float) -> ThrSummary:
    sizes = np.sort(clusters.groupby("rep")["member"].size().to_numpy())
    inst = clusters.assign(n=clusters["member"].map(instances)).groupby("rep")["n"].sum()
    n_unique = len(clusters)
    return ThrSummary(
        thr=thr,
        n_clusters=int(len(sizes)),
        singleton_clusters=int((sizes == 1).sum()),
        singleton_clusters_pct=100.0 * float((sizes == 1).sum()) / len(sizes),
        max_unique=int(sizes[-1]),
        p50=int(np.percentile(sizes, 50, method="nearest")),
        p90=int(np.percentile(sizes, 90, method="nearest")),
        p99=int(np.percentile(sizes, 99, method="nearest")),
        top10_share_pct=100.0 * float(sizes[::-1][:10].sum()) / n_unique,
        max_instances=int(inst.max()),
    )


def main() -> None:
    cfg = parse_args(__doc__)
    paths, ccfg = cfg["paths"], cfg["clustering"]
    out_dir = Path(paths["clustering_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(paths["checked"], keep_default_na=False, na_values=[""],
                     dtype={"pdb_id": str, "molA_chain_id": str, "molB_chain_id": str})

    uid_a = df["molA_seq"].map(seq_uid)
    uid_b = df["molB_seq"].map(seq_uid)
    uid_to_seq = dict(zip(pd.concat([uid_a, uid_b]), pd.concat([df["molA_seq"], df["molB_seq"]])))
    instances = pd.concat([uid_a, uid_b]).value_counts().to_dict()
    log.info("%d rows, %d unique sequences, %d chain instances", len(df), len(uid_to_seq), len(uid_a) * 2)

    fasta = out_dir / "unique_sequences.fasta"
    with fasta.open("w", encoding="utf-8") as fh:
        for uid, seq in uid_to_seq.items():
            fh.write(f">{uid}\n{seq}\n")

    summaries = []
    for thr in ccfg["thresholds"]:
        tag = thr_tag(float(thr))
        log.info("MMseqs2 %s at min-seq-id %.2f", ccfg["engine"], thr)
        tsv = mmseqs_cluster(fasta, out_dir / tag, float(thr), ccfg, int(cfg["threads"]))
        clusters = pd.read_csv(tsv, sep="\t", header=None, names=["rep", "member"], dtype=str)
        summaries.append(summarize(clusters, instances, float(thr)))

        member_to_rep = dict(zip(clusters["member"], clusters["rep"]))
        a = uid_a.map(member_to_rep)
        b = uid_b.map(member_to_rep)
        df[f"molA_cluster_{tag}"] = a
        df[f"molB_cluster_{tag}"] = b
        lo, hi = a.where(a <= b, b), b.where(a <= b, a)
        df[f"pair_cluster_{tag}"] = (lo + "||" + hi).where(a.notna() & b.notna())

    dest = ensure_parent(paths["clustered"])
    tmp = dest.with_name(dest.name + ".part")
    df.to_csv(tmp, index=False)
    tmp.replace(dest)

    generated = datetime.now().isoformat(timespec="seconds")
    (out_dir / "summary.json").write_text(json.dumps({
        "generated": generated, "engine": ccfg["engine"], "coverage": ccfg["coverage"],
        "cov_mode": ccfg["cov_mode"], "thresholds": ccfg["thresholds"],
        "unique_sequences": len(uid_to_seq), "chain_instances": int(len(uid_a) * 2),
        "per_threshold": [asdict(s) for s in summaries],
    }, indent=2), encoding="utf-8")

    header = f"{'thr':>5} {'clusters':>9} {'singletons':>11} {'max':>6} {'p50':>4} {'p90':>4} {'p99':>4} {'top10%':>7} {'max_inst':>9}"
    lines = [f"MMseqs2 chain clustering ({generated})",
             f"engine={ccfg['engine']}  coverage={ccfg['coverage']}  cov_mode={ccfg['cov_mode']}",
             f"unique sequences={len(uid_to_seq)}  chain instances={len(uid_a) * 2}", "", header]
    for s in summaries:
        lines.append(f"{s.thr:5.2f} {s.n_clusters:9d} {s.singleton_clusters:6d} ({s.singleton_clusters_pct:4.1f}%)"
                     f" {s.max_unique:5d} {s.p50:4d} {s.p90:4d} {s.p99:4d} {s.top10_share_pct:6.2f}% {s.max_instances:9d}")
    (out_dir / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    log.info("Wrote %s and %s/summary.txt", dest, out_dir)


if __name__ == "__main__":
    main()
