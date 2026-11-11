#!/usr/bin/env python3
"""Step 09 - Cluster all chain sequences with MMseqs2 at several identity thresholds.

Input : paths.prodigy_filtered (step 08b) if prodigy.filter is true, else paths.checked (step 08)
Output: paths.clustered        input rows + molA/molB/pair cluster IDs per threshold
        paths.clustering_dir/  unique_sequences.fasta, cXX/clusters.tsv, summary.txt/json

clustering.method:
  connected_components  (default) one all-vs-all MMseqs2 search; at each
        threshold, two sequences are linked when identity >= threshold and the
        coverage criterion holds, and clusters are the connected components
        (single linkage). No two sequences of different clusters are similar
        above the threshold, so cluster-based splits do not leak homologues.
  set_cover  MMseqs2 cluster/linclust as in PPI_PDB_Scraper Seq8: every member is
        similar to its representative, but members of different clusters can be
        similar above the threshold.

Unique sequences are identified by sha1(sequence)[:16]. For each threshold XX:
  molA_cluster_cXX / molB_cluster_cXX  representative seq_uid of the chain's cluster
                                       (smallest member ID for connected components)
  pair_cluster_cXX                     order-invariant "min||max" of both chain clusters
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from common import ensure_parent, get_logger, homology_linked, parse_args, seq_uid

log = get_logger("09_cluster_chains")

def required_levels(cfg: Dict) -> set:
    """Identity levels that later steps read, from the configuration."""
    req = {round(float(cfg["test_set"]["level"]), 2)} | {round(float(x), 2) for x in cfg["splits"]["schemes"]}
    hn = cfg["hard_negatives"]
    if hn.get("enabled"):
        req |= {int(str(hn[k])[1:]) / 100 for k in ("conflict_level", "family_level", "pair_cluster_key")
                if re.fullmatch(r"c\d+", str(hn[k]))}
    if cfg["pruning"].get("cluster_pruning"):
        req.add(int(re.search(r"(\d+)$", cfg["pruning"]["pair_key"]).group(1)) / 100)
    return req


def thr_tag(thr: float) -> str:
    return f"c{int(round(thr * 100)):02d}"


def run(cmd: List[str]) -> None:
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\n{p.stdout}\n{p.stderr}")


def mmseqs_version() -> str:
    p = subprocess.run(["mmseqs", "version"], capture_output=True, text=True)
    return p.stdout.strip()


def all_vs_all_hits(fasta: Path, work_dir: Path, ccfg: Dict, threads: int) -> pd.DataFrame:
    """MMseqs2 search of every sequence against all, down to the lowest threshold.
    Coverage is applied afterwards with the same rule as the leakage audit (14)."""
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True)
    hits = work_dir / "hits.m8"
    run(["mmseqs", "easy-search", str(fasta), str(fasta), str(hits), str(work_dir / "tmp"),
         "--min-seq-id", str(min(ccfg["thresholds"])), "-c", "0", "-e", str(ccfg["evalue"]),
         "-s", str(ccfg["search_sensitivity"]), "--max-seqs", str(ccfg["max_seqs"]), "--threads", str(threads),
         "--format-output", "query,target,fident,qcov,tcov"])
    df = pd.read_csv(hits, sep="\t", header=None, names=["query", "target", "fident", "qcov", "tcov"],
                     dtype={"query": str, "target": str})
    shutil.rmtree(work_dir)
    return df.loc[homology_linked(df, ccfg), ["query", "target", "fident"]]


def components_tsv(uids: List[str], hits: pd.DataFrame, thr: float, out_dir: Path) -> Path:
    """Connected components of the hits with identity >= thr; returns clusters.tsv
    (rep<TAB>member), the representative being the smallest member ID."""
    index = {u: i for i, u in enumerate(uids)}
    edges = hits[hits["fident"] >= thr]
    graph = coo_matrix((np.ones(len(edges)), (edges["query"].map(index), edges["target"].map(index))),
                       shape=(len(uids), len(uids)))
    _, labels = connected_components(graph, directed=False)
    members = pd.DataFrame({"member": uids, "label": labels})
    members["rep"] = members.groupby("label")["member"].transform("min")
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)
    tsv = out_dir / "clusters.tsv"
    members[["rep", "member"]].sort_values(["rep", "member"]).to_csv(tsv, sep="\t", header=False, index=False)
    return tsv


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

    src = paths["prodigy_filtered"] if cfg.get("prodigy", {}).get("filter") else paths["checked"]
    df = pd.read_csv(src, keep_default_na=False, na_values=[""],
                     dtype={"pdb_id": str, "molA_chain_id": str, "molB_chain_id": str, "revision": str})

    uid_a = df["molA_seq"].map(seq_uid)
    uid_b = df["molB_seq"].map(seq_uid)
    uid_to_seq = dict(zip(pd.concat([uid_a, uid_b]), pd.concat([df["molA_seq"], df["molB_seq"]])))
    instances = pd.concat([uid_a, uid_b]).value_counts().to_dict()
    log.info("%d rows, %d unique sequences, %d chain instances", len(df), len(uid_to_seq), len(uid_a) * 2)

    fasta = out_dir / "unique_sequences.fasta"
    with fasta.open("w", encoding="utf-8") as fh:
        for uid, seq in uid_to_seq.items():
            fh.write(f">{uid}\n{seq}\n")

    method = ccfg["method"]
    required = required_levels(cfg)
    missing = sorted(required - {round(float(t), 2) for t in ccfg["thresholds"]})
    if missing:
        raise SystemExit(f"clustering.thresholds must include {sorted(required)} (read by steps 10-14 with "
                         f"this configuration); missing: {missing}")
    if method == "connected_components":
        log.info("MMseqs2 all-vs-all search (min-seq-id %.2f, -s %s)", min(ccfg["thresholds"]), ccfg["search_sensitivity"])
        hits = all_vs_all_hits(fasta, out_dir / "search", ccfg, int(cfg["threads"]))
    elif method != "set_cover":
        raise SystemExit(f"Unknown clustering.method: {method}")

    summaries = []
    for thr in ccfg["thresholds"]:
        tag = thr_tag(float(thr))
        if method == "connected_components":
            tsv = components_tsv(list(uid_to_seq), hits, float(thr), out_dir / tag)
        else:
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
        "generated": generated, "mmseqs_version": mmseqs_version(), "method": method, "engine": ccfg["engine"],
        "evalue": ccfg["evalue"], "search_sensitivity": ccfg["search_sensitivity"], "max_seqs": ccfg["max_seqs"],
        "coverage": ccfg["coverage"],
        "cov_mode": ccfg["cov_mode"], "thresholds": ccfg["thresholds"],
        "unique_sequences": len(uid_to_seq), "chain_instances": int(len(uid_a) * 2),
        "per_threshold": [asdict(s) for s in summaries],
    }, indent=2), encoding="utf-8")

    header = f"{'thr':>5} {'clusters':>9} {'singletons':>11} {'max':>6} {'p50':>4} {'p90':>4} {'p99':>4} {'top10%':>7} {'max_inst':>9}"
    lines = [f"MMseqs2 chain clustering ({generated})",
             f"method={method}  coverage={ccfg['coverage']}  cov_mode={ccfg['cov_mode']}",
             f"unique sequences={len(uid_to_seq)}  chain instances={len(uid_a) * 2}", "", header]
    for s in summaries:
        lines.append(f"{s.thr:5.2f} {s.n_clusters:9d} {s.singleton_clusters:6d} ({s.singleton_clusters_pct:4.1f}%)"
                     f" {s.max_unique:5d} {s.p50:4d} {s.p90:4d} {s.p99:4d} {s.top10_share_pct:6.2f}% {s.max_instances:9d}")
    (out_dir / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    log.info("Wrote %s and %s/summary.txt", dest, out_dir)


if __name__ == "__main__":
    main()
