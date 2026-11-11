#!/usr/bin/env python3
"""Independent leakage control (optional): global-alignment identities instead of MMseqs2.

Step 14 audits leakage with MMseqs2 and the linking rule of step 09 (identity,
coverage, E-value), so it checks that the splits are consistent with the
clustering; it cannot see what that rule does not link. This script audits the
same sets against the same reference data (step 14's evaluation_sets) with
measures that the pipeline does not use to build the splits:

  nw_full          identical positions / length of the global alignment, gaps
                   included (GraphPart's default identity, Teufel et al. 2023)
  af3              identical positions / length of the evaluation chain
                   (AlphaFold 3's sequence identity for its low-homology sets,
                   Abramson et al. 2024, Supplementary Information 6.1)
  af3_significant  af3, counted only when the two chains also have a
                   significant local alignment (MMseqs2, E <= 1e-3, any
                   coverage); without it, short chains reach 40 % by chance
                   inside much longer ones

Alignments: EMBOSS needleall with GraphPart's defaults (EBLOSUM62, gap open 10,
extend 0.5, end gaps not penalised; pairs scoring below its -minscore 1.0 count
as 0 identical positions). To make the all-against-all comparison feasible,
parasail first aligns every pair with the same scoring (semi-global,
BLOSUM62 x 2, open 20, extend 1) and only the pairs that can reach identity
PRE (0.35) by either measure are realigned with needleall. On the test
configuration, aligning every pair with needleall showed that this prefilter
misses no pair at or above 0.40.

Leaks are counted with the rules of step 14 (chain or pair level, at each set's
design threshold, identity >= threshold; identical sequences have identity 1).
The script only reports; it does not change any dataset.

Runs in its own environment (about 1.5 h on 20 cores for 164,000 PDB IDs):
    conda env create -f envs/leakage_control.yml
    snakemake --cores 20 leakage_control

Input : paths.test, paths.trainval_splits, paths.hard_negatives (as step 14)
Output: reports_dir/leakage_control/summary.json     leaks per set at its design threshold
        reports_dir/leakage_control/by_threshold.tsv every set at every threshold
        reports_dir/leakage_control/leaks.tsv        each leaking chain or pair with its reference
        reports_dir/leakage_control/hits.tsv         every realigned pair (identical positions m, length L)
"""
from __future__ import annotations

import importlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
import parasail

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from common import get_logger, parse_args  # noqa: E402

evaluation_sets = importlib.import_module("14_leakage_audit").evaluation_sets
log = get_logger("leakage_control")

PRE = 0.35                   # prefilter margin below the lowest design threshold (0.40)
THRESHOLDS = [0.95, 0.90, 0.60, 0.40]
OPEN2, EXT2 = 20, 1          # EMBOSS 10 / 0.5 with the matrix doubled (parasail needs integers)
NEEDLE = ["-gapopen", "10", "-gapextend", "0.5", "-endopen", "10", "-endextend", "0.5", "-datafile", "EBLOSUM62"]
MEASURES = ("nw_full", "af3", "af3_significant")


def matrix_x2():
    m = parasail.blosum62.copy()
    for i in range(m.size):
        for j in range(m.size):
            m.set_value(i, j, 2 * int(parasail.blosum62.matrix[i][j]))
    return m


# --- stage 1: parasail over every (evaluation, reference) pair ------------------
_REF, _MAT = None, None


def _init(ref):
    global _REF, _MAT
    _REF, _MAT = ref, matrix_x2()


def prefilter(item):
    q, sq = item
    prof = parasail.profile_create_stats_sat(sq, _MAT)
    lq, out = len(sq), []
    for r, sr in _REF:
        if r == q:
            out.append(r)
            continue
        lr = len(sr)
        if min(lq, lr) / lq < PRE and min(lq, lr) / max(lq, lr) < PRE:
            continue                       # identical positions <= min(lq, lr): neither measure can reach PRE
        m = parasail.sg_stats_striped_profile_sat(prof, sr, OPEN2, EXT2).matches
        if m / lq >= PRE or m / max(lq, lr) >= PRE:
            out.append(r)
    return q, out


# --- stage 2: EMBOSS needleall on the candidates --------------------------------
def needle(item):
    q, sq, cands = item                    # cands: (uid, sequence), uid != q
    if not cands:
        return []
    with tempfile.TemporaryDirectory() as d:
        a, b = Path(d) / "a.fa", Path(d) / "b.fa"
        a.write_text(f">{q}\n{sq}\n")
        b.write_text("".join(f">{u}\n{s}\n" for u, s in cands))
        p = subprocess.run(["needleall", "-asequence", str(a), "-bsequence", str(b), *NEEDLE, "-aformat3", "pair",
                            "-outfile", "stdout", "-errfile", str(Path(d) / "err"), "-auto"],
                           capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(p.stderr)
    res, target = [], None
    for line in p.stdout.splitlines():
        if line.startswith("# 2:"):
            target = line.split(":", 1)[1].strip()
        elif line.startswith("# Length:"):
            length = int(line.split(":")[1])
        elif line.startswith("# Identity:"):
            res.append((q, target, int(line.split(":")[1].split("/")[0]), length))
    seen = {t for _, t, _, _ in res}
    res += [(q, u, 0, max(len(sq), len(s))) for u, s in cands if u not in seen]   # below -minscore
    return res


def significant_pairs(hits, seq_of, cfg, work: Path, threads: int):
    """(query, target) pairs with an MMseqs2 alignment at E <= 1e-3, any coverage."""
    c = hits[(hits["af3"] >= min(THRESHOLDS)) & (hits["query"] != hits["target"])]
    if c.empty:
        return set()
    work.mkdir(parents=True, exist_ok=True)
    (work / "q.fa").write_text("".join(f">{u}\n{seq_of[u]}\n" for u in c["query"].unique()))
    (work / "t.fa").write_text("".join(f">{u}\n{seq_of[u]}\n" for u in c["target"].unique()))
    subprocess.run(["mmseqs", "easy-search", str(work / "q.fa"), str(work / "t.fa"), str(work / "h.m8"), str(work / "tmp"),
                    "--min-seq-id", "0", "-c", "0", "-e", "1e-3", "-s", str(cfg["clustering"]["search_sensitivity"]),
                    "--max-seqs", "100000", "--threads", str(threads), "--format-output", "query,target"],
                   check=True, capture_output=True)
    m = pd.read_csv(work / "h.m8", sep="\t", header=None, names=["query", "target"], dtype=str)
    return set(zip(m["query"], m["target"])) | set(zip(m["target"], m["query"]))


def count_leaks(sets, hits, measure):
    """Step-14 rules on a table of identities; returns per-threshold rows, the
    design-threshold summary and the leaking items at the design threshold."""
    rows, summary, items = [], {}, []
    for name, pairs, ref, level, design_thr, only in sets:
        ref_set = {u for p in ref for u in p}
        pairs_of = defaultdict(set)
        for i, (a, b) in enumerate(ref):
            pairs_of[a].add(i)
            pairs_of[b].add(i)
        chains = {u for p in pairs for u in p} if only is None else set(only)
        sub = hits[hits["query"].isin(chains) & hits["target"].isin(ref_set)]
        for thr in THRESHOLDS:
            ok = sub[sub[measure] >= thr]
            hom = defaultdict(set)
            for q, t in zip(ok["query"], ok["target"]):
                hom[q].add(t)

            def ref_pair(a, b):
                ha, hb = hom.get(a, set()), hom.get(b, set())
                if not ha or not hb:
                    return None
                for i in set().union(*(pairs_of[x] for x in ha)) & set().union(*(pairs_of[y] for y in hb)):
                    x, y = ref[i]
                    if (x in ha and y in hb) or (y in ha and x in hb):
                        return x, y
                return None

            leak_chain = {u for u in chains if hom.get(u)}
            leak_pairs = {(a, b): rp for a, b in pairs if (rp := ref_pair(a, b))}
            rec = {"set": name, "measure": measure, "threshold": thr, "rows": len(pairs), "chains": len(chains),
                   "chains_with_homologue": len(leak_chain),
                   "rows_any_chain": int(sum(a in leak_chain or b in leak_chain for a, b in pairs)),
                   "rows_pair_homologue": int(sum(1 for a, b in pairs if (a, b) in leak_pairs))}
            rows.append(rec)
            if design_thr is not None and abs(thr - design_thr) < 1e-9:
                key = "chains_with_homologue" if level == "chain" else "rows_pair_homologue"
                summary[name] = {"design_level": level, "design_threshold": thr, "leaks": rec[key],
                                 "of": rec["chains"] if level == "chain" else rec["rows"]}
                if level == "chain":
                    best = ok.sort_values(measure, ascending=False).drop_duplicates("query").set_index("query")
                    items += [{"measure": measure, "set": name, "level": level, "chains": u,
                               "reference": best.at[u, "target"], "identity": round(best.at[u, measure], 3)}
                              for u in sorted(leak_chain)]
                else:
                    items += [{"measure": measure, "set": name, "level": level, "chains": f"{a}|{b}",
                               "reference": f"{x}|{y}", "identity": None}
                              for (a, b), (x, y) in sorted(leak_pairs.items())]
    return rows, summary, items


def main() -> None:
    cfg = parse_args(__doc__)
    os.chdir(REPO)
    workers = int(os.environ.get("LEAKAGE_CONTROL_WORKERS", os.cpu_count() or 1))
    out = Path(cfg["paths"]["reports_dir"]) / "leakage_control"
    out.mkdir(parents=True, exist_ok=True)
    sets, seq_of = evaluation_sets(cfg)
    Q = sorted({u for _, pairs, _, _, _, _ in sets for p in pairs for u in p})
    R = sorted({u for _, _, ref, _, _, _ in sets for p in ref for u in p})
    log.info("%d evaluation chains against %d reference chains, %d workers", len(Q), len(R), workers)

    t0, cand = time.time(), {}
    with Pool(workers, initializer=_init, initargs=([(u, seq_of[u]) for u in R],)) as pool:
        for i, (q, o) in enumerate(pool.imap_unordered(prefilter, [(u, seq_of[u]) for u in Q], chunksize=4), 1):
            cand[q] = o
            if i % 500 == 0:
                log.info("prefilter: %d of %d evaluation chains (%.0f s)", i, len(Q), time.time() - t0)
    log.info("prefilter: %d candidate pairs in %.0f s", sum(map(len, cand.values())), time.time() - t0)

    rows = []
    with Pool(workers) as pool:
        for res in pool.imap_unordered(needle, [(q, seq_of[q], [(r, seq_of[r]) for r in o if r != q])
                                                for q, o in cand.items()], chunksize=4):
            rows.extend(res)
    rows += [(q, q, len(seq_of[q]), len(seq_of[q])) for q, o in cand.items() if q in o]
    hits = pd.DataFrame(rows, columns=["query", "target", "m", "L"])
    hits["len_query"] = hits["query"].map(lambda u: len(seq_of[u]))
    hits["len_target"] = hits["target"].map(lambda u: len(seq_of[u]))
    hits["nw_full"] = hits["m"] / hits["L"]
    hits["af3"] = hits["m"] / hits["len_query"]
    sig = significant_pairs(hits, seq_of, cfg, out / "mmseqs", min(workers, 32))
    hits["af3_significant"] = [a if (q == t or (q, t) in sig) else 0.0
                               for q, t, a in zip(hits["query"], hits["target"], hits["af3"])]
    shutil.rmtree(out / "mmseqs", ignore_errors=True)
    hits.to_csv(out / "hits.tsv", sep="\t", index=False)

    all_rows, summaries, all_items = [], {}, []
    for measure in MEASURES:
        r, s, it = count_leaks(sets, hits, measure)
        all_rows += r
        summaries[measure] = s
        all_items += it
    pd.DataFrame(all_rows).to_csv(out / "by_threshold.tsv", sep="\t", index=False)
    pd.DataFrame(all_items, columns=["measure", "set", "level", "chains", "reference", "identity"]).to_csv(
        out / "leaks.tsv", sep="\t", index=False)
    (out / "summary.json").write_text(json.dumps({
        "parameters": {"prefilter": "parasail sg_stats_striped_profile_sat, BLOSUM62 x 2, open 20, extend 1",
                       "prefilter_margin": PRE, "alignment": "EMBOSS needleall " + " ".join(NEEDLE),
                       "significance": "MMseqs2 easy-search, E <= 1e-3, any coverage", "thresholds": THRESHOLDS},
        "n_evaluation_chains": len(Q), "n_reference_chains": len(R),
        "design_threshold_leaks": summaries}, indent=2))
    for measure, s in summaries.items():
        for name, v in s.items():
            log.info("%-16s %-32s %s level at %.2f: %d of %d", measure, name, v["design_level"],
                     v["design_threshold"], v["leaks"], v["of"])


if __name__ == "__main__":
    main()
