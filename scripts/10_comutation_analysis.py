#!/usr/bin/env python3
"""
10_comutation_analysis.py
Project 130 - Colorectal cancer (TCGA-COAD)  --  EXTENSION (PI request)

Motivation: single-mutation neoantigens are often ineffective (tumours escape
by losing one clone/epitope). Multi-epitope combinations are stronger targets.
This script quantifies DOUBLE and TRIPLE co-mutation among established
colorectal-cancer driver genes: how often pairs/triples of drivers occur
together in the same tumour, whether that co-occurrence exceeds chance
(enrichment) or falls below it (mutual exclusivity).

Level: gene-level (a gene is "mutated" in a sample if ANY of its filtered
missense SNVs is present in that sample), restricted to a curated CRC driver
panel to keep results interpretable and avoid gene-size artefacts (e.g. TTN).

Statistics (implemented in pure Python; no scipy dependency):
  - Pairs: Fisher's exact test via the hypergeometric distribution. Reports
    the odds ratio and one-sided p-values for co-occurrence (right tail) and
    mutual exclusivity (left tail). Co-occurrence p-values are FDR-corrected
    (Benjamini-Hochberg).
  - Triples: observed vs expected-under-independence count, an enrichment
    ratio, and a Poisson right-tail p-value (P[X >= observed] with
    lambda = expected).

Inputs:  results/03_integrated_mutation_expression.tsv
         results/neoantigen_candidates_shortlist.tsv (to flag which drivers
         carry candidate neoantigens)
Outputs: results/comutation_driver_pairs.tsv
         results/comutation_driver_triples.tsv
         results/comutation_summary.txt
"""
import itertools
import math
import os
import sys
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
INT = os.path.join(RES, "03_integrated_mutation_expression.tsv")
SHORT = os.path.join(RES, "neoantigen_candidates_shortlist.tsv")
OUT_PAIRS = os.path.join(RES, "comutation_driver_pairs.tsv")
OUT_TRIPLES = os.path.join(RES, "comutation_driver_triples.tsv")
OUT_SUM = os.path.join(RES, "comutation_summary.txt")

# ---- Curated colorectal-cancer driver genes (well established in the -------
# literature / TCGA CRC). Only those actually present in the data are used.
DRIVERS = [
    "APC", "TP53", "KRAS", "PIK3CA", "FBXW7", "SMAD4", "TCF7L2", "NRAS",
    "SMAD2", "CTNNB1", "BRAF", "SOX9", "ARID1A", "AMER1", "FAM123B", "ATM",
    "KMT2C", "KMT2D", "ERBB2", "ERBB3", "PTEN", "ACVR2A", "GNAS", "BMPR1A",
    "TGFBR2", "RNF43", "B2M", "POLE", "MSH6", "CASP8", "ELF3", "PCBP1",
    "AXIN2", "MAP2K4", "CDC27",
]

def log(m): print(f"[10]", m, flush=True)

# --- pure-python exact statistics ------------------------------------------
def hyper_pmf(x, N, K, n):
    # P(X = x) for X ~ Hypergeometric(N population, K successes, n draws)
    if x < max(0, n - (N - K)) or x > min(K, n):
        return 0.0
    return math.comb(K, x) * math.comb(N - K, n - x) / math.comb(N, n)

def fisher_tails(a, nA, nB, N):
    # a = #both; nA = #geneA; nB = #geneB; N = #samples.
    # Right tail = P(X >= a) (co-occurrence); Left tail = P(X <= a) (exclusivity).
    lo, hi = max(0, nA + nB - N), min(nA, nB)
    right = sum(hyper_pmf(x, N, nA, nB) for x in range(a, hi + 1))
    left = sum(hyper_pmf(x, N, nA, nB) for x in range(lo, a + 1))
    return min(1.0, right), min(1.0, left)

def poisson_sf(obs, lam):
    # P(X >= obs) for X ~ Poisson(lam) = 1 - P(X <= obs-1)
    if obs <= 0:
        return 1.0
    cdf = sum(math.exp(-lam) * lam**k / math.factorial(k) for k in range(obs))
    return max(0.0, min(1.0, 1.0 - cdf))

def bh_fdr(pvals):
    # Benjamini-Hochberg FDR correction. Returns q-values aligned to input.
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = m - rank + 1
        val = min(prev, pvals[i] * m / k)
        q[i] = val
        prev = val
    return q

def main():
    # ---- Build gene-level presence for driver genes only ------------------
    with open(INT) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    # sample columns start at the first TCGA-prefixed header (robust to extra
    # metadata columns such as GeneLevelTPM_SD)
    s0 = next(i for i, c in enumerate(header) if c.startswith("TCGA"))
    samples = header[s0:]
    N = len(samples)
    driverset = set(DRIVERS)
    present = {}                       # gene -> boolean numpy array over samples
    with open(INT) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            g = p[0]
            if g not in driverset:
                continue
            vec = np.fromiter((1 if v == "1" else 0 for v in p[s0:]),
                              dtype=np.int8, count=N)
            if g in present:
                present[g] = present[g] | vec        # OR across the gene's mutations
            else:
                present[g] = vec.copy()
    genes = [g for g in DRIVERS if g in present]      # keep curated order, only found
    log(f"Samples: {N}; driver genes found: {len(genes)} of {len(DRIVERS)}")
    counts = {g: int(present[g].sum()) for g in genes}

    # ---- flag drivers that carry a prioritised candidate neoantigen -------
    cand_genes = set()
    if os.path.exists(SHORT):
        with open(SHORT) as fh:
            fh.readline()
            for line in fh:
                cand_genes.add(line.split("\t", 1)[0])

    # ---- PAIRS (doubles): Fisher's exact --------------------------------
    pair_rows = []
    for gA, gB in itertools.combinations(genes, 2):
        A, B = present[gA], present[gB]
        a = int((A & B).sum())
        nA, nB = counts[gA], counts[gB]
        right, left = fisher_tails(a, nA, nB, N)
        # odds ratio from the 2x2 table (Haldane 0.5 correction to avoid /0)
        b = nA - a; c = nB - a; d = N - nA - nB + a
        orr = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        expected = nA * nB / N
        pair_rows.append(dict(GeneA=gA, GeneB=gB, nA=nA, nB=nB, nBoth=a,
                              Expected=round(expected, 2), OddsRatio=round(orr, 3),
                              p_cooccur=right, p_exclusive=left,
                              GeneA_neo=gA in cand_genes, GeneB_neo=gB in cand_genes))
    # FDR on the co-occurrence p-values
    q = bh_fdr([r["p_cooccur"] for r in pair_rows])
    qex = bh_fdr([r["p_exclusive"] for r in pair_rows])
    for r, qi, qe in zip(pair_rows, q, qex):
        r["FDR_cooccur"] = qi
        r["FDR_exclusive"] = qe
        r["Relationship"] = ("Co-occurring" if r["OddsRatio"] > 1 and qi < 0.05
                             else "Mutually exclusive" if r["OddsRatio"] < 1 and qe < 0.05
                             else "n.s.")
    pair_rows.sort(key=lambda r: r["nBoth"], reverse=True)

    cols_p = ["GeneA", "GeneB", "nA", "nB", "nBoth", "Expected", "OddsRatio",
              "p_cooccur", "FDR_cooccur", "p_exclusive", "FDR_exclusive",
              "Relationship", "GeneA_neo", "GeneB_neo"]
    with open(OUT_PAIRS, "w") as fh:
        fh.write("\t".join(cols_p) + "\n")
        for r in pair_rows:
            fh.write("\t".join(str(_fmt(r[c])) for c in cols_p) + "\n")
    log(f"Wrote {OUT_PAIRS}: {len(pair_rows)} driver pairs")

    # ---- TRIPLES: observed vs expected + Poisson enrichment --------------
    trip_rows = []
    for gA, gB, gC in itertools.combinations(genes, 3):
        both = int((present[gA] & present[gB] & present[gC]).sum())
        if both == 0:
            continue
        expected = counts[gA] * counts[gB] * counts[gC] / (N * N)
        enr = both / expected if expected > 0 else float("inf")
        pval = poisson_sf(both, expected)
        trip_rows.append(dict(GeneA=gA, GeneB=gB, GeneC=gC, nAll3=both,
                              Expected=round(expected, 3),
                              Enrichment=round(enr, 2), p_cooccur=pval,
                              AnyNeo=any(g in cand_genes for g in (gA, gB, gC))))
    if trip_rows:
        qt = bh_fdr([r["p_cooccur"] for r in trip_rows])
        for r, qi in zip(trip_rows, qt):
            r["FDR"] = qi
    trip_rows.sort(key=lambda r: r["nAll3"], reverse=True)
    cols_t = ["GeneA", "GeneB", "GeneC", "nAll3", "Expected", "Enrichment",
              "p_cooccur", "FDR", "AnyNeo"]
    with open(OUT_TRIPLES, "w") as fh:
        fh.write("\t".join(cols_t) + "\n")
        for r in trip_rows:
            fh.write("\t".join(str(_fmt(r[c])) for c in cols_t) + "\n")
    log(f"Wrote {OUT_TRIPLES}: {len(trip_rows)} co-occurring driver triples")

    # ---- Per-tumour driver-mutation load ---------------------------------
    load = np.zeros(N, dtype=int)
    for g in genes:
        load += present[g]
    dist = {k: int((load == k).sum()) for k in range(0, load.max() + 1)}
    n_ge2 = int((load >= 2).sum()); n_ge3 = int((load >= 3).sum())

    with open(OUT_SUM, "w") as fh:
        fh.write("DRIVER-GENE CO-MUTATION SUMMARY (colorectal, TCGA-COAD)\n")
        fh.write(f"Tumour samples: {N}\n")
        fh.write(f"Driver genes analysed ({len(genes)}): {', '.join(genes)}\n\n")
        fh.write("Per-gene mutation frequency (samples):\n")
        for g in sorted(genes, key=lambda x: -counts[x]):
            tag = " [has candidate neoantigen]" if g in cand_genes else ""
            fh.write(f"  {g:8s} {counts[g]:4d} ({100*counts[g]/N:4.1f}%){tag}\n")
        fh.write("\nPer-tumour number of co-mutated driver genes:\n")
        for k, v in dist.items():
            fh.write(f"  {k} drivers: {v} tumours\n")
        fh.write(f"\nTumours with >=2 co-mutated drivers (doubles+): {n_ge2} "
                 f"({100*n_ge2/N:.1f}%)\n")
        fh.write(f"Tumours with >=3 co-mutated drivers (triples+): {n_ge3} "
                 f"({100*n_ge3/N:.1f}%)\n\n")
        fh.write("Top 12 co-occurring driver pairs (by tumours carrying both):\n")
        for r in pair_rows[:12]:
            fh.write(f"  {r['GeneA']:7s}+{r['GeneB']:7s}  both={r['nBoth']:3d}  "
                     f"OR={r['OddsRatio']:.2f}  FDR_co={r['FDR_cooccur']:.2e}  "
                     f"{r['Relationship']}\n")
        fh.write("\nSignificant mutually exclusive driver pairs (FDR<0.05):\n")
        excl = [r for r in pair_rows if r["Relationship"] == "Mutually exclusive"]
        excl.sort(key=lambda r: r["FDR_exclusive"])
        for r in excl[:12]:
            fh.write(f"  {r['GeneA']:7s}/{r['GeneB']:7s}  both={r['nBoth']:3d}  "
                     f"OR={r['OddsRatio']:.2f}  FDR_excl={r['FDR_exclusive']:.2e}\n")
        fh.write("\nTop 10 co-occurring driver triples (by tumours carrying all three):\n")
        for r in trip_rows[:10]:
            fh.write(f"  {r['GeneA']:7s}+{r['GeneB']:7s}+{r['GeneC']:7s}  "
                     f"all3={r['nAll3']:3d}  enrich={r['Enrichment']:.1f}x  "
                     f"FDR={r.get('FDR',float('nan')):.2e}\n")
    log(f"Wrote {OUT_SUM}")
    log(f">=2 drivers: {n_ge2} tumours ({100*n_ge2/N:.1f}%); "
        f">=3 drivers: {n_ge3} tumours ({100*n_ge3/N:.1f}%)")

def _fmt(v):
    if isinstance(v, float):
        if v != v:  # NaN
            return "NA"
        if v == 0 or (abs(v) >= 1e-4 and abs(v) < 1e6):
            return f"{v:.4g}"
        return f"{v:.3e}"
    return v

if __name__ == "__main__":
    sys.exit(main())
