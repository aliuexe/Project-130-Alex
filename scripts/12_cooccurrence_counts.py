#!/usr/bin/env python3
r"""
12_cooccurrence_counts.py
Project 130 - Colorectal Cancer (TCGA-COAD) Neoantigen Discovery Pipeline

===============================================================================
BIOLOGICAL & COMPUTATIONAL PURPOSE
===============================================================================
This script provides exact integer patient co-occurrence count matrices and
visualizations for driver gene pairs and triplets without statistical models.
Every co-occurrence figure is reported explicitly out of the total cohort size
$N = 586$ TCGA-COAD tumours ($X$ of $N$, $Y\%$).

===============================================================================
COMPUTATIONAL METHODOLOGY
===============================================================================
- Bitwise Matrix Intersections: Evaluates patient set overlaps using boolean vector
  intersections (`present[gA] & present[gB]`).
- Matrix Construction: Outputs a symmetric $34 \times 34$ driver co-occurrence matrix
  where diagonal values represent individual gene mutation counts.

===============================================================================
INPUT & OUTPUT CONTRACTS
===============================================================================
Input:
  - `results/03_integrated_mutation_expression.tsv`

Outputs:
  - `results/cooccurrence_counts_matrix.tsv` (Symmetric count matrix)
  - `results/cooccurrence_pairs_simple.tsv` (Ranked pairs with sample counts)
  - `results/cooccurrence_triples_simple.tsv` (Ranked triplets with sample counts)
  - `figures/fig17_cooccurrence_counts_heatmap.png`
  - `figures/fig18_top_cooccurring_pairs.png`
"""

import itertools, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =============================================================================
# FILE PATHS & RESOURCE RESOLUTION
# =============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results"); FIG = os.path.join(BASE, "figures")
INT = os.path.join(RES, "03_integrated_mutation_expression.tsv")
os.makedirs(FIG, exist_ok=True)

# Curated CRC Driver Panel
DRIVERS = ["APC","TP53","KRAS","PIK3CA","FBXW7","SMAD4","TCF7L2","NRAS","SMAD2",
    "CTNNB1","BRAF","SOX9","ARID1A","AMER1","ATM","KMT2C","KMT2D","ERBB2","ERBB3",
    "PTEN","ACVR2A","GNAS","BMPR1A","TGFBR2","RNF43","B2M","POLE","MSH6","CASP8",
    "ELF3","PCBP1","AXIN2","MAP2K4","CDC27"]

def log(m):
    """Prints timestamped progress messages to stdout with line flushing."""
    print("[12]", m, flush=True)

# =============================================================================
# MAIN PIPELINE EXECUTION
# =============================================================================
def main():
    # =========================================================================
    # STEP 1: CONSTRUCT GENE-LEVEL PRESENCE VECTORS
    # =========================================================================
    with open(INT) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    s0 = next(i for i, c in enumerate(header) if c.startswith("TCGA"))
    N = len(header) - s0                                  # Total number of tumour samples
    dset = set(DRIVERS); present = {}
    with open(INT) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t"); g = p[0]
            if g not in dset: continue
            vec = np.fromiter((1 if v == "1" else 0 for v in p[s0:]),
                              dtype=np.int8, count=N)
            present[g] = present[g] | vec if g in present else vec.copy()
    genes = [g for g in DRIVERS if g in present]
    counts = {g: int(present[g].sum()) for g in genes}
    log(f"TOTAL tumour samples (N) = {N}; driver genes = {len(genes)}")

    # =========================================================================
    # STEP 2: PAIRWISE CO-OCCURRENCE COUNTS
    # =========================================================================
    pairs = []
    for gA, gB in itertools.combinations(genes, 2):
        both = int((present[gA] & present[gB]).sum())
        pairs.append((gA, gB, both))
    pairs.sort(key=lambda t: -t[2])

    with open(os.path.join(RES, "cooccurrence_pairs_simple.tsv"), "w") as fh:
        fh.write(f"# TotalSamples(N)\t{N}\n")
        fh.write("GeneA\tGeneB\tnBoth\tPctOfTotal\tnGeneA\tnGeneB\tTotalSamples\n")
        for gA, gB, both in pairs:
            fh.write(f"{gA}\t{gB}\t{both}\t{100*both/N:.1f}\t{counts[gA]}\t"
                     f"{counts[gB]}\t{N}\n")
    log("wrote cooccurrence_pairs_simple.tsv")

    # =========================================================================
    # STEP 3: TRIPLET CO-OCCURRENCE COUNTS
    # =========================================================================
    trips = []
    for gA, gB, gC in itertools.combinations(genes, 3):
        both = int((present[gA] & present[gB] & present[gC]).sum())
        if both > 0:
            trips.append((gA, gB, gC, both))
    trips.sort(key=lambda t: -t[3])
    with open(os.path.join(RES, "cooccurrence_triples_simple.tsv"), "w") as fh:
        fh.write(f"# TotalSamples(N)\t{N}\n")
        fh.write("GeneA\tGeneB\tGeneC\tnAll3\tPctOfTotal\tTotalSamples\n")
        for gA, gB, gC, both in trips:
            fh.write(f"{gA}\t{gB}\t{gC}\t{both}\t{100*both/N:.1f}\t{N}\n")
    log("wrote cooccurrence_triples_simple.tsv")

    # =========================================================================
    # STEP 4: CO-OCCURRENCE COUNT MATRIX EXPORT
    # =========================================================================
    order = sorted(genes, key=lambda g: -counts[g])
    with open(os.path.join(RES, "cooccurrence_counts_matrix.tsv"), "w") as fh:
        fh.write(f"# Values = number of tumours mutated in BOTH genes. "
                 f"Diagonal = gene's own total. TotalSamples(N) = {N}\n")
        fh.write("Gene\t" + "\t".join(order) + "\n")
        for gA in order:
            row = [gA]
            for gB in order:
                if gA == gB:
                    row.append(str(counts[gA]))
                else:
                    row.append(str(int((present[gA] & present[gB]).sum())))
            fh.write("\t".join(row) + "\n")
    log("wrote cooccurrence_counts_matrix.tsv")

    # =========================================================================
    # STEP 5: RENDER FIGURE 17 — COUNT HEATMAP (TOP 16 DRIVERS)
    # =========================================================================
    top = order[:16]; n = len(top)
    M = np.zeros((n, n), dtype=int)
    for i, gA in enumerate(top):
        for j, gB in enumerate(top):
            M[i, j] = counts[gA] if i == j else int((present[gA] & present[gB]).sum())
    fig, ax = plt.subplots(figsize=(11, 9))
    im = ax.imshow(M, cmap="Reds")
    ax.set_xticks(range(n)); ax.set_xticklabels(top, rotation=90, fontsize=10)
    ax.set_yticks(range(n)); ax.set_yticklabels(top, fontsize=10)
    thr = M.max() * 0.6
    for i in range(n):
        for j in range(n):
            ax.text(j, i, str(M[i, j]), ha="center", va="center", fontsize=8,
                    color="white" if M[i, j] > thr else "black")
    cb = fig.colorbar(im, shrink=0.7); cb.set_label("Tumours carrying both genes")
    ax.set_title(f"Driver co-mutation counts (TCGA-COAD)\n"
                 f"cell = tumours with BOTH; diagonal = gene total; N = {N} tumours")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig17_cooccurrence_counts_heatmap.png"))
    plt.close(fig); log("wrote fig17_cooccurrence_counts_heatmap.png")

    # =========================================================================
    # STEP 6: RENDER FIGURE 18 — TOP CO-OCCURRING PAIRS BAR CHART
    # =========================================================================
    toppairs = pairs[:14][::-1]
    labels = [f"{a}+{b}" for a, b, _ in toppairs]
    vals = [both for _, _, both in toppairs]
    fig, ax = plt.subplots(figsize=(10, 6.5))
    ax.barh(range(len(toppairs)), vals, color="#B4433B")
    for i, (a, b, both) in enumerate(toppairs):
        ax.text(both + N*0.005, i, f"{both} / {N}  ({100*both/N:.1f}%)",
                va="center", fontsize=10)
    ax.set_yticks(range(len(toppairs))); ax.set_yticklabels(labels, fontsize=11)
    ax.set_xlabel(f"Number of tumours carrying BOTH  (out of N = {N})")
    ax.set_xlim(0, max(vals) * 1.35)
    ax.set_title(f"Most frequent driver co-mutations (TCGA-COAD, N = {N})")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig18_top_cooccurring_pairs.png"))
    plt.close(fig); log("wrote fig18_top_cooccurring_pairs.png")

    # Console preview logging
    log(f"Top co-occurring driver pairs (of N={N}):")
    for gA, gB, both in pairs[:10]:
        print(f"    {gA:7s}+{gB:8s} {both:3d} / {N}  ({100*both/N:.1f}%)")
    log(f"Top co-occurring driver triples (of N={N}):")
    for gA, gB, gC, both in trips[:6]:
        print(f"    {gA:6s}+{gB:6s}+{gC:6s} {both:3d} / {N}  ({100*both/N:.1f}%)")

if __name__ == "__main__":
    sys.exit(main())
