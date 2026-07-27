#!/usr/bin/env python3
r"""
19_longtail_coverage.py
Project 130 - Colorectal Cancer (TCGA-COAD) Neoantigen Discovery Pipeline

===============================================================================
BIOLOGICAL & COMPUTATIONAL PURPOSE
===============================================================================
This script performs long-tail coverage stratification analysis across the TCGA-COAD
cohort using BigMHC presentation probabilities.

Motivation: Shared recurrent neoantigens ($\ge 2$ tumours) cover a significant portion
of the population, but personalized / private mutations ($= 1$ tumour) form a massive
"long tail" in cancer genomes. This script quantifies the additional patient population
yield unlocked by adding private neoantigens to off-the-shelf shared vaccines.

===============================================================================
STRATIFICATION BINS & ACCUMULATION SCHEME
===============================================================================
Stratifies all quality candidates (`BigMHC_EL >= 0.50`, `WT_EL < 0.50`, `TPM >= 10.0`, `Clonal`)
into 4 cumulative recurrence tiers:
  1. Recurrent ($\ge 10$ tumours)
  2. Moderate ($5\text{--}9$ tumours)
  3. Low ($2\text{--}4$ tumours)
  4. Private ($= 1$ tumour, patient-specific)

===============================================================================
INPUT & OUTPUT CONTRACTS
===============================================================================
Inputs:
  - `results/04_neoantigen_predictions.tsv`
  - `results/03_integrated_mutation_expression.tsv`
  - `results/mutation_clonality.tsv`

Outputs:
  - `results/longtail_coverage_summary.txt` (Summary log of cumulative yield)
  - `figures/fig27_longtail_coverage.png` (Cumulative population yield bar chart)
"""

import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =============================================================================
# FILE PATHS & RESOURCE RESOLUTION
# =============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results"); FIG = os.path.join(BASE, "figures")
NEO = os.path.join(RES, "04_neoantigen_predictions.tsv")
INT = os.path.join(RES, "03_integrated_mutation_expression.tsv")
CLON = os.path.join(RES, "mutation_clonality.tsv")

# Threshold Constants
EL_PRESENT, EL_WT_NONPRESENT, TPM_MIN = 0.50, 0.50, 10.0

def log(m):
    """Prints timestamped progress messages to stdout with line flushing."""
    print("[19]", m, flush=True)

# =============================================================================
# MAIN PIPELINE EXECUTION
# =============================================================================
def main():
    # =========================================================================
    # STEP 1: PARSE CLONAL MUTATIONS
    # =========================================================================
    clonal = set()
    with open(CLON) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(h)}
        for line in fh:
            q = line.rstrip("\n").split("\t")
            if q[ix["ClonalClass"]] == "Clonal":
                clonal.add((q[ix["GeneName"]], q[ix["ProteinChange"]]))

    # =========================================================================
    # STEP 2: EXTRACT ALL QUALITY CANDIDATES (INCLUDING PRIVATE MUTATIONS)
    # =========================================================================
    cand = set()
    with open(NEO) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[9] != "Mutant" or p[11] != "9" or p[14] == "NA":
                continue
            el = float(p[14]); delta = float(p[17]) if p[17] != "NA" else 0.0
            wt = el - delta
            try:
                tpm = float(p[7])
            except ValueError:
                continue
            if el >= EL_PRESENT and wt < EL_WT_NONPRESENT and tpm >= TPM_MIN:
                key = (p[0], p[6])
                if key in clonal:
                    cand.add(key)
    log(f"quality candidate mutations (any recurrence, BigMHC): {len(cand)}")

    # =========================================================================
    # STEP 3: MAP PATIENT SAMPLE VECTORS PER STRATUM
    # =========================================================================
    with open(INT) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    s0 = next(i for i, c in enumerate(header) if c.startswith("TCGA"))
    N = len(header) - s0
    sampsets, freq = {}, {}
    with open(INT) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            key = (p[0], p[2])
            if key not in cand:
                continue
            v = np.fromiter((c == "1" for c in p[s0:]), dtype=bool, count=N)
            sampsets[key] = (sampsets[key] | v) if key in sampsets else v
            freq[key] = int(sampsets[key].sum())
    log(f"total tumours N = {N}")

    # =========================================================================
    # STEP 4: CUMULATIVE RECURRENCE STRATIFICATION ANALYSIS
    # =========================================================================
    strata = [("recurrent (>=10)", lambda f: f >= 10),
              ("moderate (5-9)",   lambda f: 5 <= f <= 9),
              ("low (2-4)",        lambda f: 2 <= f <= 4),
              ("private (=1)",     lambda f: f == 1)]
    covered = np.zeros(N, dtype=bool)
    cum = []
    for name, test in strata:
        keys = [k for k in sampsets if test(freq[k])]
        for k in keys:
            covered |= sampsets[k]
        cum.append((name, len(keys), int(covered.sum()), 100*covered.sum()/N))

    shared = np.zeros(N, dtype=bool)
    for k in sampsets:
        if freq[k] >= 2:
            shared |= sampsets[k]
    allcov = covered
    only_private = int((allcov & ~shared).sum())
    n_private = sum(1 for k in freq if freq[k] == 1)

    # Export summary text log
    with open(os.path.join(RES, "longtail_coverage_summary.txt"), "w") as fh:
        fh.write("LONG-TAIL CUMULATIVE POPULATION COVERAGE (TCGA-COAD, BigMHC)\n")
        fh.write(f"Total tumours (N): {N}\n")
        fh.write(f"Quality candidate mutations (any recurrence): {len(cand)}\n")
        fh.write(f"  of which private (in exactly 1 tumour): {n_private}\n\n")
        for name, nk, cov, pct in cum:
            fh.write(f"  + {name:18s} ({nk:5d} mutations): {cov:3d}/{N} = {pct:5.1f}%\n")
        fh.write(f"\nShared-only (recurrence >=2) coverage: "
                 f"{int(shared.sum())}/{N} = {100*shared.sum()/N:.1f}%\n")
        fh.write(f"Including the private long tail:        "
                 f"{int(allcov.sum())}/{N} = {100*allcov.sum()/N:.1f}%\n")

    # =========================================================================
    # STEP 5: RENDER FIGURE 27 — CUMULATIVE LONG-TAIL COVERAGE BAR CHART
    # =========================================================================
    labels = ["≥10", "+ 5-9", "+ 2-4", "+ private (=1)"]
    ys = [c[3] for c in cum]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    bars = ax.bar(range(len(ys)), ys,
                  color=["#4477AA","#66CCEE","#EE7733","#B4433B"], edgecolor="white")
    for i, (lab, c) in enumerate(zip(labels, cum)):
        ax.text(i, c[3]+1, f"{c[3]:.0f}%\n({c[2]}/{N})", ha="center", fontsize=10)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels([f"{l}\n(+{c[1]} muts)" for l, c in zip(labels, cum)], fontsize=10)
    ax.set_ylabel(f"% of tumours with ≥1 candidate neoantigen (N={N})")
    ax.set_ylim(0, 105)
    ax.set_title("Long-tail cumulative coverage (BigMHC)")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig27_longtail_coverage.png"), dpi=160)
    plt.close(fig); log("wrote fig27_longtail_coverage.png")

if __name__ == "__main__":
    sys.exit(main())
