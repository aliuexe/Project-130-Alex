#!/usr/bin/env python3
r"""
14_practical_neoantigens_coverage.py
Project 130 - Colorectal Cancer (TCGA-COAD) Neoantigen Discovery Pipeline

===============================================================================
BIOLOGICAL & COMPUTATIONAL PURPOSE
===============================================================================
This script defines the 5-tier "Practical Neoantigen" clinical filtering framework
and evaluates population-level vaccine coverage using a greedy set-cover optimization
algorithm.

===============================================================================
PRACTICAL NEOANTIGEN SELECTION RULES (5-TIER GATES)
===============================================================================
A candidate neoantigen qualifies as a "Practical Neoantigen" if it simultaneously
passes five stringent therapeutic gates:
  1. Presentation Gate: `Mutant_EL >= 0.50` (predicted HLA Class I cell-surface presenter).
  2. Differential Agretopicity Gate: `WT_EL < 0.50` (wild-type peptide is NOT presented;
     ensures true mutant-specific neo-epitope creation).
  3. Abundant Expression Gate: `GeneLevelTPM >= 10.0` (mutated gene is abundantly transcribed).
  4. Truncal Clonality Gate: `ClonalClass == Clonal` (median VAF >= 0.25; present in 100% of tumour cells).
  5. Recurrence Gate: `MutationFrequency >= 2` (shared target present in >= 2 tumours).

===============================================================================
GREEDY SET-COVER OPTIMIZATION ALGORITHM
===============================================================================
To maximize population coverage with a minimum cocktail size:
  - Iteratively selects the practical neoantigen that maximizes cumulative union
    of covered tumour samples:
        \text{Select } k^* = \arg\max_k | S_k \setminus C_{\text{current}} |
  - Tracks cumulative coverage for both $\ge 1$ epitope per patient and $\ge 2$
    epitopes per patient (multi-epitope redundancy for immune escape prevention).

===============================================================================
INPUT & OUTPUT CONTRACTS
===============================================================================
Inputs:
  - `results/neoantigen_candidates_shortlist.tsv`
  - `results/mutation_clonality.tsv`
  - `results/03_integrated_mutation_expression.tsv`

Outputs:
  - `results/practical_neoantigens.tsv` (1,113 practical candidate database)
  - `results/vaccine_coverage_curve.tsv` (Cumulative population coverage trajectory)
  - `figures/fig20_vaccine_coverage.png` (Population coverage curve plot)
  - `figures/fig21_top_practical_neoantigens.png` (Top practical neoantigens bar chart)
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
SHORT = os.path.join(RES, "neoantigen_candidates_shortlist.tsv")
CLON = os.path.join(RES, "mutation_clonality.tsv")
INT = os.path.join(RES, "03_integrated_mutation_expression.tsv")
OUT_PRAC = os.path.join(RES, "practical_neoantigens.tsv")
OUT_COV = os.path.join(RES, "vaccine_coverage_curve.tsv")

# Threshold Constants
TPM_MIN = 10.0          # Abundant gene expression cutoff
EL_PRESENT = 0.50       # Mutant BigMHC_EL presentation cutoff
EL_WT_NONPRESENT = 0.50 # Wild-type BigMHC_EL non-presentation cutoff

def log(m):
    """Prints timestamped progress messages to stdout with line flushing."""
    print("[14]", m, flush=True)

# =============================================================================
# MAIN PIPELINE EXECUTION
# =============================================================================
def main():
    # =========================================================================
    # STEP 1: LOAD CLONALITY CLASSIFICATIONS (SCRIPT 13)
    # =========================================================================
    clonal = {}
    with open(CLON) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(h)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            clonal[(p[ix["GeneName"]], p[ix["ProteinChange"]])] = p[ix["ClonalClass"]]

    # =========================================================================
    # STEP 2: APPLY 5-TIER PRACTICAL NEOANTIGEN GATES
    # =========================================================================
    best = {}   # (gene, pchg) -> dict of best practical record
    with open(SHORT) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(h)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            gene = p[ix["GeneName"]]; pchg = p[ix["ProteinChange"]]
            try:
                mt_el = float(p[ix["BigMHC_EL"]])
                delta = float(p[ix["DeltaPresentation_MutMinusWT"]])
                tpm = float(p[ix["GeneLevelTPM"]])
                freq = int(p[ix["MutationFrequency"]])
            except ValueError:
                continue
            wt_el = mt_el - delta
            # Gate 1, 2, 3: Presentation, Differential Agretopicity, Expression
            if not (mt_el >= EL_PRESENT and wt_el < EL_WT_NONPRESENT and tpm >= TPM_MIN):
                continue
            # Gate 4: Truncal Clonality
            if clonal.get((gene, pchg)) != "Clonal":
                continue
            # Gate 5: Recurrence
            if freq < 2:
                continue
            rec = dict(GeneName=gene, ProteinChange=pchg,
                       Peptide=p[ix["Peptide"]], HLAAllele=p[ix["HLAAllele"]],
                       Mutant_EL=round(mt_el, 4), WT_EL=round(wt_el, 4),
                       DeltaPresentation=round(delta, 4), GeneLevelTPM=round(tpm, 1),
                       MutationFrequency=freq)
            k = (gene, pchg)
            if k not in best or mt_el > best[k]["Mutant_EL"]:
                best[k] = rec
    practical = list(best.values())
    log(f"Practical neoantigens (distinct mutations): {len(practical)}")

    # =========================================================================
    # STEP 3: MAP PATIENT SAMPLE COVERAGE FOR PRACTICAL NEOANTIGENS
    # =========================================================================
    wanted = set(best.keys())
    with open(INT) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    s0 = next(i for i, c in enumerate(header) if c.startswith("TCGA"))
    N = len(header) - s0
    sampsets = {k: np.zeros(N, dtype=bool) for k in wanted}
    with open(INT) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            k = (p[0], p[2])
            if k in wanted:
                v = np.fromiter((c == "1" for c in p[s0:]), dtype=bool, count=N)
                sampsets[k] |= v
    log(f"Total tumours (N) = {N}")

    for rec in practical:
        rec["TumoursCovered"] = int(sampsets[(rec["GeneName"], rec["ProteinChange"])].sum())
    practical.sort(key=lambda r: -r["TumoursCovered"])

    # Export Practical Neoantigen Database
    cols = ["GeneName", "ProteinChange", "Peptide", "HLAAllele", "Mutant_EL",
            "WT_EL", "DeltaPresentation", "GeneLevelTPM", "MutationFrequency",
            "TumoursCovered"]
    with open(OUT_PRAC, "w") as fh:
        fh.write(f"# TotalSamples(N)\t{N}\n")
        fh.write("\t".join(cols) + "\n")
        for r in practical:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    log(f"wrote {OUT_PRAC}")

    # =========================================================================
    # STEP 4: GREEDY SET-COVER OPTIMIZATION FOR POPULATION COVERAGE
    # =========================================================================
    keys = [(r["GeneName"], r["ProteinChange"]) for r in practical]
    covered_ge1 = np.zeros(N, dtype=bool)
    hit_count = np.zeros(N, dtype=int)
    remaining = set(keys)
    order = []
    
    # Greedy set cover loop
    while remaining:
        best_k = None; best_gain = -1
        for k in remaining:
            gain = int((sampsets[k] & ~covered_ge1).sum())
            if gain > best_gain:
                best_gain = gain; best_k = k
        if best_gain <= 0:
            break
        remaining.discard(best_k)
        covered_ge1 |= sampsets[best_k]
        hit_count += sampsets[best_k].astype(int)
        order.append((best_k, int(covered_ge1.sum()), int((hit_count >= 2).sum())))

    # Write coverage curve trajectory
    with open(OUT_COV, "w") as fh:
        fh.write(f"# TotalSamples(N)\t{N}\n")
        fh.write("Step\tGeneName\tProteinChange\tCumTumours_ge1\tPct_ge1\t"
                 "CumTumours_ge2\tPct_ge2\n")
        for i, ((g, pc), c1, c2) in enumerate(order, 1):
            fh.write(f"{i}\t{g}\t{pc}\t{c1}\t{100*c1/N:.1f}\t{c2}\t{100*c2/N:.1f}\n")
    log(f"wrote {OUT_COV}: {len(order)} neoantigens cover "
        f"{order[-1][1]} tumours ({100*order[-1][1]/N:.1f}%)")

    # =========================================================================
    # STEP 5: RENDER FIGURE 20 — CUMULATIVE VACCINE COVERAGE CURVE
    # =========================================================================
    xs = list(range(1, len(order) + 1))
    y1 = [100*o[1]/N for o in order]
    y2 = [100*o[2]/N for o in order]
    fig, ax = plt.subplots(figsize=(9.5, 6))
    ax.plot(xs, y1, "-o", color="#B4433B", ms=3, label="≥1 epitope per tumour")
    ax.plot(xs, y2, "-o", color="#2A9D8F", ms=3, label="≥2 epitopes per tumour")
    for target in (10, 20, 30):
        if target <= len(order):
            ax.annotate(f"{y1[target-1]:.0f}%", (target, y1[target-1]),
                        textcoords="offset points", xytext=(0, 8), fontsize=9, color="#B4433B")
    ax.set_xlabel("Number of shared neoantigens in the vaccine (greedy-selected)")
    ax.set_ylabel(f"% of tumours covered  (N = {N})")
    ax.set_title("Off-the-shelf vaccine coverage — practical shared neoantigens (BigMHC)\n"
                 "(clonal, expressed, mutant-specific presenters)")
    ax.legend(frameon=False); ax.set_ylim(0, max(y1)*1.1 + 1 if y1 else 100)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig20_vaccine_coverage.png"), dpi=160)
    plt.close(fig); log("wrote fig20_vaccine_coverage.png")

    # =========================================================================
    # STEP 6: RENDER FIGURE 21 — TOP PRACTICAL NEOANTIGENS BAR CHART
    # =========================================================================
    top = practical[:15][::-1]
    labels = [f"{r['GeneName']} {r['ProteinChange']}" for r in top]
    vals = [r["TumoursCovered"] for r in top]
    fig, ax = plt.subplots(figsize=(10, 6.5))
    ax.barh(range(len(top)), vals, color="#4477AA")
    for i, r in enumerate(top):
        ax.text(r["TumoursCovered"] + N*0.004, i,
                f"{r['TumoursCovered']} tum · {r['HLAAllele']} · EL={r['Mutant_EL']:.3f}",
                va="center", fontsize=8)
    ax.set_yticks(range(len(top))); ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel(f"Tumours covered (out of N = {N})")
    ax.set_xlim(0, max(vals)*1.5 if vals else 10)
    ax.set_title("Top practical shared neoantigens (TCGA-COAD, BigMHC)\n"
                 "clonal + expressed + mutant-specific presenter")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig21_top_practical_neoantigens.png"), dpi=160)
    plt.close(fig); log("wrote fig21_top_practical_neoantigens.png")

if __name__ == "__main__":
    sys.exit(main())
