#!/usr/bin/env python3
"""
08_aggregate_for_figures.py
Project 130 - Colorectal Cancer (TCGA-COAD) Neoantigen Discovery Pipeline

===============================================================================
BIOLOGICAL & COMPUTATIONAL PURPOSE
===============================================================================
This script performs a single, highly efficient streaming pass over the master
prediction table (`results/04_neoantigen_predictions.tsv`, 16.3M rows) to extract
and aggregate BigMHC summary statistics into a lightweight JSON summary file
(`results/figure_summary.json`).

===============================================================================
COMPUTATIONAL ARCHITECTURE & RESERVOIR SAMPLING
===============================================================================
- Memory Efficiency: Prevents loading the 2.4 GB prediction table into RAM during
  plotting. All subsequent figure generation scripts (Script 09) read from the
  compact `figure_summary.json` file in milliseconds.
- Reservoir Sampling: Maintains an unbiased 6,000-point random sample reservoir
  of (wild-type, mutant) presentation probability pairs `(WT_EL, Mut_EL)` for
  scatter plot visualization (Figure 11) without storing millions of data points.

===============================================================================
SUMMARY METRICS CAPTURED
===============================================================================
  1. `BigMHC_EL` Presentation Histograms (50 bins across [0.0, 1.0]).
  2. `DeltaPresentation` Histograms (51 bins across [-1.0, 1.0]).
  3. Directional Counts: Count of positive, negative, and zero delta pairs.
  4. Presenter Class Breakdown: Strong, Weak, Non-presenter counts per peptide type.
  5. Allele Breakdown: Strong presenters categorized per HLA allele.

===============================================================================
INPUT & OUTPUT CONTRACTS
===============================================================================
Input:
  - `results/04_neoantigen_predictions.tsv` (Deliverable 04)

Output:
  - `results/figure_summary.json` (Compact summary JSON)
"""

import json
import os
import random
import sys
import math

# =============================================================================
# FILE PATHS & RESOURCE RESOLUTION
# =============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
NEO = os.path.join(RES, "04_neoantigen_predictions.tsv")
OUT = os.path.join(RES, "figure_summary.json")

# Seed random number generator for reproducible reservoir sampling
random.seed(0)

# =============================================================================
# HISTOGRAM BINNING FUNCTIONS
# =============================================================================
# BigMHC_EL presentation probability bins [0.0 .. 1.0]
NB = 50
LO, HI = 0.0, 1.0
def bin_el(el):
    """Maps BigMHC_EL presentation probability to bin index [0, NB-1]."""
    if el <= 0: return 0
    i = int((el - LO) / (HI - LO) * NB)
    return min(max(i, 0), NB - 1)

# BigMHC DeltaPresentation bins [-1.0 .. 1.0]
DNB = 51
DLO, DHI = -1.0, 1.0
def deltabin(d):
    """Maps DeltaPresentation to bin index [0, DNB-1]."""
    d = max(DLO, min(DHI, d))
    i = int((d - DLO) / (DHI - DLO) * DNB)
    return min(max(i, 0), DNB - 1)

# =============================================================================
# MAIN PIPELINE EXECUTION
# =============================================================================
def main():
    idx = None
    el_hist = {"Mutant": [0]*NB, "WildType": [0]*NB}
    delta_hist = [0]*DNB
    delta_pos = 0; delta_neg = 0; delta_zero = 0
    presenter_by_type = {}       # (ptype, presentationclass) -> count
    strong_by_allele = {}        # allele -> strong mutant count
    allele_presenter = {}        # (allele, presentationclass) -> mutant count
    scatter = []                 # Reservoir sample of (wt_el, mut_el)
    seen_mut = 0
    n = 0

    # =========================================================================
    # STEP 1: STREAM DELIVERABLE 04 AND AGGREGATE METRICS
    # =========================================================================
    with open(NEO) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        iL = idx["PeptideLength"]; iT = idx["PeptideType"]
        iA = idx["HLAAllele"]; iEL = idx["BigMHC_EL"]
        iB = idx["PresentationClass"]; iD = idx["DeltaPresentation_MutMinusWT"]
        
        for line in fh:
            p = line.rstrip("\n").split("\t")
            n += 1
            if p[iL] != "9":
                continue
            ptype = p[iT]; allele = p[iA]; pcls = p[iB]
            presenter_by_type[(ptype, pcls)] = presenter_by_type.get((ptype, pcls), 0) + 1
            el_str = p[iEL]
            if el_str != "NA":
                el_val = float(el_str)
                el_hist[ptype][bin_el(el_val)] += 1
            if ptype == "Mutant":
                allele_presenter[(allele, pcls)] = allele_presenter.get((allele, pcls), 0) + 1
                if pcls == "Strong":
                    strong_by_allele[allele] = strong_by_allele.get(allele, 0) + 1
                d = p[iD]
                if d != "NA":
                    dv = float(d)
                    delta_hist[deltabin(dv)] += 1
                    if dv > 0: delta_pos += 1
                    elif dv < 0: delta_neg += 1
                    else: delta_zero += 1
                    
                    # Reservoir Sampling for Scatter Plot (Figure 11)
                    if el_str != "NA":
                        mut_el = float(el_str)
                        wt_el = mut_el - dv
                        seen_mut += 1
                        if len(scatter) < 6000:
                            scatter.append((wt_el, mut_el))
                        else:
                            j = random.randint(0, seen_mut - 1)
                            if j < 6000:
                                scatter[j] = (wt_el, mut_el)

    # =========================================================================
    # STEP 2: BUILD DICTIONARY & WRITE FIGURE SUMMARY JSON
    # =========================================================================
    summary = dict(
        NB=NB, LO=LO, HI=HI, DNB=DNB, DLO=DLO, DHI=DHI,
        el_hist=el_hist, delta_hist=delta_hist,
        delta_pos=delta_pos, delta_neg=delta_neg, delta_zero=delta_zero,
        binder_by_type={f"{k[0]}|{k[1]}": v for k, v in presenter_by_type.items()},
        strong_by_allele=strong_by_allele,
        allele_binder={f"{k[0]}|{k[1]}": v for k, v in allele_presenter.items()},
        scatter=scatter, rows_scanned=n,
    )
    with open(OUT, "w") as fh:
        json.dump(summary, fh)
    print(f"[08] scanned {n} rows; wrote {OUT}", flush=True)
    print(f"[08] BigMHC delta pos/neg/zero = {delta_pos}/{delta_neg}/{delta_zero}")
    print(f"[08] BigMHC strong-by-allele (mutant) = {strong_by_allele}")

if __name__ == "__main__":
    sys.exit(main())
