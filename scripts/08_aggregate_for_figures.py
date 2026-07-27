#!/usr/bin/env python3
"""
08_aggregate_for_figures.py
Project 130 - single streaming pass over deliverable 04 to collect BigMHC
summary statistics needed for presentation figures.
"""
import json
import os
import random
import sys
import math

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
NEO = os.path.join(RES, "04_neoantigen_predictions.tsv")
OUT = os.path.join(RES, "figure_summary.json")

random.seed(0)

# BigMHC_EL presentation probability bins [0.0 .. 1.0]
NB = 50
LO, HI = 0.0, 1.0
def bin_el(el):
    if el <= 0: return 0
    i = int((el - LO) / (HI - LO) * NB)
    return min(max(i, 0), NB - 1)

# BigMHC DeltaPresentation bins [-1.0 .. 1.0]
DNB = 51
DLO, DHI = -1.0, 1.0
def deltabin(d):
    d = max(DLO, min(DHI, d))
    i = int((d - DLO) / (DHI - DLO) * DNB)
    return min(max(i, 0), DNB - 1)

def main():
    idx = None
    el_hist = {"Mutant": [0]*NB, "WildType": [0]*NB}
    delta_hist = [0]*DNB
    delta_pos = 0; delta_neg = 0; delta_zero = 0
    presenter_by_type = {}       # (ptype, presentationclass) -> count
    strong_by_allele = {}        # allele -> strong mutant count
    allele_presenter = {}        # (allele, presentationclass) -> mutant count
    scatter = []                 # reservoir of (wt_el, mut_el)
    seen_mut = 0
    n = 0

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
