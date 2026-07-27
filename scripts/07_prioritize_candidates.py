#!/usr/bin/env python3
"""
07_prioritize_candidates.py
Project 130 - Colorectal Cancer (TCGA-COAD) Neoantigen Discovery Pipeline

===============================================================================
BIOLOGICAL & COMPUTATIONAL PURPOSE
===============================================================================
This script implements Section 14 of the project assignment. It filters and
prioritizes the 16.3M master prediction table (Deliverable 04) into a highly
focused candidate shortlist (`results/neoantigen_candidates_shortlist.tsv`).

===============================================================================
SHORTLIST PRIORITIZATION CRITERIA (§14)
===============================================================================
A mutant 9-mer peptide qualifies for the candidate shortlist if it passes four
simultaneous filtering gates:
  1. Presentation Gate: `PresentationClass` in (`Strong`, `Weak`), corresponding to
     `BigMHC_EL >= 0.50` (predicted cell-surface presentation).
  2. Differential Agretopicity Gate: `DeltaPresentation_MutMinusWT > 0`
     (mutant peptide is presented more strongly than its wild-type counterpart).
  3. Gene Expression Gate: `GeneLevelTPM > 1.0` (mutated gene is actively transcribed).
  4. Recurrence Gate: `MutationFrequency >= 2` (mutation is present in 2 or more
     tumour samples across the cohort).

===============================================================================
MULTI-KEY RANKING HIERARCHY
===============================================================================
Candidates are ranked deterministically by a 3-tier tuple:
  1. Primary: Presentation class (`Strong` binders prioritized over `Weak`).
  2. Secondary: `BigMHC_EL` presentation probability (descending, higher = better).
  3. Tertiary: Per-mutation recurrence `MutationFrequency` (descending).

===============================================================================
INPUT & OUTPUT CONTRACTS
===============================================================================
Inputs:
  - `results/04_neoantigen_predictions.tsv` (Deliverable 04)
  - `results/03_integrated_mutation_expression.tsv` (Deliverable 03)

Output:
  - `results/neoantigen_candidates_shortlist.tsv` (Prioritized candidate shortlist)
"""

import os
import sys

# =============================================================================
# FILE PATHS & RESOURCE RESOLUTION
# =============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
NEO = os.path.join(RES, "04_neoantigen_predictions.tsv")
INT = os.path.join(RES, "03_integrated_mutation_expression.tsv")
OUT = os.path.join(RES, "neoantigen_candidates_shortlist.tsv")

# Shortlist Output Column Order
OUT_COLS = ["GeneName", "ProteinChange", "Peptide", "PeptideLength", "MutPos",
            "HLAAllele", "BigMHC_EL", "BigMHC_IM", "PresentationClass",
            "DeltaPresentation_MutMinusWT", "GeneLevelTPM", "MutationFrequency",
            "TranscriptID", "Chromosome", "Position", "Ref", "Alt"]

def log(m):
    """Prints timestamped progress messages to stdout with line flushing."""
    print(f"[07] {m}", flush=True)

def per_mutation_freq():
    """
    Computes per-mutation recurrence map: (GeneName, ProteinChange) -> tumour sample count.
    Re-computed directly from Deliverable 03 for standalone script execution.
    """
    mf = {}
    with open(INT) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        s0 = next(i for i, c in enumerate(header) if c.startswith("TCGA"))
        for line in fh:
            p = line.rstrip("\n").split("\t")
            k = (p[0], p[2])
            mf[k] = mf.get(k, 0) + sum(1 for v in p[s0:] if v == "1")
    return mf

def main():
    # =========================================================================
    # STEP 1: BUILD MUTATION RECURRENCE LOOKUP MAP
    # =========================================================================
    mf = per_mutation_freq()
    log(f"Per-mutation frequency map: {len(mf)} mutations")

    # Dynamic column index lookup for Deliverable 04 header
    with open(NEO) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(hdr)}

    # =========================================================================
    # STEP 2: STREAM DELIVERABLE 04 AND APPLY SHORTLIST GATES
    # =========================================================================
    rows = []
    n = 0
    with open(NEO) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            n += 1
            # Gate 1: Class I 9-mer Mutant peptides only
            if p[idx["PeptideType"]] != "Mutant":
                continue
            if p[idx["PeptideLength"]] != "9":
                continue
            # Gate 2: Presentation class must be Strong (>=0.70) or Weak (>=0.50)
            if p[idx["PresentationClass"]] not in ("Strong", "Weak"):
                continue
            try:
                el = float(p[idx["BigMHC_EL"]])
                delta = float(p[idx["DeltaPresentation_MutMinusWT"]])
                tpm = float(p[idx["GeneLevelTPM"]])
            except ValueError:
                continue
            # Gate 3: Differential Agretopicity (delta > 0) and Gene Expression (TPM > 1.0)
            if delta <= 0 or tpm < 1.0:
                continue
            gene = p[idx["GeneName"]]
            pchg = p[idx["ProteinChange"]]
            freq = mf.get((gene, pchg), 0)
            # Gate 4: Mutation Recurrence across tumours (>= 2)
            if freq < 2:
                continue
            rec = {c: p[idx[c]] for c in hdr if c in idx}
            rec["MutationFrequency"] = freq
            rows.append((rec, el, freq, p[idx["PresentationClass"]]))

    # =========================================================================
    # STEP 3: MULTI-KEY RANKING & SORTING
    # =========================================================================
    # Rank: Strong presentation first, then highest BigMHC_EL, then highest recurrence
    rows.sort(key=lambda t: (0 if t[3] == "Strong" else 1, -t[1], -t[2]))

    # =========================================================================
    # STEP 4: WRITE CANDIDATE SHORTLIST FILE
    # =========================================================================
    with open(OUT, "w") as fh:
        fh.write("\t".join(OUT_COLS) + "\n")
        for rec, _, _, _ in rows:
            fh.write("\t".join(str(rec.get(c, "NA")) for c in OUT_COLS) + "\n")

    log(f"Scanned {n} rows; wrote {len(rows)} candidates to {OUT}")
    genes = {r[0]["GeneName"] for r in rows}
    log(f"Distinct genes in shortlist: {len(genes)}")

if __name__ == "__main__":
    sys.exit(main())
