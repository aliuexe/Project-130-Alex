#!/usr/bin/env python3
# =============================================================================
# 07_prioritize_candidates.py   (ANNOTATED teaching copy)
# =============================================================================
# WHAT THIS SCRIPT DOES:
#   The full neoantigen table (Deliverable 04) has ~16 million rows — far too
#   many to act on. This script filters it down to a short list of the most
#   PROMISING candidate neoantigens using the criteria in assignment Section 14.
#
# THE FOUR PRIORITISATION CRITERIA (why each one?):
#   1. Strong or Weak binder  -> the HLA molecule can actually display it.
#   2. Mutant binds better than wild-type (DeltaAffinity > 0)
#                              -> the mutation MAKES the peptide more visible;
#                                 the normal peptide is less well presented.
#   3. Gene is expressed (GeneLevelTPM > 1)
#                              -> the mutated protein is actually being made.
#   4. Mutation recurs in >= 2 tumours
#                              -> a shared neoantigen relevant to more patients.
#   A candidate that ticks all four is a much better bet than a random mutation.
#   (These are PRIORITISATION features, not proof of immunogenicity — Rule 8.)
#
# INPUT : results/04_neoantigen_predictions.tsv + the integrated matrix (03)
# OUTPUT: results/neoantigen_candidates_shortlist.tsv
# =============================================================================

import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
NEO = os.path.join(RES, "04_neoantigen_predictions.tsv")          # the full table
INT = os.path.join(RES, "03_integrated_mutation_expression.tsv")  # for recurrence
OUT = os.path.join(RES, "neoantigen_candidates_shortlist.tsv")

# The columns we keep in the shortlist (a readable subset of Deliverable 04).
OUT_COLS = ["GeneName", "ProteinChange", "Peptide", "PeptideLength", "MutPos",
            "HLAAllele", "BindingAffinity", "BindingRank", "BinderClass",
            "DeltaAffinity_WTminusMut", "GeneLevelTPM", "MutationFrequency",
            "TranscriptID", "Chromosome", "Position", "Ref", "Alt"]

def log(m): print(f"[07] {m}", flush=True)

def per_mutation_freq():
    # Recompute, from the integrated matrix, how many tumour samples carry each
    # SPECIFIC mutation. (We do this here too so the shortlist is correct even
    # if run on its own.)  Key = (gene, protein-change), value = sample count.
    mf = {}
    with open(INT) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # locate sample columns by their "TCGA" prefix (robust to metadata
        # columns like GeneLevelTPM_SD)
        s0 = next(i for i, c in enumerate(header) if c.startswith("TCGA"))
        for line in fh:
            p = line.rstrip("\n").split("\t")
            k = (p[0], p[2])                            # (gene, protein change)
            mf[k] = mf.get(k, 0) + sum(1 for v in p[s0:] if v == "1")
    return mf

def main():
    mf = per_mutation_freq()
    log(f"Per-mutation frequency map: {len(mf)} mutations")

    # Find which column number holds each field name in Deliverable 04.
    with open(NEO) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(hdr)}

    rows = []
    n = 0
    with open(NEO) as fh:
        fh.readline()
        for line in fh:                                 # go through all ~16M rows
            p = line.rstrip("\n").split("\t")
            n += 1
            # --- apply the four criteria, cheapest checks first ---
            if p[idx["PeptideType"]] != "Mutant":       # only mutant peptides
                continue
            if p[idx["PeptideLength"]] != "9":          # only class-I 9-mers
                continue
            if p[idx["BinderClass"]] not in ("Strong", "Weak"):   # must bind
                continue
            try:
                aff = float(p[idx["BindingAffinity"]])
                delta = float(p[idx["DeltaAffinity_WTminusMut"]])
                tpm = float(p[idx["GeneLevelTPM"]])
            except ValueError:                          # NA values -> skip
                continue
            if delta <= 0 or tpm < 1:                   # better-than-WT + expressed
                continue
            gene = p[idx["GeneName"]]
            pchg = p[idx["ProteinChange"]]
            freq = mf.get((gene, pchg), 0)
            if freq < 2:                                # recurrent in >=2 tumours
                continue
            rec = {c: p[idx[c]] for c in hdr if c in idx}
            rec["MutationFrequency"] = freq             # use the corrected recurrence
            rows.append((rec, aff, freq, p[idx["BinderClass"]]))

    # ---- Rank the survivors ---------------------------------------------------
    # Sort key: Strong binders first (0 before 1), then strongest affinity
    # (lowest nM), then most recurrent (negative so higher counts come first).
    rows.sort(key=lambda t: (0 if t[3] == "Strong" else 1, t[1], -t[2]))

    with open(OUT, "w") as fh:
        fh.write("\t".join(OUT_COLS) + "\n")
        for rec, _, _, _ in rows:
            fh.write("\t".join(str(rec.get(c, "NA")) for c in OUT_COLS) + "\n")

    log(f"Scanned {n} rows; wrote {len(rows)} candidates to {OUT}")
    genes = {r[0]["GeneName"] for r in rows}
    log(f"Distinct genes in shortlist: {len(genes)}")

if __name__ == "__main__":
    sys.exit(main())
