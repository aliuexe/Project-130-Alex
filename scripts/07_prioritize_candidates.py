#!/usr/bin/env python3
"""
07_prioritize_candidates.py
Project 130 - Colorectal cancer (TCGA-COAD)  --  ADVANCED component

Produces the prioritised neoantigen candidate shortlist (Section 14) using BigMHC.

Priority neoantigen candidates combine:
  - strong mutant peptide presentation (BigMHC_EL >= 0.50);
  - better mutant presentation than wild-type (positive DeltaPresentation);
  - expression of the mutated gene (GeneLevelTPM > 1.0);
  - occurrence in one or more tumour samples (MutationFrequency >= 2).

Input:  results/04_neoantigen_predictions.tsv
        results/03_integrated_mutation_expression.tsv
Output: results/neoantigen_candidates_shortlist.tsv
"""
import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
NEO = os.path.join(RES, "04_neoantigen_predictions.tsv")
INT = os.path.join(RES, "03_integrated_mutation_expression.tsv")
OUT = os.path.join(RES, "neoantigen_candidates_shortlist.tsv")

OUT_COLS = ["GeneName", "ProteinChange", "Peptide", "PeptideLength", "MutPos",
            "HLAAllele", "BigMHC_EL", "BigMHC_IM", "PresentationClass",
            "DeltaPresentation_MutMinusWT", "GeneLevelTPM", "MutationFrequency",
            "TranscriptID", "Chromosome", "Position", "Ref", "Alt"]

def log(m): print(f"[07] {m}", flush=True)

def per_mutation_freq():
    """(GeneName, AminoAcidChange) -> tumour samples carrying that mutation."""
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
    mf = per_mutation_freq()
    log(f"Per-mutation frequency map: {len(mf)} mutations")

    with open(NEO) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(hdr)}

    rows = []
    n = 0
    with open(NEO) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            n += 1
            if p[idx["PeptideType"]] != "Mutant":
                continue
            if p[idx["PeptideLength"]] != "9":
                continue
            if p[idx["PresentationClass"]] not in ("Strong", "Weak"):
                continue
            try:
                el = float(p[idx["BigMHC_EL"]])
                delta = float(p[idx["DeltaPresentation_MutMinusWT"]])
                tpm = float(p[idx["GeneLevelTPM"]])
            except ValueError:
                continue
            if delta <= 0 or tpm < 1.0:
                continue
            gene = p[idx["GeneName"]]
            pchg = p[idx["ProteinChange"]]
            freq = mf.get((gene, pchg), 0)
            if freq < 2:
                continue
            rec = {c: p[idx[c]] for c in hdr if c in idx}
            rec["MutationFrequency"] = freq
            rows.append((rec, el, freq, p[idx["PresentationClass"]]))

    # Rank: Strong presentation first, then highest BigMHC_EL presentation probability, then highest recurrence
    rows.sort(key=lambda t: (0 if t[3] == "Strong" else 1, -t[1], -t[2]))

    with open(OUT, "w") as fh:
        fh.write("\t".join(OUT_COLS) + "\n")
        for rec, _, _, _ in rows:
            fh.write("\t".join(str(rec.get(c, "NA")) for c in OUT_COLS) + "\n")

    log(f"Scanned {n} rows; wrote {len(rows)} candidates to {OUT}")
    genes = {r[0]["GeneName"] for r in rows}
    log(f"Distinct genes in shortlist: {len(genes)}")

if __name__ == "__main__":
    sys.exit(main())
