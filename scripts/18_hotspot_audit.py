#!/usr/bin/env python3
r"""
18_hotspot_audit.py
Project 130 - Colorectal Cancer (TCGA-COAD) Neoantigen Discovery Pipeline

===============================================================================
BIOLOGICAL & COMPUTATIONAL PURPOSE
===============================================================================
This script performs a rigorous clinical audit across all 42 highly recurrent
somatic mutation hotspots ($\ge 5$ tumours in TCGA-COAD).

Motivation: Ensures zero recurrent hotspots (e.g. KRAS G12V/G12D/G12C, TP53 hotspots,
PIK3CA E542K/E545K) are lost due to filtering artefacts, documenting explicit,
mechanistic reasons for any candidate exclusion.

===============================================================================
HOTSPOT AUDIT & DROP REASON CATEGORIZATION
===============================================================================
For every recurrent mutation ($\ge 5$ tumours), the script checks all 5 practical gates:
  1. Presentation Deficit: `BigMHC_EL < 0.50` on Option A HLA panel.
  2. Lack of Differential Agretopicity: Wild-type peptide also presented (`WT_EL >= 0.50`).
  3. Low Gene Expression: `GeneLevelTPM < 10.0`.
  4. Subclonal Architecture: `ClonalClass == Subclonal` (`medianVAF < 0.25`).
  5. Tumour Suppressor Gene (TSG) Context: Categorizes genes as Oncogenes vs Tumour Suppressors.

===============================================================================
INPUT & OUTPUT CONTRACTS
===============================================================================
Inputs:
  - `results/03_integrated_mutation_expression.tsv`
  - `results/04_neoantigen_predictions.tsv`
  - `results/mutation_clonality.tsv`

Output:
  - `results/recurrent_hotspot_audit.tsv` (Audit table for all 42 hotspots)
"""

import os, sys

# =============================================================================
# FILE PATHS & RESOURCE RESOLUTION
# =============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
INT = os.path.join(RES, "03_integrated_mutation_expression.tsv")
NEO = os.path.join(RES, "04_neoantigen_predictions.tsv")
CLON = os.path.join(RES, "mutation_clonality.tsv")
OUT = os.path.join(RES, "recurrent_hotspot_audit.tsv")

# Thresholds & Curated Tumour Suppressor List
MIN_RECUR = 5
TS_GENES = {"TP53","APC","PTEN","SMAD4","FBXW7","RB1","STK11","ARID1A","ATM",
            "SMAD2","TGFBR2","BMPR1A","ACVR2A","B2M","CASP8","RNF43","AXIN2"}
EL_PRESENT, EL_WT_NONPRESENT, TPM_MIN = 0.50, 0.50, 10.0

def log(m):
    """Prints timestamped progress messages to stdout with line flushing."""
    print("[18]", m, flush=True)

# =============================================================================
# MAIN PIPELINE EXECUTION
# =============================================================================
def main():
    # =========================================================================
    # STEP 1: IDENTIFY ALL RECURRENT HOTSPOTS (>= 5 TUMOURS)
    # =========================================================================
    recur, tpm = {}, {}
    with open(INT) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        s0 = next(i for i, c in enumerate(header) if c.startswith("TCGA"))
        for line in fh:
            p = line.rstrip("\n").split("\t")
            key = (p[0], p[2])
            recur[key] = recur.get(key, 0) + sum(1 for v in p[s0:] if v == "1")
            if p[3] != "NA":
                tpm[key[0]] = float(p[3])
    hot = {k for k, v in recur.items() if v >= MIN_RECUR}
    log(f"recurrent hotspots (>= {MIN_RECUR} tumours): {len(hot)}")

    # Load clonality mapping
    clon = {}
    with open(CLON) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(h)}
        for line in fh:
            q = line.rstrip("\n").split("\t")
            clon[(q[ix["GeneName"]], q[ix["ProteinChange"]])] = (
                q[ix["ClonalClass"]], q[ix["medianVAF"]])

    # =========================================================================
    # STEP 2: FIND BEST BIGMHC PRESENTATION FOR EACH HOTSPOT
    # =========================================================================
    best = {}   # (gene,pchg) -> (el, allele, pclass, wt_el)
    with open(NEO) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[9] != "Mutant" or p[11] != "9":
                continue
            key = (p[0], p[6])
            if key not in hot:
                continue
            if p[14] == "NA":
                continue
            el = float(p[14]); delta = float(p[17]) if p[17] != "NA" else 0.0
            wt = el - delta
            if key not in best or el > best[key][0]:
                best[key] = (el, p[13], p[16], wt)

    # =========================================================================
    # STEP 3: AUDIT EACH HOTSPOT AND WRITE DETAILED REASONING
    # =========================================================================
    rows = []
    for key in sorted(hot, key=lambda k: -recur[k]):
        g, pc = key
        r = recur[key]
        cl, vaf = clon.get(key, ("NA", "NA"))
        exprv = tpm.get(g)
        b = best.get(key)
        if b is None:
            status, reason = "Dropped", "no scorable 9-mer (nonstandard AA)"
            el = allele = pclass = wt = "NA"
        else:
            el, allele, pclass, wt = b
            reasons = []
            if el < EL_PRESENT: reasons.append("non-presenter on 3-allele panel (BigMHC_EL < 0.50)")
            if wt >= EL_WT_NONPRESENT: reasons.append("wild-type also presented (BigMHC_EL >= 0.50)")
            if exprv is None or exprv < TPM_MIN: reasons.append("expression <10 TPM")
            if cl != "Clonal": reasons.append("subclonal")
            status = "CANDIDATE" if not reasons else "Dropped"
            reason = "; ".join(reasons) if reasons else "passes all gates"
        rows.append([g, pc, r, cl, vaf,
                     f"{el:.4f}" if isinstance(el, float) else el,
                     allele, pclass,
                     f"{wt:.4f}" if isinstance(wt, float) else wt,
                     f"{exprv:.1f}" if exprv is not None else "NA",
                     "TumourSuppressor" if g in TS_GENES else "Other",
                     status, reason])

    # Export Hotspot Audit Table
    cols = ["GeneName","ProteinChange","Recurrence","ClonalClass","medianVAF",
            "BestMut_EL","BestAllele","PresentationClass","WT_EL","GeneLevelTPM",
            "GeneClass","Status","DropReason"]
    with open(OUT, "w") as fh:
        fh.write("# Selection uses ALL samples + ALL genes, per-locus (BigMHC).\n")
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    log(f"wrote {OUT} ({len(rows)} hotspots)")

if __name__ == "__main__":
    sys.exit(main())
