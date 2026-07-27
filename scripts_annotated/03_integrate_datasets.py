#!/usr/bin/env python3
# =============================================================================
# 03_integrate_datasets.py   (ANNOTATED teaching copy)
# =============================================================================
# WHAT THIS SCRIPT DOES:
#   It joins the two grids we already built — the mutation grid (01) and the
#   expression grid (02) — into one combined table, adding a per-gene
#   expression summary: the median TPM AND its spread (standard deviation).
#
# WHY WE NEED IT (biology + a data reality):
#   We want each mutation to carry information about how expressed its gene is.
#   BUT the patients whose DNA was sequenced for mutations are not exactly the
#   same patients as those whose RNA was sequenced. So we cannot line them up
#   sample-by-sample. Instead we summarise expression to ONE representative
#   value per gene for the whole cancer type: the MEDIAN TPM across all tumours.
#
# WHY THE MEDIAN (not the average/mean)?
#   The median (middle value) is "robust to outliers": a few tumours with
#   freakishly high expression won't distort it, whereas the mean would be
#   dragged up. The assignment recommends the median for exactly this reason.
#
# WHY ALSO REPORT THE SPREAD (SD / IQR)?  (added at the PI's request)
#   A single median hides how variable a gene's expression is. If a gene's TPM
#   swings wildly across tumours, the median is a less trustworthy summary. So
#   we also report:
#     GeneLevelTPM_SD  = standard deviation (the spread the PI asked for)
#     GeneLevelTPM_IQR = inter-quartile range (Q3-Q1); a robust spread measure
#                        that pairs naturally with the median
#     n_expr_samples   = how many tumours had a real (non-missing) value
#   A large SD/IQR relative to the median = "treat this gene's median with
#   caution". (Bulk TPM still cannot tell mutant- from wild-type-allele
#   expression; that would need read-level RNA data - see the methodology guide.)
#
# INPUT : results/01_mutation_by_sample.tsv  and  results/02_gene_by_sample_TPM.tsv
# OUTPUT: results/03_integrated_mutation_expression.tsv  (+ gene_level_tpm.tsv)
# =============================================================================

import os
import sys
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MUT = os.path.join(BASE, "results", "01_mutation_by_sample.tsv")   # mutation grid
TPM = os.path.join(BASE, "results", "02_gene_by_sample_TPM.tsv")   # expression grid
OUT = os.path.join(BASE, "results", "03_integrated_mutation_expression.tsv")
GLT = os.path.join(BASE, "results", "gene_level_tpm.tsv")          # gene->stats lookup

def log(m): print(f"[03] {m}", flush=True)

def main():
    log("Computing GeneLevelTPM from Deliverable 02")
    tpm = pd.read_csv(TPM, sep="\t", low_memory=False)
    tpm_sample_cols = [c for c in tpm.columns if c != "GeneName"]   # sample columns
    tpm_idx = tpm.set_index("GeneName")[tpm_sample_cols]            # genes x samples

    # ---- Summarise each gene's expression across all tumour samples ----------
    # For each gene (row) we compute, ignoring missing (NA) cells:
    med = tpm_idx.median(axis=1, skipna=True)                       # median TPM
    sd = tpm_idx.std(axis=1, skipna=True)                           # standard deviation
    iqr = (tpm_idx.quantile(0.75, axis=1)                           # Q3 - Q1 (spread)
           - tpm_idx.quantile(0.25, axis=1))
    n_expr = tpm_idx.notna().sum(axis=1)                            # #samples with a value

    # Save a small per-gene summary table with all four numbers.
    gene_level = pd.DataFrame({
        "GeneName": med.index,
        "GeneLevelTPM": med.values.round(4),
        "GeneLevelTPM_SD": sd.values.round(4),
        "GeneLevelTPM_IQR": iqr.values.round(4),
        "n_expr_samples": n_expr.values.astype(int),
    })
    gene_level.to_csv(GLT, sep="\t", index=False)
    log(f"Computed GeneLevelTPM (median + SD + IQR) for {len(gene_level)} genes")

    # Fast lookup tables: gene name -> its median, and gene name -> its SD.
    med_map = dict(zip(med.index, med.round(4)))
    sd_map = dict(zip(sd.index, sd.round(4)))

    # ---- Add GeneLevelTPM and GeneLevelTPM_SD into the mutation grid ---------
    # The mutation grid is large, so we read it line by line ("streaming") and
    # write the new file as we go, to keep memory use tiny.
    # New column order: GeneName, Mutation, AminoAcidChange, GeneLevelTPM,
    # GeneLevelTPM_SD, then the sample columns.
    log("Streaming Deliverable 01 and inserting GeneLevelTPM + SD by gene")
    n_rows = 0
    n_na = 0

    def fmt(v):
        # format a number to 4 decimals, or "NA" if it's missing.
        return "NA" if v is None or (isinstance(v, float) and np.isnan(v)) \
               else f"{v:.4f}"

    with open(MUT) as fin, open(OUT, "w") as fout:
        header = fin.readline().rstrip("\n").split("\t")
        sample_cols = header[3:]                                    # columns 3+ are samples
        out_header = (["GeneName", "Mutation", "AminoAcidChange",
                       "GeneLevelTPM", "GeneLevelTPM_SD"] + sample_cols)
        fout.write("\t".join(out_header) + "\n")
        for line in fin:                                           # every mutation row
            parts = line.rstrip("\n").split("\t")
            gene = parts[0]                                        # gene name = first field
            glt = med_map.get(gene, None)                         # look up its median TPM
            sdv = sd_map.get(gene, None)                          # look up its SD
            # If the gene has no expression value, we write "NA" (never 0 -
            # a 0 would falsely claim the gene is silent, which we don't know).
            if glt is None or (isinstance(glt, float) and np.isnan(glt)):
                n_na += 1
            out = ([parts[0], parts[1], parts[2], fmt(glt), fmt(sdv)]
                   + parts[3:])
            fout.write("\t".join(out) + "\n")
            n_rows += 1
    log(f"Mutations with no expression match (GeneLevelTPM=NA): {n_na}")
    log(f"Wrote {OUT}: {n_rows} mutations x {len(sample_cols)} "
        f"sample columns (+GeneLevelTPM, +GeneLevelTPM_SD)")

if __name__ == "__main__":
    sys.exit(main())
