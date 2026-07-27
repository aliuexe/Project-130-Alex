#!/usr/bin/env python3
# =============================================================================
# 02_build_expression_matrix.py   (ANNOTATED teaching copy)
# =============================================================================
# WHAT THIS SCRIPT DOES:
#   It turns RNA-sequencing measurements into a tidy grid of "how strongly is
#   each gene switched on" in each tumour, on the standard TPM scale.
#
# WHY WE NEED IT (biology):
#   A mutation only matters for a neoantigen if the gene is actually being READ
#   OUT (transcribed) into RNA and made into protein. RNA-seq measures how much
#   RNA each gene produces = a proxy for how "expressed" (active) the gene is.
#   Later we use this to prefer mutations in genes that are actually expressed.
#
# KEY IDEA — TPM (Transcripts Per Million):
#   Raw RNA-seq numbers can't be compared between samples directly (one sample
#   may have been sequenced more deeply than another). TPM rescales each sample
#   so its values sum to 1,000,000; then a gene's TPM is comparable across
#   samples. Our source file gives "RSEM" values whose per-sample totals are
#   ~18,000,000 (NOT already TPM), so we must convert. The assignment forbids
#   calling non-TPM numbers "TPM", so we do a real, documented conversion.
#
# INPUT : data_mrna_seq_v2_rsem.txt  (genes x samples RSEM expression, cBioPortal)
# OUTPUT: results/02_gene_by_sample_TPM.tsv  and  results/expression_metadata.txt
# =============================================================================

import os
import sys
import numpy as np      # fast maths on arrays of numbers
import pandas as pd     # the table library

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RSEM = os.path.join(BASE, "coadread_tcga_pan_can_atlas_2018",
                    "data_mrna_seq_v2_rsem.txt")               # input file
OUT = os.path.join(BASE, "results", "02_gene_by_sample_TPM.tsv")   # output grid
META = os.path.join(BASE, "results", "expression_metadata.txt")   # required metadata

def log(m): print(f"[02] {m}", flush=True)

def main():
    log(f"Reading RSEM matrix: {RSEM}")
    df = pd.read_csv(RSEM, sep="\t", low_memory=False)   # read the genes-x-samples table
    # Every column except the two identifier columns is a tumour sample:
    sample_cols = [c for c in df.columns
                   if c not in ("Hugo_Symbol", "Entrez_Gene_Id")]
    n_genes_raw, n_samples = df.shape[0], len(sample_cols)
    log(f"Raw RSEM matrix: {n_genes_raw} rows x {n_samples} samples")

    # ---- Remove rows that have no gene name (we can't use them) --------------
    df = df[df["Hugo_Symbol"].notna() & (df["Hugo_Symbol"] != "")].copy()

    # ---- Combine duplicate gene names ----------------------------------------
    # Occasionally the same gene symbol appears on more than one row. Expression
    # (abundance) is additive, so we SUM the rows for each gene, leaving exactly
    # one row per unique gene name.
    n_dup = df["Hugo_Symbol"].duplicated().sum()
    log(f"Duplicated gene symbols collapsed (summed): {n_dup}")
    # groupby = "gather all rows with the same Hugo_Symbol and combine them".
    # min_count=1 keeps a gene as NA (missing) if it had no real values at all.
    expr = df.groupby("Hugo_Symbol", as_index=True)[sample_cols].sum(min_count=1)

    # ---- The actual RSEM -> TPM conversion ------------------------------------
    # For each sample (column), add up all its gene values (col_totals). Then
    # divide every value in that column by its total and multiply by 1,000,000.
    # Result: each sample column now sums to exactly 1e6 = the definition of TPM.
    col_totals = expr.sum(axis=0, skipna=True)       # per-sample grand total
    tpm = expr.divide(col_totals, axis=1) * 1e6      # rescale to per-million

    # Double-check the conversion actually worked (each column must sum to 1e6).
    # np.allclose allows a tiny rounding tolerance. If this fails the script stops.
    check = tpm.sum(axis=0, skipna=True)
    assert np.allclose(check.dropna(), 1e6, rtol=1e-6), \
        f"TPM columns do not sum to 1e6: {check.describe()}"
    log("Verified: all sample columns sum to 1,000,000 (TPM property).")

    # ---- Count missing values (assignment Section 8 asks for this) -----------
    n_cells = tpm.shape[0] * tpm.shape[1]            # total number of cells
    n_missing = int(tpm.isna().sum().sum())          # how many are blank (NA)
    pct_missing = 100.0 * n_missing / n_cells
    log(f"Missing TPM cells: {n_missing}/{n_cells} ({pct_missing:.4f}%)")

    # ---- Save the TPM grid ---------------------------------------------------
    out = tpm.reset_index().rename(columns={"Hugo_Symbol": "GeneName"})
    out = out.round(4)                               # 4 decimals is plenty
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    log(f"Wrote {OUT}: {out.shape[0]} genes x {n_samples} samples")

    # ---- Write the metadata the assignment requires (Section 6) --------------
    # This documents exactly where the data came from and how TPM was produced,
    # so the work is transparent and reproducible.
    with open(META, "w") as fh:
        fh.write("GENE-EXPRESSION DATASET METADATA (Section 6)\n")
        fh.write("Source study: coadread_tcga_pan_can_atlas_2018 "
                 "(TCGA Colorectal, PanCancer Atlas), via cBioPortal.\n")
        fh.write("Source file: data_mrna_seq_v2_rsem.txt "
                 "(mRNA Expression, RSEM, batch-normalized from "
                 "Illumina HiSeq RNASeqV2).\n")
        fh.write("Assay: RNA-seq (Illumina HiSeq RNASeqV2), NOT microarray.\n")
        fh.write("GEO accession: not applicable - matched TCGA RNA-seq used "
                 "(same samples as mutation cohort where overlapping).\n")
        fh.write("Sequencing platform: Illumina HiSeq 2000 (RNASeqV2).\n")
        fh.write(f"Number of tumour samples (expression): {n_samples}\n")
        fh.write("Number of normal samples: provided separately in "
                 "normals/ subfolder; not used for tumour TPM matrix.\n")
        fh.write("Whether TPM was downloaded or calculated: CALCULATED. "
                 "Source values were RSEM normalized estimates "
                 "(per-sample sum ~1.8e7), converted to TPM by rescaling "
                 "each sample to sum to 1e6.\n")
        fh.write("Gene identifier type: HGNC/Hugo gene symbols "
                 "(Entrez IDs also present in source).\n")
        fh.write("Method used to map identifiers to gene symbols: source "
                 "file is already keyed on Hugo_Symbol; no external mapping "
                 "required.\n")
        fh.write("Method used to handle duplicated gene symbols: TPM values "
                 "summed across duplicate rows to one row per unique symbol; "
                 f"{n_dup} duplicate rows collapsed.\n")
        fh.write(f"Genes in final TPM matrix: {out.shape[0]}\n")
        fh.write(f"Percentage missing expression values: {pct_missing:.4f}%\n")
    log(f"Wrote metadata: {META}")

if __name__ == "__main__":
    sys.exit(main())
