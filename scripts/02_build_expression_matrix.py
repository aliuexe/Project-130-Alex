#!/usr/bin/env python3
r"""
02_build_expression_matrix.py
Project 130 - Colorectal Cancer (TCGA-COAD) Neoantigen Discovery Pipeline

===============================================================================
BIOLOGICAL & COMPUTATIONAL PURPOSE
===============================================================================
This script implements Section 6 of the project assignment and satisfies Rule 6.
It parses matched TCGA colorectal RNA-seq data (cBioPortal PanCancer Atlas:
`coadread_tcga_pan_can_atlas_2018 / data_mrna_seq_v2_rsem.txt`) and constructs
the gene-by-sample TPM (Transcripts Per Million) expression matrix.

===============================================================================
MATHEMATICAL FORMULATION — RSEM TO TPM CONVERSION (RULE 6)
===============================================================================
Source RSEM values are batch-normalized transcript abundance estimates from
Illumina HiSeq RNASeqV2 (column sums ~1.8 * 10^7). To satisfy Assignment Rule 6
(values must NOT be labeled as TPM without appropriate per-million conversion),
we convert each sample column to exact Transcripts Per Million (TPM):

    TPM(g, s) = ( RSEM(g, s) / \sum_{g'} RSEM(g', s) ) * 1,000,000

Where:
  - g is a specific gene
  - s is a specific tumour sample
  - \sum_{g'} RSEM(g', s) is the per-sample column sum over all genes

Defining Property Asserted: After conversion, \sum_g TPM(g, s) == 1,000,000.0
for every sample column s (verified with numerical tolerance rtol = 1e-6).

===============================================================================
DUPLICATED GENE SYMBOL HANDLING (§6)
===============================================================================
If a Hugo_Symbol appears on multiple rows (isoform splits in original annotation),
abundance is additive. TPM values are summed across duplicate rows, producing a
unique, non-redundant row per gene symbol. Blank/NA gene symbols are removed.

===============================================================================
INPUT & OUTPUT CONTRACTS
===============================================================================
Input:
  - `coadread_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt`

Outputs:
  - `results/02_gene_by_sample_TPM.tsv` (Deliverable 02 matrix)
  - `results/expression_metadata.txt` (Section 6 required metadata record)
"""

import os
import sys
import numpy as np
import pandas as pd

# =============================================================================
# FILE PATHS & RESOURCE RESOLUTION
# =============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RSEM = os.path.join(BASE, "coadread_tcga_pan_can_atlas_2018",
                    "data_mrna_seq_v2_rsem.txt")
OUT = os.path.join(BASE, "results", "02_gene_by_sample_TPM.tsv")
META = os.path.join(BASE, "results", "expression_metadata.txt")

def log(m):
    """Prints timestamped progress messages to stdout with line flushing."""
    print(f"[02] {m}", flush=True)

def main():
    # =========================================================================
    # STEP 1: PARSE SOURCE RSEM EXPRESSION MATRIX
    # =========================================================================
    log(f"Reading RSEM matrix: {RSEM}")
    df = pd.read_csv(RSEM, sep="\t", low_memory=False)
    
    # Isolate sample columns (excluding Hugo_Symbol and Entrez_Gene_Id)
    sample_cols = [c for c in df.columns
                   if c not in ("Hugo_Symbol", "Entrez_Gene_Id")]
    n_genes_raw, n_samples = df.shape[0], len(sample_cols)
    log(f"Raw RSEM matrix: {n_genes_raw} rows x {n_samples} samples")

    # =========================================================================
    # STEP 2: REMOVE UNANNOTATED ROWS
    # =========================================================================
    # Drop rows without a valid Hugo gene symbol
    df = df[df["Hugo_Symbol"].notna() & (df["Hugo_Symbol"] != "")].copy()

    # =========================================================================
    # STEP 3: COLLAPSE DUPLICATE GENE SYMBOLS (§6)
    # =========================================================================
    # Sum abundance across duplicate rows for identical HGNC/Hugo symbols
    n_dup = df["Hugo_Symbol"].duplicated().sum()
    log(f"Duplicated gene symbols collapsed (summed): {n_dup}")
    expr = df.groupby("Hugo_Symbol", as_index=True)[sample_cols].sum(min_count=1)

    # =========================================================================
    # STEP 4: RSEM TO TPM CONVERSION (RULE 6)
    # =========================================================================
    # Calculate per-sample column sums over all genes
    col_totals = expr.sum(axis=0, skipna=True)
    # Rescale each sample column to sum to exactly 1,000,000 (TPM)
    tpm = expr.divide(col_totals, axis=1) * 1e6

    # Verify mathematical defining property of TPM: column sum == 1,000,000
    check = tpm.sum(axis=0, skipna=True)
    assert np.allclose(check.dropna(), 1e6, rtol=1e-6), \
        f"TPM columns do not sum to 1e6: {check.describe()}"
    log("Verified: all sample columns sum to 1,000,000 (TPM property).")

    # =========================================================================
    # STEP 5: MISSING VALUE BOOKKEEPING (§8)
    # =========================================================================
    n_cells = tpm.shape[0] * tpm.shape[1]
    n_missing = int(tpm.isna().sum().sum())
    pct_missing = 100.0 * n_missing / n_cells
    log(f"Missing TPM cells: {n_missing}/{n_cells} ({pct_missing:.4f}%)")

    # =========================================================================
    # STEP 6: WRITE DELIVERABLE 02 MATRIX TO DISK
    # =========================================================================
    out = tpm.reset_index().rename(columns={"Hugo_Symbol": "GeneName"})
    out = out.round(4)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    log(f"Wrote {OUT}: {out.shape[0]} genes x {n_samples} samples")

    # =========================================================================
    # STEP 7: WRITE REQUIRED METADATA RECORD (§6)
    # =========================================================================
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
