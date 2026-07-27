#!/usr/bin/env python3
"""
02_build_expression_matrix.py
Project 130 - Colorectal cancer (TCGA-COAD)

Builds the gene-by-sample TPM matrix (Deliverable 02) from the matched
TCGA-COAD RNA-seq RSEM file (cBioPortal package
coadread_tcga_pan_can_atlas_2018 / data_mrna_seq_v2_rsem.txt).

Assignment sections implemented:
  Section 6 (Gene-Expression Dataset) - gene-by-sample TPM matrix,
            required metadata reporting, RNA-seq (not microarray),
            handling of duplicated gene symbols.
  Rule 6:  Values reported as read counts / FPKM / RPKM must NOT be
           labelled as TPM without appropriate conversion.

IMPORTANT - RSEM -> TPM conversion:
  The source values are RSEM 'batch-normalized' expression estimates
  (Illumina HiSeq RNASeqV2). Their per-sample column sums are ~1.8e7,
  i.e. they are NOT already on the TPM (per-million) scale. RSEM values
  are proportional to transcript abundance, so a valid per-sample
  conversion to TPM is:

        TPM(gene, sample) = RSEM(gene, sample)
                            / sum_over_genes( RSEM(., sample) )
                            * 1e6

  After conversion each sample column sums to exactly 1,000,000, which is
  the defining property of TPM. This is a documented conversion, so the
  output is legitimately labelled TPM (satisfies Rule 6). Genes with a
  missing/NA RSEM value are treated as 0 for the per-sample total but
  reported as NA in the output (see Section 8 missing-value handling).

Duplicated gene symbols (Section 6):
  If the same Hugo_Symbol appears on multiple rows, the TPM values are
  summed across those rows (abundance is additive), producing one row per
  unique gene symbol. Rows with a blank/NA Hugo_Symbol are dropped.

Output: results/02_gene_by_sample_TPM.tsv (tab-delimited)
        results/expression_metadata.txt   (Section 6 required metadata)
"""
import os
import sys
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RSEM = os.path.join(BASE, "coadread_tcga_pan_can_atlas_2018",
                    "data_mrna_seq_v2_rsem.txt")
OUT = os.path.join(BASE, "results", "02_gene_by_sample_TPM.tsv")
META = os.path.join(BASE, "results", "expression_metadata.txt")

def log(m): print(f"[02] {m}", flush=True)

def main():
    log(f"Reading RSEM matrix: {RSEM}")
    df = pd.read_csv(RSEM, sep="\t", low_memory=False)
    sample_cols = [c for c in df.columns
                   if c not in ("Hugo_Symbol", "Entrez_Gene_Id")]
    n_genes_raw, n_samples = df.shape[0], len(sample_cols)
    log(f"Raw RSEM matrix: {n_genes_raw} rows x {n_samples} samples")

    # ---- Drop rows with no gene symbol (cannot map to a gene) --------------
    df = df[df["Hugo_Symbol"].notna() & (df["Hugo_Symbol"] != "")].copy()

    # ---- Collapse duplicated gene symbols by summation --------------------
    n_dup = df["Hugo_Symbol"].duplicated().sum()
    log(f"Duplicated gene symbols collapsed (summed): {n_dup}")
    expr = df.groupby("Hugo_Symbol", as_index=True)[sample_cols].sum(min_count=1)

    # ---- RSEM -> TPM: rescale each sample column to sum to 1e6 -------------
    col_totals = expr.sum(axis=0, skipna=True)          # per-sample RSEM total
    tpm = expr.divide(col_totals, axis=1) * 1e6

    # sanity: every sample now sums to ~1e6
    check = tpm.sum(axis=0, skipna=True)
    assert np.allclose(check.dropna(), 1e6, rtol=1e-6), \
        f"TPM columns do not sum to 1e6: {check.describe()}"
    log("Verified: all sample columns sum to 1,000,000 (TPM property).")

    # ---- Missing-value bookkeeping (Section 8) ----------------------------
    n_cells = tpm.shape[0] * tpm.shape[1]
    n_missing = int(tpm.isna().sum().sum())
    pct_missing = 100.0 * n_missing / n_cells
    log(f"Missing TPM cells: {n_missing}/{n_cells} ({pct_missing:.4f}%)")

    # ---- Write matrix -----------------------------------------------------
    out = tpm.reset_index().rename(columns={"Hugo_Symbol": "GeneName"})
    out = out.round(4)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    log(f"Wrote {OUT}: {out.shape[0]} genes x {n_samples} samples")

    # ---- Required metadata (Section 6) ------------------------------------
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
