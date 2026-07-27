#!/usr/bin/env python3
"""
03_integrate_datasets.py
Project 130 - Colorectal Cancer (TCGA-COAD) Neoantigen Discovery Pipeline

===============================================================================
BIOLOGICAL & COMPUTATIONAL PURPOSE
===============================================================================
This script implements Section 7 of the project assignment (Integration of
Mutations and Gene Expression). It merges the binary mutation-by-sample matrix
(Deliverable 01) with the gene-level expression summary derived from the TPM
matrix (Deliverable 02) to produce the integrated dataset (Deliverable 03).

===============================================================================
AGGREGATION METHODOLOGY & DISPERSION METRICS (§7)
===============================================================================
- Central Value (`GeneLevelTPM`): For each gene, we compute the MEDIAN TPM across
  all tumour expression samples. The median is the recommended central metric
  because gene expression distributions in cancer are right-skewed and prone to
  extreme outlier expression spikes.
- Dispersion Metrics (`GeneLevelTPM_SD` & `GeneLevelTPM_IQR`):
  - Standard Deviation (SD): Measures parametric variance across tumours.
  - Interquartile Range (IQR = Q3 - Q1): Non-parametric robust dispersion metric
    that pairs naturally with the median.
- Missing Value Rule: Genes present in the mutation matrix but absent from the
  expression matrix receive `GeneLevelTPM = NA` (reported as `NA`, never `0`).

===============================================================================
SAMPLE COHORT INTEGRATION NOTE (§7)
===============================================================================
The mutation cohort (586 TCGA-COAD tumours) and expression cohort (592 TCGA-COAD
tumours) represent overlapping but distinct sample sets from the same cancer study.
Per Section 7 guidance, expression is summarized to a single cancer-level
`GeneLevelTPM` value per gene, maintaining matched sample-level mutation calls.

===============================================================================
INPUT & OUTPUT CONTRACTS
===============================================================================
Inputs:
  - `results/01_mutation_by_sample.tsv` (Deliverable 01)
  - `results/02_gene_by_sample_TPM.tsv` (Deliverable 02)

Outputs:
  - `results/03_integrated_mutation_expression.tsv` (Deliverable 03)
  - `results/gene_level_tpm.tsv` (Per-gene median + SD + IQR lookup table)
"""

import os
import sys
import numpy as np
import pandas as pd

# =============================================================================
# FILE PATHS & RESOURCE RESOLUTION
# =============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MUT = os.path.join(BASE, "results", "01_mutation_by_sample.tsv")
TPM = os.path.join(BASE, "results", "02_gene_by_sample_TPM.tsv")
OUT = os.path.join(BASE, "results", "03_integrated_mutation_expression.tsv")
GLT = os.path.join(BASE, "results", "gene_level_tpm.tsv")

def log(m):
    """Prints timestamped progress messages to stdout with line flushing."""
    print(f"[03] {m}", flush=True)

def main():
    # =========================================================================
    # STEP 1: COMPUTE PER-GENE EXPRESSION SUMMARY (MEDIAN, SD, IQR)
    # =========================================================================
    log("Computing GeneLevelTPM from Deliverable 02")
    tpm = pd.read_csv(TPM, sep="\t", low_memory=False)
    tpm_sample_cols = [c for c in tpm.columns if c != "GeneName"]
    tpm_idx = tpm.set_index("GeneName")[tpm_sample_cols]

    # Calculate median (central value), SD, and IQR dispersion per gene
    med = tpm_idx.median(axis=1, skipna=True)
    sd = tpm_idx.std(axis=1, skipna=True)
    iqr = (tpm_idx.quantile(0.75, axis=1) - tpm_idx.quantile(0.25, axis=1))
    n_expr = tpm_idx.notna().sum(axis=1)

    # Construct gene-level summary lookup table
    gene_level = pd.DataFrame({
        "GeneName": med.index,
        "GeneLevelTPM": med.values.round(4),
        "GeneLevelTPM_SD": sd.values.round(4),
        "GeneLevelTPM_IQR": iqr.values.round(4),
        "n_expr_samples": n_expr.values.astype(int),
    })
    gene_level.to_csv(GLT, sep="\t", index=False)
    log(f"Computed GeneLevelTPM (median + SD + IQR) for {len(gene_level)} genes")

    # Fast dictionary mapping for streaming lookup
    med_map = dict(zip(med.index, med.round(4)))
    sd_map = dict(zip(sd.index, sd.round(4)))

    # =========================================================================
    # STEP 2: STREAM DELIVERABLE 01 AND INTEGRATE EXPRESSION METRICS
    # =========================================================================
    log("Streaming Deliverable 01 and inserting GeneLevelTPM + SD by gene")
    n_rows = 0
    n_na = 0

    def fmt(v):
        """Format floating point numbers to 4 decimals or 'NA' if missing."""
        return "NA" if v is None or (isinstance(v, float) and np.isnan(v)) \
               else f"{v:.4f}"

    # Stream Deliverable 01 line-by-line to prevent RAM bloat
    with open(MUT) as fin, open(OUT, "w") as fout:
        header = fin.readline().rstrip("\n").split("\t")
        sample_cols = header[3:]
        out_header = (["GeneName", "Mutation", "AminoAcidChange",
                       "GeneLevelTPM", "GeneLevelTPM_SD"] + sample_cols)
        fout.write("\t".join(out_header) + "\n")
        
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            gene = parts[0]
            glt = med_map.get(gene, None)
            sdv = sd_map.get(gene, None)
            if glt is None or (isinstance(glt, float) and np.isnan(glt)):
                n_na += 1
            out = ([parts[0], parts[1], parts[2], fmt(glt), fmt(sdv)]
                   + parts[3:])
            fout.write("\t".join(out) + "\n")
            n_rows += 1

    # =========================================================================
    # STEP 3: LOG INTEGRATION SUMMARY STATISTICS
    # =========================================================================
    log(f"Mutations with no expression match (GeneLevelTPM=NA): {n_na}")
    log(f"Wrote {OUT}: {n_rows} mutations x {len(sample_cols)} "
        f"sample columns (+GeneLevelTPM, +GeneLevelTPM_SD)")

if __name__ == "__main__":
    sys.exit(main())
