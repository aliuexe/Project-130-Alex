#!/usr/bin/env python3
"""
03_integrate_datasets.py
Project 130 - Colorectal cancer (TCGA-COAD)

Integrates the mutation-by-sample matrix (Deliverable 01) with the
gene-level expression summary derived from the TPM matrix (Deliverable 02),
producing the integrated mutation+expression matrix (Deliverable 03).

Assignment section implemented:
  Section 7 (Integration of Mutations and Gene Expression):
    GeneLevelTPM = median TPM across tumour samples
    Integrated matrix is tab-delimited with columns:
      GeneName, Mutation, AminoAcidChange, GeneLevelTPM, <sample columns...>

Aggregation method (stated per Section 7):
  For each gene we compute the MEDIAN TPM across all tumour expression
  samples (recommended in the assignment because the median is robust to
  extreme values). Genes present in the mutation matrix but absent from the
  expression matrix receive GeneLevelTPM = NA (reported as NA, never 0).

Note on sample sets (Section 7): the mutation cohort (586 samples) and the
expression cohort (592 samples) are both TCGA-COAD but not identical; we
therefore summarise expression to a single cancer-level GeneLevelTPM value
per gene rather than matching individual samples. The per-sample mutation
columns from Deliverable 01 are retained.

Output: results/03_integrated_mutation_expression.tsv (tab-delimited)
"""
import os
import sys
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MUT = os.path.join(BASE, "results", "01_mutation_by_sample.tsv")
TPM = os.path.join(BASE, "results", "02_gene_by_sample_TPM.tsv")
OUT = os.path.join(BASE, "results", "03_integrated_mutation_expression.tsv")
GLT = os.path.join(BASE, "results", "gene_level_tpm.tsv")

def log(m): print(f"[03] {m}", flush=True)

def main():
    log("Computing GeneLevelTPM from Deliverable 02")
    tpm = pd.read_csv(TPM, sep="\t", low_memory=False)
    tpm_sample_cols = [c for c in tpm.columns if c != "GeneName"]
    tpm_idx = tpm.set_index("GeneName")[tpm_sample_cols]

    # ---- Per-gene expression summary across tumour samples ----------------
    # GeneLevelTPM  = median (robust central value, used for integration)
    # GeneLevelTPM_SD  = standard deviation (dispersion; how trustworthy the
    #                    single median is - a large SD means expression varies
    #                    widely across tumours)
    # GeneLevelTPM_IQR = inter-quartile range (Q3-Q1); robust dispersion that
    #                    pairs naturally with the median
    # n_expr_samples   = number of tumour samples with a non-missing value
    med = tpm_idx.median(axis=1, skipna=True)
    sd = tpm_idx.std(axis=1, skipna=True)
    iqr = (tpm_idx.quantile(0.75, axis=1) - tpm_idx.quantile(0.25, axis=1))
    n_expr = tpm_idx.notna().sum(axis=1)
    gene_level = pd.DataFrame({
        "GeneName": med.index,
        "GeneLevelTPM": med.values.round(4),
        "GeneLevelTPM_SD": sd.values.round(4),
        "GeneLevelTPM_IQR": iqr.values.round(4),
        "n_expr_samples": n_expr.values.astype(int),
    })
    gene_level.to_csv(GLT, sep="\t", index=False)
    log(f"Computed GeneLevelTPM (median + SD + IQR) for {len(gene_level)} genes")

    med_map = dict(zip(med.index, med.round(4)))
    sd_map = dict(zip(sd.index, sd.round(4)))

    # ---- Stream Deliverable 01 line-by-line, insert GeneLevelTPM + SD ------
    # Avoids loading the full 154k x 586 matrix into RAM alongside the merge.
    # Column order: GeneName, Mutation, AminoAcidChange, GeneLevelTPM,
    # GeneLevelTPM_SD, then the sample columns. (Downstream scripts locate the
    # sample columns by their "TCGA" prefix, so inserting GeneLevelTPM_SD here
    # does not shift their parsing.)
    log("Streaming Deliverable 01 and inserting GeneLevelTPM + SD by gene")
    n_rows = 0
    n_na = 0

    def fmt(v):
        return "NA" if v is None or (isinstance(v, float) and np.isnan(v)) \
               else f"{v:.4f}"

    with open(MUT) as fin, open(OUT, "w") as fout:
        header = fin.readline().rstrip("\n").split("\t")
        # header[0:3] = Gene_Name, Mutation, AminoAcid_Change; rest = samples
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
    log(f"Mutations with no expression match (GeneLevelTPM=NA): {n_na}")
    log(f"Wrote {OUT}: {n_rows} mutations x {len(sample_cols)} "
        f"sample columns (+GeneLevelTPM, +GeneLevelTPM_SD)")

if __name__ == "__main__":
    sys.exit(main())
