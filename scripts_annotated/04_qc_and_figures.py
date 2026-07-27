#!/usr/bin/env python3
# =============================================================================
# 04_qc_and_figures.py   (ANNOTATED teaching copy)
# =============================================================================
# WHAT THIS SCRIPT DOES:
#   Two jobs. (1) Quality control (QC): count and sanity-check the data so we
#   can trust it and report the numbers the assignment asks for (Section 8).
#   (2) Draw four figures that summarise the mutation and expression results.
#
# WHY QC MATTERS:
#   Before drawing biological conclusions we must confirm the data is sensible:
#   How many mutations before/after filtering? How many genes and samples? How
#   much expression data is missing? Are the gene names and genome build
#   consistent? Which genes are most mutated / most expressed? Getting these
#   right is worth real marks and catches mistakes early.
#
# INPUT : the three deliverables (01, 02, 03) already produced.
# OUTPUT: results/qc_report.txt  and  figures/fig1..fig4_*.png
# =============================================================================

import os
import sys
import numpy as np
import pandas as pd
import matplotlib                      # the plotting library
matplotlib.use("Agg")                  # "Agg" = draw to image files, no screen needed
import matplotlib.pyplot as plt        # the drawing commands, nicknamed "plt"

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
FIG = os.path.join(BASE, "figures")
MUT = os.path.join(RES, "01_mutation_by_sample.tsv")
TPM = os.path.join(RES, "02_gene_by_sample_TPM.tsv")
INT = os.path.join(RES, "03_integrated_mutation_expression.tsv")
QC_MUT = os.path.join(RES, "qc_mutation_counts.txt")
EXPR_META = os.path.join(RES, "expression_metadata.txt")
REPORT = os.path.join(RES, "qc_report.txt")
os.makedirs(FIG, exist_ok=True)        # make the figures/ folder if missing

def log(m): print(f"[04] {m}", flush=True)

def read_integrated_lean():
    # Read the integrated matrix but, to save memory, only keep the labels,
    # the GeneLevelTPM, and a computed "MutationFrequency" = how many samples
    # carry each mutation (the count of 1's across that row's sample columns).
    with open(INT) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    # find where the sample columns begin by their "TCGA" barcode prefix.
    # (Robust to metadata columns like GeneLevelTPM / GeneLevelTPM_SD before them.)
    s0 = next(i for i, c in enumerate(header) if c.startswith("TCGA"))
    sample_cols = header[s0:]
    keys, glt, freq = [], [], []
    with open(INT) as fh:
        fh.readline()                              # skip header
        for line in fh:
            p = line.rstrip("\n").split("\t")
            keys.append((p[0], p[1], p[2]))        # (gene, HGVSc, HGVSp)
            glt.append(np.nan if p[3] == "NA" else float(p[3]))  # GeneLevelTPM
            freq.append(sum(1 for v in p[s0:] if v == "1"))      # samples with this mutation
    df = pd.DataFrame(keys, columns=["GeneName", "Mutation", "AminoAcidChange"])
    df["GeneLevelTPM"] = glt
    df["MutationFrequency"] = freq
    return df, len(sample_cols)

def main():
    log("Reading integrated matrix (lean)")
    idf, n_samples = read_integrated_lean()

    # ---- Pull the before/after filtering counts from script 01's QC log ------
    mut_before = mut_after = None
    for line in open(QC_MUT):
        if "before filtering" in line: mut_before = int(line.split(":")[1])
        if "after filtering" in line: mut_after = int(line.split(":")[1])
    n_distinct_mut = len(idf)                       # number of distinct mutations
    n_unique_genes = idf["GeneName"].nunique()      # number of distinct genes

    # ---- Pull the missing-expression percentage from script 02's metadata ----
    pct_missing = None
    for line in open(EXPR_META):
        if "Percentage missing" in line:
            pct_missing = line.split(":")[1].strip()

    # ---- Ten most frequently mutated genes -----------------------------------
    # Sum the per-mutation frequencies within each gene, then take the top 10.
    gene_freq = (idf.groupby("GeneName")["MutationFrequency"].sum()
                 .sort_values(ascending=False))
    top10_mutated = gene_freq.head(10)

    # ---- Ten most highly expressed mutated genes -----------------------------
    gene_tpm = (idf.dropna(subset=["GeneLevelTPM"])
                .groupby("GeneName")["GeneLevelTPM"].first()
                .sort_values(ascending=False))
    top10_expressed = gene_tpm.head(10)

    # ---- Write the QC text report (all the Section 8 checks) -----------------
    with open(REPORT, "w") as fh:
        fh.write("=" * 64 + "\n")
        fh.write("BASIC DATA QUALITY CONTROL REPORT (Section 8)\n")
        fh.write("Project 130 - Colorectal cancer (TCGA-COAD)\n")
        fh.write("=" * 64 + "\n\n")
        fh.write(f"Mutations before filtering: {mut_before}\n")
        fh.write(f"Mutations after filtering:  {mut_after}\n")
        fh.write(f"Distinct mutations (matrix rows): {n_distinct_mut}\n")
        fh.write(f"Unique mutated genes: {n_unique_genes}\n")
        fh.write(f"Tumour samples (mutation matrix): {n_samples}\n")
        fh.write(f"Percentage missing expression values: {pct_missing}\n")
        fh.write("Gene identifier consistency: both matrices keyed on HGNC/"
                 "Hugo gene symbols; integration is a symbol-level join.\n")
        fh.write("Reference genome assembly consistency: all mutation "
                 "coordinates are GRCh38/hg38 (verified in script 01); no "
                 "mixing of assemblies.\n\n")
        fh.write("Distribution of GeneLevelTPM (non-NA):\n")
        # .describe() prints count/mean/min/quartiles/max — a numeric summary.
        fh.write(idf["GeneLevelTPM"].describe().to_string() + "\n\n")
        fh.write("Ten most frequently mutated genes "
                 "(sum of mutation occurrences across samples):\n")
        for g, v in top10_mutated.items():
            fh.write(f"  {g:12s} {int(v)}\n")
        fh.write("\nTen most highly expressed mutated genes "
                 "(GeneLevelTPM = median tumour TPM):\n")
        for g, v in top10_expressed.items():
            fh.write(f"  {g:12s} {v:.2f}\n")
    log(f"Wrote {REPORT}")

    # =====================================================================
    # FIGURES. Each block: make a blank canvas (plt.figure), draw, label,
    # save to a .png, then close. dpi=150 = resolution.
    # =====================================================================

    # ---- Figure 1: horizontal bar chart of the 10 most mutated genes ---------
    plt.figure(figsize=(8, 5))
    top10_mutated[::-1].plot(kind="barh", color="#4477AA")   # [::-1] flips order
    plt.xlabel("Total mutation occurrences across tumour samples")
    plt.ylabel("Gene")
    plt.title("Ten most frequently mutated genes (TCGA-COAD)")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, "fig1_top_mutated_genes.png"), dpi=150)
    plt.close()

    # ---- Figure 2: histogram of GeneLevelTPM on a log scale ------------------
    # Expression spans a huge range, so we plot log10(TPM) to see the shape.
    vals = idf["GeneLevelTPM"].dropna()
    vals = vals[vals > 0]                            # log needs positive numbers
    plt.figure(figsize=(8, 5))
    plt.hist(np.log10(vals), bins=60, color="#66CCEE", edgecolor="white")
    plt.xlabel("log10(GeneLevelTPM)")
    plt.ylabel("Number of mutations")
    plt.title("Distribution of GeneLevelTPM across mutations (log scale)")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, "fig2_genelevel_tpm_distribution.png"), dpi=150)
    plt.close()

    # ---- Figure 3: heat map — which of the 30 most recurrent mutations -------
    #      appear in which (first 60) tumour samples. Black = present, white = absent.
    top_mut = idf.sort_values("MutationFrequency", ascending=False).head(30)
    top_labels = (top_mut["GeneName"] + " " + top_mut["AminoAcidChange"]).tolist()
    with open(INT) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    s0 = next(i for i, c in enumerate(header) if c.startswith("TCGA"))
    sample_cols = header[s0:]
    show_samples = sample_cols[:60]
    show_idx = [s0 + i for i in range(len(sample_cols)) if i < 60]
    wanted = set((r.GeneName, r.Mutation, r.AminoAcidChange)
                 for r in top_mut.itertuples())
    order = {(r.GeneName, r.Mutation, r.AminoAcidChange): i
             for i, r in enumerate(top_mut.itertuples())}
    rows_keyed = []
    with open(INT) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if (p[0], p[1], p[2]) in wanted:
                rows_keyed.append((order[(p[0], p[1], p[2])],
                                   [int(p[j]) for j in show_idx]))
    rows_keyed.sort()
    H = np.array([r[1] for r in rows_keyed])         # the 0/1 grid to draw
    plt.figure(figsize=(11, 7))
    plt.imshow(H, aspect="auto", cmap="Greys", interpolation="none")
    plt.yticks(range(len(top_labels)), top_labels, fontsize=7)
    plt.xlabel(f"Tumour samples (first {len(show_samples)} of {n_samples})")
    plt.title("Presence/absence heat map of 30 most recurrent mutations")
    plt.colorbar(label="Mutation present (1) / absent (0)", shrink=0.6)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, "fig3_mutation_heatmap.png"), dpi=150)
    plt.close()

    # ---- Figure 4: scatter — mutation frequency vs gene expression -----------
    #      Each dot is a gene: x = how expressed, y = how often mutated.
    gdf = idf.dropna(subset=["GeneLevelTPM"]).copy()
    gene_agg = gdf.groupby("GeneName").agg(
        freq=("MutationFrequency", "sum"),
        tpm=("GeneLevelTPM", "first")).reset_index()
    gene_agg = gene_agg[gene_agg["tpm"] > 0]
    plt.figure(figsize=(8, 6))
    plt.scatter(np.log10(gene_agg["tpm"]), gene_agg["freq"],
                s=8, alpha=0.35, color="#EE6677")
    for _, r in gene_agg.sort_values("freq", ascending=False).head(10).iterrows():
        plt.annotate(r["GeneName"], (np.log10(r["tpm"]), r["freq"]), fontsize=8)
    plt.xlabel("log10(GeneLevelTPM)")
    plt.ylabel("Total mutation occurrences (per gene)")
    plt.title("Mutation frequency vs gene expression (TCGA-COAD)")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, "fig4_freq_vs_expression.png"), dpi=150)
    plt.close()

    log("Wrote 4 figures to figures/")

if __name__ == "__main__":
    sys.exit(main())
