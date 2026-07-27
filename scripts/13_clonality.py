#!/usr/bin/env python3
r"""
13_clonality.py
Project 130 - Colorectal Cancer (TCGA-COAD) Neoantigen Discovery Pipeline

===============================================================================
BIOLOGICAL & COMPUTATIONAL PURPOSE
===============================================================================
This script estimates the clonality architecture of somatic missense SNVs across
the TCGA-COAD cohort to differentiate CLONAL (truncal) from SUBCLONAL (branched)
mutations.

Motivation: Clonal neoantigens are present in 100% of tumour cells and represent
optimal therapeutic targets for cancer vaccines. Subclonal neoantigens permit
immune escape through selective outgrowth of antigen-negative sub-clones.

===============================================================================
MATHEMATICAL FORMULATION — VARIANT ALLELE FRACTION (VAF)
===============================================================================
For each variant observation in a tumour sample:
    VAF = t_alt_count / t_depth

Where:
  - `t_alt_count`: Number of sequencing reads covering the mutated allele.
  - `t_depth`: Total sequencing depth at the genomic locus.

Clonality Classification Rule (§13):
For each distinct mutation `(GeneName, ProteinChange)`, the median VAF is computed
across all tumours carrying that mutation:
  - `CLONAL`: $\text{median VAF} \ge 0.25$
  - `SUBCLONAL`: $\text{median VAF} < 0.25$

Biological Rationale for Threshold 0.25:
For a heterozygous mutation in a 100% pure diploid tumour, expected $\text{VAF} = 0.50$.
Accounting for normal cell contamination (tumour purity $\sim 50\text{--}60\%$),
$\text{VAF} \ge 0.25$ provides a robust surrogate for truncal clonality.

===============================================================================
INPUT & OUTPUT CONTRACTS
===============================================================================
Input:
  - `cohortMAF.2026-07-15.maf.gz` (with `t_depth` and `t_alt_count`)

Outputs:
  - `results/mutation_clonality.tsv` (Clonality classification summary)
  - `figures/fig19_vaf_clonality.png` (Cohort VAF distribution histogram)
"""

import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =============================================================================
# FILE PATHS & RESOURCE RESOLUTION
# =============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAF = "/tmp/coad.maf" if os.path.exists("/tmp/coad.maf") else \
      os.path.join(BASE, "cohortMAF.2026-07-15.maf.gz")
RES = os.path.join(BASE, "results"); FIG = os.path.join(BASE, "figures")
OUT = os.path.join(RES, "mutation_clonality.tsv")

# Documented Clonal/Subclonal VAF Cutoff Threshold
CLONAL_VAF = 0.25

def log(m):
    """Prints timestamped progress messages to stdout with line flushing."""
    print("[13]", m, flush=True)

# =============================================================================
# MAIN PIPELINE EXECUTION
# =============================================================================
def main():
    # =========================================================================
    # STEP 1: PARSE MAF READ DEPTHS & CALCULATE VAF
    # =========================================================================
    comp = "gzip" if MAF.endswith(".gz") else None
    use = ["Hugo_Symbol", "HGVSc", "HGVSp_Short", "Variant_Classification",
           "Variant_Type", "BIOTYPE", "GDC_FILTER", "t_depth", "t_alt_count"]
    log(f"Reading MAF: {MAF}")
    df = pd.read_csv(MAF, sep="\t", comment="#", usecols=use, dtype=str,
                     compression=comp, low_memory=False)
    
    # Filter for protein-coding PASS missense SNVs
    keep = ((df["BIOTYPE"] == "protein_coding") &
            (df["GDC_FILTER"].fillna("").isin(["", "PASS"])) &
            (df["Variant_Classification"] == "Missense_Mutation") &
            (df["Variant_Type"] == "SNP"))
    df = df[keep].copy()
    
    # Cast depth and alt count columns to numeric
    df["t_depth"] = pd.to_numeric(df["t_depth"], errors="coerce")
    df["t_alt_count"] = pd.to_numeric(df["t_alt_count"], errors="coerce")
    df = df[(df["t_depth"] > 0) & df["t_alt_count"].notna()]
    df["VAF"] = df["t_alt_count"] / df["t_depth"]
    log(f"Mutation-sample observations with VAF: {len(df)}")

    # =========================================================================
    # STEP 2: SUMMARIZE CLONALITY PER DISTINCT MUTATION
    # =========================================================================
    g = df.groupby(["Hugo_Symbol", "HGVSp_Short"])["VAF"]
    summ = g.agg(n_samples="count", medianVAF="median", meanVAF="mean").reset_index()
    
    # Apply VAF threshold rule (>= 0.25 -> Clonal)
    summ["ClonalClass"] = np.where(summ["medianVAF"] >= CLONAL_VAF,
                                   "Clonal", "Subclonal")
    summ = summ.rename(columns={"Hugo_Symbol": "GeneName",
                                "HGVSp_Short": "ProteinChange"})
    summ["medianVAF"] = summ["medianVAF"].round(4)
    summ["meanVAF"] = summ["meanVAF"].round(4)
    summ = summ.sort_values(["n_samples", "medianVAF"], ascending=False)
    summ.to_csv(OUT, sep="\t", index=False)
    
    n_cl = int((summ["ClonalClass"] == "Clonal").sum())
    log(f"Distinct mutations: {len(summ)}; clonal: {n_cl} "
        f"({100*n_cl/len(summ):.1f}%)")

    # =========================================================================
    # STEP 3: RENDER FIGURE 19 — VAF DISTRIBUTION HISTOGRAM
    # =========================================================================
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.hist(df["VAF"], bins=60, color="#4477AA", edgecolor="white")
    ax.axvline(CLONAL_VAF, color="#EE6677", ls="--", lw=2,
               label=f"clonal/subclonal cut-off (VAF={CLONAL_VAF})")
    ax.set_xlabel("Variant allele fraction (VAF = mutant reads / total reads)")
    ax.set_ylabel("Mutation observations")
    ax.set_title("Tumour VAF distribution — clonality proxy (TCGA-COAD)")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig19_vaf_clonality.png"), dpi=160)
    plt.close(fig)
    log("wrote fig19_vaf_clonality.png")

    # Console preview logging for key driver mutations
    log("Clonality of notable driver mutations:")
    for gene, pc in [("KRAS","p.G12D"),("KRAS","p.G12V"),("TP53","p.R248Q"),
                     ("PIK3CA","p.E542K"),("BRAF","p.V640E"),("SMAD4","p.R361H")]:
        r = summ[(summ.GeneName==gene) & (summ.ProteinChange==pc)]
        if len(r):
            x = r.iloc[0]
            print(f"    {gene:7s}{pc:9s} n={int(x.n_samples):3d}  "
                  f"medianVAF={x.medianVAF:.2f}  {x.ClonalClass}")

if __name__ == "__main__":
    sys.exit(main())
