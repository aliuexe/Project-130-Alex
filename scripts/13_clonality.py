#!/usr/bin/env python3
"""
13_clonality.py
Project 130 - Colorectal cancer (TCGA-COAD)  --  practicality extension

Estimates whether each mutation is CLONAL (present in essentially all tumour
cells) or SUBCLONAL (present in only a fraction). Clonal neoantigens are much
better therapeutic targets: a subclonal one lets the tumour escape simply by
losing the minority clone that carries it.

Method: from the MAF read counts we compute the variant allele fraction
        VAF = t_alt_count / t_depth   (fraction of sequencing reads at the site
        that carry the mutant base) for every mutation-in-a-sample. For each
        distinct mutation (gene + protein change) we summarise the VAF across
        the tumours carrying it (median / mean / n) and classify:
          CLONAL     if median VAF >= 0.25
          SUBCLONAL  otherwise
CAVEAT: raw VAF is only a proxy for clonality. A rigorous cancer-cell fraction
        (CCF) would correct for tumour purity and local copy number (e.g. with
        ABSOLUTE); those data are not used here, so the classes are indicative,
        not definitive. A heterozygous clonal mutation in a pure diploid tumour
        gives VAF ~0.5; purity < 1 lowers it, hence the 0.25 cut-off.

Input:  MAF (filtered to protein-coding PASS missense SNV, as elsewhere)
Output: results/mutation_clonality.tsv
        figures/fig19_vaf_clonality.png
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAF = "/tmp/coad.maf" if os.path.exists("/tmp/coad.maf") else \
      os.path.join(BASE, "cohortMAF.2026-07-15.maf.gz")
RES = os.path.join(BASE, "results"); FIG = os.path.join(BASE, "figures")
OUT = os.path.join(RES, "mutation_clonality.tsv")
CLONAL_VAF = 0.25          # documented clonal/subclonal threshold (VAF proxy)

def log(m): print("[13]", m, flush=True)

def main():
    comp = "gzip" if MAF.endswith(".gz") else None
    use = ["Hugo_Symbol", "HGVSc", "HGVSp_Short", "Variant_Classification",
           "Variant_Type", "BIOTYPE", "GDC_FILTER", "t_depth", "t_alt_count"]
    log(f"Reading MAF: {MAF}")
    df = pd.read_csv(MAF, sep="\t", comment="#", usecols=use, dtype=str,
                     compression=comp, low_memory=False)
    # same biological filter as the rest of the pipeline
    keep = ((df["BIOTYPE"] == "protein_coding") &
            (df["GDC_FILTER"].fillna("").isin(["", "PASS"])) &
            (df["Variant_Classification"] == "Missense_Mutation") &
            (df["Variant_Type"] == "SNP"))
    df = df[keep].copy()
    df["t_depth"] = pd.to_numeric(df["t_depth"], errors="coerce")
    df["t_alt_count"] = pd.to_numeric(df["t_alt_count"], errors="coerce")
    df = df[(df["t_depth"] > 0) & df["t_alt_count"].notna()]
    df["VAF"] = df["t_alt_count"] / df["t_depth"]
    log(f"Mutation-sample observations with VAF: {len(df)}")

    # summarise per distinct mutation (gene + protein change)
    g = df.groupby(["Hugo_Symbol", "HGVSp_Short"])["VAF"]
    summ = g.agg(n_samples="count", medianVAF="median", meanVAF="mean").reset_index()
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

    # ---- figure: VAF distribution (per mutation-sample observation) --------
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

    # console preview: clonality of the key driver neoantigens
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
