#!/usr/bin/env python3
"""
17_composite_score.py
Project 130 - Colorectal cancer (TCGA-COAD)  --  final ranked neoantigen table (BigMHC)

Combines five NON-REDUNDANT axes:
  1. Presentation   mutant peptide is presented   (BigMHC_EL, higher = better)
  2. Immunogenicity T cell can recognise it       (BigMHC_IM / Calis score, higher = better)
  3. Expression     protein is actually made      (GeneLevelTPM, higher = better)
  4. Clonality      present in all tumour cells    (median VAF, higher = better)
  5. Recurrence     shared across patients         (tumours covered, higher = better)

Input:  results/practical_neoantigens_scored.tsv, results/mutation_clonality.tsv
Output: results/neoantigen_ranked_final.tsv
        figures/fig26_composite_top_candidates.png
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results"); FIG = os.path.join(BASE, "figures")
SCORED = os.path.join(RES, "practical_neoantigens_scored.tsv")
CLON = os.path.join(RES, "mutation_clonality.tsv")
OUT = os.path.join(RES, "neoantigen_ranked_final.tsv")

AXES = ["Presentation", "Immunogenicity", "Expression", "Clonality", "Recurrence"]
COLORS = {"Presentation":"#4477AA","Immunogenicity":"#EE6677","Expression":"#228833",
          "Clonality":"#EE7733","Recurrence":"#AA3377"}

def log(m): print("[17]", m, flush=True)

def main():
    df = pd.read_csv(SCORED, sep="\t", comment="#")
    clon = pd.read_csv(CLON, sep="\t")[["GeneName","ProteinChange","medianVAF"]]
    df = df.merge(clon, on=["GeneName","ProteinChange"], how="left")
    n0 = len(df)

    need = ["Mutant_EL","ImmunogenicityScore","GeneLevelTPM","medianVAF","TumoursCovered"]
    for c in need:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=need).reset_index(drop=True)
    log(f"candidates scored with BigMHC: {len(df)} (of {n0})")

    # Percentile ranks (higher = better)
    df["pct_Presentation"]   = df["Mutant_EL"].rank(pct=True)      # higher BigMHC_EL -> higher pct
    df["pct_Immunogenicity"] = df["ImmunogenicityScore"].rank(pct=True)
    df["pct_Expression"]     = df["GeneLevelTPM"].rank(pct=True)
    df["pct_Clonality"]      = df["medianVAF"].rank(pct=True)
    df["pct_Recurrence"]     = (np.log10(1 + df["TumoursCovered"])
                                / np.log10(1 + df["TumoursCovered"].max()))
    pcts = [f"pct_{a}" for a in AXES]

    df["CompositeScore"] = df[pcts].mean(axis=1).round(4)
    df = df.sort_values("CompositeScore", ascending=False).reset_index(drop=True)
    df["Rank"] = np.arange(1, len(df)+1)

    keep = (["Rank","GeneName","ProteinChange","Peptide","HLAAllele",
             "Mutant_EL","WT_EL","ImmunogenicityScore","MutationAtAnchor",
             "GeneLevelTPM","medianVAF","MutationFrequency","TumoursCovered",
             "CompositeScore"] + pcts)
    df[keep].to_csv(OUT, sep="\t", index=False)
    log(f"wrote {OUT}")

    log("Top 12 by BigMHC composite score:")
    for _,r in df.head(12).iterrows():
        print(f"  #{int(r.Rank):2d} {r.GeneName:8s}{r.ProteinChange:9s} {r.Peptide} "
              f"{r.HLAAllele:12s} comp={r.CompositeScore:.3f}  "
              f"(pres {r.pct_Presentation:.2f} imm {r.pct_Immunogenicity:.2f} "
              f"expr {r.pct_Expression:.2f} clon {r.pct_Clonality:.2f} rec {r.pct_Recurrence:.2f})")

    # Stacked bar plot for top 15
    top = df.head(15).iloc[::-1]
    fig, ax = plt.subplots(figsize=(11.5, 7))
    left = np.zeros(len(top))
    for a in AXES:
        contrib = top[f"pct_{a}"].values * (1/len(AXES))
        ax.barh(range(len(top)), contrib, left=left, color=COLORS[a], label=a, edgecolor="white")
        left += contrib
    labels = [f"{r.GeneName} {r.ProteinChange}" for _,r in top.iterrows()]
    ax.set_yticks(range(len(top))); ax.set_yticklabels(labels, fontsize=9)
    for i,(_,r) in enumerate(top.iterrows()):
        ax.text(left[i]+0.005, i, f"{r.CompositeScore:.2f} · cov {int(r.TumoursCovered)}",
                va="center", fontsize=8)
    ax.set_xlabel("Composite neoantigen quality score (BigMHC presentation + immunogenicity + expression + clonality + recurrence)")
    ax.set_title("Top-ranked practical neoantigens (TCGA-COAD, BigMHC)")
    ax.legend(ncol=5, frameon=False, loc="lower right", fontsize=9)
    ax.set_xlim(0,1.05); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(FIG,"fig26_composite_top_candidates.png"),dpi=160)
    plt.close(fig); log("wrote fig26_composite_top_candidates.png")

if __name__ == "__main__":
    sys.exit(main())
