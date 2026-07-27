#!/usr/bin/env python3
"""
09_make_figures.py
Project 130 - Colorectal Cancer (TCGA-COAD) Neoantigen Discovery Pipeline

===============================================================================
BIOLOGICAL & COMPUTATIONAL PURPOSE
===============================================================================
This script renders publication-quality figures (Figures 5 through 13) for the
oral presentation and final project report. It reads pre-aggregated summary
data from `results/figure_summary.json` (Script 08) and `results/neoantigen_candidates_shortlist.tsv`
(Script 07), enabling fast rendering without re-parsing multi-gigabyte matrices.

===============================================================================
FIGURE INDEX & VISUALIZATIONS PRODUCED
===============================================================================
  1. Figure 5 (`fig5_methods_funnel.png`): Quantitative analysis funnel from 310k
     raw MAF records down to 1,536 candidate neoantigens.
  2. Figure 6 (`fig6_variant_classification.png`): MAF variant class distribution
     highlighting missense SNVs.
  3. Figure 7 (`fig7_mut_vs_wt_affinity.png`): BigMHC presentation probability
     distribution comparison (Mutant vs Wild-Type).
  4. Figure 8 (`fig8_delta_affinity.png`): DeltaPresentation distribution (`Mut_EL - WT_EL`).
  5. Figure 9 (`fig9_binder_class_mut_vs_wt.png`): Presenter class counts (Strong, Weak, Non-presenter).
  6. Figure 10 (`fig10_strong_binders_by_allele.png`): Strong presenters per HLA allele.
  7. Figure 11 (`fig11_mut_vs_wt_scatter.png`): Scatter plot of Mutant vs WT BigMHC_EL.
  8. Figure 12 (`fig12_top_candidates.png`): Top 15 shortlisted neoantigen candidates.
  9. Figure 13 (`fig13_workflow_schematic.png`): 8-step pipeline workflow diagram.

===============================================================================
INPUT & OUTPUT CONTRACTS
===============================================================================
Inputs:
  - `results/figure_summary.json`
  - `results/neoantigen_candidates_shortlist.tsv`

Outputs:
  - `figures/fig5_methods_funnel.png`
  - `figures/fig6_variant_classification.png`
  - `figures/fig7_mut_vs_wt_affinity.png`
  - `figures/fig8_delta_affinity.png`
  - `figures/fig9_binder_class_mut_vs_wt.png`
  - `figures/fig10_strong_binders_by_allele.png`
  - `figures/fig11_mut_vs_wt_scatter.png`
  - `figures/fig12_top_candidates.png`
  - `figures/fig13_workflow_schematic.png`
"""

import json
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# =============================================================================
# FILE PATHS & RESOURCE RESOLUTION
# =============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
FIG = os.path.join(BASE, "figures")
SUM = os.path.join(RES, "figure_summary.json")
SHORT = os.path.join(RES, "neoantigen_candidates_shortlist.tsv")
os.makedirs(FIG, exist_ok=True)

# =============================================================================
# CONSISTENT FIGURE STYLING & COLOR PALETTE
# =============================================================================
plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 160,
    "font.size": 13, "axes.titlesize": 16, "axes.titleweight": "bold",
    "axes.labelsize": 13, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "figure.facecolor": "white",
    "axes.facecolor": "white",
})
BLUE, RED, TEAL = "#4477AA", "#EE6677", "#66CCEE"
GREEN, ORANGE, GREY, PURPLE = "#228833", "#EE7733", "#8899AA", "#AA3377"

def save(fig, name):
    """Saves figure with tight bounding box and closes canvas."""
    path = os.path.join(FIG, name)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"[09] wrote {name}", flush=True)

def log(m):
    """Prints timestamped progress messages to stdout with line flushing."""
    print(f"[09] {m}", flush=True)

# =============================================================================
# INDIVIDUAL FIGURE RENDERING FUNCTIONS
# =============================================================================
def fig_funnel():
    """Figure 5: Quantitative analysis funnel."""
    stages = [
        ("Raw MAF records", 310472),
        ("Filtered missense SNVs\n(protein-coding, PASS)", 184574),
        ("Distinct mutations", 153996),
        ("Protein-annotated\nmutations", 145612),
        ("Mutant + WT peptides", 6873140),
        ("Unique 9-mers scored\n(x3 HLA alleles)", 2460296),
        ("Strongly presented mutant\n9-mer-allele pairs", 30260),
        ("Prioritised candidates", 1536),
    ]
    labels = [s[0] for s in stages]
    vals = [s[1] for s in stages]
    fig, ax = plt.subplots(figsize=(10, 7))
    cmap = plt.cm.viridis(np.linspace(0.1, 0.85, len(stages)))
    maxv = max(vals)
    for i, (lab, v) in enumerate(stages):
        w = v / maxv
        ax.barh(i, w, color=cmap[i], edgecolor="white", height=0.72)
        ax.text(w + 0.01, i, f"{v:,}", va="center", ha="left",
                fontsize=11, fontweight="bold")
    ax.set_yticks(range(len(stages)))
    ax.set_yticklabels(labels, fontsize=11)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.18)
    ax.set_xticks([])
    ax.spines["bottom"].set_visible(False)
    ax.grid(False)
    ax.set_title("Analysis funnel: from raw mutations to neoantigen candidates\nTCGA-COAD (BigMHC)")
    save(fig, "fig5_methods_funnel.png")

def fig_variant_class():
    """Figure 6: Variant classification bar chart."""
    data = [
        ("Missense", 185538), ("Silent", 65426), ("Frameshift Del", 20070),
        ("Nonsense", 16613), ("Frameshift Ins", 5577), ("Intron", 3690),
        ("Splice Site", 3620), ("3'UTR", 2698), ("Splice Region", 2078),
        ("RNA", 1876), ("In-frame Del", 1163), ("5'UTR", 707),
        ("Other", 479 + 453 + 206 + 164 + 108 + 6),
    ]
    data.sort(key=lambda x: x[1])
    labels = [d[0] for d in data]; vals = [d[1] for d in data]
    colors = [RED if l == "Missense" else BLUE for l in labels]
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.barh(range(len(vals)), vals, color=colors, edgecolor="white")
    ax.set_yticks(range(len(vals))); ax.set_yticklabels(labels)
    for i, v in enumerate(vals):
        ax.text(v + 1500, i, f"{v:,}", va="center", fontsize=9)
    ax.set_xlabel("Number of variant records (raw MAF)")
    ax.set_title("Somatic variant classification (TCGA-COAD)\nMissense SNVs (red) drive the neoantigen analysis")
    ax.set_xlim(0, max(vals) * 1.15)
    save(fig, "fig6_variant_classification.png")

def fig_mut_vs_wt_affinity(S):
    """Figure 7: Presentation probability histogram (Mutant vs Wild-Type)."""
    NB, LO, HI = S["NB"], S["LO"], S["HI"]
    edges = np.linspace(LO, HI, NB + 1)
    centers = (edges[:-1] + edges[1:]) / 2
    mut = np.array(S["el_hist"]["Mutant"], float)
    wt = np.array(S["el_hist"]["WildType"], float)
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.bar(centers, wt, width=(HI-LO)/NB, color=GREY, alpha=0.6,
           label="Wild-type", edgecolor="none")
    ax.bar(centers, mut, width=(HI-LO)/NB, color=RED, alpha=0.55,
           label="Mutant", edgecolor="none")
    ax.axvline(0.50, color=GREEN, ls="--", lw=2,
               label="0.50 BigMHC_EL (Weak presenter cutoff)")
    ax.axvline(0.70, color=PURPLE, ls="--", lw=2,
               label="0.70 BigMHC_EL (Strong presenter cutoff)")
    ax.set_xlabel("BigMHC Eluted-Ligand Presentation Probability  BigMHC_EL")
    ax.set_ylabel("Number of 9-mer-allele predictions")
    ax.set_title("HLA class I presentation: mutant vs wild-type 9-mers\n(BigMHC Neural Predictor)")
    ax.legend(frameon=False, fontsize=11)
    save(fig, "fig7_mut_vs_wt_affinity.png")

def fig_delta(S):
    """Figure 8: DeltaPresentation distribution."""
    DNB, DLO, DHI = S["DNB"], S["DLO"], S["DHI"]
    edges = np.linspace(DLO, DHI, DNB + 1)
    centers = (edges[:-1] + edges[1:]) / 2
    h = np.array(S["delta_hist"], float)
    colors = [GREEN if c > 0 else RED for c in centers]
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.bar(centers, h, width=(DHI-DLO)/DNB, color=colors, alpha=0.8,
           edgecolor="none")
    ax.axvline(0, color="black", lw=1.2)
    pos, neg = S["delta_pos"], S["delta_neg"]
    tot = pos + neg if (pos + neg) > 0 else 1
    ax.text(0.97, 0.95,
            f"Mutant presented STRONGER: {pos:,} ({100*pos/tot:.0f}%)\n"
            f"Mutant presented weaker: {neg:,} ({100*neg/tot:.0f}%)",
            transform=ax.transAxes, ha="right", va="top", fontsize=11,
            bbox=dict(boxstyle="round", fc="white", ec=GREY))
    ax.set_xlabel("DeltaPresentation = Mut_EL - WT_EL;  >0 = stronger mutant presentation")
    ax.set_ylabel("Number of mutant 9-mer-allele pairs")
    ax.set_title("Wild-type vs mutant presentation change (BigMHC)")
    save(fig, "fig8_delta_affinity.png")

def fig_binder_class(S):
    """Figure 9: Presenter class comparison."""
    bt = S["binder_by_type"]
    cats = ["Strong", "Weak", "Non-presenter"]
    mut = [bt.get(f"Mutant|{c}", 0) for c in cats]
    wt = [bt.get(f"WildType|{c}", 0) for c in cats]
    x = np.arange(len(cats)); w = 0.38
    fig, ax = plt.subplots(figsize=(9, 6))
    b1 = ax.bar(x - w/2, mut, w, color=RED, label="Mutant")
    b2 = ax.bar(x + w/2, wt, w, color=GREY, label="Wild-type")
    for b in list(b1) + list(b2):
        ax.text(b.get_x()+b.get_width()/2, b.get_height(),
                f"{int(b.get_height()):,}", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(x); ax.set_xticklabels(cats)
    ax.set_ylabel("Number of 9-mer-allele predictions")
    ax.set_title("Presenter classification: mutant vs wild-type 9-mers (BigMHC)")
    ax.legend(frameon=False)
    save(fig, "fig9_binder_class_mut_vs_wt.png")

def fig_strong_by_allele(S):
    """Figure 10: Presenters per HLA allele."""
    ab = S["allele_binder"]
    alleles = ["HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01"]
    cats = ["Strong", "Weak"]
    strong = [ab.get(f"{a}|Strong", 0) for a in alleles]
    weak = [ab.get(f"{a}|Weak", 0) for a in alleles]
    x = np.arange(len(alleles)); w = 0.38
    fig, ax = plt.subplots(figsize=(9, 6))
    b1 = ax.bar(x - w/2, strong, w, color=PURPLE, label="Strong (BigMHC_EL >= 0.70)")
    b2 = ax.bar(x + w/2, weak, w, color=TEAL, label="Weak (BigMHC_EL 0.50-0.70)")
    for b in list(b1) + list(b2):
        ax.text(b.get_x()+b.get_width()/2, b.get_height(),
                f"{int(b.get_height()):,}", ha="center", va="bottom", fontsize=10)
    ax.set_xticks(x); ax.set_xticklabels(alleles)
    ax.set_ylabel("Mutant 9-mers")
    ax.set_title("Predicted mutant-peptide presenters per HLA class I allele (BigMHC)")
    ax.legend(frameon=False)
    save(fig, "fig10_strong_binders_by_allele.png")

def fig_scatter(S):
    """Figure 11: Scatter plot of WT vs Mutant BigMHC_EL."""
    pts = np.array(S["scatter"], float)
    wt = pts[:, 0]; mut = pts[:, 1]
    favor = mut > wt
    fig, ax = plt.subplots(figsize=(7.8, 7.2))
    ax.scatter(wt[~favor], mut[~favor], s=6, alpha=0.25, color=GREY,
               label="Mutant weaker/equal")
    ax.scatter(wt[favor], mut[favor], s=6, alpha=0.3, color=GREEN,
               label="Mutant stronger (neoantigen-favourable)")
    lim = [-0.05, 1.05]
    ax.plot(lim, lim, color="black", lw=1.2, ls="--")
    ax.axvline(0.50, color=PURPLE, ls=":", lw=1, alpha=0.7)
    ax.axhline(0.50, color=PURPLE, ls=":", lw=1, alpha=0.7)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("Wild-type presentation probability  BigMHC_EL")
    ax.set_ylabel("Mutant presentation probability  BigMHC_EL")
    ax.set_title("Mutant vs wild-type presentation (BigMHC)\nabove diagonal = stronger mutant presentation")
    ax.legend(frameon=False, fontsize=10, loc="upper left")
    save(fig, "fig11_mut_vs_wt_scatter.png")

def fig_top_candidates():
    """Figure 12: Top shortlisted candidates horizontal bar plot."""
    rows = {}
    with open(SHORT) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            key = (p[ix["GeneName"]], p[ix["ProteinChange"]])
            el = float(p[ix["BigMHC_EL"]])
            freq = int(p[ix["MutationFrequency"]])
            allele = p[ix["HLAAllele"]]
            r = rows.get(key)
            if r is None or el > r["el"]:
                rows[key] = {"freq": freq, "el": el, "allele": allele}
            rows[key]["freq"] = max(rows[key]["freq"], freq)
    items = sorted(rows.items(), key=lambda kv: kv[1]["freq"], reverse=True)[:15]
    items = items[::-1]
    labels = [f"{g} {pc}" for (g, pc), _ in items]
    freqs = [v["freq"] for _, v in items]
    els = [v["el"] for _, v in items]
    alleles = [v["allele"] for _, v in items]
    fig, ax = plt.subplots(figsize=(10, 7))
    bars = ax.barh(range(len(items)), freqs, color=ORANGE, edgecolor="white")
    for i, (v, e, al) in enumerate(zip(freqs, els, alleles)):
        ax.text(v + 0.4, i, f"{v} samples  |  EL: {e:.3f}  |  {al}",
                va="center", fontsize=9)
    ax.set_yticks(range(len(items))); ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel("Per-mutation recurrence (tumour samples)")
    ax.set_xlim(0, max(freqs) * 1.5)
    ax.set_title("Top prioritised neoantigen candidates (TCGA-COAD, BigMHC)\nrecurrent, expressed, stronger-than-WT presenters")
    save(fig, "fig12_top_candidates.png")

def fig_workflow():
    """Figure 13: 8-step pipeline workflow schematic diagram."""
    steps = [
        ("1. Somatic mutations\nGDC TCGA-COAD MAF\n(GRCh38)", BLUE),
        ("2. Filter + mutation-\nby-sample matrix\n(01)", BLUE),
        ("3. RNA-seq\nRSEM -> TPM\nmatrix (02)", TEAL),
        ("4. Integrate:\nGeneLevelTPM\nby gene (03)", TEAL),
        ("5. QC +\nfigures (04)", GREY),
        ("6. Variant->protein;\nmutant + WT\npeptides (05)", GREEN),
        ("7. BigMHC presentation\n+ WT delta (06)", ORANGE),
        ("8. Prioritise\nneoantigen\ncandidates (07)", RED),
    ]
    fig, ax = plt.subplots(figsize=(12, 6.5))
    ax.set_xlim(0, 10.4); ax.set_ylim(0, 8); ax.axis("off")
    bw, bh, gap = 2.0, 1.7, 0.45
    xs = [0.3 + i * (bw + gap) for i in range(4)]
    top_y, bot_y = 5.3, 1.4
    positions = [(xs[0], top_y), (xs[1], top_y), (xs[2], top_y), (xs[3], top_y),
                 (xs[3], bot_y), (xs[2], bot_y), (xs[1], bot_y), (xs[0], bot_y)]
    for (txt, col), (x, y) in zip(steps, positions):
        box = FancyBboxPatch((x, y), bw, bh, boxstyle="round,pad=0.04",
                             fc=col, ec="white", lw=2, alpha=0.92)
        ax.add_patch(box)
        ax.text(x + bw/2, y + bh/2, txt, ha="center", va="center",
                color="white", fontsize=10, fontweight="bold")

    def arrow(p1, p2):
        ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=20,
                                     color="#444444", lw=2.0))
    for i in range(3):
        x, y = positions[i]
        arrow((x + bw + 0.02, y + bh/2), (positions[i+1][0] - 0.02, y + bh/2))
    x4, y4 = positions[3]
    arrow((x4 + bw/2, y4 - 0.02), (positions[4][0] + bw/2, positions[4][1] + bh + 0.02))
    for i in range(4, 7):
        x, y = positions[i]
        arrow((x - 0.02, y + bh/2), (positions[i+1][0] + bw + 0.02, y + bh/2))
    ax.set_title("Computational workflow (BigMHC pipeline)",
                 fontsize=16, fontweight="bold", pad=16)
    save(fig, "fig13_workflow_schematic.png")

# =============================================================================
# MAIN PIPELINE EXECUTION
# =============================================================================
def main():
    with open(SUM) as fh:
        S = json.load(fh)
    fig_funnel()
    fig_variant_class()
    fig_mut_vs_wt_affinity(S)
    fig_delta(S)
    fig_binder_class(S)
    fig_strong_by_allele(S)
    fig_scatter(S)
    fig_top_candidates()
    fig_workflow()
    log("all presentation figures written to figures/ (BigMHC)")

if __name__ == "__main__":
    sys.exit(main())
