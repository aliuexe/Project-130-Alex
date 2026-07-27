#!/usr/bin/env python3
r"""
16_immunogenicity.py
Project 130 - Colorectal Cancer (TCGA-COAD) Neoantigen Discovery Pipeline

===============================================================================
BIOLOGICAL & COMPUTATIONAL PURPOSE
===============================================================================
This script adds a quantitative T-cell IMMUNOGENICITY score (Calis et al. 2013 / IEDB
class-I immunogenicity model) to evaluate whether a presented peptide is capable of
triggering a TCR (T-Cell Receptor) recognition response.

Motivation: HLA presentation (`BigMHC_EL`) measures cell-surface display, but display alone
does not guarantee T-cell activation. Immunogenicity evaluates amino acid physicochemical
properties at central TCR-contact residue positions.

===============================================================================
MATHEMATICAL FORMULATION (CALIS ET AL. 2013 / IEDB MODEL)
===============================================================================
For a 9-mer peptide $P = (a_1, a_2, \dots, a_9)$:

    \text{ImmunogenicityScore} = \sum_{i=1}^{9} \text{Weight}[i] \times \text{ImmunoScale}[a_i]

Where:
  - `ImmunoScale[a]`: Empirical amino acid immunogenicity weight (large/aromatic residues
    $W, I, F, E$ enrich immunogenicity; $K, S, M$ deplete it).
  - `Weight[i]`: Position-specific TCR contact weighting for 9-mers:
      $\text{Weight} = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]$
  - Masking Rule: Positions P1, P2, and P9 are MASKED (weights $= 0.00$) because they serve
    as HLA binding anchors rather than TCR contact positions. Positions P4–P6 dominate.

===============================================================================
INPUT & OUTPUT CONTRACTS
===============================================================================
Inputs:
  - `results/practical_neoantigens.tsv`
  - `results/neoantigen_candidates_shortlist.tsv`

Outputs:
  - `results/practical_neoantigens_scored.tsv` (Immunogenicity-scored candidates)
  - `figures/fig24_immunogenicity_mut_vs_wt.png` (Mutant vs WT immunogenicity histogram)
  - `figures/fig25_top_immunogenic_candidates.png` (Top immunogenic candidates plot)
"""

import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =============================================================================
# FILE PATHS & RESOURCE RESOLUTION
# =============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results"); FIG = os.path.join(BASE, "figures")
PRAC = os.path.join(RES, "practical_neoantigens.tsv")
SHORT = os.path.join(RES, "neoantigen_candidates_shortlist.tsv")
os.makedirs(FIG, exist_ok=True)

# =============================================================================
# CALIS 2013 / IEDB PHYSICOMATHEMATICAL MODEL CONSTANTS
# =============================================================================
# Per-amino-acid immunogenicity values
IMMUNOSCALE = {
    "A": 0.127, "C": -0.175, "D": 0.072, "E": 0.325, "F": 0.380,
    "G": 0.110, "H": 0.105, "I": 0.432, "K": -0.700, "L": -0.036,
    "M": -0.570, "N": -0.021, "P": -0.036, "Q": -0.376, "R": 0.168,
    "S": -0.537, "T": -0.062, "V": 0.134, "W": 0.719, "Y": -0.012,
}
# 9-mer position importance weights (1-indexed P1-P9)
POS_WEIGHT = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]
# Default IEDB masking (N-terminal anchors P1, P2 and C-terminal anchor P9)
MASK_9 = {1, 2, 9}

def log(m):
    """Prints timestamped progress messages to stdout with line flushing."""
    print("[16]", m, flush=True)

def immunogenicity(pep):
    """Calculates Calis/IEDB Class I immunogenicity score for a 9-mer peptide."""
    if len(pep) != 9 or any(a not in IMMUNOSCALE for a in pep):
        return None
    s = 0.0
    for i, aa in enumerate(pep):
        pos = i + 1
        if pos in MASK_9:
            continue
        s += POS_WEIGHT[i] * IMMUNOSCALE[aa]
    return s

def ref_aa_from_protchange(pchg):
    """Extracts reference amino acid single-letter code from HGVSp string (e.g. p.G12D -> 'G')."""
    t = pchg[2:] if pchg.startswith("p.") else pchg
    return t[0] if t and t[0].isalpha() else None

# =============================================================================
# MAIN PIPELINE EXECUTION
# =============================================================================
def main():
    # Model sanity check logging
    demo = {"WWFWWFWWF": "aromatic-rich (should be HIGH)",
            "KSKSKSKSK": "K/S-rich (should be LOW)",
            "VVGADGVGK": "KRAS G12D mutant 9-mer"}
    log("Validation (score, higher = more immunogenic):")
    for p, d in demo.items():
        print(f"    {p}  {immunogenicity(p):+.3f}   {d}")

    # Map exact 1-based mutation position in peptide from shortlist
    mutpos_map = {}
    with open(SHORT) as fh:
        sh = fh.readline().rstrip("\n").split("\t"); si = {c: i for i, c in enumerate(sh)}
        for l in fh:
            q = l.rstrip("\n").split("\t")
            mutpos_map[(q[si["GeneName"]], q[si["ProteinChange"]], q[si["Peptide"]])] = int(q[si["MutPos"]])

    # =========================================================================
    # STEP 1: SCORE MUTANT AND WILD-TYPE PRACTICAL CANDIDATES
    # =========================================================================
    with open(PRAC) as fh:
        first = fh.readline()
        header = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(header)}
        rows = [l.rstrip("\n").split("\t") for l in fh]

    out_rows = []
    for p in rows:
        pep = p[ix["Peptide"]]
        pchg = p[ix["ProteinChange"]]
        refaa = ref_aa_from_protchange(pchg)
        imm_mut = immunogenicity(pep)
        mutpos = mutpos_map.get((p[ix["GeneName"]], pchg, pep))
        wt_pep = None
        if refaa and mutpos and 1 <= mutpos <= len(pep):
            i0 = mutpos - 1
            wt_pep = pep[:i0] + refaa + pep[i0+1:]
        imm_wt = immunogenicity(wt_pep) if wt_pep else None
        d = dict(zip(header, p))
        d["ImmunogenicityScore"] = None if imm_mut is None else round(imm_mut, 4)
        d["Immunogenicity_WT"] = None if imm_wt is None else round(imm_wt, 4)
        d["Immunogenicity_delta"] = (None if (imm_mut is None or imm_wt is None)
                                     else round(imm_mut - imm_wt, 4))
        d["MutPosInPeptide"] = mutpos
        d["MutationAtAnchor"] = ("NA" if mutpos is None
                                 else ("Anchor(binding)" if mutpos in MASK_9
                                       else "TCR-contact(recognition)"))
        out_rows.append(d)

    # Sort candidates by immunogenicity score
    out_rows.sort(key=lambda d: (float(d["ImmunogenicityScore"]) if d["ImmunogenicityScore"] is not None else -9),
                  reverse=True)

    # Export scored practical candidates
    out_cols = header + ["ImmunogenicityScore", "Immunogenicity_WT",
                         "Immunogenicity_delta", "MutPosInPeptide", "MutationAtAnchor"]
    with open(os.path.join(RES, "practical_neoantigens_scored.tsv"), "w") as fh:
        fh.write(first)
        fh.write("\t".join(out_cols) + "\n")
        for d in out_rows:
            fh.write("\t".join(str(d.get(c, "NA")) for c in out_cols) + "\n")
    log(f"wrote practical_neoantigens_scored.tsv ({len(out_rows)} candidates)")

    imm_vals = [d["ImmunogenicityScore"] for d in out_rows if d["ImmunogenicityScore"] is not None]
    n_pos_delta = sum(1 for d in out_rows if d["Immunogenicity_delta"] not in (None,"NA")
                      and float(d["Immunogenicity_delta"]) > 0)
    n_tcr = sum(1 for d in out_rows if d["MutationAtAnchor"] == "TCR-contact(recognition)")
    log(f"immunogenic (score>0): {sum(1 for v in imm_vals if v>0)}/{len(imm_vals)}; "
        f"mutation increases immunogenicity (delta>0): {n_pos_delta}; "
        f"mutation at TCR-contact position: {n_tcr}")

    # =========================================================================
    # STEP 2: RENDER FIGURE 24 — IMMUNOGENICITY HISTOGRAM (MUTANT VS WT)
    # =========================================================================
    mm = [d["ImmunogenicityScore"] for d in out_rows if d["ImmunogenicityScore"] is not None]
    ww = [d["Immunogenicity_WT"] for d in out_rows if d["Immunogenicity_WT"] is not None]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    bins = np.linspace(min(mm+ww)-0.02, max(mm+ww)+0.02, 30)
    ax.hist(ww, bins=bins, color="#AAB2B8", alpha=0.7, label="Wild-type peptide")
    ax.hist(mm, bins=bins, color="#B4433B", alpha=0.6, label="Mutant peptide")
    ax.axvline(0, color="black", lw=1, ls="--")
    ax.set_xlabel("Calis/IEDB immunogenicity score (higher = more immunogenic)")
    ax.set_ylabel("Practical neoantigens")
    ax.set_title("T-cell immunogenicity: mutant vs wild-type (practical candidates)")
    ax.legend(frameon=False); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig24_immunogenicity_mut_vs_wt.png"), dpi=160)
    plt.close(fig); log("wrote fig24_immunogenicity_mut_vs_wt.png")

    # =========================================================================
    # STEP 3: RENDER FIGURE 25 — TOP IMMUNOGENIC CANDIDATES
    # =========================================================================
    top = [d for d in out_rows if d["ImmunogenicityScore"] is not None][:15][::-1]
    labels = [f"{d['GeneName']} {d['ProteinChange']} ({d['Peptide']})" for d in top]
    vals = [float(d["ImmunogenicityScore"]) for d in top]
    colors = ["#2A9D8F" if d["MutationAtAnchor"].startswith("TCR") else "#EE7733" for d in top]
    fig, ax = plt.subplots(figsize=(11, 6.5))
    ax.barh(range(len(top)), vals, color=colors)
    for i, d in enumerate(top):
        ax.text(vals[i] + 0.002, i, f"cov {d.get('TumoursCovered','?')} · {d['HLAAllele']}",
                va="center", fontsize=8)
    ax.set_yticks(range(len(top))); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Immunogenicity score (Calis/IEDB)")
    ax.set_title("Most immunogenic practical neoantigens (TCGA-COAD)\n"
                 "teal = mutation at TCR-contact position, orange = at anchor")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig25_top_immunogenic_candidates.png"), dpi=160)
    plt.close(fig); log("wrote fig25_top_immunogenic_candidates.png")

    # Console preview logging for KRAS driver neoantigens
    log("Immunogenicity of key driver neoantigens (practical set):")
    for d in out_rows:
        if (d["GeneName"], d["ProteinChange"]) in [("KRAS","p.G12V"),("KRAS","p.G12D"),
                ("KRAS","p.G12C"),("PIK3CA","p.E542K"),("SMAD4","p.R361H"),("FBXW7","p.S582L")]:
            print(f"    {d['GeneName']:7s}{d['ProteinChange']:9s} {d['Peptide']}  "
                  f"immuno={d['ImmunogenicityScore']}  WT={d['Immunogenicity_WT']}  "
                  f"delta={d['Immunogenicity_delta']}  {d['MutationAtAnchor']}")

if __name__ == "__main__":
    sys.exit(main())
