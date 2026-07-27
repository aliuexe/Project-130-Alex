#!/usr/bin/env python3
"""
16_immunogenicity.py
Project 130 - Colorectal cancer (TCGA-COAD)  --  immunogenicity extension

Adds a T-cell IMMUNOGENICITY score, distinct from HLA binding. Binding tells us
a peptide can be DISPLAYED; immunogenicity estimates whether a displayed peptide
can actually be RECOGNISED by a T cell and trigger a response. This fills the
ImmunogenicityScore column that was previously NA.

Model: Calis et al. 2013, "Properties of MHC class I presented peptides that
enhance immunogenicity", PLoS Comput Biol (the IEDB class-I immunogenicity
predictor). It is a position-weighted sum of per-amino-acid immunogenicity
values:
    score = sum over positions i ( weight[i] * immunoscale[aa_i] )
with the default masking used by IEDB: positions 1, 2 and the C-terminus are
MASKED (they are HLA anchors that affect binding, not T-cell recognition), and
the central positions 4-6 carry the most weight (main TCR-contact residues).
Higher score = more likely to be immunogenic. Validated on 9-mers only.
Source constants: Calis 2013 / IEDB tools.iedb.org/immunogenicity.

For each candidate we score the MUTANT and WILD-TYPE 9-mer and report the
difference: a good neoantigen is not only immunogenic but MORE immunogenic than
its wild-type counterpart (the mutation increases T-cell visibility). We also
flag whether the mutation sits at a masked/anchor position (then it changes
binding, not immunogenicity) or at a TCR-contact position (then it can change
what the T cell sees).

Inputs:  results/practical_neoantigens.tsv, results/neoantigen_candidates_shortlist.tsv
Outputs: results/practical_neoantigens_scored.tsv
         figures/fig24_immunogenicity_mut_vs_wt.png
         figures/fig25_top_immunogenic_candidates.png
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results"); FIG = os.path.join(BASE, "figures")
PRAC = os.path.join(RES, "practical_neoantigens.tsv")
SHORT = os.path.join(RES, "neoantigen_candidates_shortlist.tsv")
os.makedirs(FIG, exist_ok=True)

# ---- Calis 2013 / IEDB constants ------------------------------------------
# per-amino-acid immunogenicity values (aromatic/large = enriching; K/S/M = depleting)
IMMUNOSCALE = {
    "A": 0.127, "C": -0.175, "D": 0.072, "E": 0.325, "F": 0.380,
    "G": 0.110, "H": 0.105, "I": 0.432, "K": -0.700, "L": -0.036,
    "M": -0.570, "N": -0.021, "P": -0.036, "Q": -0.376, "R": 0.168,
    "S": -0.537, "T": -0.062, "V": 0.134, "W": 0.719, "Y": -0.012,
}
# position importance weights for a 9-mer (1-indexed); central P4-6 dominate.
POS_WEIGHT = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]
# default IEDB masking: first, second and C-terminal residue are masked.
MASK_9 = {1, 2, 9}

def log(m): print("[16]", m, flush=True)

def immunogenicity(pep):
    """Calis/IEDB class-I immunogenicity for a 9-mer (higher = more immunogenic)."""
    if len(pep) != 9 or any(a not in IMMUNOSCALE for a in pep):
        return None                       # only defined for standard-AA 9-mers
    s = 0.0
    for i, aa in enumerate(pep):
        pos = i + 1
        if pos in MASK_9:                 # anchors masked (binding, not recognition)
            continue
        s += POS_WEIGHT[i] * IMMUNOSCALE[aa]
    return s

def ref_aa_from_protchange(pchg):
    # p.G12D -> 'G'  (the reference amino acid, the letter right after 'p.')
    t = pchg[2:] if pchg.startswith("p.") else pchg
    return t[0] if t and t[0].isalpha() else None

def main():
    # ---- validation sanity check -----------------------------------------
    demo = {"WWFWWFWWF": "aromatic-rich (should be HIGH)",
            "KSKSKSKSK": "K/S-rich (should be LOW)",
            "VVGADGVGK": "KRAS G12D mutant 9-mer"}
    log("Validation (score, higher = more immunogenic):")
    for p, d in demo.items():
        print(f"    {p}  {immunogenicity(p):+.3f}   {d}")

    # ---- exact MutPos per (gene, proteinchange, peptide) from the shortlist -
    mutpos_map = {}
    with open(SHORT) as fh:
        sh = fh.readline().rstrip("\n").split("\t"); si = {c: i for i, c in enumerate(sh)}
        for l in fh:
            q = l.rstrip("\n").split("\t")
            mutpos_map[(q[si["GeneName"]], q[si["ProteinChange"]], q[si["Peptide"]])] = int(q[si["MutPos"]])

    # ---- score the practical neoantigens ---------------------------------
    with open(PRAC) as fh:
        first = fh.readline()
        header = fh.readline().rstrip("\n").split("\t")   # (line 1 is a # comment)
        ix = {c: i for i, c in enumerate(header)}
        rows = [l.rstrip("\n").split("\t") for l in fh]

    out_rows = []
    for p in rows:
        pep = p[ix["Peptide"]]
        pchg = p[ix["ProteinChange"]]
        refaa = ref_aa_from_protchange(pchg)          # wild-type residue, e.g. G in p.G12D
        imm_mut = immunogenicity(pep)
        # WT 9-mer = mutant peptide with the mutated position reverted to refAA,
        # using the exact MutPos recorded in the shortlist.
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

    # composite priority: strong binder + immunogenic + expressed + recurrent
    def keyf(d):
        imm = d["ImmunogenicityScore"] if d["ImmunogenicityScore"] is not None else -9
        return (imm, d.get("TumoursCovered", 0))
    out_rows.sort(key=lambda d: (float(d["ImmunogenicityScore"]) if d["ImmunogenicityScore"] is not None else -9),
                  reverse=True)

    out_cols = header + ["ImmunogenicityScore", "Immunogenicity_WT",
                         "Immunogenicity_delta", "MutPosInPeptide", "MutationAtAnchor"]
    with open(os.path.join(RES, "practical_neoantigens_scored.tsv"), "w") as fh:
        fh.write(first)                      # keep the # TotalSamples comment
        fh.write("\t".join(out_cols) + "\n")
        for d in out_rows:
            fh.write("\t".join(str(d.get(c, "NA")) for c in out_cols) + "\n")
    log(f"wrote practical_neoantigens_scored.tsv ({len(out_rows)} candidates)")

    # summary stats
    imm_vals = [d["ImmunogenicityScore"] for d in out_rows if d["ImmunogenicityScore"] is not None]
    n_pos_delta = sum(1 for d in out_rows if d["Immunogenicity_delta"] not in (None,"NA")
                      and float(d["Immunogenicity_delta"]) > 0)
    n_tcr = sum(1 for d in out_rows if d["MutationAtAnchor"] == "TCR-contact(recognition)")
    log(f"immunogenic (score>0): {sum(1 for v in imm_vals if v>0)}/{len(imm_vals)}; "
        f"mutation increases immunogenicity (delta>0): {n_pos_delta}; "
        f"mutation at TCR-contact position: {n_tcr}")

    # ---- FIGURE 24: mutant vs WT immunogenicity of practical candidates ----
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

    # ---- FIGURE 25: top candidates by immunogenicity ----------------------
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

    # ---- console: KRAS drivers immunogenicity -----------------------------
    log("Immunogenicity of key driver neoantigens (practical set):")
    for d in out_rows:
        if (d["GeneName"], d["ProteinChange"]) in [("KRAS","p.G12V"),("KRAS","p.G12D"),
                ("KRAS","p.G12C"),("PIK3CA","p.E542K"),("SMAD4","p.R361H"),("FBXW7","p.S582L")]:
            print(f"    {d['GeneName']:7s}{d['ProteinChange']:9s} {d['Peptide']}  "
                  f"immuno={d['ImmunogenicityScore']}  WT={d['Immunogenicity_WT']}  "
                  f"delta={d['Immunogenicity_delta']}  {d['MutationAtAnchor']}")

if __name__ == "__main__":
    sys.exit(main())
