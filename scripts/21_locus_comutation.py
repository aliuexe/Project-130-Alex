#!/usr/bin/env python3
"""
21_locus_comutation.py
Project 130 - Colorectal cancer (TCGA-COAD)  --  EXTENSION (PI request)

LOCUS-LEVEL co-mutation, ALL SAMPLES, NO hypermutator filtering.

Unlike scripts 10/12/15 (which work at GENE level - a gene is "mutated" if ANY
of its missense SNVs is present), this script tests co-occurrence between
SPECIFIC MUTATIONS / SPECIFIC LOCI, i.e. exact protein changes such as
KRAS p.G12D + TP53 p.R175H. Amino-acid resolution is preserved, so a strong,
tumour-specific hotspot is never merged into a gene-level "any KRAS" bucket.

Design choices (this run):
  * Level      : specific locus = (GeneName, AminoAcidChange). Rows that share a
                 protein change (e.g. two HGVSc for K117N) are OR-ed together.
  * Gene scope : curated CRC driver panel (34 genes) - keeps results
                 interpretable and avoids gene-size artefacts (TTN, MUC16).
  * Loci kept  : recurrence >= MIN_REC (3 tumours). Loci below this have
                 essentially no statistical power for co-occurrence.
  * Samples    : ALL 586 tumours. NO hypermutator exclusion ANYWHERE. (The
                 >=200-SNV cut in script 15 is a co-mutation *robustness* control
                 only; it is deliberately NOT applied to this analysis.)
  * Reporting  : to answer "is a specific-locus co-mutation just a hypermutator
                 artefact?", every co-occurring pair ALSO reports how many of the
                 co-mutated tumours are hypermutated vs not - this is annotation
                 only, it never filters anything.

Statistics (pure Python, no scipy): Fisher's exact via the hypergeometric
distribution (odds ratio + one-sided co-occurrence / mutual-exclusivity
p-values, Benjamini-Hochberg FDR); triples via observed-vs-expected + Poisson.

Inputs:  results/03_integrated_mutation_expression.tsv
Outputs: results/locus_comutation_driver_pairs.tsv
         results/locus_comutation_driver_triples.tsv
         results/locus_comutation_summary.txt
         figures/fig28_locus_comutation_pairs.png
"""
import itertools
import math
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
FIG = os.path.join(BASE, "figures")
INT = os.path.join(RES, "03_integrated_mutation_expression.tsv")
OUT_PAIRS = os.path.join(RES, "locus_comutation_driver_pairs.tsv")
OUT_TRIPLES = os.path.join(RES, "locus_comutation_driver_triples.tsv")
OUT_SUM = os.path.join(RES, "locus_comutation_summary.txt")
os.makedirs(FIG, exist_ok=True)

MIN_REC = 3        # a specific locus must occur in >= this many tumours to test
HYPER_CUT = 200    # missense-SNV burden flag - REPORTING ONLY, never a filter

# ---- Curated colorectal-cancer driver genes (same panel as scripts 10/12) ---
DRIVERS = [
    "APC", "TP53", "KRAS", "PIK3CA", "FBXW7", "SMAD4", "TCF7L2", "NRAS",
    "SMAD2", "CTNNB1", "BRAF", "SOX9", "ARID1A", "AMER1", "FAM123B", "ATM",
    "KMT2C", "KMT2D", "ERBB2", "ERBB3", "PTEN", "ACVR2A", "GNAS", "BMPR1A",
    "TGFBR2", "RNF43", "B2M", "POLE", "MSH6", "CASP8", "ELF3", "PCBP1",
    "AXIN2", "MAP2K4", "CDC27",
]

def log(m): print("[21]", m, flush=True)

# --- pure-python exact statistics (identical to script 10) ------------------
def hyper_pmf(x, N, K, n):
    if x < max(0, n - (N - K)) or x > min(K, n):
        return 0.0
    return math.comb(K, x) * math.comb(N - K, n - x) / math.comb(N, n)

def fisher_tails(a, nA, nB, N):
    lo, hi = max(0, nA + nB - N), min(nA, nB)
    right = sum(hyper_pmf(x, N, nA, nB) for x in range(a, hi + 1))
    left = sum(hyper_pmf(x, N, nA, nB) for x in range(lo, a + 1))
    return min(1.0, right), min(1.0, left)

def poisson_sf(obs, lam):
    if obs <= 0:
        return 1.0
    cdf = sum(math.exp(-lam) * lam**k / math.factorial(k) for k in range(obs))
    return max(0.0, min(1.0, 1.0 - cdf))

def bh_fdr(pvals):
    m = len(pvals)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = m - rank + 1
        val = min(prev, pvals[i] * m / k)
        q[i] = val
        prev = val
    return q

def _fmt(v):
    if isinstance(v, float):
        if v != v:
            return "NA"
        if v == 0 or (abs(v) >= 1e-4 and abs(v) < 1e6):
            return f"{v:.4g}"
        return f"{v:.3e}"
    return v

def main():
    # ---- read header, locate sample columns by TCGA prefix ----------------
    with open(INT) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    s0 = next(i for i, c in enumerate(header) if c.startswith("TCGA"))
    samples = header[s0:]
    N = len(samples)
    driverset = set(DRIVERS)

    # ---- single pass: per-sample burden (ALL rows) + driver locus vectors --
    burden = np.zeros(N, dtype=np.int32)      # missense-SNV count per tumour
    loci = {}                                  # (gene, aachange) -> bool vector
    with open(INT) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            vec = np.fromiter((1 if v == "1" else 0 for v in p[s0:]),
                              dtype=np.int8, count=N)
            burden += vec                      # burden uses EVERY mutation row
            g = p[0]
            if g not in driverset:
                continue
            key = (g, p[2])                    # (gene, protein change) = a locus
            if key in loci:
                loci[key] |= vec               # OR rows sharing this locus
            else:
                loci[key] = vec.copy()

    # ---- hypermutator flag: REPORTING ONLY, applied to nothing ------------
    hyper = burden >= HYPER_CUT
    n_hyper = int(hyper.sum())
    log(f"Samples: {N}  |  hypermutators (>= {HYPER_CUT} SNVs, report-only): "
        f"{n_hyper} ({100*n_hyper/N:.1f}%)  |  NO samples excluded")

    # ---- keep only recurrent loci (>= MIN_REC tumours) --------------------
    rec = {k: int(v.sum()) for k, v in loci.items()}
    kept = [k for k in loci if rec[k] >= MIN_REC]
    # order by gene (curated order) then descending recurrence
    gorder = {g: i for i, g in enumerate(DRIVERS)}
    kept.sort(key=lambda k: (gorder.get(k[0], 999), -rec[k], k[1]))
    log(f"Driver loci with recurrence >= {MIN_REC}: {len(kept)} "
        f"(from {len(loci)} distinct driver loci)")

    # ---- PAIRS: Fisher's exact at locus level -----------------------------
    rows = []
    for kA, kB in itertools.combinations(kept, 2):
        A, B = loci[kA], loci[kB]
        both = A & B
        a = int(both.sum())
        nA, nB = rec[kA], rec[kB]
        right, left = fisher_tails(a, nA, nB, N)
        b = nA - a; c = nB - a; d = N - nA - nB + a
        orr = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        a_hyper = int((both & hyper).sum())
        rows.append(dict(
            GeneA=kA[0], LocusA=kA[1], recA=nA,
            GeneB=kB[0], LocusB=kB[1], recB=nB,
            nBoth=a, nBoth_hyper=a_hyper, nBoth_nonhyper=a - a_hyper,
            PctHyper=(round(100 * a_hyper / a, 1) if a else float("nan")),
            Expected=round(nA * nB / N, 3), OddsRatio=round(orr, 3),
            p_cooccur=right, p_exclusive=left,
            SameGene=(kA[0] == kB[0])))
    # FDR over ALL tested pairs
    for r, q in zip(rows, bh_fdr([r["p_cooccur"] for r in rows])):
        r["FDR_cooccur"] = q
    for r, q in zip(rows, bh_fdr([r["p_exclusive"] for r in rows])):
        r["FDR_exclusive"] = q
    for r in rows:
        r["Relationship"] = (
            "Co-occurring" if r["OddsRatio"] > 1 and r["FDR_cooccur"] < 0.05
            else "Mutually exclusive" if r["OddsRatio"] < 1 and r["FDR_exclusive"] < 0.05
            else "n.s.")

    # write only pairs that ACTUALLY co-occur (nBoth >= 1), most co-mutated first
    co_rows = [r for r in rows if r["nBoth"] >= 1]
    co_rows.sort(key=lambda r: (-r["nBoth"], -r["OddsRatio"]))
    cols = ["GeneA", "LocusA", "recA", "GeneB", "LocusB", "recB", "nBoth",
            "nBoth_hyper", "nBoth_nonhyper", "PctHyper", "Expected", "OddsRatio",
            "p_cooccur", "FDR_cooccur", "p_exclusive", "FDR_exclusive",
            "SameGene", "Relationship"]
    with open(OUT_PAIRS, "w") as fh:
        fh.write(f"# Locus-level co-mutation. Level = (gene, protein change). "
                 f"ALL {N} samples; NO hypermutator filtering. "
                 f"Loci with recurrence >= {MIN_REC}. Hyper columns are "
                 f"report-only (burden >= {HYPER_CUT} SNVs).\n")
        fh.write("\t".join(cols) + "\n")
        for r in co_rows:
            fh.write("\t".join(str(_fmt(r[c])) for c in cols) + "\n")
    log(f"Wrote {OUT_PAIRS}: {len(co_rows)} co-occurring locus pairs "
        f"({len(rows)} tested)")

    # ---- TRIPLES: observed vs expected + Poisson --------------------------
    trip = []
    for kA, kB, kC in itertools.combinations(kept, 3):
        both = int((loci[kA] & loci[kB] & loci[kC]).sum())
        if both == 0:
            continue
        expd = rec[kA] * rec[kB] * rec[kC] / (N * N)
        h = int((loci[kA] & loci[kB] & loci[kC] & hyper).sum())
        trip.append(dict(
            GeneA=kA[0], LocusA=kA[1], GeneB=kB[0], LocusB=kB[1],
            GeneC=kC[0], LocusC=kC[1], nAll3=both, nAll3_hyper=h,
            Expected=round(expd, 4),
            Enrichment=round(both / expd, 2) if expd > 0 else float("inf"),
            p_cooccur=poisson_sf(both, expd)))
    trip.sort(key=lambda r: -r["nAll3"])
    tcols = ["GeneA", "LocusA", "GeneB", "LocusB", "GeneC", "LocusC",
             "nAll3", "nAll3_hyper", "Expected", "Enrichment", "p_cooccur"]
    with open(OUT_TRIPLES, "w") as fh:
        fh.write(f"# Locus-level co-mutation TRIPLES. ALL {N} samples; no "
                 f"hypermutator filtering. Loci recurrence >= {MIN_REC}.\n")
        fh.write("\t".join(tcols) + "\n")
        for r in trip:
            fh.write("\t".join(str(_fmt(r[c])) for c in tcols) + "\n")
    log(f"Wrote {OUT_TRIPLES}: {len(trip)} co-occurring locus triples")

    # ---- SUMMARY ----------------------------------------------------------
    sig_co = sorted([r for r in co_rows if r["Relationship"] == "Co-occurring"],
                    key=lambda r: r["FDR_cooccur"])
    with open(OUT_SUM, "w") as fh:
        fh.write("LOCUS-LEVEL DRIVER CO-MUTATION (colorectal, TCGA-COAD)\n")
        fh.write("Level = specific mutation (gene, protein change); "
                 "e.g. KRAS p.G12D + TP53 p.R175H.\n")
        fh.write(f"Samples: ALL {N} tumours. NO hypermutator filtering applied "
                 f"anywhere.\n")
        fh.write(f"Hypermutators (burden >= {HYPER_CUT} SNVs) present in the "
                 f"cohort: {n_hyper} ({100*n_hyper/N:.1f}%) - REPORTED, not "
                 f"removed.\n")
        fh.write(f"Driver loci tested (recurrence >= {MIN_REC}): {len(kept)}.\n\n")
        fh.write("Top 20 co-occurring specific-locus pairs "
                 "(by tumours carrying BOTH exact mutations):\n")
        fh.write(f"  {'LocusA':>16} {'LocusB':>16}  both  hyper nonhyper  "
                 f"%hyp    OR   FDR_co  rel\n")
        for r in co_rows[:20]:
            la = f"{r['GeneA']} {r['LocusA']}"
            lb = f"{r['GeneB']} {r['LocusB']}"
            ph = "NA" if r["PctHyper"] != r["PctHyper"] else f"{r['PctHyper']:.0f}"
            fh.write(f"  {la:>16} {lb:>16}  {r['nBoth']:4d}  {r['nBoth_hyper']:4d} "
                     f"{r['nBoth_nonhyper']:7d}  {ph:>4}  {r['OddsRatio']:5.2f} "
                     f"{r['FDR_cooccur']:.1e}  {r['Relationship']}\n")
        fh.write(f"\nStatistically significant co-occurring locus pairs "
                 f"(FDR<0.05): {len(sig_co)}\n")
        for r in sig_co[:20]:
            la = f"{r['GeneA']} {r['LocusA']}"
            lb = f"{r['GeneB']} {r['LocusB']}"
            fh.write(f"  {la} + {lb}: both={r['nBoth']} "
                     f"(hyper={r['nBoth_hyper']}, non-hyper={r['nBoth_nonhyper']}), "
                     f"OR={r['OddsRatio']:.2f}, FDR={r['FDR_cooccur']:.2e}\n")
        # how much of locus co-mutation is carried by hypermutators?
        tot_both = sum(r["nBoth"] for r in co_rows)
        tot_hyper = sum(r["nBoth_hyper"] for r in co_rows)
        fh.write(f"\nAcross all co-occurring locus pairs, "
                 f"{tot_hyper}/{tot_both} co-mutation events "
                 f"({100*tot_hyper/tot_both:.0f}%) occur in hypermutated tumours "
                 f"(which are {100*n_hyper/N:.0f}% of the cohort).\n")
        fh.write("Interpretation: a specific-locus co-mutation concentrated in "
                 "non-hypermutated tumours is the stronger candidate for a real, "
                 "generalisable association; one seen only in hypermutators may "
                 "reflect burden. Either way NOTHING is filtered - all are kept.\n")
    log(f"Wrote {OUT_SUM}")

    # ---- FIGURE 28: top co-occurring locus pairs, hyper vs non-hyper split -
    top = co_rows[:14][::-1]
    if top:
        labels = [f"{r['GeneA']} {r['LocusA']} + {r['GeneB']} {r['LocusB']}"
                  for r in top]
        nh = [r["nBoth_nonhyper"] for r in top]
        hy = [r["nBoth_hyper"] for r in top]
        y = range(len(top))
        fig, ax = plt.subplots(figsize=(11, 7))
        ax.barh(y, nh, color="#2E6E9E", label="non-hypermutated tumours")
        ax.barh(y, hy, left=nh, color="#C0603A", label="hypermutated tumours")
        for i, r in enumerate(top):
            ax.text(r["nBoth"] + 0.05, i, str(r["nBoth"]), va="center", fontsize=9)
        ax.set_yticks(list(y)); ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel(f"Tumours carrying BOTH specific mutations (of N = {N})")
        ax.set_title("Locus-level driver co-mutation (TCGA-COAD, all samples)\n"
                     "specific protein change + protein change; "
                     "hypermutators retained, not filtered")
        ax.legend(loc="lower right", fontsize=9)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "fig28_locus_comutation_pairs.png"), dpi=110)
        plt.close(fig)
        log("Wrote figures/fig28_locus_comutation_pairs.png")
    else:
        log("No co-occurring locus pairs to plot.")

    # ---- console preview --------------------------------------------------
    log("Top locus co-mutations (both = tumours with BOTH exact mutations):")
    for r in co_rows[:10]:
        print(f"    {r['GeneA']} {r['LocusA']:8s} + {r['GeneB']} {r['LocusB']:8s}  "
              f"both={r['nBoth']:3d}  (hyper={r['nBoth_hyper']}, "
              f"non-hyper={r['nBoth_nonhyper']})  OR={r['OddsRatio']:.2f}  "
              f"{r['Relationship']}")

if __name__ == "__main__":
    sys.exit(main())
