#!/usr/bin/env python3
r"""
15_tmb_stratified_comutation.py
Project 130 - Colorectal Cancer (TCGA-COAD) Neoantigen Discovery Pipeline

===============================================================================
BIOLOGICAL & COMPUTATIONAL PURPOSE
===============================================================================
This script performs Tumour Mutational Burden (TMB) stratified co-mutation analysis
to eliminate TMB-driven statistical confounding (Ciriello et al., Nature Genetics).

Motivation: Hypermutated tumours (MSI-High / POLE-mutant CRC subset, ~16% of TCGA-COAD)
harbour hundreds of passenger mutations. In unstratified analyses, hypermutated
samples artificially inflate co-occurrence statistics, creating false-positive driver
co-mutation artefacts.

===============================================================================
TMB THRESHOLDING & STRATIFICATION METHODOLOGY
===============================================================================
1. Burden Calculation: Computes total missense SNV burden per tumour.
2. Hypermutator Thresholding: Cutoff $\ge 200$ missense mutations (located at the
   bimodal inflection point of the cohort TMB distribution, capturing 15.5% hypermutators).
3. Comparative Fisher's Exact Tests: Re-evaluates all 595 driver gene pairs across
   non-hypermutated tumours ($N = 495$) vs unstratified cohort ($N = 586$).
4. Categorization Verdicts:
     - `Robust co-occurrence`: Significant ($q < 0.05$) in both unstratified and TMB-controlled.
     - `TMB-artefact (lost)`: Significant in unstratified but collapses ($q \ge 0.05$) when hypermutators are excluded.
     - `Robust exclusivity`: Significant mutual exclusivity ($q < 0.05$).

===============================================================================
INPUT & OUTPUT CONTRACTS
===============================================================================
Input:
  - `results/03_integrated_mutation_expression.tsv`

Outputs:
  - `results/hypermutator_flags.tsv` (Sample TMB audit table)
  - `results/comutation_pairs_TMBcontrolled.tsv` (Comparative pair verdicts)
  - `figures/fig22_tmb_distribution.png` (Bimodal log10 TMB distribution)
  - `figures/fig23_comutation_all_vs_nonhyper.png` (All vs Non-hyper Odds Ratio scatter plot)
"""

import itertools, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =============================================================================
# FILE PATHS & RESOURCE RESOLUTION
# =============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results"); FIG = os.path.join(BASE, "figures")
INT = os.path.join(RES, "03_integrated_mutation_expression.tsv")

# Hypermutator Missense SNV Burden Threshold
HYPER_CUT = 200

# Curated Driver Panel
DRIVERS = ["APC","TP53","KRAS","PIK3CA","FBXW7","SMAD4","TCF7L2","NRAS","SMAD2",
    "CTNNB1","BRAF","SOX9","ARID1A","AMER1","ATM","KMT2C","KMT2D","ERBB2","ERBB3",
    "PTEN","ACVR2A","GNAS","BMPR1A","TGFBR2","RNF43","B2M","POLE","MSH6","CASP8",
    "ELF3","PCBP1","AXIN2","MAP2K4","CDC27"]

def log(m):
    """Prints timestamped progress messages to stdout with line flushing."""
    print("[15]", m, flush=True)

# =============================================================================
# PURE-PYTHON STATISTICAL ENGINE
# =============================================================================
def hyper_pmf(x,N,K,n):
    """Calculates Hypergeometric PMF."""
    if x<max(0,n-(N-K)) or x>min(K,n): return 0.0
    return math.comb(K,x)*math.comb(N-K,n-x)/math.comb(N,n)

def fisher_tails(a,nA,nB,N):
    """Calculates right-tail and left-tail Fisher's exact p-values."""
    lo,hi=max(0,nA+nB-N),min(nA,nB)
    right=sum(hyper_pmf(x,N,nA,nB) for x in range(a,hi+1))
    left=sum(hyper_pmf(x,N,nA,nB) for x in range(lo,a+1))
    return min(1,right),min(1,left)

def bh(pv):
    """Applies Benjamini-Hochberg FDR correction."""
    m=len(pv); order=sorted(range(m),key=lambda i:pv[i]); q=[0]*m; prev=1
    for rank,i in enumerate(reversed(order),1):
        k=m-rank+1; prev=min(prev,pv[i]*m/k); q[i]=prev
    return q

def relationship(a,nA,nB,N):
    """Computes Haldane-corrected odds ratio and Fisher p-values."""
    right,left=fisher_tails(a,nA,nB,N)
    b=nA-a;c=nB-a;d=N-nA-nB+a
    orr=((a+.5)*(d+.5))/((b+.5)*(c+.5))
    return a,orr,right,left

def run_cohort(present, genes, idx):
    """Executes pairwise Fisher's exact tests for a specified sample subset (boolean index)."""
    sub={g:present[g][idx] for g in genes}
    Ns=int(idx.sum()); cnt={g:int(sub[g].sum()) for g in genes}
    rows=[]
    for gA,gB in itertools.combinations(genes,2):
        a=int((sub[gA]&sub[gB]).sum())
        _,orr,right,left=relationship(a,cnt[gA],cnt[gB],Ns)
        rows.append([gA,gB,a,cnt[gA],cnt[gB],orr,right,left])
    q=bh([r[6] for r in rows]); qe=bh([r[7] for r in rows])
    for r,qi,qei in zip(rows,q,qe):
        r.append(qi); r.append(qei)
        r.append("Co-occurring" if r[5]>1 and qi<0.05
                 else "Mutually exclusive" if r[5]<1 and qei<0.05 else "n.s.")
    return rows, Ns, cnt

# =============================================================================
# MAIN PIPELINE EXECUTION
# =============================================================================
def main():
    # =========================================================================
    # STEP 1: CALCULATE PER-TUMOUR MISSENSE SNV BURDEN
    # =========================================================================
    with open(INT) as fh: header=fh.readline().rstrip("\n").split("\t")
    s0=next(i for i,c in enumerate(header) if c.startswith("TCGA"))
    samples=header[s0:]; N=len(samples)
    dset=set(DRIVERS); present={}; burden=np.zeros(N,dtype=int)
    with open(INT) as fh:
        fh.readline()
        for line in fh:
            p=line.rstrip("\n").split("\t")
            vec=np.fromiter((1 if v=="1" else 0 for v in p[s0:]),dtype=np.int8,count=N)
            burden+=vec
            g=p[0]
            if g in dset:
                present[g]=present[g]|vec if g in present else vec.copy()
    genes=[g for g in DRIVERS if g in present]

    # Flag hypermutated samples (burden >= 200)
    hyper = burden >= HYPER_CUT
    nonhyper = ~hyper
    log(f"N={N}; hypermutators (>= {HYPER_CUT} mut) = {int(hyper.sum())} "
        f"({100*hyper.sum()/N:.1f}%); non-hyper = {int(nonhyper.sum())}")

    # Export hypermutator audit flags
    with open(os.path.join(RES,"hypermutator_flags.tsv"),"w") as fh:
        fh.write("Sample\tMissenseSNV_burden\tClass\n")
        for s,b in sorted(zip(samples,burden), key=lambda x:-x[1]):
            fh.write(f"{s}\t{int(b)}\t{'Hypermutated' if b>=HYPER_CUT else 'Standard'}\n")

    # =========================================================================
    # STEP 2: RUN COHORT PAIRWISE STATISTICAL COMPARISON
    # =========================================================================
    idx_all=np.ones(N,dtype=bool)
    rows_all,Na,cnt_all = run_cohort(present,genes,idx_all)
    rows_nh,Nn,cnt_nh   = run_cohort(present,genes,nonhyper)
    
    da={(r[0],r[1]):r for r in rows_all}
    dn={(r[0],r[1]):r for r in rows_nh}

    cols=["GeneA","GeneB","nBoth_all","OR_all","rel_all",
          "nBoth_nonhyper","OR_nonhyper","rel_nonhyper","Verdict"]
    out=[]
    for k in da:
        ra,rn=da[k],dn[k]
        # Assign comparative verdict
        if ra[10]=="Co-occurring" and rn[10]=="Co-occurring": v="Robust co-occurrence"
        elif ra[10]=="Co-occurring" and rn[10]!="Co-occurring": v="TMB-artefact (lost)"
        elif ra[10]!="Co-occurring" and rn[10]=="Co-occurring": v="Revealed after control"
        elif ra[10]=="Mutually exclusive" and rn[10]=="Mutually exclusive": v="Robust exclusivity"
        else: v="n.s."
        out.append([k[0],k[1],ra[2],round(ra[5],2),ra[10],
                    rn[2],round(rn[5],2),rn[10],v])
    out.sort(key=lambda r:-r[2])
    with open(os.path.join(RES,"comutation_pairs_TMBcontrolled.tsv"),"w") as fh:
        fh.write(f"# TotalSamples={N}  Hypermutators={int(hyper.sum())}  "
                 f"NonHyper={int(nonhyper.sum())}  Cutoff={HYPER_CUT}\n")
        fh.write("\t".join(cols)+"\n")
        for r in out: fh.write("\t".join(str(x) for x in r)+"\n")
    log("wrote comutation_pairs_TMBcontrolled.tsv")

    # Console preview logging for key driver pairs
    log("Key pairs (all-samples  ->  non-hypermutated):")
    for k in [("KRAS","PIK3CA"),("PIK3CA","SMAD4"),("KRAS","BRAF"),
              ("APC","POLE"),("ATM","POLE"),("PIK3CA","KMT2D"),("TP53","KRAS")]:
        key=k if k in da else (k[1],k[0])
        if key in da:
            ra,rn=da[key],dn[key]
            print(f"    {key[0]:6s}+{key[1]:6s}: all OR={ra[5]:.2f} ({ra[10]:18s}) "
                  f"-> nonhyper OR={rn[5]:.2f} ({rn[10]})")
    
    def trip(g3, idx):
        v=present[g3[0]][idx]&present[g3[1]][idx]&present[g3[2]][idx]
        return int(v.sum())
    log("APC+ATM+POLE triple: all=%d  non-hyper=%d"
        % (trip(("APC","ATM","POLE"),idx_all), trip(("APC","ATM","POLE"),nonhyper)))

    # =========================================================================
    # STEP 3: RENDER FIGURE 22 — TMB DISTRIBUTION & CUTOFF INFLECTION
    # =========================================================================
    fig,ax=plt.subplots(figsize=(9,5.5))
    ax.hist(np.log10(burden+1), bins=50, color="#4477AA", edgecolor="white")
    ax.axvline(np.log10(HYPER_CUT+1), color="#EE6677", ls="--", lw=2,
               label=f"hypermutator cut-off = {HYPER_CUT} ( {int(hyper.sum())} tumours, {100*hyper.sum()/N:.0f}% )")
    ax.set_xlabel("log10(missense-SNV burden per tumour + 1)")
    ax.set_ylabel("Number of tumours")
    ax.set_title("Tumour mutational burden — MSI/POLE hypermutators form the tail (TCGA-COAD)")
    ax.legend(frameon=False); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(FIG,"fig22_tmb_distribution.png"),dpi=160)
    plt.close(fig); log("wrote fig22_tmb_distribution.png")

    # =========================================================================
    # STEP 4: RENDER FIGURE 23 — ALL VS NON-HYPER ODDS RATIO SCATTER
    # =========================================================================
    xs=[math.log2(max(da[k][5],1e-3)) for k in da]
    ys=[math.log2(max(dn[k][5],1e-3)) for k in da]
    col=[]
    for k in da:
        v=[r for r in out if (r[0],r[1])==k][0][8]
        col.append({"Robust co-occurrence":"#228833","TMB-artefact (lost)":"#EE6677",
                    "Robust exclusivity":"#4477AA","Revealed after control":"#EE7733"}.get(v,"#CCCCCC"))
    fig,ax=plt.subplots(figsize=(8,8))
    ax.scatter(xs,ys,c=col,s=45,edgecolor="white",zorder=3)
    lim=[min(xs+ys)-0.3,max(xs+ys)+0.3]
    ax.plot(lim,lim,"--",color="black",lw=1); ax.axhline(0,color="#999",lw=.8); ax.axvline(0,color="#999",lw=.8)
    
    # Annotate notable driver pairs on scatter plot
    for k in [("KRAS","PIK3CA"),("KRAS","BRAF"),("PIK3CA","SMAD4"),("APC","POLE"),("ATM","POLE"),("PIK3CA","POLE")]:
        key=k if k in da else (k[1],k[0])
        if key in da:
            ax.annotate(f"{key[0]}+{key[1]}",(math.log2(max(da[key][5],1e-3)),
                        math.log2(max(dn[key][5],1e-3))),fontsize=8,
                        textcoords="offset points",xytext=(5,3))
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("log2 odds ratio — ALL samples")
    ax.set_ylabel("log2 odds ratio — hypermutators EXCLUDED")
    ax.set_title("Driver co-mutation before vs after removing hypermutators\n"
                 "green=robust co-occur, red=TMB artefact, blue=robust exclusive")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(FIG,"fig23_comutation_all_vs_nonhyper.png"),dpi=160)
    plt.close(fig); log("wrote fig23_comutation_all_vs_nonhyper.png")

if __name__=="__main__":
    sys.exit(main())
