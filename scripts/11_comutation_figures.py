#!/usr/bin/env python3
"""
11_comutation_figures.py
Project 130 - Colorectal Cancer (TCGA-COAD) Neoantigen Discovery Pipeline

===============================================================================
BIOLOGICAL & COMPUTATIONAL PURPOSE
===============================================================================
This script renders the core co-mutation extension figures (Figures 14, 15, and 16).
It parses pairwise driver co-mutation statistics generated in Script 10 (`results/comutation_driver_pairs.tsv`)
and Deliverable 03 to visualize driver co-occurrence, mutual exclusivity, and
per-tumour driver mutation burden.

===============================================================================
FIGURE INDEX & VISUALIZATIONS PRODUCED
===============================================================================
  1. Figure 14 (`fig14_comutation_heatmap.png`): Symmetric $16 \times 16$ log2 odds ratio
     heatmap of top driver pairs. Statistically significant co-occurrence ($q < 0.05$)
     is annotated with `+` symbols; mutual exclusivity ($q < 0.05$) with `–` symbols.
  2. Figure 15 (`fig15_comutation_pairs.png`): Side-by-side bar plots detailing top
     co-occurring driver pairs (by sample overlap) and mutually exclusive pairs.
  3. Figure 16 (`fig16_driver_load.png`): Per-tumour co-mutation burden histogram
     displaying the frequency of tumours carrying 0, 1, 2, 3+ driver mutations.

===============================================================================
INPUT & OUTPUT CONTRACTS
===============================================================================
Inputs:
  - `results/comutation_driver_pairs.tsv`
  - `results/03_integrated_mutation_expression.tsv`

Outputs:
  - `figures/fig14_comutation_heatmap.png`
  - `figures/fig15_comutation_pairs.png`
  - `figures/fig16_driver_load.png`
"""

import os, sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =============================================================================
# FILE PATHS & RESOURCE RESOLUTION
# =============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results"); FIG = os.path.join(BASE, "figures")
PAIRS = os.path.join(RES, "comutation_driver_pairs.tsv")
INT = os.path.join(RES, "03_integrated_mutation_expression.tsv")
os.makedirs(FIG, exist_ok=True)

# Consistent Styling Options
plt.rcParams.update({"figure.dpi":150,"savefig.dpi":160,"font.size":12,
    "axes.titlesize":15,"axes.titleweight":"bold","axes.spines.top":False,
    "axes.spines.right":False})
BLUE, RED, TEAL, GREEN, ORANGE, PURPLE = "#4477AA","#EE6677","#66CCEE","#228833","#EE7733","#AA3377"

# =============================================================================
# HELPER DATA LOADING
# =============================================================================
def load_pairs():
    """Parses pairwise co-mutation statistics TSV into dictionary records."""
    rows=[]
    with open(PAIRS) as fh:
        hdr=fh.readline().rstrip("\n").split("\t"); ix={c:i for i,c in enumerate(hdr)}
        for line in fh:
            p=line.rstrip("\n").split("\t")
            rows.append(dict(A=p[ix["GeneA"]],B=p[ix["GeneB"]],nA=int(p[ix["nA"]]),
                nB=int(p[ix["nB"]]),nBoth=int(p[ix["nBoth"]]),OR=float(p[ix["OddsRatio"]]),
                fco=float(p[ix["FDR_cooccur"]]),fex=float(p[ix["FDR_exclusive"]]),
                rel=p[ix["Relationship"]]))
    return rows

# =============================================================================
# FIGURE RENDERING FUNCTIONS
# =============================================================================
def fig_heatmap(rows):
    """Figure 14: Log2 Odds Ratio Co-mutation Heatmap."""
    freq={}
    for r in rows:
        freq[r["A"]]=r["nA"]; freq[r["B"]]=r["nB"]
    top=[g for g,_ in sorted(freq.items(), key=lambda kv:-kv[1])][:16]
    idx={g:i for i,g in enumerate(top)}; n=len(top)
    M=np.full((n,n), np.nan); sig=np.zeros((n,n))
    for r in rows:
        if r["A"] in idx and r["B"] in idx:
            i,j=idx[r["A"]],idx[r["B"]]
            v=math.log2(max(r["OR"],1e-3))
            M[i,j]=M[j,i]=v
            s=1 if (r["rel"]=="Co-occurring") else (-1 if r["rel"]=="Mutually exclusive" else 0)
            sig[i,j]=sig[j,i]=s
    fig,ax=plt.subplots(figsize=(9.5,8))
    im=ax.imshow(M, cmap="RdBu_r", vmin=-3, vmax=3)
    ax.set_xticks(range(n)); ax.set_xticklabels(top, rotation=90, fontsize=10)
    ax.set_yticks(range(n)); ax.set_yticklabels(top, fontsize=10)
    for i in range(n):
        for j in range(n):
            if sig[i,j]==1: ax.text(j,i,"+",ha="center",va="center",color="black",fontsize=11,fontweight="bold")
            elif sig[i,j]==-1: ax.text(j,i,"–",ha="center",va="center",color="black",fontsize=12,fontweight="bold")
    cb=fig.colorbar(im,shrink=0.7); cb.set_label("log2 odds ratio  (red=co-occur, blue=exclusive)")
    ax.set_title("Driver co-mutation in colorectal cancer (TCGA-COAD)\n+ = co-occurring, – = mutually exclusive (FDR<0.05)")
    fig.tight_layout(); fig.savefig(os.path.join(FIG,"fig14_comutation_heatmap.png")); plt.close(fig)
    print("[11] wrote fig14_comutation_heatmap.png")

def fig_pairs(rows):
    """Figure 15: Top Co-occurring and Mutually Exclusive Pairs."""
    co=sorted([r for r in rows if r["rel"]=="Co-occurring"], key=lambda r:-r["nBoth"])[:8]
    ex=sorted([r for r in rows if r["rel"]=="Mutually exclusive"], key=lambda r:r["fex"])[:6]
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(13,5.5))
    lab=[f"{r['A']}+{r['B']}" for r in co][::-1]; val=[r["nBoth"] for r in co][::-1]
    ax1.barh(range(len(co)), val, color=RED)
    for i,r in enumerate(co[::-1]): ax1.text(r["nBoth"]+0.5,i,f"OR {r['OR']:.1f}",va="center",fontsize=9)
    ax1.set_yticks(range(len(co))); ax1.set_yticklabels(lab); ax1.set_xlabel("Tumours carrying BOTH")
    ax1.set_title("Top co-occurring driver pairs")
    lab2=[f"{r['A']}/{r['B']}" for r in ex][::-1]; val2=[r["OR"] for r in ex][::-1]
    ax2.barh(range(len(ex)), val2, color=BLUE)
    ax2.axvline(1, color="black", lw=1, ls="--")
    for i,r in enumerate(ex[::-1]): ax2.text(r["OR"]+0.01,i,f"FDR {r['fex']:.0e}",va="center",fontsize=9)
    ax2.set_yticks(range(len(ex))); ax2.set_yticklabels(lab2); ax2.set_xlabel("Odds ratio (<1 = exclusive)")
    ax2.set_title("Mutually exclusive driver pairs")
    fig.suptitle("Driver co-mutation vs mutual exclusivity (TCGA-COAD)", fontsize=15, fontweight="bold")
    fig.tight_layout(); fig.savefig(os.path.join(FIG,"fig15_comutation_pairs.png")); plt.close(fig)
    print("[11] wrote fig15_comutation_pairs.png")

def fig_load():
    """Figure 16: Per-tumour Driver Burden Histogram."""
    with open(INT) as fh: header=fh.readline().rstrip("\n").split("\t")
    s0=next(i for i,c in enumerate(header) if c.startswith("TCGA"))
    samples=header[s0:]; N=len(samples)
    DRIVERS=set(["APC","TP53","KRAS","PIK3CA","FBXW7","SMAD4","TCF7L2","NRAS","SMAD2",
        "CTNNB1","BRAF","SOX9","ARID1A","AMER1","ATM","KMT2C","KMT2D","ERBB2","ERBB3",
        "PTEN","ACVR2A","GNAS","BMPR1A","TGFBR2","RNF43","B2M","POLE","MSH6","CASP8",
        "ELF3","PCBP1","AXIN2","MAP2K4","CDC27"])
    pres={}
    with open(INT) as fh:
        fh.readline()
        for line in fh:
            p=line.rstrip("\n").split("\t"); g=p[0]
            if g not in DRIVERS: continue
            v=np.fromiter((1 if x=="1" else 0 for x in p[s0:]),dtype=np.int8,count=N)
            pres[g]=pres[g]|v if g in pres else v.copy()
    load=np.zeros(N,dtype=int)
    for g in pres: load+=pres[g]
    maxk=int(load.max()); counts=[int((load==k).sum()) for k in range(maxk+1)]
    colors=["#BBBBBB" if k<2 else (ORANGE if k==2 else RED) for k in range(maxk+1)]
    fig,ax=plt.subplots(figsize=(11,5.5))
    ax.bar(range(maxk+1), counts, color=colors, edgecolor="white")
    n2=int((load>=2).sum()); n3=int((load>=3).sum())
    ax.set_xlabel("Number of co-mutated driver genes in a tumour")
    ax.set_ylabel("Number of tumours")
    ax.set_title("Per-tumour driver co-mutation load (TCGA-COAD)")
    ax.text(0.98,0.9,f"≥2 drivers: {n2} tumours ({100*n2/N:.0f}%)\n≥3 drivers: {n3} tumours ({100*n3/N:.0f}%)",
            transform=ax.transAxes, ha="right", fontsize=12,
            bbox=dict(boxstyle="round",fc="white",ec="#999999"))
    ax.text(1.4, max(counts)*0.5, "singles / none", color="#888888", fontsize=10, rotation=90, va="center")
    fig.tight_layout(); fig.savefig(os.path.join(FIG,"fig16_driver_load.png")); plt.close(fig)
    print(f"[11] wrote fig16_driver_load.png  (>=2:{n2}, >=3:{n3})")

# =============================================================================
# MAIN PIPELINE EXECUTION
# =============================================================================
def main():
    rows=load_pairs(); fig_heatmap(rows); fig_pairs(rows); fig_load()

if __name__=="__main__":
    sys.exit(main())
