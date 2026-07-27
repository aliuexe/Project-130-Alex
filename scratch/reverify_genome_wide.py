#!/usr/bin/env python3
"""
scratch/reverify_genome_wide.py
Unbiased Genome-Wide Co-Mutation Re-Verification:
1. What are the TRUE top pairs by Statistical Significance (lowest p-value) across all 17,585 genes in N=586?
2. What are the TRUE top pairs by Statistical Significance in Standard Tumours (N=495, non-hypermutators)?
3. Why did TTN/MUC16/SYNE1 show OR ~ 2.5-4.5 in N=586? Is it TMB confounding (hypermutator co-occurrence)?
"""
import os
import sys
import time
import numpy as np
import pandas as pd
import scipy.stats as st

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MUT = os.path.join(BASE, "results", "01_mutation_by_sample.tsv")
OUT_REV = os.path.join(BASE, "results", "reverify_genome_wide_summary.txt")

def log(msg):
    print(f"[ReVerify] {msg}", flush=True)

def analyze_cohort(M, genes, sample_mask, label):
    N = int(sample_mask.sum())
    M_sub = M[:, sample_mask]
    counts = M_sub.sum(axis=1)
    n_genes = len(genes)
    
    # We only care about genes mutated in at least 2 samples in this sub-cohort
    valid_genes_mask = counts >= 2
    g_valid = genes[valid_genes_mask]
    c_valid = counts[valid_genes_mask]
    M_valid = M_sub[valid_genes_mask, :]
    n_v = len(g_valid)
    
    log(f"[{label}] N={N} samples. Valid genes (>=2 mut): {n_v:,}")
    t0 = time.time()
    C = M_valid @ M_valid.T
    i_idx, j_idx = np.triu_indices(n_v, k=1)
    a_all = C[i_idx, j_idx]
    
    # Filter pairs with a >= 3 to focus on meaningful overlaps
    mask = a_all >= 3
    i_sub = i_idx[mask]
    j_sub = j_idx[mask]
    a = a_all[mask]
    n_active = len(a)
    log(f"[{label}] Active pairs (a >= 3): {n_active:,} out of {n_v*(n_v-1)//2:,}")
    
    gA = g_valid[i_sub]
    gB = g_valid[j_sub]
    nA = c_valid[i_sub]
    nB = c_valid[j_sub]
    
    b_val = nA - a
    c_val = nB - a
    d_val = N - nA - nB + a
    
    expected = (nA.astype(float) * nB.astype(float)) / float(N)
    odds_ratio = ((a.astype(float) + 0.5) * (d_val.astype(float) + 0.5)) / \
                 ((b_val.astype(float) + 0.5) * (c_val.astype(float) + 0.5))
    
    p_val = st.hypergeom.sf(a - 1, N, nA, nB)
    
    # FDR across total possible pairs in this valid gene set
    total_pairs = n_v * (n_v - 1) // 2
    order = np.argsort(p_val)
    sorted_p = p_val[order]
    ranks = np.arange(1, n_active + 1, dtype=float)
    q_sorted = sorted_p * (float(total_pairs) / ranks)
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0.0, 1.0)
    
    q_val = np.zeros(n_active, dtype=float)
    q_val[order] = q_sorted
    
    df_res = pd.DataFrame({
        "GeneA": gA,
        "nA": nA,
        "GeneB": gB,
        "nB": nB,
        "nBoth": a,
        "Expected": expected,
        "OddsRatio": odds_ratio,
        "p_value": p_val,
        "FDR_q": q_val
    })
    return df_res

def main():
    log(f"Loading mutation matrix: {MUT}")
    df = pd.read_csv(MUT, sep="\t", low_memory=False)
    sample_cols = [c for c in df.columns if c.startswith("TCGA")]
    N_total = len(sample_cols)
    
    gene_presence = df.groupby("Gene_Name")[sample_cols].apply(lambda g: (g == 1).any(axis=0).astype(np.int8))
    genes = np.array(gene_presence.index.tolist())
    M = gene_presence.values # shape: (17585, 586)
    
    # Calculate per-sample TMB (missense SNVs)
    sample_burdens = M.sum(axis=0)
    hyper_mask = sample_burdens >= 200
    nonhyper_mask = ~hyper_mask
    log(f"Total samples: {N_total}. Hypermutators (>=200): {hyper_mask.sum()}. Standard (<200): {nonhyper_mask.sum()}")
    
    # 1. Analyze All Tumours (N=586)
    df_all = analyze_cohort(M, genes, np.ones(N_total, dtype=bool), "All_Tumours_N586")
    
    # 2. Analyze Standard Tumours (N=495)
    df_std = analyze_cohort(M, genes, nonhyper_mask, "Standard_Tumours_N495")
    
    with open(OUT_REV, "w") as fh:
        fh.write("=" * 80 + "\n")
        fh.write("UNBIASED GENOME-WIDE CO-MUTATION RE-VERIFICATION REPORT\n")
        fh.write("=" * 80 + "\n\n")
        
        # A. Top 20 by Statistical Significance in ALL Tumours (N=586)
        fh.write("-" * 80 + "\n")
        fh.write("A. TOP 20 MOST STATISTICALLY SIGNIFICANT PAIRS IN ALL TUMOURS (N=586, LOWEST P-VALUE)\n")
        fh.write("-" * 80 + "\n")
        fh.write(f"{'Gene A':10s} {'nA':5s} {'Gene B':10s} {'nB':5s} {'nBoth':6s} {'Expected':8s} {'OddsRatio':10s} {'p-value':11s} {'FDR q':11s}\n")
        top_p_all = df_all.sort_values(by=["p_value", "OddsRatio"], ascending=[True, False]).head(20)
        for _, r in top_p_all.iterrows():
            fh.write(f"{r['GeneA']:10s} {r['nA']:5d} {r['GeneB']:10s} {r['nB']:5d} {r['nBoth']:6d} {r['Expected']:8.2f} {r['OddsRatio']:10.2f} {r['p_value']:11.2e} {r['FDR_q']:11.2e}\n")
        
        # B. What happens to TTN + MUC16, TTN + SYNE1, TTN + OBSCN in N=586 vs N=495?
        fh.write("\n" + "-" * 80 + "\n")
        fh.write("B. TTN / MUC16 / SYNE1 / FAT4 / OBSCN BEHAVIOUR: ALL (N=586) vs STANDARD (N=495)\n")
        fh.write("-" * 80 + "\n")
        test_pairs = [
            ("TTN", "MUC16"), ("TTN", "SYNE1"), ("TTN", "FAT4"), ("TTN", "OBSCN"), ("TTN", "RYR2"),
            ("MUC16", "SYNE1"), ("MUC16", "FAT4"), ("SYNE1", "FAT4"), ("KRAS", "PIK3CA"), ("PIK3CA", "SMAD4")
        ]
        fh.write(f"{'Pair':20s} | {'--- N=586 (All) ---':30s} | {'--- N=495 (Standard) ---':30s}\n")
        fh.write(f"{'GeneA + GeneB':20s} | {'nBoth':6s} {'OR':7s} {'p-val':9s} {'FDR q':8s} | {'nBoth':6s} {'OR':7s} {'p-val':9s} {'FDR q':8s}\n")
        for g1, g2 in test_pairs:
            m_all = df_all[((df_all["GeneA"]==g1)&(df_all["GeneB"]==g2)) | ((df_all["GeneA"]==g2)&(df_all["GeneB"]==g1))]
            m_std = df_std[((df_std["GeneA"]==g1)&(df_std["GeneB"]==g2)) | ((df_std["GeneA"]==g2)&(df_std["GeneB"]==g1))]
            
            str_all = f"{int(m_all.iloc[0]['nBoth']):5d} {m_all.iloc[0]['OddsRatio']:7.2f} {m_all.iloc[0]['p_value']:9.1e} {m_all.iloc[0]['FDR_q']:8.1e}" if len(m_all) else "  N/A    N/A       N/A      N/A"
            str_std = f"{int(m_std.iloc[0]['nBoth']):5d} {m_std.iloc[0]['OddsRatio']:7.2f} {m_std.iloc[0]['p_value']:9.1e} {m_std.iloc[0]['FDR_q']:8.1e}" if len(m_std) else "  N/A    N/A       N/A      N/A"
            fh.write(f"{g1+' + '+g2:20s} | {str_all:30s} | {str_std:30s}\n")
            
        # C. Top 20 by Statistical Significance in STANDARD Tumours (N=495, lowest p-value)
        fh.write("\n" + "-" * 80 + "\n")
        fh.write("C. TOP 20 MOST STATISTICALLY SIGNIFICANT PAIRS IN STANDARD TUMOURS (N=495, NO HYPERMUTATORS)\n")
        fh.write("-" * 80 + "\n")
        fh.write(f"{'Gene A':10s} {'nA':5s} {'Gene B':10s} {'nB':5s} {'nBoth':6s} {'Expected':8s} {'OddsRatio':10s} {'p-value':11s} {'FDR q':11s}\n")
        top_p_std = df_std.sort_values(by=["p_value", "OddsRatio"], ascending=[True, False]).head(20)
        for _, r in top_p_std.iterrows():
            fh.write(f"{r['GeneA']:10s} {r['nA']:5d} {r['GeneB']:10s} {r['nB']:5d} {r['nBoth']:6d} {r['Expected']:8.2f} {r['OddsRatio']:10.2f} {r['p_value']:11.2e} {r['FDR_q']:11.2e}\n")

    log(f"Re-verification report written to {OUT_REV}")

if __name__ == "__main__":
    main()
