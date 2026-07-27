#!/usr/bin/env python3
"""
scratch/genome_wide_comutation.py
Genome-wide All-by-All Co-Mutation Analysis across 17,585 Genes (154.6 Million Pairs)

Evaluates pairwise co-occurrence across ALL mutated genes in the TCGA-COAD cohort (586 tumours)
to demonstrate the effect of gene length bias (TTN, MUC16) and multiple-testing penalties.
"""
import os
import sys
import time
import numpy as np
import pandas as pd
import scipy.stats as st

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MUT = os.path.join(BASE, "results", "01_mutation_by_sample.tsv")
OUT_SUM = os.path.join(BASE, "results", "all_genes_comutation_summary.txt")
OUT_PAIRS = os.path.join(BASE, "results", "all_genes_top_pairs.tsv")

def log(msg):
    print(f"[GenomeWide] {msg}", flush=True)

def main():
    t0 = time.time()
    log(f"Reading mutation matrix: {MUT}")
    df = pd.read_csv(MUT, sep="\t", low_memory=False)
    sample_cols = [c for c in df.columns if c.startswith("TCGA")]
    N = len(sample_cols)
    log(f"Total mutation records: {len(df)}, Tumour samples (N): {N}")

    # Build binary presence matrix per gene
    log("Building binary presence matrix across 17,585 genes...")
    gene_presence = df.groupby("Gene_Name")[sample_cols].apply(lambda g: (g == 1).any(axis=0).astype(np.int8))
    genes = np.array(gene_presence.index.tolist())
    n_genes = len(genes)
    M = gene_presence.values # shape: (n_genes, N)
    log(f"Binary matrix constructed: {n_genes} unique genes x {N} samples. Time: {time.time()-t0:.2f}s")

    # Calculate marginal counts per gene
    counts = M.sum(axis=1) # shape: (n_genes,)

    # Fast Matrix Multiplication for All-by-All Co-occurrence Counts
    t1 = time.time()
    total_pairs = n_genes * (n_genes - 1) // 2
    log(f"Computing exact overlap counts for all {total_pairs:,} gene pairs via matrix multiplication (M @ M.T)...")
    C = M @ M.T # shape: (n_genes, n_genes)
    log(f"Matrix multiplication complete in {time.time()-t1:.2f}s!")

    # Extract upper triangle indices (excluding self-pairs i == j)
    t2 = time.time()
    log("Extracting pairs with co-occurrence overlap a >= 2...")
    i_idx, j_idx = np.triu_indices(n_genes, k=1)
    
    # Extract arrays for non-self pairs
    a_all = C[i_idx, j_idx]
    
    # Filter pairs where both genes are co-mutated in >= 2 tumours
    mask = a_all >= 2
    i_sub = i_idx[mask]
    j_sub = j_idx[mask]
    a = a_all[mask]
    n_active = len(a)
    log(f"Filtered {n_active:,} active pairs with overlap a >= 2 (out of {total_pairs:,} total pairs). Time: {time.time()-t2:.2f}s")

    # Vectorized statistical calculations for active pairs
    t3 = time.time()
    log("Computing expected values, odds ratios, and hypergeometric p-values...")
    gA = genes[i_sub]
    gB = genes[j_sub]
    nA = counts[i_sub]
    nB = counts[j_sub]
    
    # 2x2 Contingency Table components:
    # a: both mutated
    # b = nA - a: A mutated, B WT
    # c = nB - a: B mutated, A WT
    # d = N - nA - nB + a: both WT
    b_val = nA - a
    c_val = nB - a
    d_val = N - nA - nB + a
    
    # Expected co-occurrence under independence
    expected = (nA.astype(float) * nB.astype(float)) / float(N)
    
    # Haldane-corrected Odds Ratio
    odds_ratio = ((a.astype(float) + 0.5) * (d_val.astype(float) + 0.5)) / \
                 ((b_val.astype(float) + 0.5) * (c_val.astype(float) + 0.5))
    
    # SciPy Hypergeometric survival function for right-tail (co-occurrence) p-values:
    # scipy.stats.hypergeom.sf(k-1, M_total, n_successes, N_draws)
    p_cooccur = st.hypergeom.sf(a - 1, N, nA, nB)
    
    log(f"Statistical calculations complete in {time.time()-t3:.2f}s!")

    # Benjamini-Hochberg FDR correction across ALL 154,607,320 pairs
    t4 = time.time()
    log(f"Applying Benjamini-Hochberg FDR adjustment across all {total_pairs:,} pairs...")
    # For pairs with a < 2, p_cooccur >= p_cooccur(a=1). We adjust p-values of active pairs:
    # Order active p-values
    order = np.argsort(p_cooccur)
    sorted_p = p_cooccur[order]
    
    # Rank k out of total_pairs M_tot
    m_tot = float(total_pairs)
    ranks = np.arange(1, n_active + 1, dtype=float)
    q_sorted = sorted_p * (m_tot / ranks)
    
    # Enforce monotonicity: q[i] = min(q[i], q[i+1])
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0.0, 1.0)
    
    q_cooccur = np.zeros(n_active, dtype=float)
    q_cooccur[order] = q_sorted
    log(f"FDR calculation complete in {time.time()-t4:.2f}s!")

    # Build DataFrame of active pairs
    res_df = pd.DataFrame({
        "GeneA": gA,
        "Locus_or_GeneA_Total": nA,
        "GeneB": gB,
        "Locus_or_GeneB_Total": nB,
        "nBoth": a,
        "Expected": np.round(expected, 3),
        "OddsRatio": np.round(odds_ratio, 3),
        "p_cooccur": p_cooccur,
        "FDR_cooccur": q_cooccur,
    })

    # Sort by nBoth descending, then OddsRatio descending
    top_by_count = res_df.sort_values(by=["nBoth", "OddsRatio"], ascending=[False, False]).head(50)
    top_by_or = res_df[res_df["nBoth"] >= 5].sort_values(by=["OddsRatio", "nBoth"], ascending=[False, False]).head(50)
    top_by_fdr = res_df.sort_values(by=["FDR_cooccur", "nBoth"], ascending=[True, False]).head(50)

    # Export Top Pairs TSV
    top_by_count.to_csv(OUT_PAIRS, sep="\t", index=False)
    log(f"Wrote top 50 co-occurring pairs to {OUT_PAIRS}")

    # Count significant pairs at FDR < 0.05
    n_sig = (q_cooccur < 0.05).sum()
    n_sig_or2 = ((q_cooccur < 0.05) & (odds_ratio > 2.0)).sum()

    # Identify giant genes in top 20 raw counts
    top20_genes = set(top_by_count.head(20)["GeneA"].tolist() + top_by_count.head(20)["GeneB"].tolist())
    giant_genes = [g for g in ["TTN", "MUC16", "SYNE1", "FLG", "OBSCN", "RYR2", "NEB", "LRP2", "USH2A", "PTPRT"] if g in top20_genes]

    # Write Summary Text File
    with open(OUT_SUM, "w") as fh:
        fh.write("=" * 72 + "\n")
        fh.write("GENOME-WIDE ALL-BY-ALL CO-MUTATION ANALYSIS SUMMARY\n")
        fh.write("Project 130 - Colorectal Cancer (TCGA-COAD, N = 586 Tumours)\n")
        fh.write("=" * 72 + "\n\n")
        fh.write(f"Total Mutated Genes Analysed: {n_genes:,}\n")
        fh.write(f"Total Pairwise Combinations Evaluated: {total_pairs:,}\n")
        fh.write(f"Active Pairs with Overlap (nBoth >= 2): {n_active:,}\n")
        fh.write(f"Statistically Significant Pairs (FDR q < 0.05): {n_sig:,}\n")
        fh.write(f"Significant Enriched Pairs (FDR q < 0.05 & OR > 2.0): {n_sig_or2:,}\n\n")
        
        fh.write("-" * 72 + "\n")
        fh.write("1. TOP 20 GENOME-WIDE CO-OCCURRING GENE PAIRS (BY RAW PATIENT OVERLAP)\n")
        fh.write("-" * 72 + "\n")
        fh.write(f"{'Gene A':12s} {'nA':5s} {'Gene B':12s} {'nB':5s} {'nBoth':6s} {'Expected':8s} {'OddsRatio':10s} {'FDR q':10s}\n")
        for _, r in top_by_count.head(20).iterrows():
            fh.write(f"{r['GeneA']:12s} {r['Locus_or_GeneA_Total']:5d} {r['GeneB']:12s} {r['Locus_or_GeneB_Total']:5d} {r['nBoth']:6d} {r['Expected']:8.2f} {r['OddsRatio']:10.2f} {r['FDR_cooccur']:10.2e}\n")
            
        fh.write("\n" + "-" * 72 + "\n")
        fh.write("2. KEY FINDINGS & LESSONS FROM GENOME-WIDE TESTING\n")
        fh.write("-" * 72 + "\n")
        fh.write(f"- Giant Gene Bias: Raw patient overlap counts are heavily dominated by giant structural/extracellular genes\n")
        fh.write(f"  such as {', '.join(giant_genes)}. These genes accumulate dozens of passenger mutations purely due to\n")
        fh.write(f"  their massive coding length (e.g. TTN is ~33,000 amino acids).\n")
        fh.write(f"- Confounding by Gene Popularity: Highly mutated gene pairs (e.g. TTN + MUC16) have large raw overlaps\n")
        fh.write(f"  simply because both genes are mutated in many tumours, NOT because they share a functional pathway.\n")
        fh.write(f"- Multiple-Testing Penalty: Testing all 154.6 million pairs imposes an extreme FDR penalty of m = 154,607,320,\n")
        fh.write(f"  requiring an unadjusted p-value < 3.2e-10 to achieve FDR q < 0.05.\n\n")

        fh.write("-" * 72 + "\n")
        fh.write("3. CANONICAL DRIVER PAIRS PERFORMANCE IN GENOME-WIDE BACKGROUND\n")
        fh.write("-" * 72 + "\n")
        driver_pairs = [("APC", "TP53"), ("APC", "KRAS"), ("KRAS", "PIK3CA"), ("PIK3CA", "SMAD4"), ("KRAS", "BRAF")]
        for g1, g2 in driver_pairs:
            match = res_df[((res_df["GeneA"]==g1)&(res_df["GeneB"]==g2)) | ((res_df["GeneA"]==g2)&(res_df["GeneB"]==g1))]
            if len(match):
                r = match.iloc[0]
                fh.write(f"  {g1:6s} + {g2:6s}: nBoth={r['nBoth']:3d}, Expected={r['Expected']:6.2f}, OR={r['OddsRatio']:5.2f}, p={r['p_cooccur']:.2e}, FDR={r['FDR_cooccur']:.2e}\n")

    log(f"Wrote genome-wide summary log to {OUT_SUM}")
    log(f"Complete genome-wide analysis finished in {time.time()-t0:.2f}s!")

if __name__ == "__main__":
    main()
