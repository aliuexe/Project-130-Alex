#!/usr/bin/env python3
"""
scratch/variant_level_comutation.py
Variant-Level (Locus-Specific) Co-Mutation & Hotspot Analysis on TTN and Genome-Wide

1. Are there specific recurrent amino-acid mutations on TTN (e.g. TTN p.X123Y) that occur across multiple tumours?
2. Do any specific TTN variants significantly co-occur with canonical driver variants (KRAS p.G12D, TP53 p.R248Q, PIK3CA p.E542K, etc.)?
3. What are the top co-occurring VARIANT PAIRS across the entire TCGA-COAD cohort (N=586 and N=495)?
"""
import os
import sys
import time
import numpy as np
import pandas as pd
import scipy.stats as st

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MUT = os.path.join(BASE, "results", "01_mutation_by_sample.tsv")
OUT_VAR = os.path.join(BASE, "results", "variant_level_comutation_summary.txt")

def log(msg):
    print(f"[VariantLevel] {msg}", flush=True)

def main():
    t0 = time.time()
    log(f"Reading mutation matrix: {MUT}")
    df = pd.read_csv(MUT, sep="\t", low_memory=False)
    sample_cols = [c for c in df.columns if c.startswith("TCGA")]
    N = len(sample_cols)
    log(f"Total distinct variant records: {len(df):,}, Tumour samples (N): {N}")

    # Create distinct Variant ID: "Gene_Name p.AminoAcid_Change"
    df["Variant_ID"] = df["Gene_Name"] + " " + df["AminoAcid_Change"]
    
    # Calculate recurrence count per variant across the cohort
    sample_mat = df[sample_cols].values.astype(np.int8)
    variant_counts = sample_mat.sum(axis=1)
    df["Recurrence"] = variant_counts

    # 1. Inspect TTN Variants specifically
    ttn_df = df[df["Gene_Name"] == "TTN"].copy()
    log(f"Total distinct missense SNVs on TTN: {len(ttn_df):,} across {N} samples.")
    ttn_sorted = ttn_df.sort_values(by="Recurrence", ascending=False)
    
    # How many TTN variants occur in >= 2, >= 3, >= 4, >= 5 tumours?
    ttn_ge2 = (ttn_sorted["Recurrence"] >= 2).sum()
    ttn_ge3 = (ttn_sorted["Recurrence"] >= 3).sum()
    ttn_ge4 = (ttn_sorted["Recurrence"] >= 4).sum()
    ttn_ge5 = (ttn_sorted["Recurrence"] >= 5).sum()
    log(f"TTN recurrence counts: >=2={ttn_ge2}, >=3={ttn_ge3}, >=4={ttn_ge4}, >=5={ttn_ge5}")

    # 2. Compare against Top Driver Hotspots across the whole cohort
    top_overall = df.sort_values(by="Recurrence", ascending=False).head(20)

    # 3. Variant-Level All-by-All Co-Mutation on Recurrent Variants (Recurrence >= 3)
    rec_mask = variant_counts >= 3
    df_rec = df[rec_mask].reset_index(drop=True)
    n_rec = len(df_rec)
    log(f"Total recurrent variants (>= 3 tumours across entire cohort): {n_rec:,}")
    
    V_ids = df_rec["Variant_ID"].values
    V_counts = df_rec["Recurrence"].values
    M_rec = df_rec[sample_cols].values.astype(np.int8) # shape: (n_rec, 586)

    # Calculate TMB mask for standard tumours (N=495)
    sample_burdens = sample_mat.sum(axis=0)
    nonhyper_mask = sample_burdens < 200
    N_std = int(nonhyper_mask.sum())
    M_rec_std = M_rec[:, nonhyper_mask]
    V_counts_std = M_rec_std.sum(axis=1)

    # Matrix multiplication for variant co-occurrence
    C_all = M_rec @ M_rec.T # shape: (n_rec, n_rec)
    C_std = M_rec_std @ M_rec_std.T

    i_idx, j_idx = np.triu_indices(n_rec, k=1)
    a_all = C_all[i_idx, j_idx]
    a_std = C_std[i_idx, j_idx]

    # Filter for pairs where both variants co-occur in >= 2 tumours in All Tumours
    active_mask = a_all >= 2
    i_sub = i_idx[active_mask]
    j_sub = j_idx[active_mask]
    a_sub_all = a_all[active_mask]
    a_sub_std = a_std[active_mask]

    # Compute p-values and Odds Ratios for active variant pairs in N=586
    nA_all = V_counts[i_sub]
    nB_all = V_counts[j_sub]
    b_all = nA_all - a_sub_all
    c_all = nB_all - a_sub_all
    d_all = N - nA_all - nB_all + a_sub_all
    
    expected_all = (nA_all.astype(float) * nB_all.astype(float)) / float(N)
    or_all = ((a_sub_all.astype(float) + 0.5) * (d_all.astype(float) + 0.5)) / \
             ((b_all.astype(float) + 0.5) * (c_all.astype(float) + 0.5))
    p_all = st.hypergeom.sf(a_sub_all - 1, N, nA_all, nB_all)

    # Compute for Standard Tumours N=495
    nA_std = V_counts_std[i_sub]
    nB_std = V_counts_std[j_sub]
    b_std = nA_std - a_sub_std
    c_std = nB_std - a_sub_std
    d_std = N_std - nA_std - nB_std + a_sub_std
    
    expected_std = (nA_std.astype(float) * nB_std.astype(float)) / float(N_std)
    or_std = ((a_sub_std.astype(float) + 0.5) * (d_std.astype(float) + 0.5)) / \
             ((b_std.astype(float) + 0.5) * (c_std.astype(float) + 0.5))
    p_std = st.hypergeom.sf(np.maximum(0, a_sub_std - 1), N_std, nA_std, nB_std)

    var_df = pd.DataFrame({
        "VariantA": V_ids[i_sub],
        "nA_586": nA_all,
        "VariantB": V_ids[j_sub],
        "nB_586": nB_all,
        "nBoth_586": a_sub_all,
        "Expected_586": expected_all,
        "OddsRatio_586": or_all,
        "p_val_586": p_all,
        "nA_495": nA_std,
        "nB_495": nB_std,
        "nBoth_495": a_sub_std,
        "OddsRatio_495": or_std,
        "p_val_495": p_std
    })

    # Sort top co-occurring variant pairs in N=586 by p-value
    top_var_586 = var_df.sort_values(by=["p_val_586", "nBoth_586"], ascending=[True, False]).head(25)
    # Sort top co-occurring variant pairs in N=495 by p-value
    top_var_495 = var_df[var_df["nBoth_495"] >= 2].sort_values(by=["p_val_495", "nBoth_495"], ascending=[True, False]).head(25)

    # Check if ANY TTN variant is in top co-occurring variant pairs
    ttn_var_pairs = var_df[var_df["VariantA"].str.startswith("TTN ") | var_df["VariantB"].str.startswith("TTN ")].sort_values(by="nBoth_586", ascending=False)

    with open(OUT_VAR, "w") as fh:
        fh.write("=" * 80 + "\n")
        fh.write("VARIANT-LEVEL (LOCUS-SPECIFIC) CO-MUTATION & TTN HOTSPOT ANALYSIS REPORT\n")
        fh.write("Project 130 - Colorectal Cancer (TCGA-COAD, N = 586 Tumours)\n")
        fh.write("=" * 80 + "\n\n")

        fh.write("-" * 80 + "\n")
        fh.write("1. TOP 15 MOST RECURRENT INDIVIDUAL MISSENSE VARIANTS IN ENTIRE COHORT\n")
        fh.write("-" * 80 + "\n")
        fh.write(f"{'Rank':5s} {'Variant ID':20s} {'Tumour Count (N=586)':22s} {'Percentage':15s}\n")
        for idx, (_, r) in enumerate(top_overall.head(15).iterrows(), 1):
            fh.write(f"{idx:<5d} {r['Variant_ID']:20s} {r['Recurrence']:<22d} {100*r['Recurrence']/N:.2f}%\n")
        fh.write("\n")

        fh.write("-" * 80 + "\n")
        fh.write("2. TOP 15 MOST RECURRENT INDIVIDUAL MISSENSE VARIANTS ON TTN\n")
        fh.write("-" * 80 + "\n")
        fh.write(f"{'Rank':5s} {'TTN Variant ID':25s} {'Tumour Count (N=586)':22s} {'Percentage':15s}\n")
        for idx, (_, r) in enumerate(ttn_sorted.head(15).iterrows(), 1):
            fh.write(f"{idx:<5d} {r['Variant_ID']:25s} {r['Recurrence']:<22d} {100*r['Recurrence']/N:.2f}%\n")
        fh.write("\n")

        fh.write("-" * 80 + "\n")
        fh.write("3. TOP 15 CO-OCCURRING VARIANT PAIRS IN ALL TUMOURS (N=586)\n")
        fh.write("-" * 80 + "\n")
        fh.write(f"{'Variant A':18s} {'nA':4s} {'Variant B':18s} {'nB':4s} {'nBoth':6s} {'Exp':6s} {'OR_586':8s} {'p_val_586':10s} {'nBoth_495':10s}\n")
        for _, r in top_var_586.head(15).iterrows():
            fh.write(f"{r['VariantA']:18s} {r['nA_586']:4d} {r['VariantB']:18s} {r['nB_586']:4d} {r['nBoth_586']:6d} {r['Expected_586']:6.2f} {r['OddsRatio_586']:8.2f} {r['p_val_586']:10.2e} {r['nBoth_495']:10d}\n")
        fh.write("\n")

        fh.write("-" * 80 + "\n")
        fh.write("4. TOP 15 CO-OCCURRING VARIANT PAIRS IN STANDARD TUMOURS (N=495, TMB < 200)\n")
        fh.write("-" * 80 + "\n")
        fh.write(f"{'Variant A':18s} {'nA':4s} {'Variant B':18s} {'nB':4s} {'nBoth_495':10s} {'OR_495':8s} {'p_val_495':10s}\n")
        for _, r in top_var_495.head(15).iterrows():
            fh.write(f"{r['VariantA']:18s} {r['nA_495']:4d} {r['VariantB']:18s} {r['nB_495']:4d} {r['nBoth_495']:10d} {r['OddsRatio_495']:8.2f} {r['p_val_495']:10.2e}\n")
        fh.write("\n")

        fh.write("-" * 80 + "\n")
        fh.write("5. MOST CO-OCCURRING VARIANT PAIRS INVOLVING A TTN VARIANT\n")
        fh.write("-" * 80 + "\n")
        if len(ttn_var_pairs) == 0:
            fh.write("No TTN variant pairs found with co-occurrence nBoth >= 2 among recurrent variants.\n")
        else:
            fh.write(f"{'TTN Variant Pair':38s} {'nBoth_586':10s} {'OR_586':8s} {'p_val_586':10s} {'nBoth_495':10s} {'p_val_495':10s}\n")
            for _, r in ttn_var_pairs.head(15).iterrows():
                pair_str = f"{r['VariantA']} + {r['VariantB']}"
                fh.write(f"{pair_str:38s} {r['nBoth_586']:10d} {r['OddsRatio_586']:8.2f} {r['p_val_586']:10.2e} {r['nBoth_495']:10d} {r['p_val_495']:10.2e}\n")

    log(f"Variant-level analysis summary written to {OUT_VAR}")

if __name__ == "__main__":
    main()
