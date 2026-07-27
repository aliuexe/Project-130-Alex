#!/usr/bin/env python3
"""
scripts/22_class2_binding_and_coverage.py
Project 130 - Colorectal Cancer (TCGA-COAD, N = 586 Tumours)
Extension: MHC Class II (15-mer) Binding Prediction & Dual Class I + Class II Vaccine Cocktail Optimization

===============================================================================
BIOLOGICAL & COMPUTATIONAL PURPOSE
===============================================================================
CD4+ helper T-cells recognize 15-mer peptides presented on MHC Class II molecules
(HLA-DRB1*15:01, *07:01). CD4+ help is essential for priming and sustaining
CD8+ cytotoxic T-cell responses against cancer neoantigens.

This script:
1. Predicts Class II binding affinity (IC50 in nM) and presentation probability
   for all 15-mer mutant and wild-type peptides in `results/peptides_all.tsv`
   against HLA-DRB1*15:01 and HLA-DRB1*07:01 using a simplified anchor-position
   PSSM (positions P1, P4, P6, P9 of the 9-mer binding core).
   LIMITATION: This is a 4-anchor-position model, not a full SMM-align/NetMHCIIpan
   predictor. It captures the dominant MHC-II anchor preferences but does not score
   non-anchor positions (P2, P3, P5, P7, P8). This is a stated limitation.
2. Identifies Practical Class II Neoantigens (Mutant IC50 < 500 nM, differential
   agretopicity or TCR contact loop novelty, TPM >= 10, Clonal, Recurrence >= 2).
3. Expands the Greedy Set-Cover algorithm to design a Dual Class I + Class II
   Vaccine Cocktail with **class-aware coverage tracking**, measuring:
   - Overall coverage (>= 1 epitope of any class per tumour)
   - Genuine CD4+/CD8+ synergy (>= 1 Class I AND >= 1 Class II epitope per tumour)
"""

import os
import sys
import time
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
FIG_DIR = os.path.join(BASE, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

PEP_ALL = os.path.join(RES, "peptides_all.tsv")
INT_FILE = os.path.join(RES, "03_integrated_mutation_expression.tsv")
CLONAL_FILE = os.path.join(RES, "mutation_clonality.tsv")
PRAC_CLASS1 = os.path.join(RES, "practical_neoantigens.tsv")

OUT_CLASS2_PRAC = os.path.join(RES, "practical_class2_neoantigens.tsv")
OUT_DUAL_CURVE = os.path.join(RES, "dual_vaccine_coverage_curve.tsv")
OUT_DUAL_SUMM = os.path.join(RES, "dual_vaccine_coverage_summary.txt")
OUT_FIG = os.path.join(FIG_DIR, "22_dual_class1_class2_vaccine_coverage.png")

def log(msg):
    print(f"[22 Class II] {msg}", flush=True)

# =============================================================================
# ANCHOR-POSITION PSSM FOR HLA-DRB1*15:01 & *07:01
# =============================================================================
# Simplified 4-anchor model scoring P1, P4, P6, P9 of the 9-mer binding core.
# Values represent log-energy contributions; negative = favorable binding.
# LIMITATION: Positions P2, P3, P5, P7, P8 are not scored.
AAS = "ACDEFGHIKLMNPQRSTVWY"

# HLA-DRB1*15:01 prefers aromatic/large hydrophobic at P1, aliphatic at P4, small/basic at P6, P9
DRB1_1501_WEIGHTS = {
    1: {"F": -2.1, "Y": -2.0, "W": -2.2, "L": -1.8, "I": -1.7, "V": -1.4, "M": -1.6, "A": -0.5},
    4: {"L": -1.5, "I": -1.4, "V": -1.3, "M": -1.4, "A": -1.0, "F": -1.2, "Y": -1.1},
    6: {"S": -1.2, "T": -1.1, "N": -1.0, "Q": -1.0, "K": -1.3, "R": -1.4, "A": -0.9, "G": -0.7},
    9: {"A": -1.4, "S": -1.3, "T": -1.2, "G": -1.1, "V": -1.2, "L": -1.0, "I": -0.8}
}

# HLA-DRB1*07:01 prefers aromatic/hydrophobic at P1, aromatic/aliphatic at P4, small/polar at P6, P9
DRB1_0701_WEIGHTS = {
    1: {"F": -2.2, "Y": -2.1, "W": -2.0, "L": -1.7, "I": -1.6, "V": -1.3, "M": -1.5},
    4: {"Y": -1.6, "F": -1.5, "W": -1.4, "L": -1.4, "I": -1.3, "V": -1.2, "M": -1.3},
    6: {"S": -1.3, "T": -1.2, "A": -1.1, "G": -1.0, "N": -0.9},
    9: {"A": -1.5, "S": -1.4, "T": -1.3, "V": -1.2, "L": -1.1, "I": -1.0, "G": -1.0}
}

def score_9mer_core(core, weights_dict, base_ic50=8000.0):
    """
    Scores a 9-mer core sequence against an HLA-DRB1 anchor-position weight matrix.
    Returns predicted IC50 in nM.
    NOTE: Only anchor positions (P1, P4, P6, P9) contribute to the score.
    """
    log_affinity = math.log10(base_ic50)
    for pos, aa_weights in weights_dict.items():
        aa = core[pos - 1]
        w = aa_weights.get(aa, 0.2)  # default slight penalty for non-preferred residues
        log_affinity += (w * 0.45)
    ic50 = 10.0 ** max(0.5, min(4.5, log_affinity))
    return ic50

def predict_class2_15mer(peptide, allele):
    """
    Evaluates all 7 possible 9-mer cores in a 15-mer peptide.
    Returns (best_ic50_nM, best_core, presentation_prob, best_offset).
    """
    if len(peptide) != 15 or any(c not in AAS for c in peptide):
        return 10000.0, peptide[:9], 0.0, 0
    
    weights = DRB1_1501_WEIGHTS if "15:01" in allele else DRB1_0701_WEIGHTS
    best_ic50 = 10000.0
    best_core = peptide[:9]
    best_offset = 0
    for i in range(7):
        core = peptide[i:i+9]
        ic50 = score_9mer_core(core, weights)
        if ic50 < best_ic50:
            best_ic50 = ic50
            best_core = core
            best_offset = i
            
    # Logistic presentation mapping where IC50 = 500 nM -> 0.50 presentation prob
    pres_prob = 1.0 / (1.0 + (best_ic50 / 500.0) ** 2)
    return best_ic50, best_core, pres_prob, best_offset

def predict_class2_15mer_constrained(peptide, allele, mutpos_1based):
    """
    Like predict_class2_15mer, but ONLY considers 9-mer cores where the mutation
    (at mutpos_1based in the 15-mer) falls within the 9-mer core (positions 1-9).
    Returns (best_ic50_nM, best_core, presentation_prob, best_offset, core_mut_pos).
    Returns None if no valid core contains the mutation.
    """
    if len(peptide) != 15 or any(c not in AAS for c in peptide):
        return None

    weights = DRB1_1501_WEIGHTS if "15:01" in allele else DRB1_0701_WEIGHTS
    best_ic50 = 10000.0
    best_core = None
    best_offset = 0
    best_cmp = 0

    for i in range(7):
        core_mut_pos = mutpos_1based - i  # 1-based position in core
        # Only consider cores where the mutation is within the 9-mer
        if core_mut_pos < 1 or core_mut_pos > 9:
            continue
        core = peptide[i:i+9]
        ic50 = score_9mer_core(core, weights)
        if ic50 < best_ic50:
            best_ic50 = ic50
            best_core = core
            best_offset = i
            best_cmp = core_mut_pos

    if best_core is None:
        return None

    pres_prob = 1.0 / (1.0 + (best_ic50 / 500.0) ** 2)
    return best_ic50, best_core, pres_prob, best_offset, best_cmp

def load_cohort_metadata():
    """Loads TPM, Clonality, and Sample-presence arrays for all mutations."""
    log(f"Loading clonality data: {CLONAL_FILE}")
    clonal_map = {}
    with open(CLONAL_FILE) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 6:
                clonal_map[(p[0], p[1])] = p[5]
                
    log(f"Loading integrated mutation & TPM data: {INT_FILE}")
    tpm_map = {}
    sampsets = {}
    with open(INT_FILE) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        s0 = next(i for i, c in enumerate(header) if c.startswith("TCGA"))
        N = len(header) - s0
        for line in fh:
            p = line.rstrip("\n").split("\t")
            gene = p[0]
            change = p[2]
            key = (gene, change)
            try:
                tpm = float(p[3]) if p[3] != "NA" else 0.0
            except:
                tpm = 0.0
            tpm_map[key] = max(tpm_map.get(key, 0.0), tpm)
            v = np.fromiter((c == "1" for c in p[s0:]), dtype=bool, count=N)
            if key in sampsets:
                sampsets[key] |= v
            else:
                sampsets[key] = v
    log(f"Total cohort tumours N = {N}. Loaded {len(sampsets):,} mutation sample vectors.")
    return N, tpm_map, clonal_map, sampsets

def run_greedy_cover(candidates, sampsets, N, max_steps=35):
    """
    Runs greedy set-cover optimization on a candidate pool (single MHC class).
    Returns trajectory of (gene, change, cov_ge1, pct_ge1, cov_ge2, pct_ge2).
    """
    keys = [(r["GeneName"], r["ProteinChange"]) for r in candidates]
    covered_ge1 = np.zeros(N, dtype=bool)
    hit_count = np.zeros(N, dtype=int)
    remaining = set(keys)
    order = []
    
    while remaining and len(order) < max_steps:
        best_k = None
        best_gain = -1
        for k in remaining:
            gain = int((sampsets[k] & ~covered_ge1).sum())
            if gain > best_gain:
                best_gain = gain
                best_k = k
        if best_gain <= 0:
            # If all ge1 covered, pick candidate that adds most ge2 hits
            for k in remaining:
                gain2 = int((sampsets[k] & (hit_count == 1)).sum())
                if gain2 > best_gain:
                    best_gain = gain2
                    best_k = k
            if best_gain <= 0:
                break
        remaining.discard(best_k)
        covered_ge1 |= sampsets[best_k]
        hit_count += sampsets[best_k].astype(int)
        c1 = int(covered_ge1.sum())
        c2 = int((hit_count >= 2).sum())
        order.append((best_k[0], best_k[1], c1, 100.0 * c1 / N, c2, 100.0 * c2 / N))
    return order

def run_dual_greedy_cover(class1_candidates, class2_candidates, sampsets, N, max_steps=30):
    """
    Class-aware greedy set-cover for the dual Class I + Class II vaccine cocktail.
    
    Tracks Class I and Class II hits SEPARATELY per tumour so that genuine
    CD4+/CD8+ synergy (>= 1 Class I AND >= 1 Class II epitope in a tumour)
    can be measured.
    
    At each step, selects the candidate from the combined pool that maximises
    marginal gain in overall coverage (>= 1 epitope of any class).
    
    Returns trajectory of tuples:
        (gene, change, mhc_class, cov_any, pct_any, cov_ge2, pct_ge2,
         cov_synergy, pct_synergy, n_c1_selected, n_c2_selected)
    where:
        cov_any = tumours with >= 1 epitope (any class)
        cov_ge2 = tumours with >= 2 distinct mutations (any class)
        cov_synergy = tumours with >= 1 Class I AND >= 1 Class II (genuine CD4+/CD8+)
    """
    # Build keyed candidate pools, tagging each with its MHC class
    c1_keys = {(r["GeneName"], r["ProteinChange"]) for r in class1_candidates}
    c2_keys = {(r["GeneName"], r["ProteinChange"]) for r in class2_candidates}
    
    # Class II-only mutations (not targetable by Class I)
    c2_exclusive = c2_keys - c1_keys
    # Class I-only mutations (not targetable by Class II)
    c1_exclusive = c1_keys - c2_keys
    # Mutations targetable by both classes
    both = c1_keys & c2_keys
    
    log(f"  Class I-exclusive keys: {len(c1_exclusive)}, Class II-exclusive keys: {len(c2_exclusive)}, Both: {len(both)}")
    
    # Build candidate entries: each entry is (gene, change, mhc_class)
    # For mutations available in both classes, create TWO entries
    remaining = set()
    for k in c1_exclusive:
        remaining.add((k[0], k[1], "ClassI"))
    for k in c2_exclusive:
        remaining.add((k[0], k[1], "ClassII"))
    for k in both:
        remaining.add((k[0], k[1], "ClassI"))
        remaining.add((k[0], k[1], "ClassII"))
    
    # Per-tumour tracking arrays
    c1_hits = np.zeros(N, dtype=int)   # Class I epitope hits per tumour
    c2_hits = np.zeros(N, dtype=int)   # Class II epitope hits per tumour
    covered_any = np.zeros(N, dtype=bool)
    hit_count = np.zeros(N, dtype=int)
    
    order = []
    n_c1_selected = 0
    n_c2_selected = 0
    
    while remaining and len(order) < max_steps:
        best_entry = None
        best_gain = -1
        
        for entry in remaining:
            gene, change, mhc_class = entry
            k = (gene, change)
            if k not in sampsets:
                continue
            gain = int((sampsets[k] & ~covered_any).sum())
            if gain > best_gain:
                best_gain = gain
                best_entry = entry
        
        if best_gain <= 0:
            # If no new tumours to cover, maximise ge2 (secondary objective)
            for entry in remaining:
                gene, change, mhc_class = entry
                k = (gene, change)
                if k not in sampsets:
                    continue
                gain2 = int((sampsets[k] & (hit_count == 1)).sum())
                if gain2 > best_gain:
                    best_gain = gain2
                    best_entry = entry
            if best_gain <= 0:
                break
        
        if best_entry is None:
            break
        
        gene, change, mhc_class = best_entry
        k = (gene, change)
        
        # Remove this entry and (if it exists) the other-class entry for the same mutation
        remaining.discard(best_entry)
        other_class = "ClassII" if mhc_class == "ClassI" else "ClassI"
        remaining.discard((gene, change, other_class))
        
        # Update coverage arrays
        svec = sampsets[k]
        covered_any |= svec
        hit_count += svec.astype(int)
        if mhc_class == "ClassI":
            c1_hits += svec.astype(int)
            n_c1_selected += 1
        else:
            c2_hits += svec.astype(int)
            n_c2_selected += 1
        
        # Compute metrics
        cov_any = int(covered_any.sum())
        cov_ge2 = int((hit_count >= 2).sum())
        # Genuine CD4+/CD8+ synergy: tumour has >= 1 Class I AND >= 1 Class II
        cov_synergy = int(((c1_hits >= 1) & (c2_hits >= 1)).sum())
        
        order.append((
            gene, change, mhc_class,
            cov_any, 100.0 * cov_any / N,
            cov_ge2, 100.0 * cov_ge2 / N,
            cov_synergy, 100.0 * cov_synergy / N,
            n_c1_selected, n_c2_selected
        ))
    
    return order

def main():
    t0 = time.time()
    N, tpm_map, clonal_map, sampsets = load_cohort_metadata()
    
    # =========================================================================
    # STEP 1: PREDICT MHC CLASS II (15-MER) BINDING & PRESENTATION
    # =========================================================================
    log(f"Scanning 15-mers in {PEP_ALL}...")
    alleles = ["HLA-DRB1*15:01", "HLA-DRB1*07:01"]
    
    # Store mut and wt 15-mers per (gene, change), grouped by sliding window index
    # Key: (gene, change), Value: list of (peptide, mutpos_1based)
    mut_15mers = {}
    wt_15mers = {}
    with open(PEP_ALL) as fh:
        header = fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 14 and p[10] == "15":
                gene, change, mutpos, pep, ptype = p[0], p[8], int(p[11]), p[12], p[13]
                key = (gene, change)
                if ptype == "Mutant":
                    mut_15mers.setdefault(key, []).append((pep, mutpos))
                elif ptype == "WildType":
                    wt_15mers.setdefault(key, []).append((pep, mutpos))
                    
    log(f"Extracted 15-mers for {len(mut_15mers):,} distinct mutations.")
    
    class2_candidates = []
    n_outside_core = 0
    n_evaluated = 0
    
    for key, peplist in mut_15mers.items():
        gene, change = key
        # Check recurrence >= 2 in cohort
        svec = sampsets.get(key)
        if svec is None or svec.sum() < 2:
            continue
        # Check Clonal & TPM >= 10
        if clonal_map.get(key) != "Clonal":
            continue
        tpm = tpm_map.get(key, 0.0)
        if tpm < 10.0:
            continue
            
        freq = int(svec.sum())
        # Build WT lookup: map (peptide_sequence) -> mutpos for aligned comparison
        wt_by_mutpos = {}
        for wpep, wmutpos in wt_15mers.get(key, []):
            wt_by_mutpos.setdefault(wmutpos, []).append(wpep)
        
        # Evaluate against DRB1*15:01 and DRB1*07:01
        for pep, mutpos in peplist:
            for allele in alleles:
                n_evaluated += 1
                
                # Use constrained prediction: only consider cores where mutation is IN the core
                result = predict_class2_15mer_constrained(pep, allele, mutpos)
                if result is None:
                    n_outside_core += 1
                    continue
                
                mut_ic50, mut_core, mut_el, offset, core_mut_pos = result
                if mut_ic50 >= 500.0:
                    continue
                
                # Compare against the ALIGNED WT 15-mer (same sliding window)
                # The WT 15-mer with the same mutpos corresponds to the same protein window
                wt_ic50_aligned = 10000.0
                wt_el_aligned = 0.0
                for wpep in wt_by_mutpos.get(mutpos, []):
                    wt_result = predict_class2_15mer_constrained(wpep, allele, mutpos)
                    if wt_result is not None:
                        wic50, wcore, wel, woff, wcmp = wt_result
                        if wic50 < wt_ic50_aligned:
                            wt_ic50_aligned = wic50
                            wt_el_aligned = wel
                
                # Classify neoantigenicity mechanism based on where the mutation sits in the core
                # P1, P4, P6, P9 = anchor positions (affect MHC binding)
                # P2, P3, P5, P7, P8 = TCR-facing positions (affect T-cell recognition)
                is_tcr_contact = (core_mut_pos in {2, 3, 5, 7, 8})
                is_anchor = (core_mut_pos in {1, 4, 6, 9})
                
                # Agretopicity: mutation creates de novo binding or substantially improves it
                is_agretopic = (wt_ic50_aligned >= 500.0 or mut_ic50 < wt_ic50_aligned * 0.5)
                
                if is_anchor and is_agretopic:
                    mech = "Agretopic_Anchor"
                elif is_tcr_contact:
                    mech = "TCR_Contact_Loop"
                else:
                    # Mutation is at an anchor position but didn't change binding enough
                    continue
                    
                class2_candidates.append({
                    "GeneName": gene,
                    "ProteinChange": change,
                    "Peptide": pep,
                    "Core9mer": mut_core,
                    "CoreMutPos": core_mut_pos,
                    "HLAAllele": allele,
                    "Mutant_IC50": round(mut_ic50, 1),
                    "WT_IC50": round(wt_ic50_aligned, 1),
                    "Mutant_EL": round(mut_el, 3),
                    "WT_EL": round(wt_el_aligned, 3),
                    "ClassIIMechanism": mech,
                    "GeneLevelTPM": round(tpm, 1),
                    "MutationFrequency": freq,
                    "TumoursCovered": freq
                })
    
    log(f"Evaluated {n_evaluated:,} (peptide, allele) pairs; {n_outside_core:,} skipped (mutation outside 9-mer core).")
                    
    # Deduplicate keeping best binder per (GeneName, ProteinChange, HLAAllele)
    best_c2 = {}
    for c in class2_candidates:
        k = (c["GeneName"], c["ProteinChange"], c["HLAAllele"])
        if k not in best_c2 or c["Mutant_IC50"] < best_c2[k]["Mutant_IC50"]:
            best_c2[k] = c
            
    prac_class2 = sorted(list(best_c2.values()), key=lambda r: (-r["TumoursCovered"], r["Mutant_IC50"]))
    log(f"Discovered {len(prac_class2):,} Practical Class II (15-mer) Neoantigens!")
    
    # Count mechanisms
    mech_counts = {}
    for r in prac_class2:
        mech_counts[r["ClassIIMechanism"]] = mech_counts.get(r["ClassIIMechanism"], 0) + 1
    for mech, count in sorted(mech_counts.items()):
        log(f"  Mechanism {mech}: {count} candidates")
    
    # Write Practical Class II Neoantigen TSV
    cols_c2 = ["GeneName", "ProteinChange", "Peptide", "Core9mer", "CoreMutPos", "HLAAllele",
               "Mutant_IC50", "WT_IC50", "Mutant_EL", "WT_EL", "ClassIIMechanism",
               "GeneLevelTPM", "MutationFrequency", "TumoursCovered"]
    with open(OUT_CLASS2_PRAC, "w") as fh:
        fh.write(f"# TotalSamples(N)\t{N}\n")
        fh.write("\t".join(cols_c2) + "\n")
        for r in prac_class2:
            fh.write("\t".join(str(r[c]) for c in cols_c2) + "\n")
    log(f"Practical Class II database written to {OUT_CLASS2_PRAC}")

    # =========================================================================
    # STEP 2: LOAD CLASS I PRACTICAL NEOANTIGENS
    # =========================================================================
    prac_class1 = []
    with open(PRAC_CLASS1) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("GeneName"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 10:
                prac_class1.append({
                    "GeneName": p[0],
                    "ProteinChange": p[1],
                    "Peptide": p[2],
                    "HLAAllele": p[3],
                    "Mutant_EL": float(p[4]),
                    "WT_EL": float(p[5]),
                    "GeneLevelTPM": float(p[7]),
                    "MutationFrequency": int(p[8]),
                    "TumoursCovered": int(p[9])
                })
    log(f"Loaded {len(prac_class1):,} Practical Class I (9-mer) Neoantigens.")

    # =========================================================================
    # STEP 3: RUN GREEDY SET-COVER FOR CLASS I ONLY, CLASS II ONLY, AND DUAL
    # =========================================================================
    log("Running Greedy Set-Cover for Class I Only (9-mers)...")
    traj_c1 = run_greedy_cover(prac_class1, sampsets, N, max_steps=30)
    
    log("Running Greedy Set-Cover for Class II Only (15-mers)...")
    traj_c2 = run_greedy_cover(prac_class2, sampsets, N, max_steps=30)
    
    log("Running Class-Aware Dual Greedy Set-Cover (Class I + II with synergy tracking)...")
    traj_dual = run_dual_greedy_cover(prac_class1, prac_class2, sampsets, N, max_steps=30)
    
    # Write Dual Vaccine Coverage Curve TSV
    with open(OUT_DUAL_CURVE, "w") as fh:
        fh.write(f"# TotalSamples(N)\t{N}\n")
        fh.write("Strategy\tStep\tGeneName\tProteinChange\tTumours_ge1\tPct_ge1\tTumours_ge2\tPct_ge2\n")
        for idx, (g, pc, c1, p1, c2, p2) in enumerate(traj_c1, 1):
            fh.write(f"Class1_Only\t{idx}\t{g}\t{pc}\t{c1}\t{p1:.1f}\t{c2}\t{p2:.1f}\n")
        for idx, (g, pc, c1, p1, c2, p2) in enumerate(traj_c2, 1):
            fh.write(f"Class2_Only\t{idx}\t{g}\t{pc}\t{c1}\t{p1:.1f}\t{c2}\t{p2:.1f}\n")
        # Dual trajectory with additional columns
        fh.write("\n# DUAL COCKTAIL TRAJECTORY (class-aware)\n")
        fh.write("Strategy\tStep\tGeneName\tProteinChange\tMHC_Class\tTumours_ge1\tPct_ge1\t"
                 "Tumours_ge2\tPct_ge2\tCD4_CD8_Synergy\tPct_Synergy\tN_ClassI\tN_ClassII\n")
        for idx, entry in enumerate(traj_dual, 1):
            g, pc, mc, ca, pa, c2, p2, cs, ps, nc1, nc2 = entry
            fh.write(f"Dual_ClassAware\t{idx}\t{g}\t{pc}\t{mc}\t{ca}\t{pa:.1f}\t"
                     f"{c2}\t{p2:.1f}\t{cs}\t{ps:.1f}\t{nc1}\t{nc2}\n")
    log(f"Dual vaccine coverage curves exported to {OUT_DUAL_CURVE}")

    # Write detailed comparative summary report
    with open(OUT_DUAL_SUMM, "w") as fh:
        fh.write("=" * 90 + "\n")
        fh.write("DUAL CLASS I + CLASS II VACCINE COCKTAIL COMPARATIVE EFFICACY REPORT\n")
        fh.write(f"Cohort Size: N = {N} Colorectal Cancer Tumours\n")
        fh.write("=" * 90 + "\n\n")
        
        fh.write("1. TOP 15 DISCOVERED PRACTICAL CLASS II (15-MER) NEOANTIGENS\n")
        fh.write("-" * 90 + "\n")
        fh.write(f"{'Gene':10s} {'Change':12s} {'15-mer Peptide':17s} {'Core9mer':10s} {'CMP':4s} "
                 f"{'HLA-DRB1':14s} {'IC50(nM)':9s} {'WT_IC50':8s} {'Mech':18s} {'Freq':5s}\n")
        for r in prac_class2[:15]:
            fh.write(f"{r['GeneName']:10s} {r['ProteinChange']:12s} {r['Peptide']:17s} "
                     f"{r['Core9mer']:10s} P{r['CoreMutPos']:<3d} {r['HLAAllele']:14s} "
                     f"{r['Mutant_IC50']:<9.1f} {r['WT_IC50']:<8.1f} "
                     f"{r['ClassIIMechanism']:18s} {r['MutationFrequency']:5d}\n")
        fh.write("\n")
        
        fh.write("2. VACCINE COCKTAIL POPULATION COVERAGE AT 10, 20, AND 30 EPITOPES\n")
        fh.write("-" * 90 + "\n")
        fh.write(f"{'Strategy':25s} | {'--- 10 Epitopes ---':20s} | {'--- 20 Epitopes ---':20s} | {'--- 30 Epitopes ---':20s}\n")
        fh.write(f"{'':25s} | {'ge1 (%)':9s} {'ge2 (%)':10s} | {'ge1 (%)':9s} {'ge2 (%)':10s} | {'ge1 (%)':9s} {'ge2 (%)':10s}\n")
        
        for name, traj in [("Class I Only (9-mer)", traj_c1), ("Class II Only (15-mer)", traj_c2)]:
            s10 = traj[min(9, len(traj)-1)]
            s20 = traj[min(19, len(traj)-1)]
            s30 = traj[min(29, len(traj)-1)]
            fh.write(f"{name:25s} | {s10[3]:6.1f}%   {s10[5]:6.1f}%    | {s20[3]:6.1f}%   {s20[5]:6.1f}%    | {s30[3]:6.1f}%   {s30[5]:6.1f}%\n")
        
        # Dual trajectory uses different tuple format
        if traj_dual:
            s10d = traj_dual[min(9, len(traj_dual)-1)]
            s20d = traj_dual[min(19, len(traj_dual)-1)]
            s30d = traj_dual[min(29, len(traj_dual)-1)]
            fh.write(f"{'Dual Optimal Cocktail':25s} | {s10d[4]:6.1f}%   {s10d[6]:6.1f}%    | {s20d[4]:6.1f}%   {s20d[6]:6.1f}%    | {s30d[4]:6.1f}%   {s30d[6]:6.1f}%\n")
        fh.write("\n")
        
        fh.write("3. GENUINE CD4+ / CD8+ SYNERGY (>= 1 Class I AND >= 1 Class II per tumour)\n")
        fh.write("-" * 90 + "\n")
        if traj_dual:
            fh.write(f"{'Step':5s} {'Gene':12s} {'Change':12s} {'Class':8s} | {'Any':6s} | {'ge2':6s} | {'Synergy':8s} | {'#C1':4s} {'#C2':4s}\n")
            for idx, entry in enumerate(traj_dual[:30], 1):
                g, pc, mc, ca, pa, c2, p2, cs, ps, nc1, nc2 = entry
                fh.write(f"{idx:5d} {g:12s} {pc:12s} {mc:8s} | {pa:5.1f}% | {p2:5.1f}% | {ps:6.1f}%  | {nc1:4d} {nc2:4d}\n")
        fh.write("\n")
        
        fh.write("4. WHY CLASS-AWARE TRACKING MATTERS\n")
        fh.write("-" * 90 + "\n")
        fh.write("a) Immunological Synergy: CD4+ T-cell help (via MHC Class II 15-mers) is required\n")
        fh.write("   to prime, sustain, and prevent exhaustion of CD8+ cytotoxic T-cells (via MHC Class I 9-mers).\n")
        fh.write("b) Class-Aware vs Naive Tracking: The naive 'ge2' metric counts tumours hit by >= 2 distinct\n")
        fh.write("   mutations, but does not distinguish whether both hits are Class I (no CD4+ help) or one\n")
        fh.write("   is Class I and one is Class II (genuine CD4+/CD8+ synergy). The 'Synergy' metric\n")
        fh.write("   specifically counts tumours with >= 1 Class I AND >= 1 Class II epitope.\n")
        fh.write("c) Stated Limitation: The MHC-II binding predictions use a simplified 4-anchor-position\n")
        fh.write("   PSSM (P1, P4, P6, P9). Non-anchor positions (P2, P3, P5, P7, P8) are not scored.\n")
        fh.write("   A full SMM-align or NetMHCIIpan predictor would provide more accurate binding estimates.\n")
    log(f"Summary report written to {OUT_DUAL_SUMM}")

    # =========================================================================
    # STEP 4: PLOT FIGURE 22 — COMPARATIVE DUAL COCKTAIL COVERAGE CURVE
    # =========================================================================
    log("Rendering comparative coverage figure...")
    fig, axes = plt.subplots(1, 3, figsize=(22, 6.5))
    
    xs1 = list(range(1, len(traj_c1)+1))
    xs2 = list(range(1, len(traj_c2)+1))
    xsd = list(range(1, len(traj_dual)+1))
    
    # Panel A: >= 1 epitope per tumour (any class)
    ax = axes[0]
    ax.plot(xs1, [x[3] for x in traj_c1], "-o", color="#457B9D", ms=4, label="Class I Only (9-mer, CD8+)")
    ax.plot(xs2, [x[3] for x in traj_c2], "-s", color="#F4A261", ms=4, label="Class II Only (15-mer, CD4+)")
    ax.plot(xsd, [x[4] for x in traj_dual], "-^", color="#E63946", ms=5, linewidth=2.5, label="Dual Class I + II")
    ax.set_xlabel("Number of Neoantigens in Vaccine Cocktail", fontsize=11)
    ax.set_ylabel(f"% of Tumours Covered (≥ 1 Epitope, N={N})", fontsize=11)
    ax.set_title("A. Population Coverage (≥ 1 Epitope)", fontsize=13, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="lower right", fontsize=9)
    ax.set_ylim([0, 100])

    # Panel B: >= 2 distinct mutations per tumour (any class)
    ax = axes[1]
    ax.plot(xs1, [x[5] for x in traj_c1], "-o", color="#457B9D", ms=4, label="Class I Only (9-mer, CD8+)")
    ax.plot(xs2, [x[5] for x in traj_c2], "-s", color="#F4A261", ms=4, label="Class II Only (15-mer, CD4+)")
    ax.plot(xsd, [x[6] for x in traj_dual], "-^", color="#E63946", ms=5, linewidth=2.5, label="Dual Class I + II")
    ax.set_xlabel("Number of Neoantigens in Vaccine Cocktail", fontsize=11)
    ax.set_ylabel(f"% of Tumours Covered (≥ 2 Epitopes, N={N})", fontsize=11)
    ax.set_title("B. Multi-Epitope Coverage (≥ 2 Mutations)", fontsize=13, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="lower right", fontsize=9)
    ax.set_ylim([0, 100])
    
    # Panel C: Genuine CD4+/CD8+ Synergy (>= 1 Class I AND >= 1 Class II)
    ax = axes[2]
    ax.plot(xsd, [x[8] for x in traj_dual], "-^", color="#2A9D8F", ms=5, linewidth=2.5,
            label="CD4+/CD8+ Synergy\n(≥1 Class I AND ≥1 Class II)")
    ax.plot(xsd, [x[6] for x in traj_dual], "--", color="#8D99AE", ms=3, alpha=0.7,
            label="≥2 Mutations (any class)")
    ax.set_xlabel("Number of Neoantigens in Vaccine Cocktail", fontsize=11)
    ax.set_ylabel(f"% of Tumours with CD4+/CD8+ Synergy (N={N})", fontsize=11)
    ax.set_title("C. Genuine CD4+/CD8+ Synergy", fontsize=13, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="lower right", fontsize=9)
    ax.set_ylim([0, 100])
    
    # Annotate Panel C with Class I/II selection counts
    if traj_dual:
        last = traj_dual[-1]
        ax.annotate(f"{last[9]} Class I + {last[10]} Class II selected",
                    xy=(len(traj_dual), last[8]), fontsize=9,
                    xytext=(-80, 20), textcoords="offset points",
                    arrowprops=dict(arrowstyle="->", color="#333"))
    
    plt.suptitle("Off-the-Shelf Colorectal Cancer Vaccine: Class I + Class II Synergy\n"
                 f"(Class-Aware CD4+ / CD8+ Coverage Tracking across N={N} Tumours)",
                 fontsize=14, fontweight="bold", y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.90])
    plt.savefig(OUT_FIG, dpi=300)
    plt.close()
    log(f"Comparative figure saved to {OUT_FIG}")
    log(f"All tasks completed in {time.time()-t0:.2f} seconds.")

if __name__ == "__main__":
    main()
