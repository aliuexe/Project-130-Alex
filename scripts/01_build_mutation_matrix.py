#!/usr/bin/env python3
"""
01_build_mutation_matrix.py
Project 130 - Colorectal Cancer (TCGA-COAD) Neoantigen Discovery Pipeline

===============================================================================
BIOLOGICAL & COMPUTATIONAL PURPOSE
===============================================================================
This script implements Sections 4 and 5 of the project assignment. It parses
the NCI Genomic Data Commons (GDC) project-level Masked Somatic Mutation MAF
(Mutation Annotation Format) file for the TCGA-COAD cohort on the GRCh38/hg38
reference assembly and builds a clean binary (0/1) mutation-by-sample matrix.

===============================================================================
SELECTION & FILTERING RULES (§4)
===============================================================================
1. Genome Assembly Assertion: Verifies that every variant is anchored to GRCh38.
2. Target Capture Scoping: Restricts analysis strictly to protein-coding genes
   (`BIOTYPE == protein_coding`).
3. Quality Assurance: Filters for high-confidence somatic variants by requiring
   `GDC_FILTER` to be empty or explicitly marked as `PASS`.
4. Variant Class Scoping: Selects nonsynonymous missense mutations
   (`Variant_Classification == Missense_Mutation`), excluding silent SNVs,
   intronic variants, and stop-gains.
5. Mutation Type Scoping: Selects single-nucleotide variants (`Variant_Type == SNP`),
   ensuring unambiguous 1-to-1 amino acid substitutions for 9-mer generation.

===============================================================================
DATA STRUCTURE & HGVS NOTATION (§5)
===============================================================================
- Distinct Mutation Key: Defined by the tuple `(Gene_Name, HGVSc, HGVSp_Short)`.
- Barcode Trimming: Trims 28-character TCGA barcodes (e.g. `TCGA-3L-AA1B-01A-...`)
  to 15-character sample-level barcodes (e.g. `TCGA-3L-AA1B-01`).
- Matrix Format: Binary presence/absence matrix where rows represent distinct
  somatic mutations, columns represent tumour samples, and values are:
    1 = mutation present in that tumour sample
    0 = mutation absent in that tumour sample

===============================================================================
INPUT & OUTPUT CONTRACTS
===============================================================================
Input:
  - `cohortMAF.2026-07-15.maf.gz` (or `/tmp/coad.maf`)

Outputs:
  - `results/01_mutation_by_sample.tsv` (Deliverable 01 matrix)
  - `results/qc_mutation_counts.txt` (QC audit summary log)
"""

import gzip
import os
import sys
import pandas as pd
import numpy as np

# =============================================================================
# FILE PATHS & RESOURCE RESOLUTION
# =============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Use uncompressed MAF if cached in /tmp for fast I/O, otherwise read repo archive
MAF = "/tmp/coad.maf" if os.path.exists("/tmp/coad.maf") else \
      os.path.join(BASE, "cohortMAF.2026-07-15.maf.gz")
OUT = os.path.join(BASE, "results", "01_mutation_by_sample.tsv")
QC = os.path.join(BASE, "results", "qc_mutation_counts.txt")

# Mandatory MAF columns required for filtering, HGVS annotation, and sample mapping
USECOLS = [
    "Hugo_Symbol", "NCBI_Build", "Chromosome", "Start_Position",
    "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification",
    "Variant_Type", "Transcript_ID", "HGVSc", "HGVSp_Short",
    "Tumor_Sample_Barcode", "BIOTYPE", "GDC_FILTER",
]

def log(msg):
    """Prints timestamped progress messages to stdout with line flushing."""
    print(f"[01] {msg}", flush=True)

def main():
    # =========================================================================
    # STEP 1: PARSE SOMATIC MUTATION MAF FILE
    # =========================================================================
    log(f"Reading MAF: {MAF}")
    comp = "gzip" if MAF.endswith(".gz") else None
    
    # Read specified columns with C engine and string dtypes to prevent memory bloat
    df = pd.read_csv(
        MAF, sep="\t", comment="#", dtype=str, usecols=USECOLS,
        compression=comp, low_memory=False, engine="c",
    )
    n_total = len(df)
    log(f"Total mutation records in MAF: {n_total}")

    # =========================================================================
    # STEP 2: VERIFY REFERENCE GENOME ASSEMBLY (RULE 5)
    # =========================================================================
    # Assert zero assembly coordinate mixing (must be strictly GRCh38 / hg38)
    builds = df["NCBI_Build"].dropna().unique().tolist()
    log(f"Reference assemblies present: {builds}")
    assert builds == ["GRCh38"], f"Unexpected assembly mix: {builds}"

    # =========================================================================
    # STEP 3: NORMALIZE TUMOUR SAMPLE BARCODES
    # =========================================================================
    # Trim 28-char aliquot barcode (TCGA-3L-AA1B-01A-11D-...) to 15-char sample ID
    df["Sample"] = df["Tumor_Sample_Barcode"].str.slice(0, 15)

    # =========================================================================
    # STEP 4: APPLY SECTION 4 SOMATIC FILTERING RULES
    # =========================================================================
    # 1. Protein-coding genes only
    f_coding = df["BIOTYPE"] == "protein_coding"
    # 2. PASS / High-confidence GDC filter calls
    f_pass = df["GDC_FILTER"].fillna("").isin(["", "PASS"])
    # 3. Nonsynonymous missense mutations only
    f_missense = df["Variant_Classification"] == "Missense_Mutation"
    # 4. Single-nucleotide variants (SNVs / SNPs)
    f_snp = df["Variant_Type"] == "SNP"

    log(f"protein_coding: {f_coding.sum()}")
    log(f"PASS/high-confidence: {f_pass.sum()}")
    log(f"Missense_Mutation: {f_missense.sum()}")
    log(f"SNP: {f_snp.sum()}")

    # Combine all 4 filters with bitwise AND
    keep = f_coding & f_pass & f_missense & f_snp
    filt = df[keep].copy()
    n_filt = len(filt)
    log(f"Records after all filters: {n_filt}")

    # Require complete HGVS coding (cDNA) and protein short notation (§5 requirement)
    filt = filt[filt["HGVSc"].notna() & filt["HGVSp_Short"].notna()]
    filt = filt[(filt["HGVSc"] != "") & (filt["HGVSp_Short"] != "")]
    log(f"Records with complete HGVS notation: {len(filt)}")

    # =========================================================================
    # STEP 5: DEFINE DISTINCT MUTATION IDENTIFIERS (§5)
    # =========================================================================
    filt["Gene_Name"] = filt["Hugo_Symbol"]
    filt["Mutation"] = filt["HGVSc"]
    filt["AminoAcid_Change"] = filt["HGVSp_Short"]

    key_cols = ["Gene_Name", "Mutation", "AminoAcid_Change"]

    # =========================================================================
    # STEP 6: BUILD MEMORY-EFFICIENT BINARY MATRIX
    # =========================================================================
    # Create composite key string for ultra-fast mapping without memory bloat
    filt["_mutkey"] = (filt["Gene_Name"] + "\t" + filt["Mutation"]
                       + "\t" + filt["AminoAcid_Change"])
    
    # Deduplicate distinct mutations to build row index
    mut_index = (filt[key_cols + ["_mutkey"]]
                 .drop_duplicates("_mutkey")
                 .reset_index(drop=True))
    mut_index["_row"] = np.arange(len(mut_index))
    row_map = dict(zip(mut_index["_mutkey"], mut_index["_row"]))

    # Sort unique tumour sample IDs to build column index
    sample_cols = sorted(filt["Sample"].unique().tolist())
    col_map = {s: i for i, s in enumerate(sample_cols)}

    # Pre-allocate uint8 matrix (0/1) for zero-copy numpy indexing
    M = np.zeros((len(mut_index), len(sample_cols)), dtype=np.uint8)
    rows = filt["_mutkey"].map(row_map).to_numpy()
    cols = filt["Sample"].map(col_map).to_numpy()
    M[rows, cols] = 1

    # Concatenate mutation metadata columns with binary sample matrix
    mat = pd.concat(
        [mut_index[key_cols].reset_index(drop=True),
         pd.DataFrame(M, columns=sample_cols)],
        axis=1,
    )

    # Deterministic sorting: most recurrent mutations first, then gene name and HGVSc
    mat["_freq"] = mat[sample_cols].sum(axis=1)
    mat = mat.sort_values(["_freq", "Gene_Name", "Mutation"],
                          ascending=[False, True, True]).drop(columns="_freq")

    # =========================================================================
    # STEP 7: WRITE DELIVERABLE 01 MATRIX TO DISK
    # =========================================================================
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    mat.to_csv(OUT, sep="\t", index=False)
    log(f"Wrote {OUT}: {mat.shape[0]} mutations x {len(sample_cols)} samples")

    # =========================================================================
    # STEP 8: WRITE QC AUDIT RECORD (§8)
    # =========================================================================
    with open(QC, "w") as fh:
        fh.write("MUTATION PROCESSING QC (Section 8)\n")
        fh.write(f"Reference assembly: GRCh38/hg38\n")
        fh.write(f"Mutation records before filtering: {n_total}\n")
        fh.write(f"Mutation records after filtering:  {n_filt}\n")
        fh.write(f"Distinct mutations (rows in matrix): {mat.shape[0]}\n")
        fh.write(f"Unique mutated genes: {mat['Gene_Name'].nunique()}\n")
        fh.write(f"Tumour samples: {len(sample_cols)}\n")
    log(f"Wrote QC summary: {QC}")

if __name__ == "__main__":
    sys.exit(main())
