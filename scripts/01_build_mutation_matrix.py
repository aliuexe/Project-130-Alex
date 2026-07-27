#!/usr/bin/env python3
"""
01_build_mutation_matrix.py
Project 130 - Colorectal cancer (TCGA-COAD)

Reads the GDC project-level Masked Somatic Mutation MAF (GRCh38) and builds the
binary mutation-by-sample matrix (Deliverable 01).

Assignment sections implemented:
  Section 4 (Somatic Mutation Dataset) - filtering rules
  Section 5 (Mutation-by-Sample Matrix) - binary matrix, HGVS notation,
            one row per distinct mutation.

Filtering (Section 4, "core assignment" restrictions), all on GRCh38/hg38:
  - Protein-coding genes only              (BIOTYPE == protein_coding)
  - PASS / high-confidence variants only   (GDC_FILTER empty OR 'PASS')
  - Nonsynonymous somatic mutations        (Variant_Classification == Missense_Mutation)
  - Missense single-nucleotide variants    (Variant_Type == SNP)

Output: results/01_mutation_by_sample.tsv  (tab-delimited)
Columns: Gene_Name, Mutation (HGVSc), AminoAcid_Change (HGVSp_short), <one column per sample>
Values:  1 = mutation present in that sample, 0 = absent
"""
import gzip
import os
import sys
import pandas as pd

# ---- Paths -----------------------------------------------------------------
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Use decompressed copy if available (much faster reads), else the .gz.
MAF = "/tmp/coad.maf" if os.path.exists("/tmp/coad.maf") else \
      os.path.join(BASE, "cohortMAF.2026-07-15.maf.gz")
OUT = os.path.join(BASE, "results", "01_mutation_by_sample.tsv")
QC = os.path.join(BASE, "results", "qc_mutation_counts.txt")

# ---- Columns we need from the MAF -----------------------------------------
USECOLS = [
    "Hugo_Symbol", "NCBI_Build", "Chromosome", "Start_Position",
    "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification",
    "Variant_Type", "Transcript_ID", "HGVSc", "HGVSp_Short",
    "Tumor_Sample_Barcode", "BIOTYPE", "GDC_FILTER",
]

def log(msg):
    print(f"[01] {msg}", flush=True)

def main():
    log(f"Reading MAF: {MAF}")
    # MAF has one leading version comment? GDC cohortMAF starts directly with header.
    comp = "gzip" if MAF.endswith(".gz") else None
    df = pd.read_csv(
        MAF, sep="\t", comment="#", dtype=str, usecols=USECOLS,
        compression=comp, low_memory=False, engine="c",
    )
    n_total = len(df)
    log(f"Total mutation records in MAF: {n_total}")

    # ---- Reference-genome assembly check (Rule 5: no coordinate mixing) ----
    builds = df["NCBI_Build"].dropna().unique().tolist()
    log(f"Reference assemblies present: {builds}")
    assert builds == ["GRCh38"], f"Unexpected assembly mix: {builds}"

    # ---- Collapse tumor barcode to TCGA sample-level (first 4 fields) -------
    # Barcode e.g. TCGA-3L-AA1B-01A-11D-A36X-10  -> sample TCGA-3L-AA1B-01
    df["Sample"] = df["Tumor_Sample_Barcode"].str.slice(0, 15)

    # ---- Section 4 filters -------------------------------------------------
    f_coding = df["BIOTYPE"] == "protein_coding"
    f_pass = df["GDC_FILTER"].fillna("").isin(["", "PASS"])
    f_missense = df["Variant_Classification"] == "Missense_Mutation"
    f_snp = df["Variant_Type"] == "SNP"

    log(f"protein_coding: {f_coding.sum()}")
    log(f"PASS/high-confidence: {f_pass.sum()}")
    log(f"Missense_Mutation: {f_missense.sum()}")
    log(f"SNP: {f_snp.sum()}")

    keep = f_coding & f_pass & f_missense & f_snp
    filt = df[keep].copy()
    n_filt = len(filt)
    log(f"Records after all filters: {n_filt}")

    # Require HGVS annotation to be present (Section 5 mandates HGVS notation)
    filt = filt[filt["HGVSc"].notna() & filt["HGVSp_Short"].notna()]
    filt = filt[(filt["HGVSc"] != "") & (filt["HGVSp_Short"] != "")]
    log(f"Records with complete HGVS notation: {len(filt)}")

    # ---- Define a distinct mutation (Section 5: one row per mutation) ------
    # A mutation is a unique (Gene, HGVSc, HGVSp_Short) combination.
    filt["Gene_Name"] = filt["Hugo_Symbol"]
    filt["Mutation"] = filt["HGVSc"]
    filt["AminoAcid_Change"] = filt["HGVSp_Short"]

    key_cols = ["Gene_Name", "Mutation", "AminoAcid_Change"]

    # ---- Build binary presence/absence matrix (memory-efficient) -----------
    # Build a categorical mutation key and sample key, then fill a numpy
    # uint8 matrix by integer index. Avoids pandas pivot_table blow-up.
    import numpy as np
    filt["_mutkey"] = (filt["Gene_Name"] + "\t" + filt["Mutation"]
                       + "\t" + filt["AminoAcid_Change"])
    mut_index = (filt[key_cols + ["_mutkey"]]
                 .drop_duplicates("_mutkey")
                 .reset_index(drop=True))
    mut_index["_row"] = np.arange(len(mut_index))
    row_map = dict(zip(mut_index["_mutkey"], mut_index["_row"]))

    sample_cols = sorted(filt["Sample"].unique().tolist())
    col_map = {s: i for i, s in enumerate(sample_cols)}

    M = np.zeros((len(mut_index), len(sample_cols)), dtype=np.uint8)
    rows = filt["_mutkey"].map(row_map).to_numpy()
    cols = filt["Sample"].map(col_map).to_numpy()
    M[rows, cols] = 1

    mat = pd.concat(
        [mut_index[key_cols].reset_index(drop=True),
         pd.DataFrame(M, columns=sample_cols)],
        axis=1,
    )

    # Deterministic row ordering: most frequently mutated first, then gene name
    mat["_freq"] = mat[sample_cols].sum(axis=1)
    mat = mat.sort_values(["_freq", "Gene_Name", "Mutation"],
                          ascending=[False, True, True]).drop(columns="_freq")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    mat.to_csv(OUT, sep="\t", index=False)
    log(f"Wrote {OUT}: {mat.shape[0]} mutations x {len(sample_cols)} samples")

    # ---- QC record (Section 8) --------------------------------------------
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
