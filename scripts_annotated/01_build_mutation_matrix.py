#!/usr/bin/env python3
# =============================================================================
# 01_build_mutation_matrix.py   (ANNOTATED teaching copy)
# =============================================================================
# WHAT THIS SCRIPT DOES, IN ONE SENTENCE:
#   It reads the big table of DNA mutations found in colorectal tumours and
#   turns it into a tidy grid ("matrix") that says, for each mutation, which
#   patients' tumours carry it (1 = yes, 0 = no).
#
# WHY WE NEED IT (biology):
#   Cancer arises when cells accumulate "somatic mutations" (DNA changes that
#   happen during a person's life, only in the tumour, not inherited). We want
#   a clean, computer-readable summary of which mutations occur in which
#   tumours. That grid is "Deliverable 01" and everything downstream builds on it.
#
# INPUT : a MAF file (Mutation Annotation Format) = a giant tab-separated table,
#         one row per (mutation, tumour-sample) observation, downloaded from the
#         NCI Genomic Data Commons (GDC) for the TCGA-COAD colorectal cohort.
# OUTPUT: results/01_mutation_by_sample.tsv  (the grid described above)
#
# A NOTE ON THE LANGUAGE:
#   This is Python. Lines starting with "#" are comments (notes for humans;
#   the computer ignores them). "pandas" is a Python library for working with
#   tables (think: a programmable spreadsheet). A "DataFrame" (df) is one table.
# =============================================================================

import gzip     # lets us read ".gz" compressed files without unzipping first
import os       # tools for building file paths that work on any computer
import sys      # lets the script exit cleanly
import pandas as pd   # the table/spreadsheet library; we nickname it "pd"

# ---- Where are the files? ---------------------------------------------------
# BASE = the project folder (one level up from this "scripts" folder).
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# The MAF is large. If a fast decompressed copy exists we use it; otherwise we
# read the compressed ".maf.gz" directly. Both contain identical data.
MAF = "/tmp/coad.maf" if os.path.exists("/tmp/coad.maf") else \
      os.path.join(BASE, "cohortMAF.2026-07-15.maf.gz")
OUT = os.path.join(BASE, "results", "01_mutation_by_sample.tsv")   # main output
QC = os.path.join(BASE, "results", "qc_mutation_counts.txt")       # QC log

# ---- Which columns of the MAF do we actually need? --------------------------
# The MAF has ~130 columns; reading only these keeps memory use low. Each name
# is a real column in the GDC MAF:
USECOLS = [
    "Hugo_Symbol",             # the gene's name, e.g. KRAS
    "NCBI_Build",              # which reference genome (should all be GRCh38)
    "Chromosome",              # which chromosome, e.g. chr12
    "Start_Position",          # where on the chromosome the mutation sits
    "Reference_Allele",        # the normal DNA letter(s) at that spot
    "Tumor_Seq_Allele2",       # the mutated DNA letter(s)
    "Variant_Classification",  # the type of change, e.g. Missense_Mutation
    "Variant_Type",            # SNP (single-letter change), insertion, deletion…
    "Transcript_ID",           # which version of the gene (transcript) was used
    "HGVSc",                   # standard notation for the DNA change, e.g. c.35G>A
    "HGVSp_Short",             # standard notation for the protein change, e.g. p.G12D
    "Tumor_Sample_Barcode",    # the tumour-sample ID (which patient/sample)
    "BIOTYPE",                 # is this a protein-coding gene? (we require yes)
    "GDC_FILTER",              # quality flags from GDC (empty = passed all filters)
]

def log(msg):
    # A tiny helper that prints progress messages tagged with "[01]".
    print(f"[01] {msg}", flush=True)

def main():
    log(f"Reading MAF: {MAF}")
    # If the file ends in ".gz" tell pandas to decompress it on the fly.
    comp = "gzip" if MAF.endswith(".gz") else None
    # Read the table. dtype=str means "treat every value as text for now"
    # (safer: it stops pandas guessing types wrongly). comment="#" skips any
    # header lines that begin with "#".
    df = pd.read_csv(
        MAF, sep="\t", comment="#", dtype=str, usecols=USECOLS,
        compression=comp, low_memory=False, engine="c",
    )
    n_total = len(df)                       # how many mutation rows we started with
    log(f"Total mutation records in MAF: {n_total}")

    # ---- Safety check: make sure every mutation uses the SAME reference genome.
    # Genomic positions only make sense relative to one reference genome build.
    # Mixing GRCh37 and GRCh38 coordinates would be a serious error (Rule 5),
    # so we assert (demand) that the only build present is GRCh38. If not, the
    # script stops immediately.
    builds = df["NCBI_Build"].dropna().unique().tolist()
    log(f"Reference assemblies present: {builds}")
    assert builds == ["GRCh38"], f"Unexpected assembly mix: {builds}"

    # ---- Turn the long tumour barcode into a shorter sample ID ---------------
    # A TCGA barcode like "TCGA-3L-AA1B-01A-11D-A36X-10" encodes the patient and
    # sample. The first 15 characters ("TCGA-3L-AA1B-01") identify the tumour
    # sample; "-01" means "primary tumour". We trim to those 15 characters so
    # that different lab aliquots of the same tumour collapse to one column.
    df["Sample"] = df["Tumor_Sample_Barcode"].str.slice(0, 15)

    # ---- The four biological filters (assignment Section 4) ------------------
    # We keep only the mutations that are meaningful for making neoantigens:
    f_coding   = df["BIOTYPE"] == "protein_coding"          # gene makes a protein
    f_pass     = df["GDC_FILTER"].fillna("").isin(["", "PASS"])  # high-confidence
    f_missense = df["Variant_Classification"] == "Missense_Mutation"  # changes 1 amino acid
    f_snp      = df["Variant_Type"] == "SNP"                # single-DNA-letter change

    # WHY these four?
    #  - protein_coding: only protein-coding genes can produce peptides.
    #  - PASS/high-confidence: throw out likely sequencing artefacts.
    #  - Missense_Mutation: it swaps one amino acid for another -> a NEW protein
    #    sequence that the immune system might recognise (a potential neoantigen).
    #    (Silent mutations don't change the protein; nonsense/frameshift are a
    #     separate, more complex case left as an extension.)
    #  - SNP: a clean single-letter substitution, the simplest case.
    log(f"protein_coding: {f_coding.sum()}")
    log(f"PASS/high-confidence: {f_pass.sum()}")
    log(f"Missense_Mutation: {f_missense.sum()}")
    log(f"SNP: {f_snp.sum()}")

    # Combine the four filters with "&" (AND): a row is kept only if ALL are true.
    keep = f_coding & f_pass & f_missense & f_snp
    filt = df[keep].copy()
    n_filt = len(filt)
    log(f"Records after all filters: {n_filt}")

    # Require that each mutation has both notations filled in (Section 5 asks for
    # HGVS notation). Rows missing them are dropped.
    filt = filt[filt["HGVSc"].notna() & filt["HGVSp_Short"].notna()]
    filt = filt[(filt["HGVSc"] != "") & (filt["HGVSp_Short"] != "")]
    log(f"Records with complete HGVS notation: {len(filt)}")

    # ---- Define what counts as ONE distinct mutation -------------------------
    # We label a mutation by three things together: the gene, the DNA change
    # (HGVSc) and the protein change (HGVSp_Short). Two rows with the same three
    # values are "the same mutation" seen in two different tumours.
    filt["Gene_Name"] = filt["Hugo_Symbol"]
    filt["Mutation"] = filt["HGVSc"]
    filt["AminoAcid_Change"] = filt["HGVSp_Short"]
    key_cols = ["Gene_Name", "Mutation", "AminoAcid_Change"]

    # ---- Build the 0/1 grid efficiently --------------------------------------
    # A tidy way to build a big grid without running out of memory: make a
    # numbered list of the distinct mutations (rows) and of the samples
    # (columns), then flip a "1" into the right cell for each observation.
    import numpy as np      # numpy = fast numerical arrays
    # "_mutkey" = a single text label combining the three key columns.
    filt["_mutkey"] = (filt["Gene_Name"] + "\t" + filt["Mutation"]
                       + "\t" + filt["AminoAcid_Change"])
    # mut_index = the list of unique mutations, each given a row number "_row".
    mut_index = (filt[key_cols + ["_mutkey"]]
                 .drop_duplicates("_mutkey")
                 .reset_index(drop=True))
    mut_index["_row"] = np.arange(len(mut_index))
    row_map = dict(zip(mut_index["_mutkey"], mut_index["_row"]))  # mutation -> row #

    sample_cols = sorted(filt["Sample"].unique().tolist())        # the samples
    col_map = {s: i for i, s in enumerate(sample_cols)}           # sample -> column #

    # Start with an all-zeros grid (uint8 = tiny 0-255 integers to save memory).
    M = np.zeros((len(mut_index), len(sample_cols)), dtype=np.uint8)
    rows = filt["_mutkey"].map(row_map).to_numpy()   # row number for each observation
    cols = filt["Sample"].map(col_map).to_numpy()    # column number for each observation
    M[rows, cols] = 1                                # put a 1 wherever a mutation was seen

    # Glue the three label columns to the left of the 0/1 grid to form the table.
    mat = pd.concat(
        [mut_index[key_cols].reset_index(drop=True),
         pd.DataFrame(M, columns=sample_cols)],
        axis=1,
    )

    # Sort so the most common mutations appear at the top (nice for reading).
    mat["_freq"] = mat[sample_cols].sum(axis=1)      # how many samples have each mutation
    mat = mat.sort_values(["_freq", "Gene_Name", "Mutation"],
                          ascending=[False, True, True]).drop(columns="_freq")

    # Save as a tab-separated file (Rule 1: outputs are tab-delimited .tsv).
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    mat.to_csv(OUT, sep="\t", index=False)
    log(f"Wrote {OUT}: {mat.shape[0]} mutations x {len(sample_cols)} samples")

    # ---- Write a small quality-control summary (assignment Section 8) --------
    with open(QC, "w") as fh:
        fh.write("MUTATION PROCESSING QC (Section 8)\n")
        fh.write(f"Reference assembly: GRCh38/hg38\n")
        fh.write(f"Mutation records before filtering: {n_total}\n")
        fh.write(f"Mutation records after filtering:  {n_filt}\n")
        fh.write(f"Distinct mutations (rows in matrix): {mat.shape[0]}\n")
        fh.write(f"Unique mutated genes: {mat['Gene_Name'].nunique()}\n")
        fh.write(f"Tumour samples: {len(sample_cols)}\n")
    log(f"Wrote QC summary: {QC}")

# This standard Python idiom means: "if this file is run directly, call main()."
if __name__ == "__main__":
    sys.exit(main())
