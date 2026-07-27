#!/usr/bin/env python3
"""
05_annotate_and_generate_peptides.py
Project 130 - Colorectal Cancer (TCGA-COAD) Neoantigen Discovery Pipeline

===============================================================================
BIOLOGICAL & COMPUTATIONAL PURPOSE
===============================================================================
This script implements Sections 9 and 10 of the project assignment (Advanced
Neoantigen Component). It translates somatic missense SNVs into mutant and
wild-type peptide sequences for downstream HLA binding and immunogenicity analysis.

===============================================================================
VARIANT-TO-PROTEIN ANNOTATION (§9)
===============================================================================
- Transcript Selection Rule: Maps each missense mutation to its affected
  transcript using VEP (Variant Effect Predictor) canonical annotations
  precomputed in the GDC MAF (Ensembl CANONICAL / MANE Select transcripts).
- Reference Assembly (Rule 5): Coordinates and transcripts strictly on GRCh38.

===============================================================================
MUTANT & WILD-TYPE PEPTIDE GENERATION (§10)
===============================================================================
For every eligible missense mutation, using the UniProt reference proteome:
  - Generates all mutation-containing 9-mers (MHC Class I / CD8+ T-cell epitopes).
  - Generates all mutation-containing 15-mers (MHC Class II / CD4+ T-cell epitopes).
  - Generates exact matching wild-type peptide controls for differential scoring.

===============================================================================
PEPTIDE VERIFICATION & AUDITING (§10 REQUIREMENTS)
===============================================================================
For every generated peptide, 4 strict verification rules are enforced:
  1. Reference AA Agreement: Asserts that reference amino acid matches the UniProt
     FASTA sequence at `Protein_position` (1-based index). Mismatches are audited
     and excluded (`results/peptide_audit_failed.tsv`), never silently mutated.
  2. Mutant AA Substitution: Placed cleanly at the specified 1-based `MutPos`.
  3. Boundary Verification: Sliding windows do not cross protein N- or C-termini.
  4. Pairwise Difference: Asserts that mutant and WT peptides differ exclusively
     at the single mutated residue index.

===============================================================================
INPUT & OUTPUT CONTRACTS
===============================================================================
Inputs:
  - `cohortMAF.2026-07-15.maf.gz` (VEP-annotated GDC MAF, GRCh38)
  - `uniprotkb_proteome_UP000005640_*.fasta` (Reviewed UniProt human proteome)

Outputs:
  - `results/peptides_all.tsv` (Deliverable peptide database)
  - `results/peptide_annotation.tsv` (Section 9 variant annotation record)
  - `results/peptide_audit_failed.tsv` (Audit log of excluded AA mismatches)
"""

import gzip
import glob
import os
import sys
import pandas as pd

# =============================================================================
# FILE PATHS & RESOURCE RESOLUTION
# =============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAF = "/tmp/coad.maf" if os.path.exists("/tmp/coad.maf") else \
      os.path.join(BASE, "cohortMAF.2026-07-15.maf.gz")
RES = os.path.join(BASE, "results")
OUT_PEP = os.path.join(RES, "peptides_all.tsv")
OUT_ANN = os.path.join(RES, "peptide_annotation.tsv")
OUT_AUDIT = os.path.join(RES, "peptide_audit_failed.tsv")

# Peptide window lengths: 9-mer (Class I CD8+) and 15-mer (Class II CD4+)
LENGTHS = [9, 15]

def log(m):
    """Prints timestamped progress messages to stdout with line flushing."""
    print(f"[05] {m}", flush=True)

# =============================================================================
# HELPER FUNCTIONS: FASTA LOADING & SLIDING WINDOW GENERATION
# =============================================================================
def find_fasta():
    """Locate UniProt human proteome FASTA in the project root."""
    cands = []
    for pat in ["*.fasta", "*.fasta.gz", "*.fa", "*.fa.gz",
                "uniprot*", "UP000005640*"]:
        cands += glob.glob(os.path.join(BASE, pat))
    cands = [c for c in cands if os.path.isfile(c)]
    if not cands:
        sys.exit("[05] ERROR: no protein FASTA found in project root. "
                 "Please drop the UniProt human proteome FASTA into "
                 f"{BASE}")
    # Prefer largest file (complete proteome)
    cands.sort(key=os.path.getsize, reverse=True)
    return cands[0]

def load_fasta(path):
    """Parses UniProt FASTA header (>sp|ACCESSION|NAME) into dictionary {ACC: sequence}."""
    op = gzip.open if path.endswith(".gz") else open
    seqs = {}
    acc = None
    buf = []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if acc is not None:
                    seqs[acc] = "".join(buf)
                buf = []
                hdr = line[1:].strip()
                parts = hdr.split("|")
                if len(parts) >= 2:
                    acc = parts[1]
                else:
                    acc = hdr.split()[0]
            else:
                buf.append(line.strip())
        if acc is not None:
            seqs[acc] = "".join(buf)
    return seqs

def windows(seq, m, L):
    """
    Generates all 1-based sliding windows of length L that overlap position m.
    
    Yields tuples of (peptide_sequence, 1-based_mutation_position_in_peptide, 1-based_protein_start).
    """
    n = len(seq)
    out = []
    for start in range(max(1, m - L + 1), min(n - L + 1, m) + 1):
        pep = seq[start - 1:start - 1 + L]
        if len(pep) == L:
            mutpos = m - start + 1  # 1-based mutation position in peptide
            out.append((pep, mutpos, start))
    return out

# =============================================================================
# MAIN PIPELINE EXECUTION
# =============================================================================
def main():
    # =========================================================================
    # STEP 1: LOAD UNIPROT REFERENCE PROTEOME
    # =========================================================================
    fasta_path = find_fasta()
    log(f"Using protein FASTA: {os.path.basename(fasta_path)}")
    seqs = load_fasta(fasta_path)
    log(f"Loaded {len(seqs)} protein sequences")

    # =========================================================================
    # STEP 2: PARSE VEP-ANNOTATED MAF FILE (§9)
    # =========================================================================
    log(f"Reading VEP-annotated MAF: {MAF}")
    comp = "gzip" if MAF.endswith(".gz") else None
    usecols = ["Hugo_Symbol", "NCBI_Build", "Chromosome", "Start_Position",
               "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification",
               "Variant_Type", "Transcript_ID", "HGVSc", "HGVSp_Short",
               "Protein_position", "Amino_acids", "Consequence", "ENSP",
               "SWISSPROT", "CANONICAL", "MANE", "BIOTYPE", "GDC_FILTER"]
    df = pd.read_csv(MAF, sep="\t", comment="#", dtype=str, usecols=usecols,
                     compression=comp, low_memory=False)

    # Filter for eligible protein-coding PASS missense SNVs
    keep = ((df["BIOTYPE"] == "protein_coding") &
            (df["GDC_FILTER"].fillna("").isin(["", "PASS"])) &
            (df["Variant_Classification"] == "Missense_Mutation") &
            (df["Variant_Type"] == "SNP"))
    df = df[keep].copy()
    log(f"Eligible missense SNV rows: {len(df)}")

    # Parse amino acid substitutions (RefAA/AltAA) and protein positions
    aa = df["Amino_acids"].fillna("")
    df["RefAA"] = aa.str.split("/").str[0]
    df["AltAA"] = aa.str.split("/").str[-1]
    pp = df["Protein_position"].fillna("").str.split("/").str[0]
    df["ProtPos"] = pd.to_numeric(pp, errors="coerce")
    df["UniProt"] = df["SWISSPROT"].fillna("").str.split(".").str[0]

    df = df[(df["RefAA"].str.len() == 1) & (df["AltAA"].str.len() == 1) &
            df["ProtPos"].notna() & (df["UniProt"] != "")]
    df["ProtPos"] = df["ProtPos"].astype(int)
    log(f"Rows with clean single-AA change + UniProt + position: {len(df)}")

    # Deduplicate variants across transcripts to unique protein variants
    var_cols = ["Hugo_Symbol", "Transcript_ID", "ENSP", "UniProt",
                "CANONICAL", "MANE", "Chromosome", "Start_Position",
                "Reference_Allele", "Tumor_Seq_Allele2", "HGVSc",
                "HGVSp_Short", "Consequence", "ProtPos", "RefAA", "AltAA"]
    variants = df[var_cols].drop_duplicates().reset_index(drop=True)
    log(f"Distinct protein variants: {len(variants)}")

    # =========================================================================
    # STEP 3: PEPTIDE GENERATION & REFERENCE AA AUDITING (§10)
    # =========================================================================
    ann_rows = []
    audit_rows = []
    n_ok = n_nofasta = n_mismatch = 0
    n_pep = 0

    # Stream peptide records directly to TSV to bound RAM consumption
    pep_cols = ["GeneName", "TranscriptID", "ProteinID", "UniProt",
                "Chromosome", "Position", "Ref", "Alt", "ProteinChange",
                "ProteinPosition", "PeptideLength", "MutPos", "Peptide", "Type"]
    pep_fh = open(OUT_PEP, "w")
    pep_fh.write("\t".join(pep_cols) + "\n")

    def emit(d):
        pep_fh.write("\t".join(str(d[c]) for c in pep_cols) + "\n")

    for r in variants.itertuples(index=False):
        seq = seqs.get(r.UniProt)
        base_ann = dict(
            GeneName=r.Hugo_Symbol, TranscriptID=r.Transcript_ID,
            ProteinID=r.ENSP, UniProt=r.UniProt,
            Chromosome=r.Chromosome, Position=r.Start_Position,
            Ref=r.Reference_Allele, Alt=r.Tumor_Seq_Allele2,
            HGVSc=r.HGVSc, ProteinChange=r.HGVSp_Short,
            Consequence=r.Consequence, ProteinPosition=r.ProtPos,
            RefAA=r.RefAA, AltAA=r.AltAA,
            CANONICAL=r.CANONICAL, MANE=r.MANE,
            ReferenceGenome="GRCh38",
            TranscriptSelection="Ensembl canonical / MANE (VEP pick in GDC MAF)",
        )
        if seq is None:
            n_nofasta += 1
            audit_rows.append({**base_ann, "AuditReason": "no_fasta_for_uniprot"})
            continue
        if r.ProtPos > len(seq):
            n_mismatch += 1
            audit_rows.append({**base_ann, "AuditReason": "position_out_of_range",
                               "FastaLen": len(seq)})
            continue
        # SECTION 10 CHECK: Reference AA must match UniProt FASTA sequence
        if seq[r.ProtPos - 1] != r.RefAA:
            n_mismatch += 1
            audit_rows.append({**base_ann, "AuditReason": "refAA_mismatch",
                               "FastaAA": seq[r.ProtPos - 1], "FastaLen": len(seq)})
            continue

        n_ok += 1
        ann_rows.append(base_ann)
        
        # Construct mutant full protein sequence by single AA substitution
        mut_seq = seq[:r.ProtPos - 1] + r.AltAA + seq[r.ProtPos:]
        
        # Extract sliding windows for 9-mers (Class I) and 15-mers (Class II)
        for L in LENGTHS:
            for pep_wt, mutpos, start in windows(seq, r.ProtPos, L):
                pep_mut = mut_seq[start - 1:start - 1 + L]
                # Assert mutant and WT peptides differ exclusively at the mutated index
                diffs = [i for i in range(L) if pep_wt[i] != pep_mut[i]]
                assert diffs == [mutpos - 1], (r.Hugo_Symbol, L, diffs, mutpos)
                common = dict(
                    GeneName=r.Hugo_Symbol, TranscriptID=r.Transcript_ID,
                    ProteinID=r.ENSP, UniProt=r.UniProt,
                    Chromosome=r.Chromosome, Position=r.Start_Position,
                    Ref=r.Reference_Allele, Alt=r.Tumor_Seq_Allele2,
                    ProteinChange=r.HGVSp_Short, ProteinPosition=r.ProtPos,
                    PeptideLength=L, MutPos=mutpos)
                emit({**common, "Peptide": pep_mut, "Type": "Mutant"})
                emit({**common, "Peptide": pep_wt, "Type": "WildType"})
                n_pep += 2

    pep_fh.close()
    
    # =========================================================================
    # STEP 4: WRITE ANNOTATION & AUDIT FILES
    # =========================================================================
    log(f"Variants OK: {n_ok}; no-FASTA: {n_nofasta}; ref mismatch: {n_mismatch}")
    pd.DataFrame(ann_rows).to_csv(OUT_ANN, sep="\t", index=False)
    pd.DataFrame(audit_rows).to_csv(OUT_AUDIT, sep="\t", index=False)
    log(f"Wrote {OUT_PEP}: {n_pep} peptide rows")
    log(f"Wrote {OUT_ANN}: {len(ann_rows)} annotated mutations")
    log(f"Wrote {OUT_AUDIT}: {len(audit_rows)} excluded (audit)")

if __name__ == "__main__":
    sys.exit(main())
