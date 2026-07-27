#!/usr/bin/env python3
"""
05_annotate_and_generate_peptides.py
Project 130 - Colorectal cancer (TCGA-COAD)  --  ADVANCED component

Implements assignment Sections 9 and 10:

  Section 9 (Variant-to-Protein Annotation):
    Each selected coding (missense) mutation is mapped to its affected
    transcript and protein using the VEP annotation ALREADY PRESENT in the
    GDC MAF (the MAF is VEP-annotated: Transcript_ID/ENSP/SWISSPROT,
    Protein_position, Amino_acids, Consequence, CANONICAL, MANE). No
    coordinates are recomputed, so hg38 consistency is preserved (Rule 5).

    Transcript-selection rule (must be consistent, Section 9): we use the
    single VEP-chosen transcript that GDC reports on each MAF row, which is
    the Ensembl CANONICAL / MANE Select transcript for the variant. We
    record CANONICAL and MANE flags for transparency.

  Section 10 (Mutant Peptide Generation):
    For every eligible missense mutation we generate, using the reference
    protein sequence:
      - all mutation-containing 9-mers   (MHC-I / CD8+)
      - all mutation-containing 15-mers  (MHC-II / CD4+)
      - the corresponding wild-type peptides
    Every peptide contains the altered residue. Mutation position within
    each peptide is reported with 1-based indexing (MutPos).

    Verification performed for every peptide (Section 10 requirements):
      - reference amino acid matches the reference protein sequence at the
        protein position;
      - the mutant amino acid is placed at the correct position;
      - the window does not cross the protein boundary;
      - WT and mutant peptides differ only at the mutated position.

Protein sequences: reference proteome FASTA supplied locally
  (UniProt reviewed human proteome UP000005640). Mutations are matched to a
  protein sequence by their SWISSPROT (UniProt) accession from the MAF.
  Mutations whose reference AA does not agree with the FASTA (isoform /
  version mismatch) are written to an audit file and excluded from peptide
  output (never silently mutated).

Inputs:
  cohortMAF.2026-07-15.maf.gz            (VEP-annotated MAF, GRCh38)
  <uniprot proteome fasta>.(fasta|fasta.gz)   dropped into project root
Outputs:
  results/peptides_all.tsv               (one row per peptide)
  results/peptide_annotation.tsv         (one row per mutation, Section 9)
  results/peptide_audit_failed.tsv       (ref-AA mismatches, excluded)

Peptide windows: for a mutation at 1-based protein position m and window
length L, every window [start, start+L-1] with start in
[max(1, m-L+1), min(len-L+1, m)] contains position m and stays in-bounds.
"""
import gzip
import glob
import os
import sys
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAF = "/tmp/coad.maf" if os.path.exists("/tmp/coad.maf") else \
      os.path.join(BASE, "cohortMAF.2026-07-15.maf.gz")
RES = os.path.join(BASE, "results")
OUT_PEP = os.path.join(RES, "peptides_all.tsv")
OUT_ANN = os.path.join(RES, "peptide_annotation.tsv")
OUT_AUDIT = os.path.join(RES, "peptide_audit_failed.tsv")

# three-letter -> one-letter for parsing HGVSp_Short if ever needed
LENGTHS = [9, 15]

def log(m): print(f"[05] {m}", flush=True)

# ---------------------------------------------------------------------------
def find_fasta():
    cands = []
    for pat in ["*.fasta", "*.fasta.gz", "*.fa", "*.fa.gz",
                "uniprot*", "UP000005640*"]:
        cands += glob.glob(os.path.join(BASE, pat))
    cands = [c for c in cands if os.path.isfile(c)]
    if not cands:
        sys.exit("[05] ERROR: no protein FASTA found in project root. "
                 "Please drop the UniProt human proteome FASTA into "
                 f"{BASE}")
    # prefer the largest (the proteome), not a stray small file
    cands.sort(key=os.path.getsize, reverse=True)
    return cands[0]

def load_fasta(path):
    """Return dict {uniprot_accession: sequence}. UniProt headers look like
    >sp|P01116|RASK_HUMAN ...  -> key P01116."""
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
                # parse sp|ACC|NAME
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
    """All 1-based windows of length L that include position m (1-based)."""
    n = len(seq)
    out = []
    for start in range(max(1, m - L + 1), min(n - L + 1, m) + 1):
        pep = seq[start - 1:start - 1 + L]
        if len(pep) == L:
            mutpos = m - start + 1  # 1-based position of mutation in peptide
            out.append((pep, mutpos, start))
    return out

# ---------------------------------------------------------------------------
def main():
    fasta_path = find_fasta()
    log(f"Using protein FASTA: {os.path.basename(fasta_path)}")
    seqs = load_fasta(fasta_path)
    log(f"Loaded {len(seqs)} protein sequences")

    log(f"Reading VEP-annotated MAF: {MAF}")
    comp = "gzip" if MAF.endswith(".gz") else None
    usecols = ["Hugo_Symbol", "NCBI_Build", "Chromosome", "Start_Position",
               "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification",
               "Variant_Type", "Transcript_ID", "HGVSc", "HGVSp_Short",
               "Protein_position", "Amino_acids", "Consequence", "ENSP",
               "SWISSPROT", "CANONICAL", "MANE", "BIOTYPE", "GDC_FILTER"]
    df = pd.read_csv(MAF, sep="\t", comment="#", dtype=str, usecols=usecols,
                     compression=comp, low_memory=False)

    # ---- Eligibility: protein-coding PASS missense SNV (same as core) -----
    keep = ((df["BIOTYPE"] == "protein_coding") &
            (df["GDC_FILTER"].fillna("").isin(["", "PASS"])) &
            (df["Variant_Classification"] == "Missense_Mutation") &
            (df["Variant_Type"] == "SNP"))
    df = df[keep].copy()
    log(f"Eligible missense SNV rows: {len(df)}")

    # ---- Parse Amino_acids (e.g. 'G/D') and Protein_position (e.g. '12/189')
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

    # ---- Deduplicate distinct protein variants ----------------------------
    # a distinct variant = (UniProt, ProtPos, RefAA, AltAA, Transcript_ID)
    var_cols = ["Hugo_Symbol", "Transcript_ID", "ENSP", "UniProt",
                "CANONICAL", "MANE", "Chromosome", "Start_Position",
                "Reference_Allele", "Tumor_Seq_Allele2", "HGVSc",
                "HGVSp_Short", "Consequence", "ProtPos", "RefAA", "AltAA"]
    variants = df[var_cols].drop_duplicates().reset_index(drop=True)
    log(f"Distinct protein variants: {len(variants)}")

    ann_rows = []
    audit_rows = []
    n_ok = n_nofasta = n_mismatch = 0
    n_pep = 0

    # Stream peptide rows straight to disk (millions of rows; avoids RAM blow-up)
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
        # Section 10 check: reference AA agrees with reference protein
        if seq[r.ProtPos - 1] != r.RefAA:
            n_mismatch += 1
            audit_rows.append({**base_ann, "AuditReason": "refAA_mismatch",
                               "FastaAA": seq[r.ProtPos - 1], "FastaLen": len(seq)})
            continue

        n_ok += 1
        ann_rows.append(base_ann)
        mut_seq = seq[:r.ProtPos - 1] + r.AltAA + seq[r.ProtPos:]
        for L in LENGTHS:
            for pep_wt, mutpos, start in windows(seq, r.ProtPos, L):
                pep_mut = mut_seq[start - 1:start - 1 + L]
                # sanity: differ only at mutated position
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
    log(f"Variants OK: {n_ok}; no-FASTA: {n_nofasta}; ref mismatch: {n_mismatch}")

    pd.DataFrame(ann_rows).to_csv(OUT_ANN, sep="\t", index=False)
    pd.DataFrame(audit_rows).to_csv(OUT_AUDIT, sep="\t", index=False)
    log(f"Wrote {OUT_PEP}: {n_pep} peptide rows")
    log(f"Wrote {OUT_ANN}: {len(ann_rows)} annotated mutations")
    log(f"Wrote {OUT_AUDIT}: {len(audit_rows)} excluded (audit)")

if __name__ == "__main__":
    sys.exit(main())
