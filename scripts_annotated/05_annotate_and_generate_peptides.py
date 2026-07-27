#!/usr/bin/env python3
# =============================================================================
# 05_annotate_and_generate_peptides.py   (ANNOTATED teaching copy)
# =============================================================================
# THE BIOLOGY THIS SCRIPT AUTOMATES (read this first):
#   Central dogma: DNA -> (transcription) -> RNA -> (translation) -> PROTEIN.
#   A protein is a chain of "amino acids" (20 kinds, each written as one letter,
#   e.g. G = glycine, D = aspartate). A "missense" mutation changes one DNA
#   letter, which changes ONE amino acid in the protein.
#
#   The immune system doesn't see whole proteins; it sees short PIECES of them
#   called "peptides", displayed on the cell surface by MHC/HLA molecules.
#   - MHC class I displays ~9-amino-acid peptides (9-mers) to CD8+ "killer" T cells.
#   - MHC class II displays ~15-amino-acid peptides (15-mers) to CD4+ "helper" T cells.
#   A peptide that contains a mutated amino acid, and looks "foreign", can be a
#   NEOANTIGEN (a tumour-specific flag the immune system might attack).
#
# WHAT THIS SCRIPT DOES:
#   For every missense mutation, it (1) finds the protein it affects, (2) cuts
#   out every short peptide window that contains the mutated amino acid, in both
#   the MUTANT and the normal ("wild-type") version, for lengths 9 and 15.
#
# WHY WILD-TYPE TOO?  So we can later compare mutant vs normal: a good neoantigen
#   is one the mutation makes MORE visible to the immune system than the normal
#   peptide.
#
# INPUT : the VEP-annotated MAF (already tells us the protein position & amino
#         acid change) + the UniProt human protein sequences (FASTA file).
# OUTPUT: results/peptides_all.tsv (every peptide), peptide_annotation.tsv,
#         and peptide_audit_failed.tsv (mutations we had to skip, with reasons).
# =============================================================================

import gzip
import glob     # lets us search for files by pattern (e.g. "*.fasta")
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

LENGTHS = [9, 15]      # peptide lengths to generate: 9-mers (MHC-I), 15-mers (MHC-II)

def log(m): print(f"[05] {m}", flush=True)

# ---------------------------------------------------------------------------
def find_fasta():
    # A FASTA file stores protein sequences. This finds it in the project folder
    # (it may be named a few different ways), preferring the largest file (the
    # full proteome, not some small stray file).
    cands = []
    for pat in ["*.fasta", "*.fasta.gz", "*.fa", "*.fa.gz",
                "uniprot*", "UP000005640*"]:
        cands += glob.glob(os.path.join(BASE, pat))
    cands = [c for c in cands if os.path.isfile(c)]
    if not cands:
        sys.exit("[05] ERROR: no protein FASTA found in project root. "
                 "Please drop the UniProt human proteome FASTA into "
                 f"{BASE}")
    cands.sort(key=os.path.getsize, reverse=True)
    return cands[0]

def load_fasta(path):
    # Read the FASTA into a dictionary: {protein_ID: full_amino_acid_sequence}.
    # A FASTA looks like:
    #   >sp|P01116|RASK_HUMAN GTPase KRas ...      <- a header line (starts ">")
    #   MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIED...   <- the sequence, maybe many lines
    # We use the middle part of the header (P01116) as the key.
    op = gzip.open if path.endswith(".gz") else open
    seqs = {}
    acc = None      # the current protein's ID
    buf = []        # collects the sequence lines for the current protein
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):                 # a new protein starts here
                if acc is not None:
                    seqs[acc] = "".join(buf)          # save the previous protein
                buf = []
                hdr = line[1:].strip()
                parts = hdr.split("|")                # split "sp|P01116|RASK_HUMAN"
                if len(parts) >= 2:
                    acc = parts[1]                    # -> "P01116"
                else:
                    acc = hdr.split()[0]
            else:
                buf.append(line.strip())              # add this line to the sequence
        if acc is not None:
            seqs[acc] = "".join(buf)                  # save the last protein
    return seqs

def windows(seq, m, L):
    # Given a protein sequence, a 1-based mutation position m, and a window
    # length L, return every length-L slice that CONTAINS position m and stays
    # inside the protein. Example for L=9 at position 12: the mutation can sit
    # anywhere from slot 1 to slot 9 of the 9-mer, giving up to 9 windows.
    n = len(seq)
    out = []
    for start in range(max(1, m - L + 1), min(n - L + 1, m) + 1):
        pep = seq[start - 1:start - 1 + L]           # Python counts from 0, hence -1
        if len(pep) == L:
            mutpos = m - start + 1                    # where the mutation sits (1-based)
            out.append((pep, mutpos, start))
    return out

# ---------------------------------------------------------------------------
def main():
    fasta_path = find_fasta()
    log(f"Using protein FASTA: {os.path.basename(fasta_path)}")
    seqs = load_fasta(fasta_path)                     # dictionary of protein sequences
    log(f"Loaded {len(seqs)} protein sequences")

    log(f"Reading VEP-annotated MAF: {MAF}")
    comp = "gzip" if MAF.endswith(".gz") else None
    # We read extra columns this time because we need the PROTEIN details that
    # a tool called VEP already computed and stored in the MAF:
    usecols = ["Hugo_Symbol", "NCBI_Build", "Chromosome", "Start_Position",
               "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification",
               "Variant_Type", "Transcript_ID", "HGVSc", "HGVSp_Short",
               "Protein_position", "Amino_acids", "Consequence", "ENSP",
               "SWISSPROT", "CANONICAL", "MANE", "BIOTYPE", "GDC_FILTER"]
    df = pd.read_csv(MAF, sep="\t", comment="#", dtype=str, usecols=usecols,
                     compression=comp, low_memory=False)

    # ---- Same biological filter as script 01 (protein-coding PASS missense SNV)
    keep = ((df["BIOTYPE"] == "protein_coding") &
            (df["GDC_FILTER"].fillna("").isin(["", "PASS"])) &
            (df["Variant_Classification"] == "Missense_Mutation") &
            (df["Variant_Type"] == "SNP"))
    df = df[keep].copy()
    log(f"Eligible missense SNV rows: {len(df)}")

    # ---- Pull out the protein change details ---------------------------------
    # "Amino_acids" looks like "G/D" (normal G becomes mutant D).
    # "Protein_position" looks like "12/189" (position 12 in a 189-aa protein).
    aa = df["Amino_acids"].fillna("")
    df["RefAA"] = aa.str.split("/").str[0]            # "G"  (reference amino acid)
    df["AltAA"] = aa.str.split("/").str[-1]           # "D"  (mutant amino acid)
    pp = df["Protein_position"].fillna("").str.split("/").str[0]
    df["ProtPos"] = pd.to_numeric(pp, errors="coerce")   # 12  (as a number)
    df["UniProt"] = df["SWISSPROT"].fillna("").str.split(".").str[0]  # protein ID

    # Keep only clean single-amino-acid changes with a known position & protein.
    df = df[(df["RefAA"].str.len() == 1) & (df["AltAA"].str.len() == 1) &
            df["ProtPos"].notna() & (df["UniProt"] != "")]
    df["ProtPos"] = df["ProtPos"].astype(int)
    log(f"Rows with clean single-AA change + UniProt + position: {len(df)}")

    # ---- Reduce to the list of DISTINCT protein variants ---------------------
    # (the same mutation can appear in many tumours; we only need to build its
    #  peptides once).
    var_cols = ["Hugo_Symbol", "Transcript_ID", "ENSP", "UniProt",
                "CANONICAL", "MANE", "Chromosome", "Start_Position",
                "Reference_Allele", "Tumor_Seq_Allele2", "HGVSc",
                "HGVSp_Short", "Consequence", "ProtPos", "RefAA", "AltAA"]
    variants = df[var_cols].drop_duplicates().reset_index(drop=True)
    log(f"Distinct protein variants: {len(variants)}")

    ann_rows = []       # one summary row per successfully annotated mutation
    audit_rows = []     # mutations we had to skip, with the reason
    n_ok = n_nofasta = n_mismatch = 0
    n_pep = 0

    # We write peptides straight to disk as we go (there are millions of them,
    # too many to hold in memory at once).
    pep_cols = ["GeneName", "TranscriptID", "ProteinID", "UniProt",
                "Chromosome", "Position", "Ref", "Alt", "ProteinChange",
                "ProteinPosition", "PeptideLength", "MutPos", "Peptide", "Type"]
    pep_fh = open(OUT_PEP, "w")
    pep_fh.write("\t".join(pep_cols) + "\n")
    def emit(d):                                       # helper: write one peptide row
        pep_fh.write("\t".join(str(d[c]) for c in pep_cols) + "\n")

    # ---- Go through every distinct mutation ----------------------------------
    for r in variants.itertuples(index=False):
        seq = seqs.get(r.UniProt)                      # the protein's amino-acid string
        base_ann = dict(                               # summary info for this mutation
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
        # ---- Three safety checks before we trust this mutation ---------------
        if seq is None:                                # we don't have that protein's sequence
            n_nofasta += 1
            audit_rows.append({**base_ann, "AuditReason": "no_fasta_for_uniprot"})
            continue
        if r.ProtPos > len(seq):                       # position is beyond the protein length
            n_mismatch += 1
            audit_rows.append({**base_ann, "AuditReason": "position_out_of_range",
                               "FastaLen": len(seq)})
            continue
        # CRUCIAL CHECK (Section 10): the amino acid the MAF says should be at
        # this position must actually match the reference protein. If it doesn't
        # (usually because of protein-isoform/version differences), we DON'T
        # guess — we record it in the audit file and skip it.
        if seq[r.ProtPos - 1] != r.RefAA:
            n_mismatch += 1
            audit_rows.append({**base_ann, "AuditReason": "refAA_mismatch",
                               "FastaAA": seq[r.ProtPos - 1], "FastaLen": len(seq)})
            continue

        n_ok += 1
        ann_rows.append(base_ann)
        # Make the mutant protein: copy the normal sequence but swap the one
        # amino acid at the mutated position for the mutant amino acid.
        mut_seq = seq[:r.ProtPos - 1] + r.AltAA + seq[r.ProtPos:]
        for L in LENGTHS:                              # for 9-mers, then 15-mers
            for pep_wt, mutpos, start in windows(seq, r.ProtPos, L):
                pep_mut = mut_seq[start - 1:start - 1 + L]   # the mutant version of this window
                # Sanity check: mutant and normal peptide must differ at EXACTLY
                # one spot — the mutated position. If not, something is wrong and
                # the script stops (this guards against silent bugs).
                diffs = [i for i in range(L) if pep_wt[i] != pep_mut[i]]
                assert diffs == [mutpos - 1], (r.Hugo_Symbol, L, diffs, mutpos)
                common = dict(
                    GeneName=r.Hugo_Symbol, TranscriptID=r.Transcript_ID,
                    ProteinID=r.ENSP, UniProt=r.UniProt,
                    Chromosome=r.Chromosome, Position=r.Start_Position,
                    Ref=r.Reference_Allele, Alt=r.Tumor_Seq_Allele2,
                    ProteinChange=r.HGVSp_Short, ProteinPosition=r.ProtPos,
                    PeptideLength=L, MutPos=mutpos)
                emit({**common, "Peptide": pep_mut, "Type": "Mutant"})   # mutant peptide
                emit({**common, "Peptide": pep_wt, "Type": "WildType"})  # normal peptide
                n_pep += 2

    pep_fh.close()
    log(f"Variants OK: {n_ok}; no-FASTA: {n_nofasta}; ref mismatch: {n_mismatch}")

    # Save the per-mutation summary and the audit (skipped) list.
    pd.DataFrame(ann_rows).to_csv(OUT_ANN, sep="\t", index=False)
    pd.DataFrame(audit_rows).to_csv(OUT_AUDIT, sep="\t", index=False)
    log(f"Wrote {OUT_PEP}: {n_pep} peptide rows")
    log(f"Wrote {OUT_ANN}: {len(ann_rows)} annotated mutations")
    log(f"Wrote {OUT_AUDIT}: {len(audit_rows)} excluded (audit)")

if __name__ == "__main__":
    sys.exit(main())
