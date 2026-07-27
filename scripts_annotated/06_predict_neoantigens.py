#!/usr/bin/env python3
# =============================================================================
# 06_predict_neoantigens.py   (ANNOTATED teaching copy)
# =============================================================================
# THE BIOLOGY (read first):
#   HLA (in humans) = MHC molecules. They grab short peptides inside a cell and
#   display them on the surface. T cells inspect these displayed peptides. If a
#   MUTANT peptide is displayed and looks foreign, T cells can kill that cell.
#   So a peptide is only a useful neoantigen if an HLA molecule can actually
#   BIND and display it.
#
#   Different people have different HLA "alleles" (versions), and each allele
#   binds a different set of peptides. Because we don't know each patient's HLA
#   type, we use a fixed panel of 3 common HLA class I alleles (Section 11).
#
# WHAT THIS SCRIPT DOES:
#   For every 9-mer peptide from script 05, it uses a machine-learning tool
#   called MHCflurry to PREDICT how tightly each of the 3 HLA alleles would bind
#   it. It records the strength, compares mutant vs wild-type, and writes the
#   final neoantigen table (Deliverable 04).
#
# HOW BINDING STRENGTH IS MEASURED — "IC50" in nanomolar (nM):
#   IC50 is a concentration. LOWER IC50 = STRONGER binding (counter-intuitive!).
#   Rough convention: <50 nM = strong binder, <500 nM = weak binder, else non-binder.
#   "Percentile rank" is another view: a low rank (e.g. 0.5%) means the peptide
#   binds better than 99.5% of random peptides for that allele.
#
# IMPORTANT DISTINCTION (Section 13): being DISPLAYED (binding) is not the same
#   as being IMMUNOGENIC (actually provoking a T-cell response). A strong binder
#   is not automatically immunogenic. We keep those ideas in separate columns
#   and never invent an immunogenicity number we didn't compute.
#
# INPUT : results/peptides_all.tsv (from script 05) + the integrated matrix (03)
# OUTPUT: results/04_neoantigen_predictions.tsv (Deliverable 04)
# =============================================================================

import os
import shutil    # used to check whether optional tools are installed
import sys
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
PEP = os.path.join(RES, "peptides_all.tsv")                       # peptides from 05
INT = os.path.join(RES, "03_integrated_mutation_expression.tsv")  # for TPM + frequency
SCORES9 = os.path.join(RES, "mhcflurry_scores_9mer.tsv")          # cached predictions
OUT = os.path.join(RES, "04_neoantigen_predictions.tsv")          # the final table

# The fixed HLA panels (Section 11, "Option A").
HLA_I = ["HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01"]   # class I -> 9-mers
HLA_II = ["HLA-DRB1*15:01", "HLA-DRB1*07:01"]           # class II -> 15-mers

def log(m): print(f"[06] {m}", flush=True)

# ---------------------------------------------------------------------------
# Two "predictor" classes. A class here is just a reusable object with a method.
# _MHCflurryPredictor is the real thing; _StubPredictor is a fake used only to
# test the plumbing quickly (never for the real submission).
# ---------------------------------------------------------------------------
class _StubPredictor:
    """Fake predictor for pipeline testing only (makes reproducible pretend
    numbers so we can check the code runs, without installing MHCflurry)."""
    version = "STUB-validation"
    def predict_to_dataframe(self, peptides, allele):
        import hashlib
        rows = []
        for p in peptides:
            # hash = a repeatable pseudo-random number derived from the peptide.
            h = int(hashlib.md5((p + allele).encode()).hexdigest(), 16)
            rows.append((p, float(20 + h % 5000),
                         round((h % 10000) / 100.0, 3),
                         round(((h >> 8) % 1000) / 1000.0, 4)))
        return pd.DataFrame(rows, columns=["peptide", "affinity",
                                           "affinity_percentile",
                                           "presentation_score"])

class _MHCflurryPredictor:
    """The REAL predictor: loads MHCflurry's trained neural-network models."""
    def __init__(self):
        from mhcflurry import Class1PresentationPredictor
        import mhcflurry
        self.pred = Class1PresentationPredictor.load()   # load the trained models
        self.version = mhcflurry.__version__
    def predict_to_dataframe(self, peptides, allele):
        # Ask MHCflurry to score a list of peptides against ONE HLA allele.
        r = self.pred.predict(
            peptides=list(peptides), alleles=[allele],
            include_affinity_percentile=True, verbose=0)
        # Keep just the columns we need from MHCflurry's output.
        keep = {}
        keep["peptide"] = r["peptide"]
        keep["affinity"] = r.get("affinity")                    # IC50 nM (lower=stronger)
        keep["affinity_percentile"] = r.get("affinity_percentile")  # % rank (lower=better)
        keep["presentation_score"] = r.get("presentation_score")    # 0-1 (higher=better)
        return pd.DataFrame(keep)

def load_predictor():
    # Decide which predictor to use. Normally the real MHCflurry; if it isn't
    # installed we return None and scores become NA (never faked).
    if os.environ.get("P130_STUB") == "1":
        log("Using STUB predictor (validation only, not for submission)")
        return _StubPredictor()
    try:
        p = _MHCflurryPredictor()
        log(f"Loaded MHCflurry {p.version}")
        return p
    except Exception as e:
        log(f"MHCflurry unavailable: {e}")
        return None

# ---------------------------------------------------------------------------
def build_maps():
    # Reads the integrated matrix (03) and builds two quick lookup tables:
    #  1) mutfreq: how many tumour samples carry each SPECIFIC mutation
    #     (this is the "recurrence" of the neoantigen — a key priority feature).
    #  2) tpm: each gene's expression level (GeneLevelTPM).
    mutfreq, tpm = {}, {}
    with open(INT) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # sample columns start at the first "TCGA" header (robust to metadata
        # columns such as GeneLevelTPM_SD sitting before them)
        s0 = next(i for i, c in enumerate(header) if c.startswith("TCGA"))
        for line in fh:
            p = line.rstrip("\n").split("\t")
            gene, aachange = p[0], p[2]                 # gene + protein change (e.g. p.G12D)
            n_present = sum(1 for v in p[s0:] if v == "1")  # count 1's = samples with it
            key = (gene, aachange)
            mutfreq[key] = mutfreq.get(key, 0) + n_present
            if p[3] != "NA":                            # column 3 is GeneLevelTPM
                tpm[gene] = float(p[3])
    return mutfreq, tpm

def unique_9mers():
    # Collect the SET of distinct 9-mer peptide sequences (a "set" auto-removes
    # duplicates). Scoring each unique sequence once — instead of millions of
    # repeats — saves enormous time.
    s = set()
    with open(PEP) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[10] == "9":                            # column 10 = PeptideLength
                s.add(p[12])                            # column 12 = Peptide sequence
    return sorted(s)

def score_9mers(pred):
    # Produce a lookup: (peptide, allele) -> (affinity, percentile, presentation).
    # We CACHE the scores to a file so a re-run doesn't have to recompute the
    # slow neural-network predictions.
    if os.path.exists(SCORES9):                         # already computed before?
        log(f"Loading cached scores: {SCORES9}")
        sc = pd.read_csv(SCORES9, sep="\t")
    elif pred is None:                                  # no tool and no cache
        log("No MHCflurry predictor and no score cache; "
            "class-I scores will be NA.")
        return {}
    else:
        uniq = unique_9mers()
        # MHCflurry only understands the 20 standard amino acids. A few peptides
        # from the proteome contain rare letters (e.g. 'U' = selenocysteine);
        # we skip those and report them as NA rather than guess.
        STD = set("ACDEFGHIKLMNPQRSTVWY")
        scorable = [p for p in uniq if set(p) <= STD]
        skipped = len(uniq) - len(scorable)
        log(f"Unique 9-mers: {len(uniq)}; scorable (standard AA): "
            f"{len(scorable)}; skipped nonstandard (reported NA): {skipped}")
        frames = []
        for a in HLA_I:                                 # score once per allele
            log(f"  scoring allele {a} ...")
            df = pred.predict_to_dataframe(scorable, a)
            df["HLAAllele"] = a
            frames.append(df)
        sc = pd.concat(frames, ignore_index=True)       # stack the 3 alleles' results
        sc.to_csv(SCORES9, sep="\t", index=False)       # save the cache
        # Record which MHCflurry version made these scores (for the report).
        with open(SCORES9 + ".version", "w") as vh:
            vh.write(str(getattr(pred, "version", "NA")))
        log(f"Wrote cache: {SCORES9} ({len(sc)} rows)")
    # Turn the score table into a fast dictionary keyed by (peptide, allele).
    keys = list(zip(sc["peptide"].tolist(), sc["HLAAllele"].tolist()))
    vals = list(zip(sc["affinity"].astype(float).tolist(),
                    sc["affinity_percentile"].astype(float).tolist(),
                    sc["presentation_score"].astype(float).tolist()))
    return dict(zip(keys, vals))

def binder(aff):
    # Turn an IC50 affinity into a human label. Remember: lower nM = stronger.
    if aff is None or (isinstance(aff, float) and np.isnan(aff)):
        return "NA"
    if aff < 50: return "Strong"        # <50 nM  = strong binder
    if aff < 500: return "Weak"         # <500 nM = weak binder
    return "Non-binder"

# The exact columns of the output table, in the order the assignment specifies.
OUT_COLS = ["GeneName", "Chromosome", "Position", "Ref", "Alt", "TranscriptID",
            "ProteinChange", "GeneLevelTPM", "MutationFrequency", "PeptideType",
            "Peptide", "PeptideLength", "MutPos", "HLAAllele", "BindingAffinity",
            "BindingRank", "BindingScore", "BinderClass",
            "DeltaAffinity_WTminusMut", "ImmunogenicityScore",
            "PredictionTool", "ToolVersion", "PredictionMode"]

def main():
    if not os.path.exists(PEP):
        sys.exit(f"[06] ERROR: {PEP} not found. Run script 05 first.")

    mutfreq, tpm = build_maps()        # frequency + expression lookups
    pred = load_predictor()            # the MHCflurry model (or None)
    scores = score_9mers(pred)         # (peptide, allele) -> predicted numbers

    # Record which tool/version produced the scores (for transparency, Rule 4).
    if pred is not None:
        tool, ver, mode = "MHCflurry", pred.version, "presentation"
    elif scores:                       # regenerated from the cache, no live tool
        vfile = SCORES9 + ".version"
        ver = (open(vfile).read().strip()
               if os.path.exists(vfile) else "NA")
        tool, mode = "MHCflurry", "presentation"
    else:
        tool, ver, mode = "NA", "NA", "NA"

    # Optional extra tools. If they aren't installed, related columns stay NA.
    prime = shutil.which("PRIME") or shutil.which("PRIME.sh")   # immunogenicity tool
    netmhc2 = shutil.which("netMHCIIpan")                       # class-II binding tool
    imm_tool = "PRIME" if prime else "NA"
    log(f"PRIME on PATH: {bool(prime)}; NetMHCIIpan on PATH: {bool(netmhc2)}")

    # ---- PASS 1: gather mutant & wild-type affinities so we can compare them --
    # For each mutation window and allele we remember both the Mutant and the
    # WildType binding affinity; the difference is the Section-14 "delta".
    log("Pass 1/2: collecting 9-mer affinities for WT-vs-Mut delta")
    aff_by_key = {}   # (gene, proteinChange, protPos, mutPos, allele) -> {Mutant:..., WildType:...}
    with open(PEP) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[10] != "9":                            # only 9-mers get class-I scores
                continue
            gene, pchg, ppos, mutpos, pep_seq, ptype = (
                p[0], p[8], p[9], p[11], p[12], p[13])
            for a in HLA_I:
                aff = scores.get((pep_seq, a), (None, None, None))[0]
                k = (gene, pchg, ppos, mutpos, a)
                aff_by_key.setdefault(k, {})[ptype] = aff

    # ---- PASS 2: write the final table, one row per (peptide, allele) ---------
    log("Pass 2/2: writing neoantigen table")
    n_out = 0
    with open(PEP) as fin, open(OUT, "w") as fout:
        fout.write("\t".join(OUT_COLS) + "\n")
        fin.readline()
        for line in fin:
            p = line.rstrip("\n").split("\t")
            (gene, tid, ensp, uni, chrom, pos, ref, alt, pchg, ppos,
             length, mutpos, pep_seq, ptype) = p
            glt = tpm.get(gene, None)
            glt_s = "NA" if glt is None else f"{glt:.4f}"       # gene expression
            mfreq = mutfreq.get((gene, pchg), 0)                # per-mutation recurrence

            if length == "9":                                   # class I (9-mers)
                for a in HLA_I:
                    aff, rank, pres = scores.get((pep_seq, a),
                                                 (None, None, None))
                    # Look up the paired Mutant/WildType affinities for the delta.
                    k = (gene, pchg, ppos, mutpos, a)
                    pair = aff_by_key.get(k, {})
                    wt, mut = pair.get("WildType"), pair.get("Mutant")
                    if wt is not None and mut is not None:
                        # delta = WT - Mut. Positive => mutant binds STRONGER.
                        delta = f"{(wt - mut):.3f}"
                    else:
                        delta = "NA"
                    imm = "NA"   # immunogenicity: NA unless PRIME is run (never faked)
                    row = [gene, chrom, pos, ref, alt, tid, pchg, glt_s,
                           str(mfreq), ptype, pep_seq, length, mutpos, a,
                           "NA" if aff is None else f"{aff:.3f}",     # BindingAffinity
                           "NA" if rank is None else f"{rank:.4f}",   # BindingRank
                           "NA" if pres is None else f"{pres:.4f}",   # BindingScore
                           binder(aff), delta, imm, tool, ver, mode]
                    fout.write("\t".join(row) + "\n")
                    n_out += 1
            else:                                               # class II (15-mers)
                for a in HLA_II:
                    # We report the 15-mers but leave their binding NA unless a
                    # class-II tool (NetMHCIIpan) is installed — honest, not zero.
                    row = [gene, chrom, pos, ref, alt, tid, pchg, glt_s,
                           str(mfreq), ptype, pep_seq, length, mutpos, a,
                           "NA", "NA", "NA", "NA", "NA", "NA",
                           "NetMHCIIpan" if netmhc2 else "NA",
                           "NA", "class_II" if netmhc2 else "NA"]
                    fout.write("\t".join(row) + "\n")
                    n_out += 1

    log(f"Wrote {OUT}: {n_out} rows")
    log("Section 14 delta = WT_affinity - Mut_affinity (nM); "
        "positive => mutant binds more strongly than wild-type.")

if __name__ == "__main__":
    sys.exit(main())
