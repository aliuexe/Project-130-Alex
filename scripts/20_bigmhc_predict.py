#!/usr/bin/env python3
"""
20_bigmhc_predict.py
Project 130 - Colorectal Cancer (TCGA-COAD) Neoantigen Discovery Pipeline

===============================================================================
BIOLOGICAL & COMPUTATIONAL PURPOSE
===============================================================================
This script serves as the standalone BigMHC deep-learning prediction runner
(Albert et al., Nature Machine Intelligence 2023). It interface directly with
the BigMHC Python package via CLI subprocess execution.

===============================================================================
EVALUATION MODES
===============================================================================
  1. Eluted-Ligand Presentation Mode (`-m=el`): Predicts HLA Class I presentation
     probability `BigMHC_EL` in $[0, 1]$.
  2. Immunogenicity Mode (`-m=im`): Predicts transfer-learned T-cell immunogenicity
     log-odds `BigMHC_IM`.

===============================================================================
INPUT & OUTPUT CONTRACTS
===============================================================================
Inputs:
  - `results/peptides_all.tsv` (Extracts unique 9-mer peptides)
  - `BIGMHC_DIR` environment variable pointing to cloned BigMHC repository

Output:
  - `results/bigmhc_scores_9mer.tsv` (Cached BigMHC presentation and immunogenicity database)
"""

import csv, os, subprocess, sys

# =============================================================================
# FILE PATHS & RESOURCE RESOLUTION
# =============================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
PEP = os.path.join(RES, "peptides_all.tsv")
OUT = os.path.join(RES, "bigmhc_scores_9mer.tsv")
WORK = os.path.join(RES, "_bigmhc_tmp"); os.makedirs(WORK, exist_ok=True)

# Hardware & Subprocess Settings
HLA_I = ["HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01"]
PRESENT_CUTOFF = 0.5
BIGMHC_DIR = os.environ.get("BIGMHC_DIR", "")     # Path to cloned BigMHC repo
DEVICE = os.environ.get("BIGMHC_DEVICE", "cpu")   # "cpu" or "cuda"

def log(m):
    """Prints timestamped progress messages to stdout with line flushing."""
    print("[20]", m, flush=True)

# =============================================================================
# HELPER FUNCTIONS: FASTA EXTRACTION & SUBPROCESS WRAPPING
# =============================================================================
def unique_9mers():
    """Extracts unique 9-mer peptides composed of standard amino acids."""
    s = set()
    with open(PEP) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[10] == "9":
                s.add(p[12])
    STD = set("ACDEFGHIKLMNPQRSTVWY")
    return sorted(x for x in s if set(x) <= STD)

def run_bigmhc(infile, mode):
    """Invokes BigMHC `predict.py` via CLI subprocess for presentation (`el`) or immunogenicity (`im`)."""
    if not BIGMHC_DIR:
        sys.exit("[20] ERROR: set BIGMHC_DIR to your cloned bigmhc repo, e.g.\n"
                 "  BIGMHC_DIR=/path/to/bigmhc python3 scripts/20_bigmhc_predict.py")
    predict = os.path.join(BIGMHC_DIR, "src", "predict.py")
    cmd = [sys.executable, predict, f"-i={infile}", f"-m={mode}",
           "-c=0", "-p=1", f"-d={DEVICE}"]
    log("running: " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    return infile + ".prd"

def parse_prd(prd_path):
    """Parses BigMHC `.prd` output CSV into dictionary `{(allele, peptide): score}`."""
    out = {}
    with open(prd_path) as fh:
        rdr = csv.reader(fh)
        header = next(rdr)
        score_idx = len(header) - 1
        for row in rdr:
            if len(row) <= score_idx:
                continue
            allele, pep = row[0], row[1]
            try:
                out[(allele, pep)] = float(row[score_idx])
            except ValueError:
                pass
    return out

# =============================================================================
# MAIN PIPELINE EXECUTION
# =============================================================================
def main():
    peps = unique_9mers()
    log(f"unique standard-AA 9-mers: {len(peps)}")

    el_scores, im_scores = {}, {}
    for allele in HLA_I:
        infile = os.path.join(WORK, f"in_{allele.replace('*','').replace(':','')}.csv")
        with open(infile, "w", newline="") as fh:
            w = csv.writer(fh); w.writerow(["allele", "peptide"])
            for pep in peps:
                w.writerow([allele, pep])
        # Run presentation (el) and immunogenicity (im)
        el_scores.update(parse_prd(run_bigmhc(infile, "el")))
        im_scores.update(parse_prd(run_bigmhc(infile, "im")))

    # Export combined cache file
    with open(OUT, "w") as fh:
        fh.write("peptide\tallele\tBigMHC_EL\tBigMHC_IM\tPresented\n")
        for allele in HLA_I:
            for pep in peps:
                el = el_scores.get((allele, pep))
                im = im_scores.get((allele, pep))
                pres = "1" if (el is not None and el >= PRESENT_CUTOFF) else "0"
                fh.write(f"{pep}\t{allele}\t"
                         f"{'NA' if el is None else round(el,4)}\t"
                         f"{'NA' if im is None else round(im,4)}\t{pres}\n")
    log(f"wrote {OUT} ({len(peps)*len(HLA_I)} peptide-allele rows)")
    n_pres = sum(1 for v in el_scores.values() if v >= PRESENT_CUTOFF)
    log(f"presented (BigMHC_EL >= {PRESENT_CUTOFF}): {n_pres} peptide-allele pairs")

if __name__ == "__main__":
    sys.exit(main())
