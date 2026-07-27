#!/usr/bin/env python3
"""
20_bigmhc_predict.py
Project 130 - Standalone BigMHC prediction runner.

BigMHC (Albert et al., Nat Mach Intell 2023; github.com/KarchinLab/bigmhc) is a
deep-learning ensemble that predicts MHC-I **presentation** (mode `el`, an
eluted-ligand probability in [0,1]) and, via transfer learning, **immunogenicity**
(mode `im`).

WHAT THIS SCRIPT DOES
  1. Collects unique 9-mer peptides (mutant AND wild-type) from peptides_all.tsv.
  2. Writes a BigMHC input CSV (allele, peptide) for each class-I allele.
  3. Runs BigMHC in mode `el` (presentation) and `im` (immunogenicity) via subprocess.
  4. Parses BigMHC's `.prd` output and writes one tidy cache:
        results/bigmhc_scores_9mer.tsv
        columns: peptide, allele, BigMHC_EL, BigMHC_IM
"""
import csv, os, subprocess, sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
PEP = os.path.join(RES, "peptides_all.tsv")
OUT = os.path.join(RES, "bigmhc_scores_9mer.tsv")
WORK = os.path.join(RES, "_bigmhc_tmp"); os.makedirs(WORK, exist_ok=True)

HLA_I = ["HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01"]
PRESENT_CUTOFF = 0.5
BIGMHC_DIR = os.environ.get("BIGMHC_DIR", "")     # path to cloned bigmhc repo
DEVICE = os.environ.get("BIGMHC_DEVICE", "cpu")   # "cpu" or "cuda"

def log(m): print("[20]", m, flush=True)

def unique_9mers():
    s = set()
    with open(PEP) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[10] == "9":
                s.add(p[12])
    # BigMHC needs the 20 standard amino acids
    STD = set("ACDEFGHIKLMNPQRSTVWY")
    return sorted(x for x in s if set(x) <= STD)

def run_bigmhc(infile, mode):
    """Call BigMHC predict.py; returns path to the .prd output."""
    if not BIGMHC_DIR:
        sys.exit("[20] ERROR: set BIGMHC_DIR to your cloned bigmhc repo, e.g.\n"
                 "  BIGMHC_DIR=/path/to/bigmhc python3 scripts/20_bigmhc_predict.py")
    predict = os.path.join(BIGMHC_DIR, "src", "predict.py")
    # input CSV columns: 0=allele, 1=peptide  -> -c=0 -p=1
    cmd = [sys.executable, predict, f"-i={infile}", f"-m={mode}",
           "-c=0", "-p=1", f"-d={DEVICE}"]
    log("running: " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    return infile + ".prd"

def parse_prd(prd_path):
    """BigMHC writes the input CSV plus an appended score column. Return
    dict {(allele, peptide): score} using the LAST numeric column."""
    out = {}
    with open(prd_path) as fh:
        rdr = csv.reader(fh)
        header = next(rdr)
        # score column = the last column BigMHC added (named BigMHC_EL / BigMHC_IM)
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
        # presentation (EL) and immunogenicity (IM)
        el_scores.update(parse_prd(run_bigmhc(infile, "el")))
        im_scores.update(parse_prd(run_bigmhc(infile, "im")))

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
