#!/usr/bin/env bash
###############################################################################
# run_full_predictions.sh
# Project 130 - Colorectal cancer (TCGA-COAD)  --  ADVANCED component runner
#
# Runs the full-scale advanced neoantigen pipeline on YOUR machine using BigMHC
# (Albert et al., Nature Machine Intelligence 2023).
#
# Order of operations (matches the assignment):
#   1. 05_annotate_and_generate_peptides.py  (Sections 9-10) -> peptides_all.tsv
#   2. 06_predict_neoantigens.py             (Sections 11-15) -> 04_neoantigen_predictions.tsv
#
# PREREQUISITES:
#   pip install pandas numpy torch
#   # optional, to populate class-II (15-mer) columns:
#   #   install NetMHCIIpan 4.x and put 'netMHCIIpan' on PATH
#
# Also required in the project root:
#   cohortMAF.2026-07-15.maf.gz                 (the GDC VEP MAF)
#   the UniProt human proteome FASTA            (UP000005640, reviewed)
#
# USAGE:
#   cd "<project root>"
#   bash scripts/run_full_predictions.sh
###############################################################################
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
echo "[runner] project root: $ROOT"

# 1. Peptide generation (fast; runs on ALL eligible missense mutations)
echo "[runner] Step 1/2: generating peptides (Sections 9-10)"
python3 scripts/05_annotate_and_generate_peptides.py

# 2. BigMHC presentation + immunogenicity + Section 14/15 assembly
echo "[runner] Step 2/2: BigMHC presentation + immunogenicity (Sections 11-15)"
echo "[runner] NOTE: scoring unique 9-mers with BigMHC..."
python3 scripts/06_predict_neoantigens.py

echo "[runner] DONE."
echo "[runner] Deliverable written: results/04_neoantigen_predictions.tsv"
