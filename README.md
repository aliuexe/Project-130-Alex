# Project 130 — Integrating Cancer Mutations, Gene Expression, and Neoantigen Prediction
## Colorectal cancer (TCGA-COAD)

---

### Selected cancer type
Colorectal cancer, analysed through **TCGA-COAD** (Colon Adenocarcinoma) somatic
mutation data and matched TCGA colorectal RNA-seq expression data.

---

### Dataset accession numbers and sources
| Dataset | Source | Identifier / file | Access |
|---|---|---|---|
| Somatic mutations (MAF) | NCI Genomic Data Commons (GDC), TCGA-COAD | Project-level **Masked Somatic Mutation** MAF, downloaded via the GDC Data Transfer Tool (`gdc-client`) using a project-filtered manifest (`gdc_manifest.2026-07-20.150553.txt`). Local file: `cohortMAF.2026-07-15.maf.gz` | Open access |
| RNA-seq expression | TCGA Colorectal PanCancer Atlas, via cBioPortal | Study `coadread_tcga_pan_can_atlas_2018`, file `data_mrna_seq_v2_rsem.txt` (mRNA Expression, RSEM, batch-normalised, Illumina HiSeq RNASeqV2) | Open access |
| Reference proteome | UniProt | Reviewed human proteome **UP000005640** (`uniprotkb_proteome_UP000005640_*.fasta`) — used for peptide extraction / reference-AA verification | Open access |

---

### Download dates
- GDC somatic mutation MAF: **2026-07-20** (manifest generated 2026-07-20; file `cohortMAF.2026-07-15.maf.gz`).
- cBioPortal RNA-seq (RSEM): **2026-07-20**.
- UniProt human proteome FASTA: **2026-07-20**.

---

### Reference genome assembly
- **GRCh38 / hg38** throughout. Every mutation record in the MAF reports
  `NCBI_Build = GRCh38`; this is asserted in `scripts/01_build_mutation_matrix.py`.

---

### Software dependencies and versions
| Tool | Version | Used in |
|---|---|---|
| Python | 3.10+ | all `.py` scripts |
| pandas | 2.x | data processing |
| numpy | 1.2x / 2.x | numeric ops |
| matplotlib | 3.x | QC & analysis figures |
| PyTorch | 2.x | BigMHC deep-learning inference |
| BigMHC | 1.0.0 (Albert et al. 2023) | HLA class I 9-mer presentation (`BigMHC_EL`) and immunogenicity (`BigMHC_IM`) (script 06) |
| NetMHCIIpan | 4.x *(optional)* | HLA class II 15-mer binding (script 06); if absent, class-II columns are `NA` |

BigMHC stamps its version into every row of `04_neoantigen_predictions.tsv`
(`ToolVersion` column), per Rule 4.

---

### Pipeline — scripts, in execution order

| # | Script | Purpose | Assignment sections |
|---|---|---|---|
| 01 | `scripts/01_build_mutation_matrix.py` | Filter MAF (protein-coding, PASS, nonsynonymous missense SNV) → binary mutation-by-sample matrix | §4, §5 |
| 02 | `scripts/02_build_expression_matrix.py` | RSEM → TPM conversion → gene-by-sample TPM matrix | §6 |
| 03 | `scripts/03_integrate_datasets.py` | GeneLevelTPM = median TPM; merge into mutation matrix by gene | §7 |
| 04 | `scripts/04_qc_and_figures.py` | QC checks + 4 figures | §8 |
| 05 | `scripts/05_annotate_and_generate_peptides.py` | Variant→protein annotation; mutant + WT 9-mer & 15-mer generation | §9, §10 |
| 06 | `scripts/06_predict_neoantigens.py` | BigMHC presentation (`BigMHC_EL`), immunogenicity (`BigMHC_IM`), WT-vs-mutant delta → neoantigen table | §11–§15 |
| 07 | `scripts/07_prioritize_candidates.py` | Prioritised neoantigen candidate shortlist | §14 |
| 08 | `scripts/08_aggregate_for_figures.py` | Aggregate BigMHC metrics into `figure_summary.json` | — |
| 09 | `scripts/09_make_figures.py` | Render figures 5–13 with BigMHC presentation probabilities | — |
| 14 | `scripts/14_practical_neoantigens_coverage.py` | Practical neoantigens & greedy set-cover coverage curve | — |
| 16 | `scripts/16_immunogenicity.py` | Calis + BigMHC immunogenicity evaluation | — |
| 17 | `scripts/17_composite_score.py` | 5-axis composite quality scoring & final ranking | — |
| 18 | `scripts/18_hotspot_audit.py` | Recurrent hotspot audit | — |
| 19 | `scripts/19_longtail_coverage.py` | Long-tail coverage analysis | — |
| 20 | `scripts/20_bigmhc_predict.py` | Standalone BigMHC prediction runner | — |

---

### Explanation of Key Columns in `04_neoantigen_predictions.tsv`
- `GeneName, Chromosome, Position, Ref, Alt, TranscriptID, ProteinChange` — variant identity (GRCh38).
- `GeneLevelTPM` — median tumour TPM of the mutated gene (`NA` if unavailable).
- `MutationFrequency` — per-mutation recurrence: number of tumour samples carrying that mutation.
- `PeptideType` — `Mutant` or `WildType`.
- `Peptide, PeptideLength, MutPos` — peptide sequence, length (9 or 15), 1-based position of mutated residue.
- `HLAAllele` — HLA allele the prediction refers to (`HLA-A*02:01`, `HLA-A*01:01`, `HLA-A*03:01`).
- `BigMHC_EL` — predicted eluted-ligand presentation probability ($0\text{--}1$); higher = stronger presentation.
- `BigMHC_IM` — transfer-learned T-cell immunogenicity score.
- `PresentationClass` — `Strong` ($\ge 0.70$) / `Weak` ($\ge 0.50$) / `Non-presenter` ($< 0.50$).
- `DeltaPresentation_MutMinusWT` — `Mut_EL - WT_EL`. **Positive ⇒ mutant peptide is more strongly presented than wild-type**.
- `PredictionTool, ToolVersion, PredictionMode` — provenance stamp (`BigMHC`).
