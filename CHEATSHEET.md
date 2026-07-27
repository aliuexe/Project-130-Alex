# Project 130 — One-Page Cheat Sheet

**Colorectal cancer (TCGA-COAD): mutations → expression → peptides → BigMHC presentation → immunogenicity → ranked neoantigens.** All GRCh38. Predictions, not experimental validation.

## Pipeline in one line each
`01` filter MAF (protein-coding, PASS, missense SNV) → mutation×sample matrix · `02` RNA-seq RSEM→TPM matrix · `03` integrate: median GeneLevelTPM (+SD/IQR) · `04` QC + figures · `05` peptides (9- & 15-mers, mut+WT) · `06` BigMHC presentation & immunogenicity (3 HLA-I) · `07` shortlist · **ext:** `10-12` co-mutation · `15` TMB control · `13` clonality · `14` practical filter + coverage · `16` immunogenicity · `17` composite rank.

## Key numbers (with the one-line "why")
| Item | Value | Why |
|---|---|---|
| Raw → filtered mutations | 310,472 → 184,574 | protein-coding, PASS, missense, SNV |
| Distinct mutations × samples | 153,996 × 586 | one row per (gene, HGVSc, HGVSp) |
| TPM matrix | 20,511 genes × 592 | RSEM rescaled so each sample sums to 1e6 |
| GeneLevelTPM | **median** across tumours | robust to outliers (+SD/IQR reported) |
| Peptides | 6,873,140 (145,612 muts; 2,964 audited out) | all mut-containing 9/15-mers + WT |
| Unique 9-mers scored | 2,460,296 × 3 alleles | HLA-A\*02:01, A\*01:01, A\*03:01 |
| BigMHC EL Presenter cut-offs | Strong $\ge 0.70$, Weak $\ge 0.50$ | BigMHC presentation probability $[0, 1]$ |
| Strong mutant presenters | 30,260 pairs | — |
| Shortlist | 13,715 (1,827 genes) | BigMHC_EL $\ge 0.50$ + delta>0 + TPM>1 + recur≥2 |
| Clonal fraction | 60% (VAF≥0.25) | clonal = better target; all top drivers clonal |
| Practical neoantigens | **1,113** | + WT non-presenter + TPM≥10 + clonal |
| Vaccine coverage | 82.3% (≥1 epitope) | practical shared neoantigen set |
| Hypermutator cut-off | ≥200 SNVs = 15.5% | bimodal valley ≈ known ~16% MSI/POLE |
| Top ranked neoantigen | **KRAS G12V** (#1) | composite of 5 axes, recurrence log-absolute |

## Validation (our evidence it works, without wet-lab)
- KRAS G12D reproduces the **textbook peptide exactly** (VVGADGVGK, MutPos 5).
- Top mutated genes recover **TP53, KRAS, PIK3CA**; top expressed includes **CEACAM5 (CEA)**.
- **KRAS/BRAF mutually exclusive** (OR 0.07); **KRAS+PIK3CA co-occur** (real).
- Prioritised candidates = the neoantigens the field already pursues (KRAS G12x on HLA-A\*03:01).

## Likely questions → 1-line answers
- **Why missense only?** clean new peptide; frameshift/indel = extension (esp. MSI).
- **Is it real TPM?** documented per-million rescaling of RSEM (proxy); native GDC TPM = strict version.
- **Does expression prove the mutant is made?** No — gene-level, not allele-specific (needs RNA VAF).
- **Median vs mean?** robust to outliers; SD/IQR reported to judge trust.
- **Why 3 HLA alleles?** Option A common alleles; no patient HLA typing (a limitation).
- **Biggest co-mutation TP53+KRAS — real?** No, 105 obs ≈ 103 expected (chance). Real one: KRAS+PIK3CA (84 vs 61).
- **Not a burden artefact?** Excluded hypermutators → APC+ATM+POLE triple 10→**0**; KRAS+PIK3CA survives & strengthens.
- **Why APC only 7% (vs ~75%)?** missense-only filter; APC is knocked out by truncation — a co-mutation-view limitation, not an error.
- **Presentation vs immunogenicity?** separate columns; BigMHC_EL (presentation) vs BigMHC_IM / Calis (immunogenicity).
- **#1 candidate?** KRAS G12V (composite score 0.851).

## Honest limitations (say first)
Gene-level (not allele-specific) expression · RSEM→TPM proxy · missense only · fixed HLA panel · VAF clonality proxy · immunogenicity is a sequence heuristic · composite is weight-sensitive · **all predictions, not validation.**

## One-sentence take-home
*"We built a reproducible neoantigen-discovery pipeline using BigMHC that recovers the known colorectal landscape as validation, ranks practical shared candidates (KRAS G12V on top), and — using a Nature Genetics method — controls for the tumour-burden confounder that inflates naïve co-mutation counts."*
