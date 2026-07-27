# Project 130 — Every Choice We Made, and Why

### A decision log for defending the colorectal-cancer neoantigen pipeline (TCGA-COAD)

This document records **every deliberate choice** in the project — what we picked, what the alternatives were, why we chose it, and the honest caveat. It is organized in pipeline order. Two quick-reference sections are at the end: a table of **every threshold/number**, and **questions a reviewer is likely to ask with crisp answers**.

Guiding principle throughout: *make defensible, transparent choices, report every number, and never hide a limitation.*

---

## 0. Framing — what the project is actually trying to do

**Choice:** Treat this as a *neoantigen-discovery* pipeline, not just a binding predictor: mutation → peptide → HLA presentation → practical filtering → immunogenicity → a ranked candidate list, plus cohort-level co-mutation structure.

**Why:** The assignment's advanced component is about neoantigen potential, and a neoantigen is only useful if it is presented, tumour-specific, expressed, in all tumour cells, and (for a shared vaccine) recurrent. Every downstream choice serves that end goal.

**Caveat:** Everything is *computational prediction*, not experimental validation (Assignment Rule 8). We treat outputs as prioritisation, never proof.

---

## 1. Dataset choices

### 1.1 Cancer type = colorectal (TCGA-COAD)
- **Alternatives:** breast, lung, bladder, gastric, pancreatic, etc.
- **Why colorectal:** (a) large, well-characterised TCGA cohort; (b) a clear set of driver genes (APC, TP53, KRAS, PIK3CA, BRAF, SMAD4) that lets us *validate* the pipeline — if we recover them, the method works; (c) a hypermutated MSI/POLE subset that gives a high neoantigen burden and, usefully, a built-in confounder to demonstrate rigour; (d) KRAS G12x on HLA-A\*03:01 is a real, clinically pursued shared neoantigen, giving a positive control.
- **Caveat:** Results are colorectal-specific; the pipeline generalises but the biology does not.

### 1.2 Somatic mutations = GDC project-level Masked Somatic Mutation MAF (WXS, GRCh38)
- **Alternatives:** cBioPortal MAF; per-sample controlled-access MAFs; other callers.
- **Why:** It is the assignment's named canonical source; open-access (no dbGaP); already **VEP-annotated** (so protein consequences are precomputed and standardised); harmonised to a single reference genome (GRCh38).
- **Caveat:** "Masked" = germline-filtered, de-identified; a single consensus call set, so caller-specific sensitivity differences are not explored.

### 1.3 Expression = matched TCGA-COAD RNA-seq (cBioPortal RSEM)
- **Alternatives:** a separate GEO series; GDC STAR-Counts (native TPM).
- **Why:** Same cancer type and largely the same cohort as the mutations (the assignment's "matched TCGA" option), one download, fully reproducible.
- **Caveat:** RSEM is not TPM (see 2.2); and mutation vs expression samples are overlapping but not identical (see 3.2).

### 1.4 Protein sequences = UniProt reviewed human proteome (UP000005640)
- **Why:** Canonical, curated, one file; its accessions match the MAF's `SWISSPROT` column, giving a clean key to fetch each protein.
- **Caveat:** Isoform/sequence-version differences vs the MAF's transcript cause a small number of reference-AA mismatches, which we audit and exclude (see 5.4).

### 1.5 Reference genome = GRCh38/hg38 throughout
- **Why:** GDC data is all GRCh38; mixing assemblies would put mutations at wrong coordinates (Assignment Rule 5). Script 01 *asserts* only GRCh38 is present and stops otherwise.

---

## 2. Core pipeline choices

### Script 01 — mutation-by-sample matrix

**2.1 The four filters (protein-coding, PASS, missense, SNV).**
- *Protein-coding (`BIOTYPE == protein_coding`):* only protein-coding genes can produce peptides.
- *High-confidence (`GDC_FILTER` empty or PASS):* removes panel-of-normals and non-exonic flags — likely artefacts.
- *Missense (`Variant_Classification == Missense_Mutation`):* swaps one amino acid → a genuinely new protein sequence the immune system might see. Silent mutations don't change protein; nonsense/frameshift are a distinct, more complex class (see Limitations).
- *SNV (`Variant_Type == SNP`):* the simplest, cleanest change; keeps the peptide-generation logic exact.
- **Result:** 310,472 raw → 184,574 filtered records.
- **Caveat:** Restricting to missense excludes frameshift neoantigens, which are important in the MSI subset (flagged as an extension).

**2.2 "Distinct mutation" = (gene, HGVSc, HGVSp).**
- **Why:** A mutation is defined by its gene, DNA change, and protein change together; two tumours sharing all three share "the same" mutation. HGVS is the universal standard (Assignment Rule requiring HGVS).
- **Result:** 153,996 distinct mutations.

**2.3 Sample ID = first 15 characters of the TCGA barcode.**
- **Why:** Collapses lab aliquots of the same tumour to one sample-level column (`TCGA-XX-XXXX-01`), and the `-01` code confirms primary tumour. This matches the expression matrix's sample IDs so the two datasets share identifiers.

**2.4 Binary 0/1 matrix (present/absent), not counts or VAF.**
- **Why:** The assignment specifies a binary mutation-by-sample matrix; it's the simplest faithful "who has what," and feeds recurrence, heat maps, and co-mutation directly. (VAF is added later, separately, for clonality.)
- **Result:** 586 tumour samples.

### Script 02 — expression TPM matrix

**2.5 RSEM → TPM by rescaling each sample column to sum to 1,000,000.**
- **Why:** The source RSEM values sum to ~1.8×10⁷ per sample — they are *not* TPM. The assignment forbids labelling non-TPM as TPM (Rule 6). Rescaling each sample to sum to 1e6 gives the defining property of TPM and makes samples comparable. We *assert* the columns sum to 1e6 as a check.
- **Alternatives:** Use GDC STAR-Counts native `tpm_unstranded` (true, length-normalised TPM).
- **Caveat (state this proactively):** Strict TPM also involves gene-length normalisation; our rescaling is a documented per-million normalisation of RSEM, a common and defensible proxy but not fully length-normalised TPM. Native GDC TPM would be the rigorous alternative.

**2.6 Duplicate gene symbols summed to one row.**
- **Why:** Expression (abundance) is additive; 7 duplicate rows were collapsed.

**2.7 Missing values = NA, never 0.**
- **Why:** 0 would falsely assert "gene silent"; NA correctly says "unknown." 5.64% of cells are NA.
- **Result:** 20,511 genes × 592 samples.

### Script 03 — integration

**2.8 GeneLevelTPM = MEDIAN TPM across tumour samples (not mean).**
- **Why:** The median is robust to outliers — a few tumours with extreme expression won't distort it; the mean would be pulled up. The assignment explicitly recommends the median.
- **Caveat:** A single value hides variability — addressed by 2.9.

**2.9 Also report SD, IQR, and n (added at the PI's request).**
- **Why:** A large spread means the median is a less trustworthy summary. SD is what the PI asked for; IQR pairs more naturally with a median (robust); n tells how many samples backed the estimate.

**2.10 Gene-level join (by symbol), not sample-matched integration.**
- **Why:** The mutation and expression cohorts are overlapping but not identical, so we cannot line them up per patient; we summarise expression to one cancer-level value per gene (Assignment §7 permits this).
- **Caveat (important, PI-flagged):** This is *gene-level* — it tells us the gene is expressed, not that the *mutant allele* is transcribed (allele-specific expression). Bulk RNA-seq pools both alleles; confirming mutant-allele expression needs read-level RNA (RNA VAF), which our data doesn't support.

### Script 04 — QC and figures

**2.11 The QC checks (counts before/after filtering, unique genes, samples, % missing, identifier and assembly consistency, TPM distribution, top-10 mutated / expressed).**
- **Why:** Confirms the pipeline behaved and satisfies Assignment §8. The built-in sanity check: top mutated genes recover TP53/KRAS/PIK3CA and top expressed includes CEACAM5 (the CEA marker) — expected biology, so we trust the rest.

---

## 3. Advanced pipeline choices

### Script 05 — variant → peptide

**3.1 Use the MAF's existing VEP annotation instead of re-running VEP.**
- **Why:** The GDC MAF is already VEP-annotated (transcript, protein position, amino-acid change), all on GRCh38. Re-running VEP risks introducing a second coordinate system; reading the MAF's own annotation keeps everything on one assembly.

**3.2 Transcript selection = the single VEP-canonical / MANE transcript GDC records per row.**
- **Why:** The assignment demands a *consistent* transcript-selection rule. Using GDC's canonical pick is consistent and reproducible; we retain the CANONICAL/MANE flags for transparency.
- **Caveat:** Occasionally the canonical transcript uses non-textbook numbering (e.g. BRAF **V640E** on the long isoform, not the familiar V600E) — correct, just a different transcript.

**3.3 Generate 9-mers AND 15-mers.**
- **Why:** They match the two immune arms — HLA class I presents ~9-mers to CD8⁺ killer T cells, class II presents ~15-mers to CD4⁺ helper T cells. Covering both is complete.

**3.4 Every window that contains the mutated residue; mutant + wild-type; 1-based MutPos.**
- **Why:** The assignment requires all mutation-containing windows and the WT counterpart (needed for the mutant-vs-WT comparison). We verify, for every peptide, that the reference AA matches the reference protein, the mutant is placed correctly, the window stays in-bounds, and mutant/WT differ at exactly one position.
- **Result:** 145,612 annotated mutations → 6,873,140 peptides.

**3.5 Exclude (audit) mutations whose reference AA doesn't match the reference protein.**
- **Why:** A mismatch means the wrong protein/isoform; we do NOT guess — we log it and drop it. 2,964 excluded.

### Script 06 — HLA binding

**3.6 Predictor = MHCflurry.**
- **Alternatives:** NetMHCpan, MixMHCpred, NetMHCIIpan.
- **Why:** MHCflurry is open, pip-installable, scriptable, and needs no academic licence (NetMHCpan requires one). The assignment lists it as an acceptable class-I predictor.
- **Caveat:** A single predictor; consensus across tools would be more robust.

**3.7 HLA panel = 3 common class-I alleles (HLA-A\*02:01, A\*01:01, A\*03:01) + 2 class-II (DRB1\*15:01, \*07:01).**
- **Why:** The assignment's "Option A" fixed panel of common alleles, used because we do not have patient-specific HLA typing. Every score is reported *with* its allele (Rule 7).
- **Caveat:** Fixed panel ≠ patient HLA; and only 3 class-I alleles limits population coverage.

**3.8 Binder thresholds: Strong < 50 nM, Weak < 500 nM (IC50).**
- **Why:** Standard immunology cut-offs used across the field.

**3.9 DeltaAffinity = WT_IC50 − Mutant_IC50 (positive ⇒ mutant binds stronger).**
- **Why:** Lower IC50 = stronger binding, so this direction makes "positive = the mutation improves presentation," which is the neoantigen-favourable case. Direction is documented on every row.

**3.10 Immunogenicity and class-II binding = NA (tools not run at this stage).**
- **Why:** Honesty — we didn't run PRIME/NetMHCIIpan here, so we write NA, never a fabricated number or 0 (Rule 8). (Immunogenicity is added later in Script 16.)

**3.11 Score UNIQUE 9-mers once per allele and cache.**
- **Why:** 2.46M unique sequences instead of ~7M repeats; caching means re-runs don't recompute the slow neural network. Provenance (tool + version) is stamped on every row (Rule 4).

### Script 07 — prioritised shortlist

**3.12 Criteria: Strong/Weak binder + DeltaAffinity > 0 + GeneLevelTPM > 1 + recurrent ≥ 2 tumours.**
- **Why each:** binder = presentable; delta>0 = the mutation improves presentation vs WT; TPM>1 = gene is on; recurrent≥2 = shared, not a one-off.
- **Result:** 1,536 candidates / 966 genes.

**3.13 MutationFrequency = per-mutation (samples with THAT variant), not per-gene.**
- **Why (a bug we caught and fixed):** "How recurrent is KRAS G12D?" must count tumours with G12D specifically (65), not all KRAS mutations (250). §14 means the specific mutation. Fixing this shrank the shortlist from an inflated 51,212 to a meaningful 1,536 and surfaced the real driver neoantigens.

---

## 4. Extension choices (added iteratively with the PI)

### 4A. Driver co-mutation (Scripts 10–12)

**4.1 Gene level, restricted to a curated CRC driver panel (34 genes).**
- **Why:** The PI asked to "count double/triple co-mutation." Doing it on all 17,585 genes is combinatorially impossible and dominated by size artefacts (TTN, MUC16 accrue mutations by length, not selection). A curated driver list keeps it interpretable and biologically meaningful.
- **Caveat:** Inherits the missense-only filter, so tumour suppressors inactivated by truncation are undercounted — clearest tell: **APC shows ~7% here vs its true ~75%** in colorectal cancer.

**4.2 A gene is "mutated" in a sample if ANY of its missense SNVs is present (OR across the gene's rows).**
- **Why:** Co-mutation is a gene-level concept; we don't want to split a gene by which specific residue.

**4.3 Statistics: Fisher's exact (pairs) and Poisson expected-vs-observed (triples), with Benjamini–Hochberg FDR.**
- **Why:** Fisher's exact via the hypergeometric distribution is the exact test for a 2×2 co-occurrence table and controls for each gene's marginal frequency (so a big count from two common genes isn't mistaken for real association). BH-FDR corrects for testing hundreds of pairs. Poisson approximates the rare triple case.
- **Validation:** KRAS/BRAF come out mutually exclusive (same MAPK pathway — textbook), KRAS+PIK3CA co-occurring — both expected.

**4.4 A separate "straightforward" counts view out of N=586 (Script 12).**
- **Why:** The PI wanted the plain "X of 586 tumours" version alongside the odds ratios, because raw counts are more intuitive.
- **Key teaching point we surfaced:** raw counts are confounded by popularity — TP53+KRAS (105) is the biggest count but exactly chance; KRAS+PIK3CA (84) is fewer but genuinely enriched.

### 4B. TMB control (Script 15) — motivated by the SelectSim paper

**4.5 Flag hypermutators at ≥ 200 missense SNVs and re-run co-mutation without them.**
- **Why the method:** The Ciriello-lab SelectSim paper shows co-occurrence is confounded by tumour mutational burden — hypermutated (MSI/POLE) tumours make every gene pair look co-occurrent. Excluding them is the simplest robust version of their per-sample-weighting idea.
- **Why 200:** The burden distribution is **bimodal** — the bulk of tumours sit at 75–124 missense SNVs, then there's a sharp jump (80th pct = 124 → 90th = 708). 200 sits in that valley. It flags 15.5% of tumours, matching the known ~16% hypermutated fraction of TCGA-COAD. The result is robust: any cut-off from 200–500 flags essentially the same tumours (≥300 → 14.7%, ≥500 → 13.8%).
- **Result:** The 71× APC+ATM+POLE triple → **0 tumours** (pure TMB artefact — 87% of POLE tumours are hypermutated, so removing them removes POLE); KRAS+PIK3CA survives and strengthens (OR 2.4→2.9).
- **Caveat:** Removing 15.5% of samples also cuts statistical power for genes concentrated in hypermutators; the definitive artefact is the triple hitting zero.

### 4C. Clonality (Script 13)

**4.6 VAF = t_alt_count / t_depth; clonal if median VAF ≥ 0.25.**
- **Why:** Clonal neoantigens (in all tumour cells) are far better targets — a subclonal one lets the tumour escape by dropping the minority clone. VAF (fraction of reads carrying the mutation) is the standard proxy. A heterozygous clonal mutation in a pure diploid tumour gives VAF ~0.5; purity <1 lowers it, hence 0.25.
- **Result:** 60% of mutations clonal; every top driver neoantigen is clonal (as expected — drivers are early events).
- **Caveat:** Raw VAF is a proxy; rigorous cancer-cell fraction needs tumour purity + copy number (e.g. ABSOLUTE), which we don't use.

### 4D. Practical filter + patient coverage (Script 14)

**4.7 "Practical neoantigen" = mutant binder (<500 nM) + wild-type NON-binder (≥500 nM) + TPM ≥ 10 + clonal + recurrent ≥ 2.**
- **Why each, especially the new ones:**
  - *Wild-type non-binder (differential agretopicity):* the single most important addition. If the WT peptide is also presented, the immune system is *tolerant* to it (thymic selection), so the near-identical mutant won't provoke a response. We require the mutant presented and the WT not.
  - *TPM ≥ 10 (not >1):* "expressed" should mean *abundantly* made, not barely detectable.
- **Result:** 1,536 → 301 practical neoantigens.

**4.8 Patient coverage via greedy set-cover; track ≥1 and ≥2 epitopes per tumour.**
- **Why:** For a shared/off-the-shelf vaccine you want the fewest neoantigens covering the most patients, ideally hitting each with ≥2 epitopes (so escape from one still leaves another). Greedy set-cover is the standard approximation.
- **Key finding:** shared missense neoantigens cover ~43% of tumours with ≥1 epitope but very few with ≥2 — because the strongest shared targets (KRAS G12x) are mutually exclusive. This quantifies *why* single neoantigens are weak (the PI's original point) and argues for combining shared + personalized, or adding frameshift neoantigens.

### 4E. Immunogenicity (Script 16)

**4.9 Model = Calis et al. 2013 / IEDB class-I immunogenicity (pure-Python).**
- **Why:** It's the standard sequence-based immunogenicity predictor, needs no external tool, and fills the previously-NA immunogenicity column — turning a *binding* predictor into a *neoantigen* predictor. Higher score = more immunogenic.
- **Implementation choices (faithful to IEDB):** amino-acid immunogenicity scale × position weights; **mask positions 1, 2 and the C-terminus** (the IEDB default — these are HLA anchors that affect binding, not T-cell recognition); central positions 4–6 carry most weight (main TCR contacts).
- **Validation:** aromatic-rich `WWFWWFWWF` → +0.90; K/S-rich `KSKSKSKSK` → −0.88 — matches the paper's stated enrichments.
- **Added interpretation:** mutant-vs-WT immunogenicity delta, and whether the mutation sits at an **anchor** (affects binding) or a **TCR-contact** position (affects recognition).
- **Key finding:** 179 of 301 practical candidates carry the mutation at an anchor — they work via *agretopicity* (creating binding), not by raising immunogenicity; only 29 have a mutation that increases immunogenicity above WT. Both mechanisms are legitimate.
- **Caveat:** Validated on 9-mers only; a modest sequence heuristic, one input among several (exact values checkable at the IEDB web tool).

### 4F. Composite quality score (Script 17)

**4.10 Five non-redundant axes: binding, immunogenicity, expression, clonality, recurrence. Tumour-specificity stays a gate, not a scored axis.**
- **Why non-redundant:** each captures a distinct requirement; keeping WT-non-binder as a prerequisite gate avoids double-counting binding.

**4.11 Percentile-rank normalization + equal weights + mean.**
- **Why:** Percentiles put all axes on a comparable 0–1 scale robust to units/outliers; equal weights are the least arbitrary default; every axis percentile is reported so the reader can re-weight. A sensitivity check (binding+immuno-heavy weighting) shows the very top is stable (6/10 top-10 overlap).
- **A FLAW WE FOUND AND FIXED (tell this as a strength):** our first version used percentile ranking on *all* axes, including recurrence — which is the *wrong* transform there, because recurrence is floor-dominated (257 of 301 candidates sit at exactly 2 tumours). Percentile ranking inflated a 3-tumour candidate to the 91st percentile and compressed the real gap between KRAS G12V (56 tumours) and a 3-tumour candidate. The symptom: **NRAS G13R (only 3 tumours) initially ranked #1**. We diagnosed this and **corrected it** by scoring recurrence on a **log-absolute** scale instead. After the fix, **KRAS G12V ranks #1** (the genuine best shared target, 56 tumours), KRAS G12A #6, and NRAS G13R drops to #4 — still a strong *per-patient* neoantigen (clonal, expressed, creates the epitope de novo: WT 26,000 nM → mutant 106 nM) but correctly no longer a top *shared* target. The other four axes keep percentile normalization (their distributions are well-spread). For a *personalized*-vaccine goal you would simply drop the recurrence axis.

### 4G. BigMHC Deep-Learning Predictor Transition (Script 06, 20)

**4.12 Transition from MHCflurry IC50 (nM) to BigMHC Eluted-Ligand Presentation (`BigMHC_EL`) and Immunogenicity (`BigMHC_IM`).**
- **Why BigMHC:** BigMHC (Albert et al., *Nature Machine Intelligence* 2023) is a modern deep-learning neural network trained on mass-spectrometry eluted-ligand datasets and transfer-learned for T-cell immunogenicity. It directly reports presentation probability `BigMHC_EL` in $[0, 1]$ (where $\ge 0.50$ indicates presentation) and immunogenicity score `BigMHC_IM`.
- **Direction & Thresholds:**
  - *Presentable:* Mutant `BigMHC_EL >= 0.50` (higher probability is better).
  - *Differential Agretopicity:* Mutant `BigMHC_EL >= 0.50` AND Wild-type `BigMHC_EL < 0.50` ($\Delta\text{Presentation} = \text{Mut\_EL} - \text{WT\_EL} > 0$).
- **Result:** Fully refactored pipeline scripts (06 to 19) operate directly on BigMHC presentation probabilities with zero reliance on MHCflurry.

### 4H. Genome-Wide All-by-All Co-Mutation vs. Curated Driver Panel (Script 10 vs. `scratch/reverify_genome_wide.py`)

**4.13 Why restrict driver co-mutation to 35 curated drivers instead of all 17,585 genes?**
- **Why:** In an unfiltered all-by-all analysis across 17,585 genes ($154.6 \times 10^6$ pairs), giant genes (*TTN*, *MUC16*, *SYNE1*, *OBSCN*) dominate raw overlap counts ($n_{\text{both}} = 68\text{--}88$) and show apparent statistical significance ($\text{OR} = 2.5\text{--}4.4, p < 10^{-6}$) in the full cohort ($N=586$).
- **The TMB Confounding Discovery:** Re-verification in `scratch/reverify_genome_wide.py` revealed that these giant gene co-occurrences are **surrogate markers of Tumour Mutational Burden (TMB) / the hypermutator phenotype** ($N=90$, $\ge 200$ SNVs). In hypermutators, almost every long gene is mutated in 100% of tumours, creating strong statistical co-occurrence across $N=586$.
- **When Hypermutators are Removed ($N=495$ Standard Tumours):** Every giant passenger pair (`TTN + MUC16`, `TTN + SYNE1`, `TTN + OBSCN`, `MUC16 + FAT4`) collapses to an Odds Ratio $\approx 1.1\text{--}1.3$ with $p > 0.10$ (completely non-significant under independent chance). Conversely, genuine oncogenic driver partnerships like **`KRAS + PIK3CA`** strengthen from $\text{OR} = 2.39$ to $\text{OR} = 2.94$ ($p = 9.36 \times 10^{-7}$), emerging empirically as the **#1 most recurrent, statistically significant co-occurring gene pair in the entire human genome** in standard colorectal cancer.

### 4I. Variant-Level (Locus-Specific) Co-Mutation & Hotspot Architecture (`scratch/variant_level_comutation.py`)

**4.14 Why does *TTN* lack recurrent hotspots compared to driver genes?**
- **Why:** To test whether any specific locus on *TTN* acts as an oncogenic hotspot or is highly co-selected, we evaluated all 153,996 distinct missense mutations at the variant level (`Gene p.Change`) in `scratch/variant_level_comutation.py`.
- **Passenger vs. Driver Architecture:** 
  - **Driver Hotspots (*KRAS*, *TP53*, *PIK3CA*, *BRAF*):** Under positive evolutionary selection, independent tumours repeatedly mutate the exact same amino acid codon (`KRAS p.G12D` in 65 tumours [11.1%], `KRAS p.G12V` in 56 tumours [9.6%], `BRAF p.V640E` in 51 tumours [8.7%]).
  - **Passenger Dispersion (*TTN*):** Out of 678 distinct missense mutations on *TTN*, exactly **1 variant** (`p.E26065K`) occurs in 3 tumours (0.51%), **14 variants** occur in 2 tumours, and **663 variants (97.8%) are completely private** (occurring in exactly 1 tumour).
- **Variant-Level Co-Occurrence:** In standard tumours ($N=495$), `TTN p.E26065K` has **zero co-occurrences** ($n_{\text{both}} = 0$). In contrast, specific driver variant pairs show genuine locus-specific co-selection (`KRAS p.G12V + TP53 p.R175H` in 10 standard tumours, $\text{OR} = 3.52, p = 3.79 \times 10^{-3}$; `KRAS p.G13D + PIK3CA p.E545K` in 7 standard tumours, $\text{OR} = 4.16, p = 6.33 \times 10^{-3}$).

### 4J. MHC Class II (15-mer) Binding Prediction & Class-Aware Dual CD4+/CD8+ Vaccine Cocktail Design (`scripts/22_class2_binding_and_coverage.py`)

**4.15 Why expand the greedy set-cover algorithm to include Class II (15-mer) epitopes alongside Class I (9-mer) epitopes?**
- **Biological Rationale:** Effective and durable antitumor immunity requires priming both CD8+ cytotoxic T-cells (via MHC Class I 9-mer peptides) and CD4+ helper T-cells (via MHC Class II 15-mer peptides). CD4+ T-cell help is biologically necessary to promote dendritic cell licensing, produce IFN-$\gamma$/IL-21, and prevent CD8+ T-cell exhaustion.
- **Two Mechanisms of Class II Neoantigenicity:**
  1. **TCR-Contact Loop Novelty (1,339 candidates):** Mutations at solvent-exposed TCR recognition positions (P2, P3, P5, P7, P8) within the presented 9-mer core. Both wild-type and mutant peptides bind MHC-II equivalently, but the amino acid change creates a novel, immunogenic T-cell visible surface.
  2. **Differential Agretopicity (288 candidates):** Mutations at MHC anchor positions (P1, P4, P6, P9) confer de novo or substantially improved binding to `HLA-DRB1*15:01` or `*07:01` (`Mutant IC50 < 500 nM` and `WT IC50 >= 500 nM` or `Mutant < WT × 0.5`).
- **Discovered Class II Shortlist:** Scanning all 15-mers against an anchor-position PSSM (P1, P4, P6, P9 only) identified **1,627 Practical Class II Neoantigens** (`results/practical_class2_neoantigens.tsv`).
- **Methodological Fixes Applied After Code Review:**
  1. **Mutation-in-core guard:** Only 9-mer cores where the mutation falls within positions 1–9 are evaluated. Previously, 40% of binding predictions used cores where the mutation sat in the flanking region of the 15-mer.
  2. **Aligned WT comparison:** The wild-type IC50 is now computed from the **same sliding window** (same `MutPos`) as the mutant, not the minimum across all WT 15-mers.
  3. **Class-aware coverage tracking:** The dual greedy set-cover tracks Class I and Class II hits **separately** per tumour, enabling measurement of genuine CD4+/CD8+ synergy ($\ge 1$ Class I AND $\ge 1$ Class II per tumour) rather than simply counting $\ge 2$ mutations of any class.
- **Stated Limitation:** The MHC-II binding predictions use a simplified 4-anchor-position PSSM, not a full SMM-align or NetMHCIIpan predictor. Non-anchor positions (P2, P3, P5, P7, P8) are not scored. This means TCR Contact Loop candidates have identical mut/wt binding scores by construction — the mechanism classification is correct but the binding affinity does not differentiate them.
- **Efficacy Comparison (30 Epitopes across $N=586$ Tumours):**
  - **Class I Only (9-mer):** `60.9%` ($\ge 1$ epitope) / `18.6%` ($\ge 2$ epitopes).
  - **Class II Only (15-mer):** `57.5%` ($\ge 1$ epitope) / `15.0%` ($\ge 2$ epitopes).
  - **Dual Optimal Cocktail (20 Class I + 10 Class II):** **`63.1%` ($\ge 1$ epitope) / `20.0%` ($\ge 2$ epitopes)** / **`10.4%` genuine CD4+/CD8+ synergy** ($\ge 1$ Class I AND $\ge 1$ Class II per tumour) (`figures/22_dual_class1_class2_vaccine_coverage.png`).

---

## 5. Cross-cutting choices (reproducibility & conventions)

- **Numbered scripts (01→07, then extensions 10–17)** in execution order — each reads the previous step's output; the whole pipeline reproduces from code (Assignment Rule 3).
- **All primary outputs tab-delimited `.tsv`**, no manual spreadsheet editing (Rules 1–2).
- **Missing values = NA, never 0**, everywhere (a 0 makes a false claim).
- **Software + version stamped** into outputs; access dates and sources in the README (Rule 4).
- **Sample columns detected by the `TCGA` prefix** (not a hard-coded index) so adding columns like GeneLevelTPM_SD doesn't silently break downstream scripts.
- **Sequence rationale:** mutations (01) and expression (02) are independent and come first; integration (03) needs both; QC (04) needs the integrated matrix; peptides (05) need the MAF + proteome; binding (06) needs peptides; prioritisation (07) needs binding + integration; every extension builds on these.

---

## 6. Every threshold / number, with its one-line justification

| Choice | Value | Why |
|---|---|---|
| Reference genome | GRCh38 | GDC-harmonised; no coordinate mixing (Rule 5) |
| Mutation filter | protein-coding, PASS, missense, SNV | only these can make a clean new peptide |
| Sample-barcode trim | 15 chars | TCGA sample level; matches expression IDs |
| RSEM→TPM | rescale each sample to 1e6 | defining property of TPM (Rule 6) |
| Duplicate genes | summed | expression is additive |
| Missing values | NA (never 0) | 0 would falsely claim "silent" |
| GeneLevelTPM | median across tumours | robust to outliers |
| Dispersion reported | SD + IQR + n | judge how trustworthy the median is |
| Peptide lengths | 9-mer + 15-mer | MHC-I (CD8) + MHC-II (CD4) |
| Transcript rule | VEP canonical/MANE | consistent, reproducible |
| Ref-AA mismatch | excluded (audited) | never guess (2,964 dropped) |
| HLA panel | A\*02:01, A\*01:01, A\*03:01 (+2 class II) | common alleles, Option A |
| Binder cut-offs | Strong <50 nM, Weak <500 nM | standard immunology |
| Delta direction | WT − Mut (positive = mutant stronger) | lower IC50 = stronger binding |
| Shortlist | binder + delta>0 + TPM>1 + recur≥2 | presentable, tumour-specific, on, shared |
| MutationFrequency | per specific mutation | §14 means the variant, not the gene |
| Co-mutation genes | 35 curated CRC drivers | avoid gene-size artefacts (TTN) |
| Genome-wide co-mutation | 17,585 genes across N=586 vs. N=495 | proves giant gene pairs are TMB/hypermutator artefacts; KRAS+PIK3CA is #1 in standard CRC |
| Variant-level hotspots | 153,996 distinct missense SNVs | proves TTN lacks driver hotspots (97.8% private) and has zero co-occurrence in standard CRC |
| Co-mutation stats | Fisher (pairs), Poisson (triples), BH-FDR | control for marginal frequency + multiple testing |
| Hypermutator cut-off | ≥ 200 missense SNVs | valley of bimodal dist; 15.5% ≈ known ~16% |
| Clonal cut-off | median VAF ≥ 0.25 | proxy for clonal (VAF ~0.5 het, purity-adjusted) |
| Practical filter | Mut<500, **WT≥500**, TPM≥10, clonal, recur≥2 | presentable + tumour-specific + made + truncal + shared |
| Coverage | greedy set-cover, ≥1 & ≥2 epitopes | fewest neoantigens covering most patients |
| Immunogenicity | Calis/IEDB, mask pos 1/2/C-term | anchors affect binding, not recognition |
| Composite | 5 axes, percentile-rank, equal weights | comparable scale, least-arbitrary weights |

---

## 7. Questions the professor may ask — and crisp answers

- **"Why missense only?"** Clean single-amino-acid change → a defined new peptide, and the simplest exact peptide logic. Frameshift/indel neoantigens (big in the MSI subset) are a documented extension.
- **"Is your TPM really TPM?"** It's a documented per-million rescaling of RSEM — a standard proxy, but not fully length-normalised. Native GDC STAR-Counts TPM would be the strict version; disclosed as a limitation.
- **"Does expression prove the mutant is made?"** No — it's gene-level, so it shows the gene is on, not that the mutant *allele* is transcribed. Allele-specific expression needs read-level RNA (RNA VAF), beyond our data.
- **"Why the median, not the mean?"** Robust to a few extreme tumours; and we also report SD/IQR so you can see when the median is unreliable.
- **"Why these three HLA alleles?"** Assignment Option A — common class-I alleles, because we lack patient HLA typing. Fixed panel is a stated limitation.
- **"Your biggest co-mutation count is TP53+KRAS — real?"** No — 105 observed vs ~103 expected by chance; it's just two common genes overlapping. The genuine one is KRAS+PIK3CA (84 vs ~61 expected).
- **"How do you know the co-mutations aren't burden artefacts?"** We excluded hypermutators and re-ran: the 71× APC+ATM+POLE triple drops to 0 (pure TMB), KRAS+PIK3CA survives and strengthens.
- **"Why is APC only 7% when it's ~75% in colorectal cancer?"** Because APC is inactivated by *truncating* mutations, which our missense-only filter removes — a limitation of the co-mutation view, not an error in the counts.
- **"How is your final list ranked, and did you sanity-check it?"** A composite of five equal axes; KRAS G12V ranks #1. We caught a normalization flaw during that check: percentile-ranking the floor-dominated recurrence axis had briefly put a 3-tumour candidate (NRAS G13R) at #1. We fixed it by scoring recurrence on a log-absolute scale, which restored KRAS G12V to #1. NRAS G13R is still an excellent *per-patient* target (creates the epitope de novo) but correctly no longer a top *shared* one.
- **"Binding vs immunogenicity — do you distinguish them?"** Yes, separate columns: binding = presented; immunogenicity (Calis) = T-cell-visible. A strong binder isn't automatically immunogenic, and 179/301 of our candidates work via anchor-driven *binding* (agretopicity), not raised immunogenicity.
- **"Have you validated anything?"** Not experimentally (predictions only). But internally: KRAS G12D reproduces the textbook peptide exactly, the driver genes come out at expected frequencies, KRAS/BRAF are mutually exclusive, and the prioritised candidates are the neoantigens the field already pursues.
- **"Why restrict driver co-mutation to 35 curated drivers instead of testing all 17,585 genes?"** When we test all 17,585 genes across the full cohort ($N=586$), giant genes (*TTN*, *MUC16*, *SYNE1*, *OBSCN*) dominate raw counts and statistical co-occurrence ($\text{OR} = 2.5\text{--}4.4$) because they act as surrogate detectors of the **hypermutator (MSI-H/POLE) phenotype** ($N=90$, $\ge 200$ SNVs). Once hypermutated tumours are removed ($N=495$ standard tumours), every giant passenger pair collapses to non-significance ($\text{OR} \rightarrow 1.2, p > 0.10$), while genuine oncogenic driver partnerships like **`KRAS + PIK3CA`** strengthen from $\text{OR} = 2.39$ to $\text{OR} = 2.94$ ($p = 9.36 \times 10^{-7}$), emerging empirically as the #1 most recurrent, statistically significant co-occurring gene pair in the entire human genome.
- **"Could a specific variant locus on TTN be an oncogenic driver hotspot or highly co-selected?"** We tested all 153,996 distinct missense mutations at the variant level (`Gene p.Change`) across the cohort. Out of 678 distinct missense mutations on *TTN*, 97.8% (663 variants) are completely private (occur in exactly 1 tumour), and the most recurrent variant (`TTN p.E26065K`) occurs in only 3 tumours (0.51%) with zero co-occurrences in standard non-hypermutated tumours ($N=495$). In contrast, true driver genes (*KRAS*, *TP53*, *PIK3CA*, *BRAF*) exhibit strong positive selection at specific catalytic codons (`KRAS p.G12D` in 65 tumours [11.1%], `TP53 p.R175H` in 39 tumours [6.7%]) and demonstrate genuine locus-specific co-selection in standard tumours.
- **"Why expand the greedy set-cover algorithm to include Class II (15-mer) epitopes alongside Class I (9-mer) epitopes?"** A clinical cancer vaccine cocktail must prime both CD4+ helper T-cells (via MHC Class II 15-mers) and CD8+ cytotoxic T-cells (via MHC Class I 9-mers). CD4+ T-cell help is biologically essential to prevent CD8+ exhaustion and sustain antitumor immunity. We use a **class-aware greedy set-cover** (`scripts/22_class2_binding_and_coverage.py`) that tracks Class I and Class II hits **separately** per tumour. This enables measurement of genuine CD4+/CD8+ synergy ($\ge 1$ Class I AND $\ge 1$ Class II per tumour) rather than simply counting $\ge 2$ mutations of any class. The dual cocktail achieves **63.1%** overall coverage and **10.4% genuine CD4+/CD8+ synergy** (61 tumours receiving both CD4+ helper and CD8+ cytotoxic targeting). **Stated limitations:** MHC-II binding uses a simplified 4-anchor-position PSSM (P1, P4, P6, P9 only), not a full NetMHCIIpan predictor; TCR Contact Loop candidates have identical mut/wt binding scores by construction since non-anchor positions are not scored.

---

## 8. Honest limitations (say these before you're asked)

1. Gene-level, not allele-specific, expression.
2. RSEM→TPM is a per-million rescaling, not length-normalised TPM.
3. Missense SNVs only — no frameshift/indel neoantigens.
4. Fixed 3-allele HLA panel, not patient-specific typing.
5. Clonality from raw VAF (no purity/copy-number correction).
6. Class-II binding uses a simplified 4-anchor-position PSSM (P1, P4, P6, P9 only), not a full SMM-align or NetMHCIIpan predictor. Non-anchor positions (P2, P3, P5, P7, P8) are not scored, so TCR Contact Loop candidates have identical mut/wt binding scores by construction. Immunogenicity (Calis heuristic, Script 16) is a sequence-based predictor without experimental calibration.
7. Composite ranking is multi-criteria and weight-sensitive (top is stable; mid-ranks shuffle) — we found and corrected a recurrence-normalization flaw, and report all axes so it can be re-weighted.
8. Everything is computational prediction, not experimental validation.

---

## 9. HOW each step is computed — rigorous methods

**Environment.** Python 3 with pandas / numpy / matplotlib; MHCflurry for class-I binding; all statistics implemented in **pure Python** (`math.comb`, `math.exp`) — no scipy dependency, so every test is transparent. Outputs are tab-delimited. Large tables (the 154k×586 matrix, the 16M-row neoantigen table) are **streamed line-by-line** rather than loaded whole, to bound memory. Downstream scripts locate sample columns by the `TCGA` barcode prefix, so inserting metadata columns never shifts the parsing.

### 01 — mutation-by-sample matrix
- Read MAF (`pandas.read_csv`, tab-sep, `dtype=str`, only needed columns).
- Four filters as boolean masks combined with `&`: `BIOTYPE=="protein_coding"`, `GDC_FILTER ∈ {"","PASS"}`, `Variant_Classification=="Missense_Mutation"`, `Variant_Type=="SNP"`.
- Assembly guard: `assert NCBI_Build.unique()==["GRCh38"]`.
- Sample ID: `Tumor_Sample_Barcode.str[:15]`. Distinct-mutation key: concatenate `gene⇥HGVSc⇥HGVSp`.
- **Matrix build (memory-safe, avoids the pandas pivot blow-up):** map each unique mutation→row index and each sample→column index (dicts); allocate `M = np.zeros((n_mut, n_samples), dtype=uint8)`; compute integer index arrays `rows`, `cols` via `.map(dict)`; set `M[rows, cols] = 1` in one vectorised assignment; concat key columns + `DataFrame(M)`; sort by row-sum; write.

### 02 — RSEM → TPM
- Read genes×samples; drop rows with empty `Hugo_Symbol`; collapse duplicates with `groupby("Hugo_Symbol").sum(min_count=1)`.
- **Conversion:** `col_totals = expr.sum(axis=0, skipna=True)`; `tpm = expr.divide(col_totals, axis=1) * 1e6`. **Check:** `assert np.allclose(tpm.sum(axis=0), 1e6)`. Missing cells stay NA.

### 03 — integration
- Per gene from the TPM matrix: `median(axis=1, skipna=True)`, `std(axis=1)`, `quantile(.75,axis=1) − quantile(.25,axis=1)` (IQR), `notna().sum(axis=1)` (n).
- Build dicts gene→median, gene→SD; **stream** deliverable 01 line-by-line, inserting `f"{median:.4f}"` and `f"{SD:.4f}"` (or `"NA"`) as columns 4–5, then the sample fields; write.

### 05 — variant → peptide
- **FASTA parse** → dict `{UniProt_accession: sequence}` (header `>sp|ACC|NAME` → key = ACC).
- Parse each variant: `RefAA,AltAA = Amino_acids.split("/")`; `ProtPos = int(Protein_position.split("/")[0])`; `UniProt = SWISSPROT.split(".")[0]`.
- **Window enumeration (core algorithm)** for length L, 1-based mutation position m:
  `for start in range(max(1, m−L+1), min(len(seq)−L+1, m)+1): pep = seq[start−1 : start−1+L]; MutPos = m − start + 1`.
- **Mutant protein:** `mut_seq = seq[:m−1] + AltAA + seq[m:]`; mutant window = `mut_seq[start−1:start−1+L]`.
- **Verification (guards/asserts):** `seq[m−1]==RefAA` else write to audit and skip; `m ≤ len(seq)`; `len(pep)==L`; `[i for i in range(L) if wt[i]≠mut[i]] == [MutPos−1]`. Peptides are streamed to disk as generated (millions of rows).

### 06 — MHC binding
- `Class1PresentationPredictor.load()`.
- Unique 9-mers = `set(...)`, restricted to the 20 standard amino acids (`set(pep) ⊆ AAs`; nonstandard → reported NA).
- Per allele: `pred.predict(peptides, alleles=[a], include_affinity_percentile=True)` → `affinity` (IC50 nM), `affinity_percentile`, `presentation_score`; concatenate the three alleles; **cache** to TSV + a `.version` sidecar (records MHCflurry version). Build lookup `{(peptide,allele): (aff,rank,pres)}`.
- **Two-pass streaming assembly:** pass 1 records, per key `(gene, proteinChange, protPos, MutPos, allele)`, the Mutant and WildType affinities; pass 2 writes each row, computing `DeltaAffinity = WT_aff − Mut_aff` and `BinderClass` (<50 Strong / <500 Weak / else Non-binder). Any missing value → `"NA"`.

### 07 — shortlist
- Per-mutation frequency map from deliverable 03: `(gene, AminoAcidChange) → Σ("1"s across sample columns)`.
- Stream deliverable 04; keep rows where `PeptideType=="Mutant"`, `PeptideLength=="9"`, `BinderClass ∈ {Strong,Weak}`, `DeltaAffinity>0`, `GeneLevelTPM≥1`, and per-mutation `frequency≥2`; sort (Strong first, then lowest IC50, then highest recurrence); write.

### Co-mutation (10–12)
- **Gene presence:** for each driver gene, `present[g] = OR over its mutation rows` of the 0/1 sample vector (numpy int8, `|=`). Pair co-occurrence `a = (present[A] & present[B]).sum()`.
- **Fisher's exact via hypergeometric** with a, nA, nB, N: right tail `P(X≥a)=Σ_{x=a}^{min(nA,nB)} C(nA,x)·C(N−nA, nB−x)/C(N, nB)` (`math.comb`); left tail symmetric; **odds ratio** `(a+.5)(d+.5)/((b+.5)(c+.5))` (Haldane correction).
- **BH-FDR:** sort p-values; `q_i = min over descending rank of (p_i·m/k)`.
- **Triples:** expected `= N·(nA/N)(nB/N)(nC/N)`; enrichment = observed/expected; **Poisson tail** `P(X≥obs)=1 − Σ_{k<obs} e^{−λ}λ^k/k!`, `λ = expected`.
- Counts view (12): same presence arrays, reported as `both / N` and `%`.

### TMB stratification (15)
- **Burden** per sample = column-sum over **all** mutation rows (every missense SNV). `hyper = burden ≥ 200`.
- Re-run the identical Fisher pair analysis on the boolean subset `present[g][~hyper]` with `N = (~hyper).sum()`; compare odds ratios and verdicts all-samples vs non-hyper.

### Clonality (13)
- Per MAF row `VAF = t_alt_count / t_depth` (drop depth 0 / NA). `groupby(gene, HGVSp_Short).VAF → median, mean, n`. `Clonal if median VAF ≥ 0.25`.

### Practical filter + coverage (14)
- Recover `WT_IC50 = Mut_IC50 + DeltaAffinity` from the shortlist.
- **Filter:** `Mut_IC50<500 ∧ WT_IC50≥500 ∧ GeneLevelTPM≥10`; join clonality; require `Clonal ∧ frequency≥2`; keep the lowest-IC50 peptide per mutation.
- **Sample sets:** stream deliverable 03 → `(gene, aachange) → boolean sample array` (OR of matching rows).
- **Greedy set-cover:** while candidates remain, pick `argmax_k (sampset_k ∧ ¬covered).sum()`; update `covered |= sampset_k` and `hit_count += sampset_k` (to count tumours with ≥2 selected epitopes); record cumulative % at each step.

### Immunogenicity (16)
- `score(pep) = Σ_{i ∉ {1,2,9}} POS_WEIGHT[i]·IMMUNOSCALE[pep_i]` (Calis constants; positions 1, 2 and C-terminus masked; central 4–6 weighted highest).
- WT peptide via the exact `MutPos` from the shortlist: `wt = mut[:MutPos−1] + RefAA + mut[MutPos:]`; `Δimmuno = score(mut) − score(wt)`; anchor flag = `MutPos ∈ {1,2,9}`.

### Composite (17)
- Four axes normalised by percentile rank: `pandas.rank(pct=True)` (binding on `−IC50` so lower IC50 → higher).
- Recurrence normalised **log-absolute**: `log10(1+cov) / log10(1+max_cov)` (deliberately not percentile — see 4.11).
- `CompositeScore = mean(5 axes)`; sort descending. **Sensitivity:** recompute with weights 0.30/0.30/0.133×3 (binding+immuno-heavy); report top-10 overlap.
