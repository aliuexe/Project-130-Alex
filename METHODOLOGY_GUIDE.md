# Project 130 — Methodology Guide

### Integrating cancer mutations, gene expression, and neoantigen prediction in colorectal cancer (TCGA-COAD)

*A step-by-step explanation written for a second-year undergraduate with a biology background and little computer-science experience. It explains **what** each step does, **why** we do it, and **how** it works — including the biology and the computing.*

---

## How to use this guide

Read Parts A–B first for the big picture and vocabulary. Then Part D walks through the pipeline one script at a time; each step tells you the biological purpose, the method, and the reasoning. Part E is a gentle introduction to the computing ideas if the code is unfamiliar. Parts F–H help you interpret results, understand the limitations, and explain the project out loud.

Alongside this guide, the folder `scripts_annotated/` contains the same code as `scripts/` but with a plain-language comment on almost every line. Read the guide and the annotated code side by side.

---

## Part A — The biological question

### What is a neoantigen, and why do we care?

Your immune system is constantly checking your cells for signs of "foreignness". Inside every cell, proteins are chopped into short fragments called **peptides**. Special molecules called **MHC** (in humans, **HLA**) grab some of these peptides and hold them up on the cell surface, like showing ID cards. Immune cells called **T cells** inspect these displayed peptides. If a peptide looks foreign, T cells can destroy the cell.

Cancer cells carry **mutations** — changes in their DNA that healthy cells don't have. Some mutations change a protein's amino-acid sequence, producing peptides that the immune system has never seen before. A mutated peptide that is displayed by HLA and can provoke a T-cell response is called a **neoantigen** ("neo" = new). Neoantigens are extremely important in modern cancer immunotherapy: they are tumour-specific flags, so a therapy aimed at them attacks the cancer while sparing healthy cells. Personalised cancer vaccines are built by predicting a patient's neoantigens computationally — which is exactly what this project does at the level of a whole cancer type.

### The central dogma, in one paragraph

DNA is the cell's instruction manual. A stretch of DNA that codes for a protein is a **gene**. To use a gene, the cell first copies it into **RNA** (this is **transcription**), and then reads the RNA to build a **protein** (this is **translation**). A protein is a chain of **amino acids**; there are 20 kinds, each abbreviated to one letter (for example, `G` = glycine, `D` = aspartate). The sequence of amino acids determines the protein's shape and job.

A **missense mutation** changes one DNA letter in a way that swaps one amino acid for another. For example, the famous colorectal/pancreatic cancer mutation **KRAS G12D** means "at position 12 of the KRAS protein, the normal glycine (G) has become aspartate (D)". That single change creates new peptides that were not in the healthy protein — potential neoantigens.

### Why colorectal cancer?

Colorectal cancer is one of the best-studied cancers. Its mutations are dominated by a handful of well-known **driver genes** — `APC`, `TP53`, `KRAS`, `PIK3CA`, `BRAF`, `SMAD4` — and a subset of tumours are "hypermutated", giving colorectal cancer a relatively high neoantigen load. Because the biology is well understood, it's a great test case: if our pipeline is correct, it should rediscover these known drivers, which gives us confidence in the parts we can't check by eye.

---

## Part B — Key terms (quick glossary)

| Term | Plain-language meaning |
|---|---|
| **Somatic mutation** | A DNA change acquired during life, present only in the tumour (not inherited). |
| **Germline mutation** | An inherited DNA change, present in every cell. (We are **not** studying these.) |
| **TCGA** | The Cancer Genome Atlas — a huge public database of tumour genomes. |
| **GDC** | Genomic Data Commons — the NCI portal that serves TCGA data. |
| **TCGA-COAD** | The colon-adenocarcinoma (colorectal) cohort within TCGA. |
| **MAF** | Mutation Annotation Format — a big table listing every mutation found, one row per mutation-in-a-sample. |
| **Reference genome (GRCh38/hg38)** | The standard "map" of the human genome that all positions are measured against. We use GRCh38 throughout and never mix it with the older GRCh37. |
| **RNA-seq** | A method that measures how much RNA each gene makes — a proxy for how "switched on" (expressed) a gene is. |
| **TPM** | Transcripts Per Million — a normalised expression unit that is comparable across samples (each sample's values sum to 1,000,000). |
| **HLA / MHC** | Molecules that display peptides on the cell surface. Class I shows ~9-amino-acid peptides to CD8+ "killer" T cells; class II shows ~15-amino-acid peptides to CD4+ "helper" T cells. |
| **Peptide (k-mer)** | A short chain of amino acids. A "9-mer" is 9 amino acids long. |
| **Wild-type** | The normal, non-mutated version (of a peptide or protein). |
| **Binding affinity (IC50, in nM)** | How tightly an HLA molecule holds a peptide. **Lower nM = stronger binding.** |
| **Neoantigen** | A mutated peptide that is displayed by HLA and may trigger a T-cell response. |
| **HGVS notation** | A standard way to write mutations, e.g. `c.35G>A` (DNA change) and `p.G12D` (protein change). |
| **Matrix** | A grid of numbers (rows × columns). Here, mutations × samples, or genes × samples. |

---

## Part C — The data, and why each piece

We combined three public datasets:

1. **Somatic mutations** — the GDC project-level **Masked Somatic Mutation MAF** for TCGA-COAD, on the GRCh38 reference genome. This tells us *which mutations occur in which tumours* and, helpfully, is already annotated by a tool called **VEP** (Variant Effect Predictor) that worked out each mutation's effect on the protein. Downloaded with the GDC Data Transfer Tool.

2. **Gene expression** — a colorectal RNA-seq dataset from **cBioPortal** (study `coadread_tcga_pan_can_atlas_2018`). This tells us *how strongly each gene is expressed*. It comes as **RSEM** values, which we convert to TPM (see Step 2).

3. **Reference protein sequences** — the **UniProt** reviewed human proteome. We need the full amino-acid sequence of each protein so we can cut out peptides around a mutation.

Everything is open-access, downloaded on the same date, and kept on a single reference genome (GRCh38). Keeping careful track of *where data came from and when* is part of doing reproducible science.

---

## Part D — The pipeline, step by step

The whole analysis is a chain of small, numbered programs ("scripts"). Each one reads the output of the previous step and produces the next. This modular design means any step can be re-run or checked independently — a core principle of reproducible bioinformatics.

```
01 mutations ─┐
              ├─ 03 integrate ─ 04 QC/figures
02 expression ┘
05 peptides ─ 06 HLA binding ─ 07 prioritise candidates
```

### Step 1 — Building the mutation matrix (`01_build_mutation_matrix.py`)

**Goal.** Turn the raw MAF into a clean grid: rows are individual mutations, columns are tumour samples, and each cell is `1` (mutation present) or `0` (absent). This is **Deliverable 01**.

**What we filter, and why.** The raw MAF has ~310,000 mutation records. We keep only those that matter for making neoantigens:

- **Protein-coding genes only** — a mutation can only create a peptide if its gene actually makes a protein.
- **High-confidence ("PASS") variants** — this removes likely sequencing artefacts, so we don't chase false mutations.
- **Missense mutations only** — these swap one amino acid for another, creating a genuinely new protein sequence. (Silent mutations don't change the protein; more disruptive types like frameshifts are a valid but more complex extension we set aside.)
- **Single-nucleotide variants (SNVs)** — the simplest, cleanest kind of change (one DNA letter).

After filtering, ~184,000 records remain, representing **153,996 distinct mutations across 586 tumour samples**.

**Key ideas.**
- A **distinct mutation** is defined by three things together: the gene, the DNA change (`c.35G>A`), and the protein change (`p.G12D`). Two tumours with the same three values share "the same" mutation.
- We write mutations in **HGVS notation**, the universal standard, so anyone can interpret them unambiguously.
- We check that **every mutation uses GRCh38**. Mixing genome versions would put mutations at wrong positions — a serious error — so the script refuses to continue if it finds any other build.

**Why a binary (0/1) matrix?** It's the simplest faithful summary of "who has what". It feeds directly into later steps (how often a mutation recurs, which samples to display in a heat map, etc.).

### Step 2 — Building the expression (TPM) matrix (`02_build_expression_matrix.py`)

**Goal.** Turn RNA-seq measurements into a grid of genes × samples, in **TPM** units. This is **Deliverable 02**.

**Why expression matters.** A mutation is only a plausible neoantigen if the gene is actually being made into protein. A mutation in a gene that is completely switched off in the tumour is unlikely to produce a displayed peptide. RNA-seq expression is our proxy for "is this gene on?".

**The RSEM → TPM conversion (an important, subtle point).** The source file gives **RSEM** values. When we add up all the values in one sample, they total about **18,000,000**, not 1,000,000 — so these numbers are **not** already TPM. The assignment explicitly forbids labelling non-TPM numbers as "TPM". The fix is a real, documented conversion: for each sample, divide every gene's value by that sample's total and multiply by 1,000,000. Now each sample's values sum to exactly 1,000,000, which is the defining property of TPM, and the values are comparable across samples.

> **Honest caveat (worth knowing).** Strictly speaking, true TPM also involves a "gene-length normalisation" step. Our rescaling makes the numbers sum to a million but does not re-do length normalisation, so this is best described as a documented per-million rescaling of RSEM — a very common and defensible approach for ranking gene expression, but not a fully length-normalised TPM. A stricter alternative would be to download native TPM values from GDC's newer STAR-Counts files. This is disclosed in the report's limitations.

**Handling duplicates and gaps.** If a gene name appears on more than one row, we sum the rows (expression is additive) so each gene has exactly one row. Cells with no measurement are recorded as **NA** ("not available"), never as `0` — because `0` would falsely claim the gene is silent when in fact we simply don't know.

Result: **20,511 genes × 592 tumour samples**, with about 5.6% of cells missing.

### Step 3 — Integrating mutations and expression (`03_integrate_datasets.py`)

**Goal.** Add expression information to each mutation, producing **Deliverable 03**.

**The problem this step solves.** The patients sequenced for mutations are **not exactly the same** patients sequenced for RNA. So we cannot line them up one-to-one. Instead we summarise expression to a single representative value per gene for the whole cancer type: **GeneLevelTPM = the median TPM across all tumour samples**.

**Why the median (not the mean/average)?** The median is the middle value and is **robust to outliers**: a few tumours with extreme expression won't distort it, whereas the mean would be pulled towards the extremes. This gives a stable "typical expression" for each gene.

**The result** keeps the mutation grid from Step 1 and adds a `GeneLevelTPM` column. Genes with no expression record get `NA` (never `0`).

### Step 4 — Quality control and figures (`04_qc_and_figures.py`)

**Goal.** Confirm the data is sensible and summarise it visually (**Section 8**).

**Why QC is not optional.** Before drawing any biological conclusion, we must be sure the pipeline behaved. QC reports the counts before/after filtering, the number of unique genes and samples, the percentage of missing expression, and confirms that gene identifiers and the genome build are consistent. It also lists the ten most-mutated and ten most-expressed genes.

**A built-in sanity check.** The top mutated genes come out as `TTN`, `MUC16`, `SYNE1`, then `TP53`, `KRAS`, `PIK3CA`. The first few are simply enormous genes (bigger genes collect more random mutations), but seeing the real drivers **TP53, KRAS, PIK3CA** at the expected ranks tells us the filtering and matrix construction are correct. The top-expressed mutated gene is `CEACAM5`, which encodes **CEA** — the classic clinical colorectal tumour marker. These expected results are reassuring.

**The four figures:** (1) most-mutated genes, (2) the spread of GeneLevelTPM, (3) a heat map of recurrent mutations across samples, and (4) mutation frequency versus expression.

### Step 5 — From variants to peptides (`05_annotate_and_generate_peptides.py`)

**Goal.** For every missense mutation, generate all the short peptides that contain the mutated amino acid — both mutant and wild-type — at lengths 9 and 15. This begins the **advanced** component (**Sections 9–10**).

**How we find the protein.** The MAF was already annotated by VEP, so for each mutation it tells us the transcript, the protein, the protein position (e.g. `12/189`), and the amino-acid change (`G/D`). Using one **consistent transcript-selection rule** — the Ensembl-canonical / MANE transcript that VEP chose — we look up the full protein sequence from UniProt.

**Why 9-mers and 15-mers.** These match the two arms of the immune system: HLA **class I** presents ~9-mers to CD8+ killer T cells, and HLA **class II** presents ~15-mers to CD4+ helper T cells. We generate both so the analysis covers both pathways.

**How the peptides are cut.** For a mutation at protein position *m*, we slide a window of length *L* across the protein and keep every window that (a) contains position *m* and (b) stays inside the protein. For a 9-mer, the mutation can sit at any of the 9 positions, so we usually get 9 windows. For each window we make the **mutant** version (with the changed amino acid) and the **wild-type** version (the normal one). We record the mutation's position within each peptide using 1-based counting (`MutPos`).

**Verification (this is where scientific care shows).** For every peptide the script checks that:
- the reference amino acid the MAF claims is really the one in the reference protein;
- the mutant amino acid was placed at the correct spot;
- the window did not run off the end of the protein;
- the mutant and wild-type peptides differ at exactly one position.

If the reference amino acid does **not** match the reference protein (usually because of protein-isoform or version differences), we do **not** guess — we log that mutation in an audit file and exclude it. About 2,964 mutations were excluded this way; the rest produced **6,873,140 peptides**. As a validation, KRAS G12D reproduces the textbook example exactly: mutant `VVGADGVGK`, wild-type `VVGAGGVGK`, differing only at position 5.

### Step 6 — Predicting HLA binding (`06_predict_neoantigens.py`)

**Goal.** Predict how tightly each HLA allele would bind each peptide, compare mutant to wild-type, and assemble **Deliverable 04** (**Sections 11–15**).

**Which HLA alleles?** Because we don't know each patient's HLA type, we use a fixed panel of three common HLA class I alleles — **HLA-A\*02:01, A\*01:01, A\*03:01** — plus two class II alleles for the 15-mers. Every predicted score is always reported next to the allele it belongs to (you can never interpret a binding score without knowing the allele).

**The prediction tool.** We use **MHCflurry**, a machine-learning model trained on large experimental datasets of which peptides bind which HLA molecules. Given a peptide and an allele, it predicts three things: the **binding affinity** (IC50 in nanomolar — *lower means stronger*), a **percentile rank** (how this peptide compares to random peptides for that allele — *lower is better*), and a **presentation score** (0–1, *higher is better*).

**Binder classes.** From the affinity we label each peptide: **Strong** (< 50 nM), **Weak** (< 500 nM), or **Non-binder**. These are standard immunology thresholds.

**Comparing mutant vs wild-type (Section 14).** For each peptide window we compute
`DeltaAffinity = WT_affinity − Mutant_affinity`. Because lower nM = stronger binding, a **positive delta means the mutant binds more strongly than the normal peptide** — exactly the situation that makes a good neoantigen (the mutation makes the peptide *more* visible to the immune system). Across the whole dataset, **54% of mutant 9-mers bind more strongly than their wild-type counterpart**.

**Being honest about what we didn't compute.** Class II (15-mer) binding needs a different tool (NetMHCIIpan), and immunogenicity needs yet another (PRIME). We did not run those, so those columns are filled with **NA** — never with a made-up number or a `0`. This honesty is required by the assignment (Rule 8) and is good scientific practice.

**A note on scale and speed.** There are ~2.46 million *unique* 9-mer sequences. We score each unique peptide once per allele (instead of millions of repeats) and **cache** the results to a file, so the slow neural-network step never has to be repeated. The final table has ~16.3 million rows and every score is stamped with the tool name and version for transparency.

### Step 7 — Prioritising candidate neoantigens (`07_prioritize_candidates.py`)

**Goal.** Shrink the 16-million-row table to a short list of the most promising neoantigens (**Section 14**).

**The four criteria (and the reasoning).** A mutant 9-mer is shortlisted only if it is:
1. a **Strong or Weak binder** (an HLA molecule can display it);
2. a **better binder than its wild-type** (positive delta — the mutation increases visibility);
3. in an **expressed gene** (GeneLevelTPM > 1 — the protein is actually made);
4. **recurrent** (the specific mutation appears in ≥ 2 tumours — it matters to more than one patient).

**Why "recurrence" is per-mutation, not per-gene.** This is a subtle but important point (and a bug we caught and fixed during review). "How often does KRAS G12D occur?" must count the tumours with **exactly** that mutation (65 samples), not all KRAS mutations combined (250). Using the specific-mutation count is what Section 14 means by "occurrence in tumour samples", and it makes the shortlist far more meaningful.

**The result.** The shortlist is **1,536 candidate peptide–allele pairs across 966 genes**. Encouragingly, its most recurrent members are the canonical colorectal neoantigens: **KRAS G12V, G12C, G12S, G12A** (mostly on HLA-A\*03:01), **PIK3CA E542K**, and **SMAD4 R361H**. The pipeline rediscovering **KRAS G12x on HLA-A\*03:01** — a real, clinically pursued shared neoantigen — is strong evidence that the whole method works.

---

## Part E — A gentle intro to the computing concepts

If the code is unfamiliar, here are the ideas you actually need.

- **A script** is just a recipe: a list of instructions the computer follows top to bottom. Ours are written in **Python**, a readable programming language. Lines starting with `#` are **comments** — notes for humans that the computer ignores.

- **A library** is a toolbox someone else wrote. We use **pandas** (for working with tables — think of a programmable spreadsheet), **numpy** (fast maths on lists of numbers), **matplotlib** (drawing charts), and **MHCflurry** (the HLA-binding predictor).

- **A DataFrame** is one table in pandas: rows and named columns, like a spreadsheet sheet. "Filtering" a DataFrame means keeping only the rows that meet a condition (e.g. "only missense mutations").

- **A dictionary** (`{key: value}`) is a lookup table: give it a key (say, a gene name) and it instantly returns the value (say, that gene's expression). We use these to attach information quickly.

- **"Streaming" a file** means reading it one line at a time instead of loading it all into memory at once. Our biggest file is 2.4 gigabytes with 16 million rows — too big to hold in memory comfortably — so we process it line by line.

- **Caching** means saving a slow computation's result so you never have to redo it. We cached the MHCflurry predictions; re-running the final step then takes minutes instead of hours.

- **Why some steps "couldn't run in the sandbox".** The specialised prediction tool (MHCflurry) has heavy dependencies that only installed on your own computer, not in the assistant's restricted environment. So those steps were run locally — a normal part of bioinformatics, where different tools live in different places.

You do **not** need to be able to write this code from scratch. You need to be able to explain, for each step, *what it does and why* — which is exactly what this guide and the annotated scripts give you.

---

## Part F — Reading and interpreting the results

- **Deliverable 01** (`01_mutation_by_sample.tsv`): rows = mutations, columns = samples, `1`/`0` = present/absent. Scan a row to see how many tumours share a mutation; scan the top rows to see the commonest mutations.

- **Deliverable 02** (`02_gene_by_sample_TPM.tsv`): rows = genes, columns = samples, values = TPM. Higher = more expressed. Housekeeping genes (e.g. `ACTB`, `GAPDH`) sit at the top; that's expected.

- **Deliverable 03** (`03_integrated_mutation_expression.tsv`): the mutation grid plus a `GeneLevelTPM` column. Use it to see whether a mutated gene is also expressed.

- **Deliverable 04** (`04_neoantigen_predictions.tsv`): one row per (peptide, HLA allele). The columns to read are `BindingAffinity` (lower = stronger), `BinderClass`, `DeltaAffinity_WTminusMut` (positive = mutant binds better), `GeneLevelTPM`, and `MutationFrequency`. `ImmunogenicityScore` and the class-II binding columns are `NA` because those tools were not run.

- **Shortlist** (`neoantigen_candidates_shortlist.tsv`): the payoff — the mutations that combine strong-ish binding, better-than-wild-type binding, expression, and recurrence. Start reading here for the biological story.

A useful rule of thumb for a good candidate: **strong mutant binding + positive delta + expressed gene + recurrent mutation**. KRAS G12V ticks all four.

---

## Part G — Limitations and honest caveats

Good science states what it did *not* establish:

- **Different sample sets.** Mutation and expression data come from overlapping but not identical patients, so expression is summarised per gene (a median), not matched per patient.
- **Missense only.** We analysed single-amino-acid changes. Insertions, deletions, and frameshifts can also make strong neoantigens and are a natural extension.
- **RSEM→TPM is a rescaling.** As noted in Step 2, it is a documented per-million rescaling, not a fully length-normalised TPM; native GDC TPM would be stricter.
- **Predictions, not experiments.** Binding and recurrence are *computational predictions*. They prioritise candidates for follow-up; they do **not** prove that any peptide actually triggers an immune response. (Binding ≠ immunogenicity.)
- **Fixed HLA panel.** We used three common class I alleles, not each patient's real HLA type.
- **Class II and immunogenicity not scored.** Those columns are honestly left as NA.

---

## Part H — Explaining the project in your own words

**One-sentence version.** "We took every colorectal-cancer mutation from TCGA, kept the ones that change a protein, checked which of those genes are expressed, turned the mutations into short peptides, and used a machine-learning tool to predict which mutant peptides would be displayed by HLA molecules more strongly than the normal peptide — giving a ranked list of candidate neoantigens."

**The five presentation points (Part III E):**
1. *Dataset selection* — colorectal cancer via TCGA-COAD; GDC mutations + cBioPortal RNA-seq; GRCh38 throughout.
2. *Workflow* — seven numbered, reproducible scripts from raw mutations to a candidate shortlist (see the workflow diagram).
3. *Major findings* — 153,996 mutations, 30,260 strong mutant binders, and a 1,536-candidate shortlist that recovers the known colorectal neoantigen landscape.
4. *One technical challenge* — comparing mutant vs wild-type binding correctly and at scale (16 million rows), and computing per-mutation (not per-gene) recurrence.
5. *One biological conclusion* — a recurrence- and expression-aware pipeline independently rediscovers KRAS G12x on HLA-A\*03:01, a real shared neoantigen — evidence the method is sound.

**If asked "did you validate anything?"** — Yes: the KRAS G12D peptide matches the textbook example exactly; the top mutated genes are the known colorectal drivers; and the prioritised candidates are the neoantigens the field already studies. These independent checks are how we gained confidence without doing wet-lab experiments.
