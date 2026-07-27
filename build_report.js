const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
  Table, TableRow, TableCell, WidthType, ShadingType, BorderStyle,
  ImageRun, PageBreak,
} = require("docx");

const FIGDIR = process.argv[3] || "figures";
function fig(file, caption) {
  const path = `${FIGDIR}/${file}`;
  const kids = [];
  if (fs.existsSync(path)) {
    kids.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 120, after: 40 },
      children: [new ImageRun({ type: "png", data: fs.readFileSync(path), transformation: { width: 560, height: 350 } })] }));
  }
  kids.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 160 },
    children: [new TextRun({ text: caption, italics: true, size: 18, color: "666666" })] }));
  return kids;
}

const ACCENT = "1F4E79";
const H = (t, lvl) => new Paragraph({ heading: lvl, spacing: { before: 220, after: 110 }, children: [new TextRun({ text: t, color: ACCENT })] });
const P = (runs) => new Paragraph({ spacing: { after: 130 }, alignment: AlignmentType.JUSTIFIED, children: Array.isArray(runs) ? runs : [new TextRun(runs)] });
const T = (t) => new TextRun(t);
const B = (t) => new TextRun({ text: t, bold: true });

function tbl(headers, rows, widths) {
  const total = widths.reduce((a, b) => a + b, 0);
  const hdr = new TableRow({
    tableHeader: true,
    children: headers.map((h, i) => new TableCell({
      width: { size: widths[i], type: WidthType.DXA },
      shading: { type: ShadingType.CLEAR, fill: ACCENT },
      children: [new Paragraph({ children: [new TextRun({ text: h, bold: true, color: "FFFFFF", size: 18 })] })],
    })),
  });
  const body = rows.map((r) => new TableRow({
    children: r.map((c, i) => new TableCell({
      width: { size: widths[i], type: WidthType.DXA },
      children: [new Paragraph({ children: [new TextRun({ text: String(c), size: 18 })] })],
    })),
  }));
  return new Table({ columnWidths: widths, width: { size: total, type: WidthType.DXA }, rows: [hdr, ...body] });
}

const doc = new Document({
  styles: { default: { document: { run: { font: "Calibri", size: 22 } } } },
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } },
    children: [
      new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 60 }, children: [new TextRun({ text: "Integrating Cancer Mutations, Gene Expression, and Neoantigen Prediction in Colorectal Cancer (TCGA-COAD)", bold: true, size: 30, color: ACCENT })] }),
      new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 40 }, children: [new TextRun({ text: "Project 130 — Brief Report", size: 22, color: "666666" })] }),
      new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 220 }, children: [new TextRun({ text: "Group members: [Group Member 1], [Group Member 2]   |   Date: 2026-07-20", size: 20, color: "666666" })] }),

      // 1. Background
      H("1. Background and selected cancer type", HeadingLevel.HEADING_1),
      P([T("Colorectal cancer is among the most common malignancies worldwide and one of the best-characterised solid tumours at the genomic level. Its mutational landscape is dominated by a small set of recurrent driver events — inactivation of "), B("APC"), T(" and "), B("TP53"), T(", and activating mutations in "), B("KRAS"), T(", "), B("PIK3CA"), T(", and "), B("BRAF"), T(" — superimposed on a long tail of passenger mutations. A subset of tumours is hypermutated (microsatellite-unstable or "), B("POLE"), T("-mutated), giving colorectal cancer a comparatively high neoantigen burden and making it a clinically relevant model for immunotherapy and neoantigen discovery. For this project we selected colorectal cancer analysed through "), B("The Cancer Genome Atlas colon adenocarcinoma cohort (TCGA-COAD)"), T(", integrating somatic mutation calls with matched TCGA colorectal RNA-seq expression, and extending the analysis to predicted mutant peptides and their neoantigen potential.")]),

      // 2. Data sources
      H("2. Data sources", HeadingLevel.HEADING_1),
      P([B("Somatic mutations. "), T("Project-level somatic mutation data were obtained from the NCI Genomic Data Commons (GDC) for TCGA-COAD as an open-access Masked Somatic Mutation file in Mutation Annotation Format (MAF), downloaded via the GDC Data Transfer Tool using a project-filtered manifest (accessed 2026-07-20; local file cohortMAF.2026-07-15.maf.gz). The MAF is aligned and annotated on the "), B("GRCh38/hg38"), T(" reference assembly and is VEP-annotated, providing transcript-, protein-, and consequence-level information used later for peptide generation.")]),
      P([B("Gene expression. "), T("RNA-seq expression for the same cancer type was obtained from the TCGA Colorectal PanCancer Atlas cohort via cBioPortal (study coadread_tcga_pan_can_atlas_2018, file data_mrna_seq_v2_rsem.txt; RSEM batch-normalised values from Illumina HiSeq RNASeqV2; accessed 2026-07-20). RNA-seq, not microarray, was used. Because the mutation and expression cohorts are both TCGA-COAD but not identical samples, expression was summarised at the cancer level (see Methods).")]),
      P([B("Reference proteome. "), T("Protein sequences for peptide extraction and reference-amino-acid verification were taken from the UniProt reviewed human proteome UP000005640 (accessed 2026-07-20).")]),

      // 3. Methods
      H("3. Methods", HeadingLevel.HEADING_1),
      P([B("Mutation processing. "), T("From 310,472 MAF records we retained protein-coding, high-confidence (PASS) nonsynonymous missense single-nucleotide variants, yielding 184,574 filtered records. Tumour barcodes were collapsed to the TCGA sample level, and a binary mutation-by-sample matrix was constructed in which each row is a distinct mutation (unique gene + HGVS coding change + protein change) and each column a tumour sample (1 = present, 0 = absent). Mutations are expressed in HGVS notation (e.g. c.35G>A, p.G12D), with multiple mutations in the same gene represented as separate rows.")]),
      P([B("Expression processing and RSEM→TPM conversion. "), T("The source RSEM values are normalised abundance estimates whose per-sample column sums are approximately 1.8×10⁷ — they are therefore not already on the TPM (per-million) scale. To comply with the requirement that non-TPM values not be mislabelled as TPM, each sample column was rescaled to sum to exactly 1,000,000 (TPM = RSEM / Σ_sample RSEM × 10⁶). Duplicated gene symbols were collapsed by summation to one row per unique symbol. This produced a gene-by-sample TPM matrix of 20,511 genes across 592 tumour samples.")]),
      P([B("Integration. "), T("For each gene we computed GeneLevelTPM as the median TPM across tumour samples (robust to extreme values) and merged this value into the mutation matrix by gene symbol. Mutations in genes absent from the expression matrix received GeneLevelTPM = NA (never zero). The integrated matrix retains the per-sample mutation columns alongside GeneLevelTPM.")]),
      P([B("Advanced component. "), T("Each eligible missense mutation was mapped to its VEP-selected transcript and protein (Ensembl canonical / MANE Select, as recorded in the GDC MAF). Using the UniProt reference protein, we generated all mutation-containing 9-mer and 15-mer peptides plus their wild-type counterparts, verifying for every peptide that the reference amino acid matched the reference protein, that the mutant residue was placed correctly, that windows stayed within the protein, and that wild-type and mutant peptides differed only at the mutated position. Class I 9-mers were evaluated with MHCflurry against a fixed HLA panel (HLA-A*02:01, HLA-A*01:01, HLA-A*03:01); class II 15-mers were reserved for NetMHCIIpan (HLA-DRB1*15:01, HLA-DRB1*07:01). For each mutant peptide we compared predicted affinity with its wild-type counterpart (DeltaAffinity = WT − Mut; positive indicates stronger mutant binding).")]),

      // 4. QC results
      H("4. Quality-control results", HeadingLevel.HEADING_1),
      P([T("Filtering reduced 310,472 raw records to 184,574 high-confidence missense SNVs, forming 153,996 distinct mutations across 17,585 genes and 586 tumour samples. The expression matrix comprised 20,511 genes over 592 tumour samples, with 5.64% missing expression values and 7 duplicated gene symbols collapsed. Gene identifiers were consistent (both matrices keyed on Hugo symbols), and all mutation coordinates were confirmed to be on GRCh38 with no assembly mixing. GeneLevelTPM values spanned several orders of magnitude (median ≈ 9.4 TPM), consistent with typical RNA-seq dynamic range.")]),
      P([B("Ten most frequently mutated genes and ten most highly expressed mutated genes:")]),
      tbl(
        ["Most frequently mutated (occurrences)", "Most highly expressed mutated (GeneLevelTPM)"],
        [
          ["TTN (693)", "CEACAM5 (5778.9)"],
          ["MUC16 (335)", "ACTB (5584.4)"],
          ["SYNE1 (279)", "EEF1A1 (4675.8)"],
          ["TP53 (258)", "ACTG1 (4510.8)"],
          ["KRAS (250)", "GAPDH (3745.8)"],
          ["FAT4 (202)", "EEF2 (3413.5)"],
          ["RYR2 (173)", "FTL (2759.4)"],
          ["PIK3CA (171)", "KRT8 (2756.2)"],
          ["NEB (150)", "RPS6 (2510.3)"],
          ["OBSCN (146)", "RPL8 (2499.0)"],
        ],
        [4680, 4680]
      ),
      P([T("Four figures accompany this report: a bar plot of the most frequently mutated genes (Fig. 1), a histogram of GeneLevelTPM (Fig. 2), a presence/absence heat map of the most recurrent mutations across samples (Fig. 3), and a scatter plot comparing mutation frequency with gene expression (Fig. 4).")]),
      ...fig("fig1_top_mutated_genes.png", "Figure 1. Ten most frequently mutated genes (TCGA-COAD)."),
      ...fig("fig2_genelevel_tpm_distribution.png", "Figure 2. Distribution of GeneLevelTPM across mutations (log scale)."),
      ...fig("fig3_mutation_heatmap.png", "Figure 3. Presence/absence heat map of the 30 most recurrent mutations."),
      ...fig("fig4_freq_vs_expression.png", "Figure 4. Mutation frequency versus gene expression."),

      // 5. Mutation and expression results
      H("5. Mutation and expression results", HeadingLevel.HEADING_1),
      P([T("The mutation results recapitulate the known biology of colorectal cancer. The apparent frequency leaders TTN, MUC16, SYNE1, RYR2, NEB and OBSCN are very large genes whose high mutation counts reflect target size rather than driver status, and they serve as an internal sanity check. Embedded among them are the canonical colorectal drivers TP53 (258 occurrences), KRAS (250) and PIK3CA (171); BRAF and FAT4 also rank highly. The recovery of these drivers at expected relative frequencies indicates the filtering and matrix construction behaved correctly.")]),
      P([T("The expression results are equally consistent with expectation. The most highly expressed mutated genes are dominated by housekeeping and structural transcripts (ACTB, EEF1A1, GAPDH, ribosomal proteins) together with CEACAM5, which encodes carcinoembryonic antigen (CEA) — the classical colorectal tumour marker — as the single most highly expressed mutated gene. Integrating the two layers shows that clinically important drivers such as KRAS are both recurrently mutated and robustly expressed (GeneLevelTPM ≈ 67), a combination that is precisely what makes a mutation a plausible neoantigen source.")]),

      // 6. Neoantigen results
      H("6. Neoantigen results", HeadingLevel.HEADING_1),
      P([T("Peptide generation was applied to all eligible missense mutations. Of 148,576 distinct protein variants, 145,612 passed reference-amino-acid verification against the UniProt proteome; 2,964 were excluded because the reference residue in the MAF disagreed with the reference protein (isoform or sequence-version differences) and were logged rather than silently altered. This produced 6,873,140 peptide rows — 1,296,171 mutant and 1,296,171 wild-type 9-mers, and 2,140,399 mutant and 2,140,399 wild-type 15-mers — each annotated with its 1-based mutation position. As an exact validation, the KRAS G12D hotspot yielded the expected mutant 9-mer VVGADGVGK and wild-type VVGAGGVGK differing only at position 5, matching the assignment's worked example on the MANE transcript ENST00000256078.")]),
      P([T("HLA class I binding for the 9-mers was predicted with MHCflurry (v2.2.1) across the three-allele panel. The 2,460,296 unique standard-amino-acid 9-mers were each scored against the three alleles (52 peptides containing selenocysteine were reported as NA); each mutant peptide was paired with its wild-type counterpart to compute the change in predicted affinity (DeltaAffinity = WT − Mut). Because lower predicted IC50 denotes stronger binding, a positive delta identifies mutations that create or strengthen an HLA-binding peptide relative to the wild-type sequence. The complete neoantigen table (deliverable 04, 16,338,622 peptide–allele rows) reports, for every pair, the affinity, percentile rank, presentation score, binder class, wild-type comparison, gene-level expression, and mutation frequency. Across the cohort, 57,545 peptide–allele pairs were predicted strong binders (<50 nM), of which 30,260 involved mutant peptides.")]),
      P([T("Applying the Section 14 prioritisation criteria — mutant class-I peptides that are strong or weak binders, bind more strongly than their wild-type counterpart (positive delta), arise in an expressed gene (GeneLevelTPM > 1), and recur in at least two tumours (per-mutation recurrence, i.e. the number of samples carrying that specific variant) — yields a focused shortlist of 1,536 candidate peptide–allele pairs spanning 966 genes (results/neoantigen_candidates_shortlist.tsv). Reassuringly, the most recurrent candidates are the canonical colorectal driver neoantigens: KRAS p.G12V (recurrent in 56 tumours; VVGAVGVGK on HLA-A*03:01), KRAS p.G12C (17 tumours), KRAS p.G12A (11), PIK3CA p.E542K (13; on HLA-A*03:01), and SMAD4 p.R361H (11; on HLA-A*01:01) — each combining recurrent mutation, gene expression, and stronger mutant-than-wild-type binding. The strongest predicted binders overall (e.g. NR1D2 p.S406L, ~9.8 nM on HLA-A*02:01) illustrate high-affinity but low-recurrence candidates. That the data-driven shortlist independently recovers the KRAS G12x hotspots on HLA-A*03:01 — an allele documented to present KRAS neoantigens — is a strong biological sanity check on the pipeline.")]),

      // 7. Biological interpretation
      H("7. Biological interpretation", HeadingLevel.HEADING_1),
      P([T("Priority neoantigen candidates are those satisfying several features simultaneously: strong predicted mutant peptide–HLA binding, better binding than the wild-type peptide, favourable immunogenicity, expression of the mutated gene, and occurrence in multiple tumours. Colorectal driver hotspots such as KRAS G12D are attractive because they are both recurrent across patients and expressed, meaning a single predicted neoantigen could in principle be relevant to many tumours sharing the same mutation. It is essential, however, to distinguish MHC binding or presentation — the prediction that a peptide can be displayed by an HLA molecule — from immunogenicity, the prediction that a displayed peptide can actually elicit a T-cell response. A strong binder is not automatically immunogenic, which is why binding and immunogenicity are reported in separate columns and immunogenicity is left as NA unless a dedicated predictor (PRIME) is run. These predictions are prioritisation features, not proof of experimental immunogenicity.")]),

      // 8. Limitations
      H("8. Limitations", HeadingLevel.HEADING_1),
      P([T("Several limitations apply. First, the mutation and expression cohorts are both TCGA-COAD but not identical samples, so expression is summarised as a cancer-level median rather than matched per sample; this discards sample-specific expression variation. Second, the analysis is restricted to missense SNVs — in-frame insertions/deletions and frameshifts, which can generate highly immunogenic neoantigens, are not included. Third, all binding and immunogenicity values are computational predictions using a fixed HLA panel rather than patient-specific typing, and must not be interpreted as experimental validation. Fourth, class II (15-mer) binding and class-I immunogenicity are reported as NA unless NetMHCIIpan and PRIME are installed and the advanced step re-run. Finally, 2,964 mutations were excluded from peptide generation owing to reference-sequence discordance, a consequence of isoform/version differences between the MAF annotation and the UniProt proteome.")]),

      // 9. Contributions
      H("9. Contributions of each group member", HeadingLevel.HEADING_1),
      P([B("[Group Member 1]: "), T("[e.g. dataset acquisition and mutation processing (scripts 01, 05); reference-genome and QC documentation].")]),
      P([B("[Group Member 2]: "), T("[e.g. expression processing and integration (scripts 02–04); MHC binding/immunogenicity and neoantigen table (script 06); report and README].")]),
      P([T("Both members understand and can explain the complete workflow, from data download through mutation and expression matrix construction, integration, and neoantigen prediction.")]),
    ],
  }],
});

Packer.toBuffer(doc).then((buf) => {
  fs.writeFileSync(process.argv[2] || "Project130_Report.docx", buf);
  console.log("wrote report");
});
