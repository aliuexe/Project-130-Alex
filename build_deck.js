const pptxgen = require("pptxgenjs");
const p = new pptxgen();
p.layout = "LAYOUT_WIDE";                 // 13.33 x 7.5
const W = 13.33, H = 7.5;

// ---- palette (immuno-genomics: deep petrol + teal + coral + mint) ----------
const DARK = "10333F", TEAL = "1C7293", MINT = "2A9D8F",
      CORAL = "EF6461", GOLD = "E9C46A", INK = "1A2A33",
      MUTED = "5B6C74", LIGHT = "FFFFFF", CODEBG = "0E2A34",
      CODEFG = "DCEAEE", PANEL = "F1F5F7";
const FIG = "/sessions/determined-epic-cannon/mnt/Project 130/figures/";
const SERIF = "Cambria", SANS = "Calibri", MONO = "Courier New";

// ---- helpers ---------------------------------------------------------------
function bg(s, c){ s.background = { color: c }; }
function title(s, t, c){
  s.addText(t, { x:0.6, y:0.42, w:12.1, h:0.9, fontFace:SERIF, fontSize:30,
    bold:true, color: c||DARK, align:"left" });
}
function stepBadge(s, n, col){
  s.addShape(p.ShapeType.ellipse, { x:0.6, y:0.5, w:0.78, h:0.78, fill:{color:col} });
  s.addText(String(n), { x:0.6, y:0.5, w:0.78, h:0.78, align:"center", valign:"middle",
    fontFace:SERIF, fontSize:30, bold:true, color:"FFFFFF" });
}
function stepTitle(s, n, t, sub, col){
  stepBadge(s, n, col);
  s.addText(t, { x:1.55, y:0.44, w:11.1, h:0.6, fontFace:SERIF, fontSize:27, bold:true, color:DARK });
  if(sub) s.addText(sub, { x:1.57, y:1.02, w:11.1, h:0.4, fontFace:SANS, fontSize:14, italic:true, color:MUTED });
}
function codeBox(s, code, x, y, w, h, fs){
  s.addShape(p.ShapeType.roundRect, { x, y, w, h, rectRadius:0.06,
    fill:{color:CODEBG}, line:{color:TEAL, width:1} });
  s.addText(code, { x:x+0.16, y:y+0.12, w:w-0.32, h:h-0.24, fontFace:MONO,
    fontSize:fs||10.5, color:CODEFG, align:"left", valign:"top", lineSpacingMultiple:1.02 });
}
function fig(s, name, x, y, w, h){
  s.addImage({ path: FIG+name, x, y, w, h, sizing:{ type:"contain", w, h } });
}
function caption(s, t, x, y, w){
  s.addText(t, { x, y, w, h:0.3, fontFace:SANS, fontSize:10, italic:true, color:MUTED, align:"center" });
}
function bullets(s, items, x, y, w, h, fs){
  s.addText(items.map((t,i)=>({ text:t, options:{ bullet:{code:"2022"}, color:INK,
    fontSize:fs||14.5, fontFace:SANS, breakLine:true, paraSpaceAfter:8 } })),
    { x, y, w, h, valign:"top" });
}
function stat(s, x, y, w, num, lab, col, labcol){
  s.addText(num, { x, y, w, h:0.8, fontFace:SERIF, fontSize:40, bold:true, color:col, align:"center" });
  s.addText(lab, { x, y:y+0.82, w, h:0.6, fontFace:SANS, fontSize:12.5, color:labcol||INK, align:"center" });
}

// ===========================================================================
// 1. TITLE (dark)
let s = p.addSlide(); bg(s, DARK);
s.addText("Integrating Cancer Mutations, Gene Expression &", { x:0.8, y:2.1, w:11.7, h:0.7,
  fontFace:SERIF, fontSize:32, bold:true, color:LIGHT });
s.addText("Neoantigen Prediction in Colorectal Cancer", { x:0.8, y:2.75, w:11.7, h:0.7,
  fontFace:SERIF, fontSize:32, bold:true, color:GOLD });
s.addText("TCGA-COAD  |  Project 130", { x:0.8, y:3.6, w:11.7, h:0.5,
  fontFace:SANS, fontSize:18, color:MINT });
s.addShape(p.ShapeType.line, { x:0.85, y:4.25, w:3.2, h:0, line:{color:CORAL, width:2.5} });
s.addText("[Group Member 1]  •  [Group Member 2]        2026-07-20", { x:0.8, y:4.45, w:11.7, h:0.4,
  fontFace:SANS, fontSize:14, color:"BBD0D8" });
s.addText("Somatic mutations  →  expression  →  mutant peptides  →  HLA binding  →  candidate neoantigens",
  { x:0.8, y:6.55, w:11.7, h:0.4, fontFace:SANS, fontSize:13, italic:true, color:"8FB2BD" });
s.addNotes("Colorectal cancer via TCGA-COAD. Core pipeline builds mutation, expression and integrated matrices; advanced component predicts mutant peptides and their neoantigen potential.");

// 2. OVERVIEW & OBJECTIVES
s = p.addSlide(); bg(s, LIGHT); title(s, "Project overview & objectives");
const objs = [
  ["Dataset", "One cancer type (colorectal / TCGA-COAD); build an integrated somatic-mutation + expression dataset."],
  ["Core", "Mutation-by-sample matrix, gene-by-sample TPM matrix, and an integrated matrix with a cancer-level expression summary."],
  ["Advanced", "Translate coding mutations into mutant peptides; predict HLA binding and neoantigen potential."],
  ["Rigour", "Fully reproducible from numbered scripts; exact software versions and access dates; tab-delimited outputs."],
];
let yy = 1.7;
objs.forEach((o,i)=>{
  const col = [TEAL,MINT,CORAL,GOLD][i];
  s.addShape(p.ShapeType.roundRect,{x:0.6,y:yy,w:0.5,h:0.5,rectRadius:0.08,fill:{color:col}});
  s.addText(o[0], { x:1.25, y:yy-0.03, w:2.5, h:0.55, fontFace:SANS, fontSize:16, bold:true, color:DARK, valign:"middle" });
  s.addText(o[1], { x:3.7, y:yy-0.05, w:9.0, h:0.6, fontFace:SANS, fontSize:14, color:INK, valign:"middle" });
  yy += 1.05;
});
s.addNotes("State the two deliverable tiers: core (80 pts) and advanced (20 pts), and the reproducibility rules.");

// 3. DATASET SELECTION
s = p.addSlide(); bg(s, LIGHT); title(s, "Dataset selection: colorectal cancer");
bullets(s, [
  "TCGA-COAD chosen: large, well-characterised colorectal cohort with a high neoantigen burden.",
  "Somatic mutations from the GDC project-level Masked Somatic Mutation MAF (GRCh38).",
  "Core analysis filtered to protein-coding, PASS, missense single-nucleotide variants.",
  "Missense SNVs (red) are the substrate for the neoantigen analysis.",
], 0.6, 1.6, 5.9, 4.2, 15);
fig(s, "fig6_variant_classification.png", 6.7, 1.55, 6.2, 4.6);
caption(s, "Somatic variant classes in the raw TCGA-COAD MAF", 6.7, 6.15, 6.2);
s.addNotes("Why colorectal: KRAS/TP53/APC/PIK3CA drivers, MSI subset with high neoantigen load. Filtering rationale from Section 4.");

// 4. DATA SOURCES
s = p.addSlide(); bg(s, LIGHT); title(s, "Data sources & provenance");
const rows = [
  [{text:"Dataset",options:{bold:true,color:"FFFFFF",fill:{color:DARK}}},
   {text:"Source",options:{bold:true,color:"FFFFFF",fill:{color:DARK}}},
   {text:"Identifier / assembly",options:{bold:true,color:"FFFFFF",fill:{color:DARK}}},
   {text:"Access",options:{bold:true,color:"FFFFFF",fill:{color:DARK}}}],
  ["Somatic mutations (MAF)","NCI GDC — TCGA-COAD","Masked Somatic Mutation, GRCh38 (gdc-client manifest)","Open"],
  ["RNA-seq expression","cBioPortal — coadread_tcga_pan_can_atlas_2018","data_mrna_seq_v2_rsem.txt (RSEM, HiSeq RNASeqV2)","Open"],
  ["Reference proteome","UniProt","Reviewed human proteome UP000005640","Open"],
];
s.addTable(rows, { x:0.6, y:1.7, w:12.1, colW:[2.6,3.6,4.5,1.4], rowH:0.5,
  fontFace:SANS, fontSize:12.5, color:INK, valign:"middle", border:{type:"solid",color:"D5DEE2",pt:1},
  fill:{color:"FFFFFF"} });
s.addText("All downloads accessed 2026-07-20. Assembly is GRCh38/hg38 throughout — no coordinate mixing (Rule 5).",
  { x:0.6, y:4.7, w:12.1, h:0.5, fontFace:SANS, fontSize:13, italic:true, color:MUTED });
s.addText("RNA-seq (not microarray). RSEM values are converted to TPM (documented) — never mislabelled (Rule 6).",
  { x:0.6, y:5.2, w:12.1, h:0.5, fontFace:SANS, fontSize:13, italic:true, color:MUTED });
s.addNotes("Emphasise open access, exact files, download dates, single assembly.");

// 5. WORKFLOW OVERVIEW
s = p.addSlide(); bg(s, LIGHT); title(s, "Computational workflow");
fig(s, "fig13_workflow_schematic.png", 0.5, 1.4, 12.3, 5.2);
s.addText("Eight numbered, reproducible scripts (01–07) from MAF + RNA-seq to prioritised neoantigens.",
  { x:0.6, y:6.7, w:12.1, h:0.4, fontFace:SANS, fontSize:13, italic:true, color:MUTED });
s.addNotes("Walk left-to-right, top row = core, bottom row = advanced. Each box is one script.");

// 6. SECTION DIVIDER — CORE
s = p.addSlide(); bg(s, DARK);
s.addText("Part I — Core pipeline", { x:0.8, y:2.7, w:11.7, h:0.9, fontFace:SERIF, fontSize:40, bold:true, color:LIGHT });
s.addText("Mutation matrix  ·  TPM matrix  ·  integration  ·  quality control", { x:0.8, y:3.7, w:11.7, h:0.5,
  fontFace:SANS, fontSize:18, color:MINT });
s.addShape(p.ShapeType.line, { x:0.85, y:3.55, w:3.0, h:0, line:{color:CORAL, width:2.5} });

// 7. STEP 1 — MUTATION PROCESSING
s = p.addSlide(); bg(s, LIGHT);
stepTitle(s, 1, "Somatic mutation processing", "Sections 4–5  •  script 01_build_mutation_matrix.py", TEAL);
bullets(s, [
  "Filter to protein-coding, PASS, nonsynonymous missense SNVs.",
  "Collapse tumour barcodes to TCGA sample level.",
  "Binary mutation-by-sample matrix; one row per distinct mutation (HGVS).",
], 0.6, 1.7, 5.7, 2.4, 14);
codeBox(s,
`# Section 4 filters (all GRCh38)
keep = ((df.BIOTYPE == "protein_coding") &
        df.GDC_FILTER.isin(["", "PASS"]) &
        (df.Variant_Classification ==
             "Missense_Mutation") &
        (df.Variant_Type == "SNP"))
# distinct mutation = Gene + HGVSc + HGVSp
# value = 1 if present in sample else 0`,
  0.6, 4.15, 5.9, 2.7, 11);
fig(s, "fig1_top_mutated_genes.png", 6.7, 1.6, 6.2, 4.6);
caption(s, "Top mutated genes — TP53/KRAS/PIK3CA drivers recovered", 6.7, 6.2, 6.2);
s.addNotes("310,472 raw → 184,574 filtered → 153,996 distinct mutations × 586 samples.");

// 8. DELIVERABLE 01 EXAMPLE
s = p.addSlide(); bg(s, LIGHT); title(s, "Deliverable 01 — mutation-by-sample matrix");
const d1 = [
  [{text:"Gene_Name",options:{bold:true,color:"FFFFFF",fill:{color:TEAL}}},
   {text:"Mutation",options:{bold:true,color:"FFFFFF",fill:{color:TEAL}}},
   {text:"AminoAcid_Change",options:{bold:true,color:"FFFFFF",fill:{color:TEAL}}},
   {text:"TCGA-…-01",options:{bold:true,color:"FFFFFF",fill:{color:TEAL}}},
   {text:"TCGA-…-02",options:{bold:true,color:"FFFFFF",fill:{color:TEAL}}},
   {text:"TCGA-…-03",options:{bold:true,color:"FFFFFF",fill:{color:TEAL}}}],
  ["KRAS","c.35G>A","p.G12D","0","1","0"],
  ["KRAS","c.35G>T","p.G12V","0","0","1"],
  ["BRAF","c.1919T>A","p.V640E","0","0","0"],
  ["TP53","c.743G>A","p.R248Q","1","0","1"],
];
s.addTable(d1, { x:0.9, y:1.9, w:11.5, colW:[1.8,2.0,2.7,1.667,1.667,1.667], rowH:0.55,
  fontFace:MONO, fontSize:13, color:INK, align:"center", valign:"middle",
  border:{type:"solid",color:"D5DEE2",pt:1} });
s.addText("153,996 distinct mutations × 586 tumour samples. Multiple mutations in a gene are separate rows (Section 5).",
  { x:0.9, y:5.2, w:11.5, h:0.5, fontFace:SANS, fontSize:14, italic:true, color:MUTED });
s.addNotes("Values strictly 0/1; HGVS notation; tab-delimited (Rule 1).");

// 9. STEP 2 — EXPRESSION RSEM→TPM
s = p.addSlide(); bg(s, LIGHT);
stepTitle(s, 2, "Expression: RSEM → TPM", "Section 6  •  script 02_build_expression_matrix.py", MINT);
s.addShape(p.ShapeType.roundRect,{x:0.6,y:1.7,w:5.9,h:1.15,rectRadius:0.06,fill:{color:PANEL},line:{color:MINT,width:1}});
s.addText([{text:"Why convert?  ",options:{bold:true,color:DARK}},
  {text:"RSEM sample sums ≈ 1.8×10⁷, so values are NOT TPM. Rescaling each sample to 1×10⁶ yields true TPM (Rule 6).",options:{color:INK}}],
  { x:0.75, y:1.8, w:5.6, h:0.95, fontFace:SANS, fontSize:12.5, valign:"middle" });
codeBox(s,
`# RSEM -> TPM: rescale each sample to 1e6
col_totals = expr.sum(axis=0, skipna=True)
tpm = expr.divide(col_totals, axis=1) * 1e6
assert np.allclose(tpm.sum(axis=0), 1e6)
# duplicated gene symbols summed -> one row`,
  0.6, 3.05, 5.9, 2.0, 11);
s.addText("20,511 genes × 592 tumour samples  ·  5.64% missing → NA (never 0)",
  { x:0.6, y:5.2, w:5.9, h:0.5, fontFace:SANS, fontSize:12.5, italic:true, color:MUTED });
fig(s, "fig2_genelevel_tpm_distribution.png", 6.7, 1.6, 6.2, 4.6);
caption(s, "Distribution of GeneLevelTPM (log scale)", 6.7, 6.2, 6.2);
s.addNotes("Emphasise the honest RSEM→TPM conversion. Every sample column sums to exactly 1e6.");

// 10. STEP 3 — INTEGRATION
s = p.addSlide(); bg(s, LIGHT);
stepTitle(s, 3, "Integration of mutations & expression", "Section 7  •  script 03_integrate_datasets.py", TEAL);
bullets(s, [
  "Mutation and expression cohorts differ → summarise expression per gene.",
  "GeneLevelTPM = median TPM across tumour samples (robust to outliers).",
  "Merge by gene; genes absent from expression → NA (never 0).",
], 0.6, 1.7, 5.9, 2.3, 14);
codeBox(s,
`# GeneLevelTPM = median TPM across samples
gene_level = (tpm.set_index("GeneName")
                 .median(axis=1, skipna=True))
# merge into mutation matrix by gene symbol`,
  0.6, 4.05, 5.9, 1.5, 11);
fig(s, "fig4_freq_vs_expression.png", 6.7, 1.6, 6.2, 4.7);
caption(s, "Mutation frequency vs gene expression per gene", 6.7, 6.25, 6.2);
s.addNotes("Deliverable 03 keeps per-sample mutation columns + GeneLevelTPM.");

// 11. QC + FIGURES
s = p.addSlide(); bg(s, LIGHT); title(s, "Quality control (Section 8)");
stat(s, 0.6, 1.7, 2.9, "310,472", "raw MAF records", TEAL);
stat(s, 3.5, 1.7, 2.9, "153,996", "distinct mutations", MINT);
stat(s, 0.6, 3.55, 2.9, "586 / 592", "mutation / expr samples", CORAL);
stat(s, 3.5, 3.55, 2.9, "5.64%", "missing expression", GOLD);
s.addText("Consistent Hugo symbols · GRCh38 throughout · counts before/after filtering reported.",
  { x:0.6, y:5.35, w:5.9, h:0.7, fontFace:SANS, fontSize:12.5, italic:true, color:MUTED });
fig(s, "fig3_mutation_heatmap.png", 6.7, 1.55, 6.2, 4.7);
caption(s, "Recurrent mutations across tumour samples", 6.7, 6.2, 6.2);
s.addNotes("Four required QC figures produced; two shown here and on prior slides.");

// 12. SECTION DIVIDER — ADVANCED
s = p.addSlide(); bg(s, DARK);
s.addText("Part II — Neoantigens", { x:0.8, y:2.7, w:11.7, h:0.9, fontFace:SERIF, fontSize:40, bold:true, color:LIGHT });
s.addText("mutant peptides  ·  HLA binding  ·  immunogenicity  ·  prioritisation", { x:0.8, y:3.7, w:11.7, h:0.5,
  fontFace:SANS, fontSize:18, color:MINT });
s.addShape(p.ShapeType.line, { x:0.85, y:3.55, w:3.0, h:0, line:{color:CORAL, width:2.5} });

// 13. STEP 4 — VARIANT → PEPTIDE
s = p.addSlide(); bg(s, LIGHT);
stepTitle(s, 4, "Variant → protein → peptides", "Sections 9–10  •  script 05_annotate_and_generate_peptides.py", MINT);
bullets(s, [
  "Map each missense mutation to its VEP-canonical transcript & protein.",
  "Generate all mutation-containing 9-mers (MHC-I) and 15-mers (MHC-II) + wild-type.",
  "Verify reference AA vs reference protein; 1-based mutation position.",
], 0.6, 1.7, 6.2, 2.3, 13.5);
codeBox(s,
`def windows(seq, m, L):        # m = 1-based pos
  for start in range(max(1, m-L+1),
                     min(len(seq)-L+1, m)+1):
    yield seq[start-1:start-1+L], m-start+1
# KRAS p.G12D  ->  mutant  VVGADGVGK  (MutPos 5)
#              ->  wild-t  VVGAGGVGK`,
  0.6, 4.05, 6.2, 2.1, 10.5);
s.addShape(p.ShapeType.roundRect,{x:7.1,y:1.7,w:5.6,h:4.5,rectRadius:0.06,fill:{color:PANEL},line:{color:MINT,width:1}});
s.addText("Worked example — KRAS G12D", { x:7.3, y:1.85, w:5.2, h:0.5, fontFace:SERIF, fontSize:18, bold:true, color:DARK });
s.addText([
  {text:"Reference protein: …LVVV",options:{color:INK}},{text:"G",options:{bold:true,color:TEAL}},{text:"AGGVGKS…",options:{color:INK,breakLine:true}},
  {text:"Mutant   protein: …LVVV",options:{color:INK}},{text:"D",options:{bold:true,color:CORAL}},{text:"AGGVGKS…",options:{color:INK,breakLine:true}},
  {text:"\n9-mer mutant:  ",options:{color:INK}},{text:"VVGA",options:{color:INK}},{text:"D",options:{bold:true,color:CORAL}},{text:"GVGK",options:{color:INK,breakLine:true}},
  {text:"9-mer wild-t:  ",options:{color:INK}},{text:"VVGA",options:{color:INK}},{text:"G",options:{bold:true,color:TEAL}},{text:"GVGK",options:{color:INK,breakLine:true}},
  {text:"\nMatches the assignment's worked example exactly.",options:{italic:true,color:MUTED}},
], { x:7.3, y:2.45, w:5.2, h:3.6, fontFace:MONO, fontSize:12.5, valign:"top" });
s.addNotes("145,612 mutations annotated; 6,873,140 peptides; 2,964 excluded for ref-AA mismatch (audited).");

// 14. PEPTIDE FUNNEL
s = p.addSlide(); bg(s, LIGHT); title(s, "From mutations to candidates — scale");
fig(s, "fig5_methods_funnel.png", 2.6, 1.4, 8.1, 5.6);
s.addNotes("Funnel shows the peptide explosion (6.87M) narrowing to 1,536 prioritised candidates.");

// 15. STEP 5 — HLA + MHC BINDING
s = p.addSlide(); bg(s, LIGHT);
stepTitle(s, 5, "HLA selection & MHC binding", "Sections 11–12  •  script 06_predict_neoantigens.py", TEAL);
bullets(s, [
  "Option A fixed HLA-I panel: A*02:01, A*01:01, A*03:01.",
  "MHCflurry 2.2.1 presentation predictor: affinity, %rank, presentation score.",
  "Binder class: Strong <50 nM, Weak <500 nM. Class-II 15-mers reserved for NetMHCIIpan.",
], 0.6, 1.7, 6.2, 2.4, 13.5);
codeBox(s,
`r = predictor.predict(
      peptides=uniq, alleles=[allele],
      include_affinity_percentile=True)
# every score carries its HLA allele (Rule 7)
# missing tool outputs -> NA, never 0`,
  0.6, 4.15, 6.2, 1.9, 10.5);
fig(s, "fig7_mut_vs_wt_affinity.png", 7.0, 1.6, 6.0, 4.6);
caption(s, "Mutant vs wild-type 9-mer binding, with 50/500 nM cutoffs", 7.0, 6.2, 6.0);
s.addNotes("2.46M unique 9-mers × 3 alleles. Present binding vs presentation vs immunogenicity distinction next.");

// 16. WT vs MUTANT COMPARISON
s = p.addSlide(); bg(s, LIGHT); title(s, "Wild-type vs mutant comparison (Section 14)");
fig(s, "fig8_delta_affinity.png", 0.5, 1.5, 6.3, 4.7);
fig(s, "fig11_mut_vs_wt_scatter.png", 7.0, 1.4, 5.9, 5.2);
s.addText("ΔAffinity = WT − Mutant.  54% of mutant 9-mers bind more strongly than wild-type (below the diagonal, right).",
  { x:0.5, y:6.35, w:6.3, h:0.7, fontFace:SANS, fontSize:12, italic:true, color:MUTED });
s.addNotes("Positive delta / below-diagonal = neoantigen-favourable. This is the core selection signal.");

// 17. BINDERS BY ALLELE & CLASS
s = p.addSlide(); bg(s, LIGHT); title(s, "Predicted binders by class & allele");
fig(s, "fig9_binder_class_mut_vs_wt.png", 0.5, 1.5, 6.2, 4.7);
fig(s, "fig10_strong_binders_by_allele.png", 6.9, 1.5, 6.0, 4.7);
s.addText("30,260 strong-binding mutant 9-mer–allele pairs; HLA-A*02:01 yields the most, then A*03:01, then A*01:01.",
  { x:0.5, y:6.3, w:12.4, h:0.5, fontFace:SANS, fontSize:12.5, italic:true, color:MUTED });
s.addNotes("Class I and class II reported separately (different pathways).");

// 18. IMMUNOGENICITY DISTINCTION
s = p.addSlide(); bg(s, LIGHT); title(s, "Binding ≠ immunogenicity (Section 13)");
s.addShape(p.ShapeType.roundRect,{x:0.6,y:1.8,w:5.8,h:3.6,rectRadius:0.08,fill:{color:PANEL},line:{color:TEAL,width:1.5}});
s.addText("MHC binding / presentation", { x:0.9, y:2.05, w:5.2, h:0.5, fontFace:SERIF, fontSize:19, bold:true, color:TEAL });
bullets(s, [
  "Peptide is predicted to be displayed by an HLA molecule.",
  "Predicted here with MHCflurry.",
  "A strong binder is NOT automatically immunogenic.",
], 0.85, 2.6, 5.3, 2.6, 13.5);
s.addShape(p.ShapeType.roundRect,{x:6.9,y:1.8,w:5.8,h:3.6,rectRadius:0.08,fill:{color:PANEL},line:{color:CORAL,width:1.5}});
s.addText("Immunogenicity", { x:7.2, y:2.05, w:5.2, h:0.5, fontFace:SERIF, fontSize:19, bold:true, color:CORAL });
bullets(s, [
  "Presented peptide can stimulate a T-cell response.",
  "Requires a dedicated predictor (e.g. PRIME).",
  "Reported as NA unless that tool is run — never fabricated.",
], 7.15, 2.6, 5.3, 2.6, 13.5);
s.addText("Binding and immunogenicity are kept in separate columns of deliverable 04.",
  { x:0.6, y:5.7, w:12.1, h:0.5, fontFace:SANS, fontSize:13, italic:true, color:MUTED });
s.addNotes("Key conceptual point examiners look for. Predictions are prioritisation features, not experimental validation (Rule 8).");

// 19. TOP CANDIDATES (biological conclusion)
s = p.addSlide(); bg(s, LIGHT); title(s, "Prioritised neoantigen candidates (Section 14)");
fig(s, "fig12_top_candidates.png", 3.0, 1.35, 7.3, 5.4);
s.addNotes("The data-driven shortlist rediscovers the known colorectal neoantigen landscape — a strong sanity check.");

// 20. MAJOR FINDINGS
s = p.addSlide(); bg(s, DARK); title(s, "Major findings", LIGHT);
stat(s, 0.7, 1.9, 3.9, "1,536", "prioritised candidates (966 genes)", GOLD, "CFE0E6");
stat(s, 4.7, 1.9, 3.9, "30,260", "strong mutant binders", MINT, "CFE0E6");
stat(s, 8.7, 1.9, 3.9, "54%", "mutant > wild-type binding", CORAL, "CFE0E6");
s.addText([
  {text:"• ",options:{color:GOLD}},{text:"Drivers recovered at expected frequency (TP53, KRAS, PIK3CA).",options:{color:LIGHT,breakLine:true}},
  {text:"• ",options:{color:GOLD}},{text:"Top candidates are canonical CRC neoantigens: KRAS G12V/C/S/A, PIK3CA E542K, SMAD4 R361H.",options:{color:LIGHT,breakLine:true}},
  {text:"• ",options:{color:GOLD}},{text:"KRAS G12x recovered on HLA-A*03:01 — a documented, clinically pursued neoantigen.",options:{color:LIGHT,breakLine:true}},
], { x:0.9, y:4.2, w:11.6, h:2.4, fontFace:SANS, fontSize:16, paraSpaceAfter:10, valign:"top" });
s.addNotes("Tie the numbers to biology.");

// 21. CHALLENGE + CONCLUSION
s = p.addSlide(); bg(s, LIGHT); title(s, "One challenge · one conclusion");
s.addShape(p.ShapeType.roundRect,{x:0.6,y:1.8,w:5.9,h:4.3,rectRadius:0.08,fill:{color:PANEL},line:{color:TEAL,width:1.5}});
s.addText("Technical challenge", { x:0.9, y:2.0, w:5.3, h:0.5, fontFace:SERIF, fontSize:19, bold:true, color:TEAL });
bullets(s, [
  "Scale: pairing WT vs mutant across ~16M peptide–allele rows.",
  "Correct per-mutation (not per-gene) recurrence.",
  "Tooling: running MHCflurry over 2.46M unique 9-mers reproducibly.",
], 0.85, 2.55, 5.4, 3.2, 13.5);
s.addShape(p.ShapeType.roundRect,{x:6.9,y:1.8,w:5.8,h:4.3,rectRadius:0.08,fill:{color:PANEL},line:{color:CORAL,width:1.5}});
s.addText("Biological conclusion", { x:7.2, y:2.0, w:5.2, h:0.5, fontFace:SERIF, fontSize:19, bold:true, color:CORAL });
bullets(s, [
  "A recurrence- and expression-aware pipeline independently rediscovers the colorectal neoantigen landscape.",
  "KRAS G12x on HLA-A*03:01 emerges as a top shared-neoantigen candidate.",
], 7.15, 2.55, 5.3, 3.2, 13.5);
s.addNotes("These map to the two required presentation points.");

// 22. LIMITATIONS & REPRODUCIBILITY
s = p.addSlide(); bg(s, LIGHT); title(s, "Limitations & reproducibility");
bullets(s, [
  "Mutation and expression cohorts are both TCGA-COAD but not identical samples (cancer-level median used).",
  "Analysis restricted to missense SNVs; in-frame indels / frameshifts are an extension.",
  "RSEM→TPM is a documented per-million rescaling; native GDC TPM is an alternative.",
  "Class-II binding and immunogenicity are NA unless NetMHCIIpan / PRIME are run.",
  "Predictions are prioritisation features, not experimental validation (Rule 8).",
], 0.6, 1.7, 7.3, 3.6, 13.5);
s.addShape(p.ShapeType.roundRect,{x:8.2,y:1.7,w:4.5,h:4.5,rectRadius:0.08,fill:{color:DARK}});
s.addText("Reproducible by design", { x:8.45, y:1.95, w:4.0, h:0.5, fontFace:SERIF, fontSize:17, bold:true, color:GOLD });
s.addText([
  {text:"7 numbered scripts (01–07)",options:{color:LIGHT,breakLine:true}},
  {text:"4 tab-delimited deliverables",options:{color:LIGHT,breakLine:true}},
  {text:"README + software versions",options:{color:LIGHT,breakLine:true}},
  {text:"Access dates + GRCh38 only",options:{color:LIGHT,breakLine:true}},
  {text:"Missing values = NA, never 0",options:{color:LIGHT,breakLine:true}},
  {text:"No manual table edits",options:{color:LIGHT,breakLine:true}},
], { x:8.45, y:2.6, w:4.0, h:3.4, fontFace:SANS, fontSize:13.5, paraSpaceAfter:9, valign:"top" });
s.addNotes("Cover Part V rules compliance.");

// 23. CLOSING
s = p.addSlide(); bg(s, DARK);
s.addText("Thank you", { x:0.8, y:2.6, w:11.7, h:1.0, fontFace:SERIF, fontSize:44, bold:true, color:LIGHT });
s.addText("Questions & discussion", { x:0.8, y:3.7, w:11.7, h:0.6, fontFace:SANS, fontSize:20, color:MINT });
s.addShape(p.ShapeType.line, { x:0.85, y:3.6, w:3.0, h:0, line:{color:CORAL, width:2.5} });
s.addText("Project 130  ·  TCGA-COAD neoantigen pipeline  ·  mutations → expression → peptides → HLA → candidates",
  { x:0.8, y:6.5, w:11.7, h:0.4, fontFace:SANS, fontSize:13, italic:true, color:"8FB2BD" });

p.writeFile({ fileName: process.argv[2] || "Project130_Presentation.pptx" }).then(f=>console.log("wrote", f));
