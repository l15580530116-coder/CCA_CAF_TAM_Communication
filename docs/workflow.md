# Analysis Workflow

Scripts are numbered by execution order.

## Step 0 — Environment Check
`00_check_environment.R` — Verifies R version, Bioconductor, and required packages.

## Step 1 — Download TCGA-CHOL
`01_download_TCGA_CHOL.R` — Downloads TCGA-CHOL RNA-seq data from GDC using TCGAbiolinks.

## Step 2 — Clinical Data
`02_clean_TCGA_CHOL_clinical.R` — Cleans and formats clinical/survival data.

## Step 3 — Expression Preprocessing
`03_preprocess_TCGA_CHOL_expression.R` — Converts raw counts to gene symbols and log2(TPM+1).

## Step 4 — GEO Bulk Preprocessing
`04_download_preprocess_GSE107943.R` — Downloads GSE107943 from GEO; extracts raw counts; normalizes to log2(CPM+1).

## Step 5 — Gene Set Curation
`05_prepare_gene_sets.R` — Curates 842 IM/CAF/TAM-related genes from KEGG, HALLMARK, and Reactome pathways.

## Step 6 — DEG and Cross-Cohort Validation
`06_DEG_TCGA_GSE107943.R` — DESeq2 (TCGA) + edgeR paired (GSE107943); cross-cohort overlap; directional concordance.
`06_fix_gse_colnames.R` — Fixes column name formatting in large GSE CSV files.
`06_fig2d_pca.R`, `06_generate_remaining_figures.R` — DEG-related figure generation.

## Step 7 — Functional Enrichment
`07_functional_enrichment_validated_DEGs.R` — GO BP and KEGG enrichment of validated DEGs using clusterProfiler.

## Step 8 — ssGSEA and Survival
`08_GSVA_ssGSEA_scores_survival.R` — ssGSEA pathway scoring; aggressive microenvironment score construction; exploratory survival analysis.

## Step 9 — Immune/Checkpoint Analysis
`10_immune_infiltration_checkpoint_analysis.R` — Marker-based immune cell scoring (13 cell types); checkpoint correlation.

## Step 10 — Single-Cell Analysis
`11_scRNA_GSE138709_preprocess_score_mapping.R` — GSE138709 preprocessing; Seurat v5 pipeline.
`11b_GSE138709_gene_identifier_diagnosis.R` — Gene identifier diagnostics.
`11c_fix_cell_annotation_scoring.R` — Cell type annotation and module scoring.

## Step 11 — CellChat
`12_CellChat_GSE138709_CAF_TAM_Epithelial.R` — CellChat v2.2.0.9001 analysis of CAF–TAM–Epithelial communication.

## Step 12 — Hub Gene Prioritization
`13_integrated_hub_gene_axis_prioritization.R` — Multi-layer integrated evidence scoring (7 evidence layers across 37 candidate genes).

## Step 13 — External Validation
`14_validate_hub_genes_GSE26566.R` — GSE26566 microarray expression validation.

## Step 14 — Druggability and Docking
`16_clinical_relevance_hub_genes_axes.R` — Axis score construction; exploratory Cox regression.
`17_druggability_candidate_drug_screening.R` — Druggability classification; drug-target curation.
`18a_docking_feasibility_check.py` — Environment check for AutoDock Vina and Open Babel.
`18b_run_vina_docking_IDO1_TGFBR1.py` — Molecular docking for IDO1 and TGFBR1.
`18c_merge_docking_results.py` — Neighbor residue analysis; merge docking results.

## Step 15 — Figure and Table Assembly
`19_generate_figure1_workflow.py` — Figure 1 workflow schematic.
`30_assemble_composite_figures_2_7.py` — Composite Figure 2–7 assembly.
`31_*.py`, `36*.py`, `38*.py` — Figure QC, repair, and expansion.
`34*.py`, `35*.py`, `39*.py`, `40*.py`, `41*.py` — Submission package generation.
