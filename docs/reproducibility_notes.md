# Reproducibility Notes

## Environment

- **R** 4.5.3 with Bioconductor 3.22
- **Python** 3.14
- **AutoDock Vina** v1.2.7
- **Open Babel** v3.1.1

## Known Dependencies

1. **Internet access required**: Scripts 01, 04, 11, and 14 download data from GDC and GEO.
2. **Rtools on Windows**: CellChat v2.2.0.9001 installation from GitHub requires Rtools for C++ compilation.
3. **AutoDock Vina and Open Babel**: Required for molecular docking (Steps 18a–18c). Must be installed separately and available in PATH.
4. **msigdbr**: If msigdbr is unavailable (e.g., Zenodo DNS issues), the gene set curation script uses built-in KEGG/HALLMARK pathway definitions.

## Computational Caveats

- **Survival analyses**: All survival analyses in this study are exploratory. Small cohort sizes (TCGA-CHOL n=35, GSE107943 n=30) produce wide confidence intervals. Multivariate Cox and LASSO-Cox were not performed due to insufficient sample size (events-per-variable ratio < 10).
- **CellChat**: CellChat provides computational inference of intercellular communication based on ligand–receptor co-expression. Inferred axes require experimental validation and should not be interpreted as established functional signaling.
- **Molecular docking**: All docking scores are computational predictions. No molecular dynamics simulations were performed. Docking results do not represent experimental binding affinities.

## Large Files

- GSE138709 processed Seurat object (~502 MB RDS) is not included in this repository.
- Raw sequencing data (FASTQ/BAM) are available from original publications via GEO/SRA.
