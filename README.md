# CAF–TAM Communication Shapes an Aggressive Microenvironment in Cholangiocarcinoma

This repository contains analysis scripts for an integrated computational study combining bulk transcriptomics, single-cell RNA-seq, CellChat-inferred intercellular communication, hub gene prioritization, druggability assessment, and in silico molecular docking in cholangiocarcinoma (CCA).

## Public Datasets

All raw data are publicly available and must be downloaded separately:

| Dataset | Accession | Platform | Samples |
|---------|-----------|----------|---------|
| TCGA-CHOL | [GDC](https://portal.gdc.cancer.gov/projects/TCGA-CHOL) | RNA-seq (STAR Counts) | 44 (35 tumor, 9 normal) |
| GSE107943 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107943) | RNA-seq | 57 (30 tumor, 27 paired normal) |
| GSE138709 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138709) | 10X scRNA-seq | 32,626 cells (5 iCCA patients) |
| GSE26566 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26566) | Microarray (GPL6104) | 169 samples |

## Repository Structure

```
├── README.md
├── LICENSE
├── .gitignore
├── environment/
│   ├── python_requirements.txt
│   └── package_versions.md
├── scripts/
│   ├── 00_check_environment.R
│   ├── 01_download_TCGA_CHOL.R
│   ├── ... (numbered by analysis order)
│   └── 41_fix_supplementary_tables.py
├── gene_sets/
│   └── curated gene set definitions
├── docs/
│   ├── workflow.md
│   ├── data_availability.md
│   └── reproducibility_notes.md
├── tables/
│   └── README.md
└── figures/
    └── README.md
```

## Analysis Workflow

| Step | Script(s) | Description |
|------|-----------|-------------|
| 0 | `00_check_environment.R` | Environment verification |
| 1 | `01_download_TCGA_CHOL.R` | TCGA-CHOL download |
| 2 | `02_clean_TCGA_CHOL_clinical.R` | Clinical data cleaning |
| 3 | `03_preprocess_TCGA_CHOL_expression.R` | Expression preprocessing |
| 4 | `04_download_preprocess_GSE107943.R` | GSE107943 download + preprocessing |
| 5 | `05_prepare_gene_sets.R` | IM/CAF/TAM gene set curation |
| 6 | `06_DEG_TCGA_GSE107943.R` | Differential expression analysis |
| 7 | `07_functional_enrichment_validated_DEGs.R` | GO/KEGG enrichment |
| 8 | `08_GSVA_ssGSEA_scores_survival.R` | ssGSEA pathway scoring |
| 9 | `10_immune_infiltration_checkpoint_analysis.R` | Immune/checkpoint correlation |
| 10 | `11*.R`, `11c_*.R` | scRNA-seq preprocessing + annotation |
| 11 | `12_CellChat_*.R` | CellChat communication inference |
| 12 | `13_integrated_hub_gene_*.R` | Hub gene prioritization |
| 13 | `14_validate_hub_genes_GSE26566.R` | GSE26566 external validation |
| 14 | `17_druggability_*.R`, `18*.py` | Druggability + molecular docking |
| 15 | `19_*.py`, `3*_*.py` | Figure and table assembly |

## Software Requirements

- **R** 4.5.3 with packages: TCGAbiolinks, DESeq2, edgeR, limma, GSVA, Seurat (v5), CellChat (v2.2.0.9001), clusterProfiler, ConsensusClusterPlus, WGCNA
- **Python** 3.14 with packages: numpy, matplotlib, scipy, openpyxl, PyMuPDF, python-docx
- **AutoDock Vina** v1.2.7
- **Open Babel** v3.1.1

## Reproducibility

Raw data are not included due to size and database redistribution policies. Users should download data from GDC and GEO using the accession numbers above. Scripts are numbered by analysis order and include headers specifying purpose, input, output, and methods.

## Funding

This work was supported by the National Natural Science Foundation of China (No. 82260136), and the Finance Science and Technology Project of Hainan Province (Nos. ZDYF2021SHFZ053 and YSPTZX202027).

## Citation

If using this repository, please cite the associated manuscript once available.

## Contact

Corresponding authors:
- Xiang Yang, xiangyang200611@126.com
- Zhang Shufang, zsf66189665@126.com

Haikou Affiliated Hospital of Central South University Xiangya School of Medicine, Haikou 570208, China
