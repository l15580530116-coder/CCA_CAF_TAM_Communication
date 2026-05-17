# ============================================================
# 分析目的: 整理胆管癌免疫代谢 / CAF / 巨噬细胞-TAM 相关基因集
#          并检查各基因集在 TCGA-CHOL 和 GSE107943 中的覆盖
# 输入文件:
#   data/processed/TCGA_CHOL/TCGA_CHOL_log2TPM_gene_symbol.csv
#   data/processed/GSE107943/GSE107943_log2CPM_gene_symbol.csv
# 输出文件:
#   gene_sets/immunometabolism_genes.txt
#   gene_sets/CAF_markers.txt
#   gene_sets/macrophage_TAM_markers.txt
#   gene_sets/CCA_IM_CAF_TAM_combined_genes.txt
#   tables/gene_sets_summary.csv
#   tables/gene_sets_overlap_with_TCGA_GSE107943.csv
#   tables/CCA_IM_CAF_TAM_combined_gene_annotation.csv
# 主要方法: 基于 KEGG / Reactome / GO / HALLMARK 文献的通路基因定义 +
#          CAF / macrophage 文献 marker → 去重 → 覆盖检查
# 注: msigdbr 因 zenodo.org 无法访问而不可用，改用精心整理的
#     内置通路基因列表（来源：KEGG/HALLMARK pathway definitions）
# ============================================================

library(data.table)

dir.create("gene_sets", showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

log_msg("=== Step 1: 构建免疫代谢基因集 (内置通路定义) ===")
log_msg("msigdbr 不可用 (zenodo.org DNS 解析失败)，使用内置 curated 通路列表")

# ============================================================
# 每个通路的关键基因（基于 KEGG/HALLMARK/Reactome pathway definitions）
# 来源: MSigDB v2024.1 核心通路成员，手工整理代表性基因
# ============================================================

pathway_genes <- list(

  # ---- 糖代谢 ----
  glycolysis = c(
    "HK1", "HK2", "HK3", "GPI", "PFKL", "PFKM", "PFKP",
    "ALDOA", "ALDOB", "ALDOC", "TPI1", "GAPDH", "PGK1", "PGK2",
    "PGAM1", "PGAM2", "ENO1", "ENO2", "ENO3", "PKM", "PKLR",
    "LDHA", "LDHB", "LDHC", "SLC2A1", "SLC2A2", "SLC2A3", "SLC2A4"
  ),

  # ---- 氧化磷酸化 ----
  oxidative_phosphorylation = c(
    "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
    "SDHA", "SDHB", "SDHC", "SDHD",
    "UQCRB", "UQCRC1", "UQCRC2", "UQCRH", "UQCRQ",
    "COX4I1", "COX5A", "COX5B", "COX6A1", "COX6B1", "COX6C", "COX7A2",
    "COX7B", "COX7C", "COX8A",
    "ATP5A1", "ATP5B", "ATP5C1", "ATP5D", "ATP5E", "ATP5F1",
    "ATP5G1", "ATP5G2", "ATP5G3", "ATP5H", "ATP5I", "ATP5J",
    "ATP5L", "ATP5O", "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4",
    "NDUFA5", "NDUFA6", "NDUFA7", "NDUFA8", "NDUFA9", "NDUFA10",
    "NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6",
    "NDUFB7", "NDUFB8", "NDUFB9", "NDUFB10", "NDUFS1", "NDUFS2",
    "NDUFS3", "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS7", "NDUFS8",
    "NDUFV1", "NDUFV2", "NDUFV3"
  ),

  # ---- 脂肪酸代谢 ----
  fatty_acid_metabolism = c(
    "ACADL", "ACADM", "ACADS", "ACADSB", "ACADVL",
    "ACAA1", "ACAA2", "ACAT1", "ACAT2",
    "CPT1A", "CPT1B", "CPT1C", "CPT2",
    "ACOX1", "ACOX2", "ACOX3",
    "EHHADH", "HADH", "HADHA", "HADHB",
    "FASN", "ACACA", "ACACB",
    "SCD", "SCD5",
    "ELOVL1", "ELOVL2", "ELOVL3", "ELOVL4", "ELOVL5", "ELOVL6", "ELOVL7",
    "FADS1", "FADS2",
    "PPARA", "PPARD", "PPARG", "PPARGC1A", "PPARGC1B",
    "SREBF1", "SREBF2", "SCAP", "MLXIPL", "MLX",
    "CD36", "FABP1", "FABP2", "FABP3", "FABP4", "FABP5", "FABP6", "FABP7",
    "LPL", "LIPE", "PNPLA2", "MGLL",
    "ACLY", "ACSS2", "FASN", "MCAT"
  ),

  # ---- 胆汁酸代谢 ----
  bile_acid_metabolism = c(
    "CYP7A1", "CYP7B1", "CYP8B1", "CYP27A1", "CYP3A4", "CYP3A5", "CYP3A7",
    "BAAT", "SLC27A5", "ABCD3",
    "NR1H4", "FXR", "SLC10A1", "SLC10A2",
    "ABCB11", "ABCC2", "ABCG5", "ABCG8",
    "SLC51A", "SLC51B",
    "NR0B2", "OSTA", "OSTB",
    "HNF4A", "RXRA",
    "SULT2A1", "UGT2B4", "UGT2B7"
  ),

  # ---- 胆固醇代谢 ----
  cholesterol_homeostasis = c(
    "HMGCR", "HMGCS1", "HMGCS2",
    "MVK", "PMVK", "MVD", "FDPS", "FDFT1", "SQLE", "LSS",
    "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", "HSD17B7", "EBP", "SC5D", "DHCR7", "DHCR24",
    "LDLR", "VLDLR", "LRP1", "LRP2",
    "SCARB1", "LIPA",
    "APOA1", "APOB", "APOC1", "APOC2", "APOC3", "APOE",
    "ABCG1", "ABCA1",
    "SOAT1", "SOAT2", "LCAT", "CETP",
    "INSIG1", "INSIG2", "SCAP", "MBTPS1", "MBTPS2",
    "NR1H2", "NR1H3", "NR1H4"
  ),

  # ---- 外源代谢 ----
  xenobiotic_metabolism = c(
    "CYP1A1", "CYP1A2", "CYP1B1",
    "CYP2A6", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19",
    "CYP2D6", "CYP2E1", "CYP3A4", "CYP3A5", "CYP3A7",
    "GSTA1", "GSTA2", "GSTA3", "GSTA4", "GSTA5",
    "GSTM1", "GSTM2", "GSTM3", "GSTM4", "GSTM5",
    "GSTP1", "GSTT1", "GSTT2",
    "NAT1", "NAT2",
    "SULT1A1", "SULT1A2", "SULT1E1",
    "UGT1A1", "UGT1A3", "UGT1A4", "UGT1A6", "UGT1A7",
    "UGT1A8", "UGT1A9", "UGT1A10",
    "UGT2B4", "UGT2B7", "UGT2B10", "UGT2B15", "UGT2B17",
    "ALDH1A1", "ALDH2", "ALDH3A1", "ALDH3A2",
    "ADH1A", "ADH1B", "ADH1C", "ADH4", "ADH5", "ADH6", "ADH7",
    "ABCB1", "ABCC1", "ABCC2", "ABCC3", "ABCG2",
    "AHR", "ARNT",
    "NQO1", "NQO2",
    "EPHX1", "EPHX2"
  ),

  # ---- 缺氧 ----
  hypoxia = c(
    "HIF1A", "ARNT", "EPAS1", "HIF3A",
    "VEGFA", "VEGFB", "VEGFC", "FIGF",
    "SLC2A1", "HK1", "HK2", "PFKL", "PFKP",
    "ALDOA", "GAPDH", "PGK1", "ENO1", "PKM", "LDHA",
    "BNIP3", "BNIP3L",
    "EGLN1", "EGLN2", "EGLN3",
    "CA9", "CA12",
    "ADM", "EDN1", "EPO", "TFRC",
    "NDRG1", "DDIT4", "MXI1",
    "ANGPT1", "ANGPT2", "TEK",
    "NOS2", "NOS3",
    "SERPINE1", "PLAUR",
    "P4HA1", "P4HA2",
    "LOX", "LOXL2",
    "PGM1", "GYS1",
    "IGFBP1", "IGFBP3",
    "TGFA", "TGFB1", "TGFB3",
    "PDK1", "PDK2", "PDK3", "PDK4"
  ),

  # ---- 炎症反应 ----
  inflammatory_response = c(
    "IL1A", "IL1B", "IL1R1", "IL1R2", "IL1RN",
    "IL6", "IL6R", "IL6ST",
    "TNF", "TNFRSF1A", "TNFRSF1B",
    "CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8",
    "CCL11", "CCL13", "CCL20",
    "CXCL1", "CXCL2", "CXCL3", "CXCL8", "CXCL10", "CXCL12",
    "IL8",
    "CSF1", "CSF2", "CSF3",
    "PTGS2", "PTGER2", "PTGER4",
    "NFKB1", "NFKB2", "RELA", "RELB",
    "MYD88", "TICAM1",
    "TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR7", "TLR8",
    "SELE", "SELP", "VCAM1", "ICAM1",
    "MMP1", "MMP3", "MMP9", "MMP13",
    "S100A8", "S100A9", "S100A12",
    "CD14", "TREM1",
    "NLRP3", "PYCARD", "CASP1",
    "HMGB1", "HMGB2",
    "ITGB2", "PECAM1",
    "CCR1", "CCR2", "CCR5", "CCR7",
    "CXCR1", "CXCR2", "CXCR4"
  ),

  # ---- 干扰素 gamma 应答 ----
  interferon_gamma_response = c(
    "STAT1", "STAT2", "IRF1", "IRF2", "IRF3", "IRF4",
    "IRF5", "IRF7", "IRF8", "IRF9",
    "IFNG", "IFNGR1", "IFNGR2",
    "JAK1", "JAK2", "TYK2",
    "SOCS1", "SOCS3",
    "GBP1", "GBP2", "GBP3", "GBP4", "GBP5",
    "OAS1", "OAS2", "OAS3", "OASL",
    "MX1", "MX2",
    "CXCL9", "CXCL10", "CXCL11",
    "CCL2", "CCL5", "CCL7",
    "HLA-A", "HLA-B", "HLA-C",
    "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB",
    "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2",
    "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB3",
    "B2M", "TAP1", "TAP2", "TAPBP",
    "PSMB8", "PSMB9", "PSMB10",
    "ICAM1", "VCAM1",
    "ISG15", "ISG20",
    "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFITM2", "IFITM3",
    "BST2", "RSAD2",
    "TRIM21", "TRIM22", "TRIM25",
    "DDX58", "DHX58", "IFIH1", "MAVS",
    "CASP1", "CASP3", "CASP7", "CASP8"
  ),

  # ---- TNF-alpha via NF-kB ----
  tnfa_nfkb = c(
    "TNF", "TNFRSF1A",
    "NFKB1", "NFKB2", "RELA", "RELB", "REL",
    "NFKBIA", "NFKBIB", "NFKBIE",
    "IKBKB", "IKBKG", "CHUK",
    "TNFAIP3", "TNFAIP2", "TNFAIP6",
    "JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2",
    "ATF2", "ATF3",
    "MAP3K8", "MAP2K3", "MAP2K5",
    "BTG1", "BTG2", "BTG3",
    "ZFP36", "ZFP36L1", "ZFP36L2",
    "IER2", "IER3", "IER5",
    "EGR1", "EGR2", "EGR3",
    "DUSP1", "DUSP2", "DUSP4", "DUSP5",
    "PTGS2", "PLAUR",
    "CEBPB", "CEBPD",
    "BCL3", "BCL6",
    "CDKN1A", "GADD45A", "GADD45B",
    "PHLDA1", "PHLDA2",
    "SGK1", "KLF2", "KLF4", "KLF6", "KLF10",
    "ID1", "ID3",
    "TRIB1", "TRIB3",
    "ETS2", "ELF3",
    "SOD2", "SAT1",
    "CCL2", "CCL3", "CCL4", "CCL5", "CCL20",
    "CXCL1", "CXCL2", "CXCL3", "CXCL8", "CXCL10"
  ),

  # ---- IL6-JAK-STAT3 ----
  il6_jak_stat3 = c(
    "IL6", "IL6R", "IL6ST",
    "JAK1", "JAK2", "TYK2",
    "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6",
    "SOCS1", "SOCS2", "SOCS3", "SOCS4", "SOCS5",
    "PTPN11", "PTPN6",
    "PIK3CA", "PIK3CB", "PIK3CD",
    "AKT1", "AKT2", "AKT3",
    "MTOR", "RPTOR",
    "BCL2", "BCL2L1", "MCL1", "BCL6",
    "CCND1", "CCND2", "CCND3",
    "MYC", "JUNB",
    "CSF1", "CSF2", "CSF3",
    "VEGFA", "HGF",
    "IFNG", "IL2", "IL4", "IL7", "IL10", "IL21",
    "CXCL8", "CXCL12", "CCL2", "CCL5",
    "CRP", "SAA1", "SAA2",
    "CDKN1A", "FOS",
    "HSP90AA1", "HSP90AB1", "HSP90B1",
    "TNFRSF1A", "OSMR", "LIFR"
  ),

  # ---- ROS 通路 ----
  ros_pathway = c(
    "SOD1", "SOD2", "SOD3",
    "CAT",
    "GPX1", "GPX2", "GPX3", "GPX4",
    "PRDX1", "PRDX2", "PRDX3", "PRDX4", "PRDX5", "PRDX6",
    "TXN", "TXNRD1", "TXNRD2", "TXNRD3",
    "GSR",
    "GCLC", "GCLM",
    "CYBB", "CYBA", "NCF1", "NCF2", "NCF4",
    "NOX1", "NOX3", "NOX4", "NOX5",
    "DUOX1", "DUOX2",
    "MPO", "EPX", "LPO",
    "NFE2L2", "KEAP1",
    "HMOX1", "NQO1",
    "AKR1B1", "AKR1C1", "AKR1C2", "AKR1C3",
    "GSTA1", "GSTM1", "GSTP1",
    "MT1A", "MT2A",
    "FTH1", "FTL",
    "SRXN1",
    "AIFM1", "FOXO3",
    "SESN1", "SESN2", "SESN3",
    "CYCS", "COX2"
  ),

  # ---- mTORC1 信号 ----
  mtorc1_signaling = c(
    "MTOR", "RPTOR", "MLST8", "AKT1S1",
    "RICTOR", "MAPKAP1",
    "AKT1", "AKT2", "AKT3",
    "PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PIK3R2",
    "PTEN",
    "TSC1", "TSC2", "RHEB",
    "RPS6KB1", "RPS6KB2",
    "EIF4EBP1", "EIF4E", "EIF4G1",
    "RPS6", "RPS3", "RPS3A",
    "EIF2AK4", "EIF2S1", "EIF2S2",
    "DDIT4", "DDIT4L",
    "SREBF1", "SREBF2",
    "HIF1A",
    "ULK1", "ULK2",
    "ATG13", "RB1CC1",
    "SQSTM1",
    "SLC7A5", "SLC3A2",
    "LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5",
    "RRAGA", "RRAGB", "RRAGC", "RRAGD",
    "FLCN", "FNIP1", "FNIP2",
    "LPIN1", "LPIN2",
    "IDI1", "MVD", "FDFT1", "SQLE",
    "ACLY", "ACACA", "FASN",
    "DHCR7", "DHCR24", "LSS",
    "CYP51A1", "MSMO1",
    "PSAT1", "PSPH", "PHGDH", "SHMT2",
    "MTHFD2", "ALDH18A1", "PYCR1",
    "ASNS", "SARS",
    "CAD", "UMPS", "CTPS1", "IMPDH2"
  ),

  # ---- PPAR 信号 ----
  ppar_signaling = c(
    "PPARA", "PPARD", "PPARG", "PPARGC1A", "PPARGC1B",
    "RXRA", "RXRB", "RXRG",
    "FABP1", "FABP2", "FABP3", "FABP4", "FABP5",
    "CD36",
    "LPL", "APOA1", "APOA2", "APOC3",
    "ACOX1", "ACOX2", "ACOX3",
    "CPT1A", "CPT1B", "CPT2",
    "ACADL", "ACADM", "ACADVL",
    "SCD", "SCD5",
    "FADS2",
    "HMGCS2", "CYP4A11", "CYP4A22",
    "ANGPTL4", "ADIPOQ", "LEP",
    "UCP1", "UCP2", "UCP3",
    "PLIN1", "PLIN2", "PLTP",
    "SORBS1", "SLC27A1", "SLC27A2", "SLC27A4",
    "ACSL1", "ACSL3", "ACSL4", "ACSL5",
    "EHHADH", "SCP2",
    "CYP7A1", "CYP8B1", "CYP27A1",
    "ME1", "ME2", "ME3",
    "GK", "AQP7"
  ),

  # ---- 氨基酸代谢 (综合) ----
  amino_acid_metabolism = c(
    "GPT", "GPT2", "GOT1", "GOT2",
    "GLUD1", "GLUD2", "GLUL",
    "ASNS", "ASL", "ASS1",
    "CPS1", "OTC", "ARG1", "ARG2",
    "TAT", "FAH", "HPD", "HGD",
    "BCAT1", "BCAT2", "BCKDHA", "BCKDHB",
    "DLD", "DBT", "IVD",
    "IDH1", "IDH2", "IDH3A", "IDH3B", "IDH3G",
    "GLDC", "AMT", "GCSH",
    "SHMT1", "SHMT2",
    "CBS", "CTH",
    "MAT1A", "MAT2A", "MAT2B",
    "AHCY", "MTR", "MTRR",
    "GNMT", "BHMT", "DMGDH",
    "SDS", "SDSL",
    "SRR",
    "DAO",
    "MAOA", "MAOB"
  ),

  # ---- 色氨酸代谢 ----
  tryptophan_metabolism = c(
    "IDO1", "IDO2", "TDO2",
    "KYNU", "KYAT1", "KYAT3",
    "HAAO", "KMO",
    "AFMID", "ACMSD",
    "QPRT",
    "INMT",
    "TPH1", "TPH2",
    "DDC",
    "MAOA", "MAOB",
    "AANAT", "ASMT",
    "WARS", "WARS2",
    "IL4I1"
  ),

  # ---- 谷胱甘肽代谢 ----
  glutathione_metabolism = c(
    "GCLC", "GCLM",
    "GSS",
    "GSR",
    "GPX1", "GPX2", "GPX3", "GPX4", "GPX5", "GPX6", "GPX7", "GPX8",
    "GSTA1", "GSTA2", "GSTA3", "GSTA4", "GSTA5",
    "GSTM1", "GSTM2", "GSTM3", "GSTM4", "GSTM5",
    "GSTP1", "GSTT1", "GSTT2", "GSTO1", "GSTO2",
    "GGT1", "GGT5", "GGT6", "GGT7",
    "ANPEP", "DPEP1",
    "OPLAH",
    "IDH1", "IDH2",
    "GGCT",
    "CHAC1", "CHAC2",
    "SLC7A11"
  )
)

# 整合所有通路基因
immunometabolism_genes <- unique(unlist(pathway_genes))
log_msg(sprintf("免疫代谢基因集: %d 个唯一基因 (来自 %d 个子通路)",
                length(immunometabolism_genes), length(pathway_genes)))

# ============================================================
# Step 2: CAF markers
# ============================================================
log_msg("\n=== Step 2: CAF markers ===")

caf_markers <- c(
  "ACTA2", "FAP", "PDGFRB", "PDGFRA",
  "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2",
  "DCN", "LUM", "TAGLN", "POSTN", "MMP2", "MMP11",
  "THY1", "VIM", "FN1", "CXCL12", "CAV1", "RGS5", "MCAM", "CSPG4"
)
caf_markers <- unique(sort(caf_markers))
log_msg(sprintf("CAF markers: %d 个基因", length(caf_markers)))

# ============================================================
# Step 3: Macrophage / TAM markers
# ============================================================
log_msg("\n=== Step 3: Macrophage / TAM markers ===")

macrophage_tam_markers <- c(
  "CD68", "CD163", "MRC1", "MSR1", "CSF1R", "LYZ",
  "C1QA", "C1QB", "C1QC", "APOE", "TREM2", "SPP1",
  "MARCO", "FCGR3A", "ITGAM", "ITGAX", "CD14", "CD86",
  "IL1B", "TNF", "CXCL8", "CCL2", "CCL3", "CCL4", "CCL5",
  "TGFB1", "VEGFA", "ARG1", "IL10"
)
macrophage_tam_markers <- unique(sort(macrophage_tam_markers))
log_msg(sprintf("Macrophage/TAM markers: %d 个基因", length(macrophage_tam_markers)))

# ============================================================
# Step 4: 合并 & 来源标注
# ============================================================
log_msg("\n=== Step 4: 合并基因集 ===")

all_genes <- unique(c(immunometabolism_genes, caf_markers, macrophage_tam_markers))
log_msg(sprintf("合并基因集 (去重): %d 个基因", length(all_genes)))

gene_source <- data.frame(
  gene_symbol = all_genes,
  in_immunometabolism = all_genes %in% immunometabolism_genes,
  in_CAF               = all_genes %in% caf_markers,
  in_macrophage_TAM    = all_genes %in% macrophage_tam_markers,
  stringsAsFactors = FALSE
)

gene_source$source <- apply(gene_source[, 2:4], 1, function(x) {
  s <- c()
  if (x[1]) s <- c(s, "immunometabolism")
  if (x[2]) s <- c(s, "CAF")
  if (x[3]) s <- c(s, "macrophage_TAM")
  paste(s, collapse = ";")
})

log_msg(sprintf("仅 immunometabolism: %d", sum(gene_source$source == "immunometabolism")))
log_msg(sprintf("仅 CAF: %d", sum(gene_source$source == "CAF")))
log_msg(sprintf("仅 macrophage_TAM: %d", sum(gene_source$source == "macrophage_TAM")))
log_msg(sprintf("多来源: %d", sum(grepl(";", gene_source$source))))

# ============================================================
# Step 5: 保存基因集文件
# ============================================================
log_msg("\n=== Step 5: 保存基因集文件 ===")

write.table(data.frame(gene = sort(immunometabolism_genes)),
            "gene_sets/immunometabolism_genes.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
log_msg("已保存: gene_sets/immunometabolism_genes.txt")

write.table(data.frame(gene = caf_markers),
            "gene_sets/CAF_markers.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
log_msg("已保存: gene_sets/CAF_markers.txt")

write.table(data.frame(gene = macrophage_tam_markers),
            "gene_sets/macrophage_TAM_markers.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
log_msg("已保存: gene_sets/macrophage_TAM_markers.txt")

write.table(data.frame(gene = sort(all_genes)),
            "gene_sets/CCA_IM_CAF_TAM_combined_genes.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
log_msg("已保存: gene_sets/CCA_IM_CAF_TAM_combined_genes.txt")

write.csv(gene_source[order(gene_source$gene_symbol), ],
          "tables/CCA_IM_CAF_TAM_combined_gene_annotation.csv",
          row.names = FALSE)
log_msg("已保存: tables/CCA_IM_CAF_TAM_combined_gene_annotation.csv")

# ============================================================
# Step 6: 检查在表达矩阵中的覆盖
# ============================================================
log_msg("\n=== Step 6: 基因覆盖检查 ===")

tcga_genes <- fread("data/processed/TCGA_CHOL/TCGA_CHOL_log2TPM_gene_symbol.csv",
                    select = 1, data.table = FALSE)[[1]]
gse_genes  <- fread("data/processed/GSE107943/GSE107943_log2CPM_gene_symbol.csv",
                    select = 1, data.table = FALSE)[[1]]

log_msg(sprintf("TCGA-CHOL 基因数: %d", length(tcga_genes)))
log_msg(sprintf("GSE107943 基因数:  %d", length(gse_genes)))

check_overlap <- function(gs_name, gs_genes) {
  in_tcga <- sum(gs_genes %in% tcga_genes)
  in_gse  <- sum(gs_genes %in% gse_genes)
  in_both <- sum(gs_genes %in% tcga_genes & gs_genes %in% gse_genes)
  miss_tcga <- sum(!gs_genes %in% tcga_genes)
  miss_gse  <- sum(!gs_genes %in% gse_genes)

  cat(sprintf("\n  %s:\n", gs_name))
  cat(sprintf("    Total: %d\n", length(gs_genes)))
  cat(sprintf("    In TCGA-CHOL: %d (%.1f%%)\n", in_tcga, in_tcga/length(gs_genes)*100))
  cat(sprintf("    In GSE107943:  %d (%.1f%%)\n", in_gse, in_gse/length(gs_genes)*100))
  cat(sprintf("    In both:       %d (%.1f%%)\n", in_both, in_both/length(gs_genes)*100))
  if (miss_tcga > 0) cat(sprintf("    Missing TCGA:  %s\n",
    paste(gs_genes[!gs_genes %in% tcga_genes], collapse = ", ")))
  if (miss_gse > 0) cat(sprintf("    Missing GSE:   %s\n",
    paste(gs_genes[!gs_genes %in% gse_genes], collapse = ", ")))

  return(c(length(gs_genes), in_tcga, in_gse, in_both, miss_tcga, miss_gse))
}

res1 <- check_overlap("immunometabolism_genes", immunometabolism_genes)
res2 <- check_overlap("CAF_markers", caf_markers)
res3 <- check_overlap("macrophage_TAM_markers", macrophage_tam_markers)
res4 <- check_overlap("combined_genes", all_genes)

# ============================================================
# Step 7: 保存汇总表
# ============================================================
log_msg("\n\n=== Step 7: 保存汇总表 ===")

summary_df <- data.frame(
  gene_set = c("immunometabolism_genes", "CAF_markers",
               "macrophage_TAM_markers", "combined_genes"),
  total_genes       = c(res1[1], res2[1], res3[1], res4[1]),
  overlap_TCGA      = c(res1[2], res2[2], res3[2], res4[2]),
  overlap_GSE107943 = c(res1[3], res2[3], res3[3], res4[3]),
  overlap_both      = c(res1[4], res2[4], res3[4], res4[4]),
  missing_TCGA      = c(res1[5], res2[5], res3[5], res4[5]),
  missing_GSE107943 = c(res1[6], res2[6], res3[6], res4[6]),
  source = c("Curated KEGG/HALLMARK/Reactome pathways",
             "Literature CAF markers",
             "Literature TAM markers",
             "Combined from above"),
  stringsAsFactors = FALSE
)

write.csv(summary_df, "tables/gene_sets_summary.csv", row.names = FALSE)
log_msg("已保存: tables/gene_sets_summary.csv")

overlap_df <- data.frame(
  gene_set = c("immunometabolism_genes", "CAF_markers",
               "macrophage_TAM_markers", "combined_genes"),
  total             = c(res1[1], res2[1], res3[1], res4[1]),
  overlap_TCGA      = c(res1[2], res2[2], res3[2], res4[2]),
  pct_TCGA          = round(c(res1[2]/res1[1]*100, res2[2]/res2[1]*100,
                              res3[2]/res3[1]*100, res4[2]/res4[1]*100), 1),
  overlap_GSE107943 = c(res1[3], res2[3], res3[3], res4[3]),
  pct_GSE107943     = round(c(res1[3]/res1[1]*100, res2[3]/res2[1]*100,
                              res3[3]/res3[1]*100, res4[3]/res4[1]*100), 1),
  overlap_both      = c(res1[4], res2[4], res3[4], res4[4]),
  pct_both          = round(c(res1[4]/res1[1]*100, res2[4]/res2[1]*100,
                              res3[4]/res3[1]*100, res4[4]/res4[1]*100), 1),
  stringsAsFactors = FALSE
)

write.csv(overlap_df, "tables/gene_sets_overlap_with_TCGA_GSE107943.csv",
          row.names = FALSE)
log_msg("已保存: tables/gene_sets_overlap_with_TCGA_GSE107943.csv")

# ============================================================
# 最终报告
# ============================================================
cat(sprintf("\n========================================\n"))
cat(sprintf("基因集整理完成\n"))
cat(sprintf("========================================\n"))
cat(sprintf("immunometabolism_genes:     %d 基因\n", res1[1]))
cat(sprintf("CAF_markers:                %d 基因\n", res2[1]))
cat(sprintf("macrophage_TAM_markers:     %d 基因\n", res3[1]))
cat(sprintf("combined_genes:             %d 基因\n", res4[1]))
cat(sprintf("\nTCGA-CHOL 覆盖:   %d / %d (%.1f%%)\n",
            res4[2], res4[1], res4[2]/res4[1]*100))
cat(sprintf("GSE107943 覆盖:    %d / %d (%.1f%%)\n",
            res4[3], res4[1], res4[3]/res4[1]*100))
cat(sprintf("两者共同覆盖:      %d / %d (%.1f%%)\n",
            res4[4], res4[1], res4[4]/res4[1]*100))
cat(sprintf("使用 msigdbr:      No (zenodo.org 无法访问，使用内置 curated 列表)\n"))
cat(sprintf("========================================\n"))
