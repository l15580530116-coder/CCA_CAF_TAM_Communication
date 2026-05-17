# ============================================================
# 分析目的: 预处理 TCGA-CHOL 表达矩阵 —— Ensembl ID → gene symbol、
#           log2(TPM+1) 标准化、样本分组、匹配生存信息
# 输入文件:
#   data/raw/TCGA_CHOL/TCGA_CHOL_STAR_counts.rds
#   data/processed/TCGA_CHOL/TCGA_CHOL_gene_annotation.csv
#   data/clinical/TCGA_CHOL/TCGA_CHOL_survival_clean.csv
# 输出文件:
#   data/processed/TCGA_CHOL/TCGA_CHOL_counts_gene_symbol.csv
#   data/processed/TCGA_CHOL/TCGA_CHOL_log2TPM_gene_symbol.csv
#   data/processed/TCGA_CHOL/TCGA_CHOL_sample_metadata.csv
#   data/processed/TCGA_CHOL/TCGA_CHOL_tumor_expression_for_survival.csv
#   data/processed/TCGA_CHOL/TCGA_CHOL_matched_survival_metadata.csv
# 主要方法: SummarizedExperiment → TPM log2 标准化 → 基因符号映射
#          → 样本亚型鉴定 → 与临床生存数据匹配
# ============================================================

# ---- 加载包 ----
library(SummarizedExperiment)
library(data.table)

# ---- 创建目录 ----
dir.create("data/processed/TCGA_CHOL", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Step 1: 读取 SE 对象，检查 assay
# ============================================================
cat("=== Step 1: 读取数据 ===\n")
se <- readRDS("data/raw/TCGA_CHOL/TCGA_CHOL_STAR_counts.rds")
cat(sprintf("可用 assay: %s\n", paste(assayNames(se), collapse = ", ")))

# ============================================================
# Step 2: 提取 raw counts 和 TPM
# ============================================================
cat("\n=== Step 2: 提取表达矩阵 ===\n")

raw_counts <- assay(se, "unstranded")
cat(sprintf("Raw counts: %d genes x %d samples\n", nrow(raw_counts), ncol(raw_counts)))

# 判断使用 TPM 还是 CPM
if ("tpm_unstrand" %in% assayNames(se)) {
  use_tpm <- TRUE
  tpm <- assay(se, "tpm_unstrand")
  log2_expr <- log2(tpm + 1)
  norm_method <- "log2(TPM+1)"
} else {
  use_tpm <- FALSE
  lib_sizes <- colSums(raw_counts)
  cpm <- t(t(raw_counts) / lib_sizes * 1e6)
  log2_expr <- log2(cpm + 1)
  norm_method <- "log2(CPM+1)"
}
cat(sprintf("标准化方法: %s\n", norm_method))

# ============================================================
# Step 3: 基因注释 —— Ensembl ID → gene symbol
# ============================================================
cat("\n=== Step 3: 基因符号映射 ===\n")

gene_annot <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_gene_annotation.csv",
                       stringsAsFactors = FALSE)
cat(sprintf("基因注释: %d 行\n", nrow(gene_annot)))

# 剥离 Ensembl ID 版本号
gene_annot$ensembl_id <- sub("\\..*", "", gene_annot$gene_id)

# 确保顺序与表达矩阵一致
stopifnot(identical(rownames(raw_counts), gene_annot$gene_id))
rownames(raw_counts) <- gene_annot$ensembl_id
rownames(log2_expr) <- gene_annot$ensembl_id

# 添加 gene symbol
gene_annot$gene_name[is.na(gene_annot$gene_name) | gene_annot$gene_name == ""] <- NA
n_no_symbol <- sum(is.na(gene_annot$gene_name))
cat(sprintf("缺失 gene symbol: %d\n", n_no_symbol))

# ============================================================
# Step 4: 聚合多个 Ensembl ID → gene symbol
# ============================================================
cat("\n=== Step 4: 聚合重复基因符号 ===\n")

# 只保留有 gene symbol 的
keep <- !is.na(gene_annot$gene_name)
raw_counts <- raw_counts[keep, , drop = FALSE]
log2_expr <- log2_expr[keep, , drop = FALSE]
gene_annot <- gene_annot[keep, , drop = FALSE]
cat(sprintf("过滤后: %d 个 Ensembl ID\n", nrow(raw_counts)))

# raw counts: 按 gene symbol 求和
cat("聚合 raw counts (sum)...\n")
gene_symbols <- gene_annot$gene_name
raw_counts_symbol <- rowsum(raw_counts, group = gene_symbols, reorder = TRUE)

# log2 TPM/CPM: 按 gene symbol 取 mean
cat("聚合 log2 expression (mean)...\n")
log2_expr_symbol <- apply(log2_expr, 2, function(col) {
  tapply(col, gene_symbols, mean, na.rm = TRUE)
})

n_genes <- nrow(raw_counts_symbol)
cat(sprintf("Gene symbol 数量: %d\n", n_genes))

# ---- 保存 counts gene symbol ----
counts_out <- data.frame(gene_symbol = rownames(raw_counts_symbol), raw_counts_symbol,
                         check.names = FALSE)
write.csv(counts_out,
          "data/processed/TCGA_CHOL/TCGA_CHOL_counts_gene_symbol.csv",
          row.names = FALSE)
cat(sprintf("已保存: data/processed/TCGA_CHOL/TCGA_CHOL_counts_gene_symbol.csv (%d 基因 x %d 样本)\n",
            nrow(raw_counts_symbol), ncol(raw_counts_symbol)))

# ---- 保存 log2 表达矩阵 ----
log2_out <- data.frame(gene_symbol = rownames(log2_expr_symbol), log2_expr_symbol,
                       check.names = FALSE)
tpm_cpm_name <- if (use_tpm) "TCGA_CHOL_log2TPM_gene_symbol.csv" else "TCGA_CHOL_log2CPM_gene_symbol.csv"
write.csv(log2_out,
          paste0("data/processed/TCGA_CHOL/", tpm_cpm_name),
          row.names = FALSE)
cat(sprintf("已保存: data/processed/TCGA_CHOL/%s (%d 基因 x %d 样本)\n",
            tpm_cpm_name, nrow(log2_expr_symbol), ncol(log2_expr_symbol)))

# ============================================================
# Step 5: 样本分组 —— 根据 TCGA barcode 判断 Tumor/Normal
# ============================================================
cat("\n=== Step 5: 样本分组 ===\n")

# TCGA barcode: TCGA-XX-XXXX-XX-...
# 位置 14-15: 01-09 = Tumor, 10-19 = Normal, 20-29 = Metastatic
barcodes <- colnames(se)
sample_types <- substr(barcodes, 14, 15)

sample_meta <- data.frame(
  sample_id        = barcodes,
  patient_id       = substr(barcodes, 1, 12),
  sample_type_code = sample_types,
  stringsAsFactors  = FALSE
)

# 分类
sample_meta$sample_type <- ifelse(sample_meta$sample_type_code == "01",
                                  "Primary Tumor",
                                  ifelse(sample_meta$sample_type_code == "11",
                                         "Solid Tissue Normal",
                                         "Other"))
cat(sprintf("Tumor 样本: %d\n", sum(sample_meta$sample_type == "Primary Tumor")))
cat(sprintf("Normal 样本: %d\n", sum(sample_meta$sample_type == "Solid Tissue Normal")))
cat(sprintf("其他: %d\n", sum(sample_meta$sample_type == "Other")))

# ---- 保存 sample metadata ----
write.csv(sample_meta,
          "data/processed/TCGA_CHOL/TCGA_CHOL_sample_metadata.csv",
          row.names = FALSE)
cat("已保存: data/processed/TCGA_CHOL/TCGA_CHOL_sample_metadata.csv\n")

# ============================================================
# Step 6: 匹配肿瘤样本与生存信息
# ============================================================
cat("\n=== Step 6: 匹配生存信息 ===\n")

survival <- read.csv("data/clinical/TCGA_CHOL/TCGA_CHOL_survival_clean.csv",
                     stringsAsFactors = FALSE)
cat(sprintf("生存数据: %d 例\n", nrow(survival)))

# 提取肿瘤样本
tumor_meta <- sample_meta[sample_meta$sample_type == "Primary Tumor", , drop = FALSE]
cat(sprintf("Tumor 样本: %d\n", nrow(tumor_meta)))

# 匹配 patient_id
tumor_meta$survival_matched <- tumor_meta$patient_id %in% survival$patient_id
n_matched <- sum(tumor_meta$survival_matched)
cat(sprintf("匹配到生存信息的 tumor 样本: %d\n", n_matched))

# 提取匹配肿瘤的表达矩阵
tumor_samples <- tumor_meta$sample_id[tumor_meta$survival_matched]
tumor_expr <- log2_out[, c("gene_symbol", tumor_samples), drop = FALSE]
write.csv(tumor_expr,
          "data/processed/TCGA_CHOL/TCGA_CHOL_tumor_expression_for_survival.csv",
          row.names = FALSE)
cat(sprintf("已保存: data/processed/TCGA_CHOL/TCGA_CHOL_tumor_expression_for_survival.csv (%d 基因 x %d 样本)\n",
            nrow(tumor_expr), ncol(tumor_expr) - 1))

# 匹配的临床信息
matched_survival <- survival[survival$patient_id %in% tumor_meta$patient_id[tumor_meta$survival_matched], , drop = FALSE]
write.csv(matched_survival,
          "data/processed/TCGA_CHOL/TCGA_CHOL_matched_survival_metadata.csv",
          row.names = FALSE)
cat(sprintf("已保存: data/processed/TCGA_CHOL/TCGA_CHOL_matched_survival_metadata.csv (%d 例)\n",
            nrow(matched_survival)))

# ============================================================
# 最终统计报告
# ============================================================
n_tumor  <- sum(sample_meta$sample_type == "Primary Tumor")
n_normal <- sum(sample_meta$sample_type == "Solid Tissue Normal")

cat(sprintf("\n========================================\n"))
cat(sprintf("TCGA-CHOL 表达预处理完成\n"))
cat(sprintf("========================================\n"))
cat(sprintf("原始的 Ensembl ID 数量:  %d\n", nrow(raw_counts)))
cat(sprintf("标准化方法:              %s\n", norm_method))
cat(sprintf("Gene symbol 数量:        %d\n", n_genes))
cat(sprintf("Tumor 样本数:            %d\n", n_tumor))
cat(sprintf("Normal 样本数:           %d\n", n_normal))
cat(sprintf("匹配到生存的 Tumor 样本: %d\n", n_matched))
cat(sprintf("========================================\n"))
