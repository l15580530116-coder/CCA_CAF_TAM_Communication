# ============================================================
# 分析目的: 从 TCGA-GDC 下载胆管癌 (CHOL) RNA-seq STAR Counts
#          数据及临床信息，保存原始对象和提取的矩阵/注释
# 输入文件: 无（通过 TCGAbiolinks 在线获取）
# 输出文件:
#   data/raw/TCGA_CHOL/TCGA_CHOL_STAR_counts.rds     — SummarizedExperiment
#   data/processed/TCGA_CHOL/TCGA_CHOL_counts_matrix.csv  — counts 矩阵
#   data/processed/TCGA_CHOL/TCGA_CHOL_gene_annotation.csv — 基因注释
#   data/clinical/TCGA_CHOL/TCGA_CHOL_clinical.csv        — 临床信息
# 主要方法: TCGAbiolinks (GDCquery → GDCdownload → GDCprepare)
# ============================================================

# ============================================================
# ⚠️ 警告: 本脚本需要从 TCGA-GDC 下载大量数据 (~50-200 MB)
#    下载时间取决于网络状况，可能需要 10-60 分钟。
#    请确保网络连接稳定，并预留足够磁盘空间。
#    确认后设置 RUN_DOWNLOAD <- TRUE 再执行。
# ============================================================

RUN_DOWNLOAD <- TRUE  # 改为 TRUE 以确认下载

if (!RUN_DOWNLOAD) {
  stop("请先阅读上方警告，确认后将 RUN_DOWNLOAD 改为 TRUE 再执行。")
}

# ---- 网络设置 ----
options(timeout = 600)  # 增加超时时间，防止大文件下载中断

# ---- 加载包 ----
library(TCGAbiolinks)
library(SummarizedExperiment)

# ---- 日志函数 ----
log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}

# ---- 创建目录 ----
dir.create("data/raw/TCGA_CHOL", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed/TCGA_CHOL", recursive = TRUE, showWarnings = FALSE)
dir.create("data/clinical/TCGA_CHOL", recursive = TRUE, showWarnings = FALSE)

log_msg("目录创建完成")

# ============================================================
# Step 1: 查询并下载 RNA-seq STAR Counts
# ============================================================
log_msg("Step 1/4: 查询 TCGA-CHOL RNA-seq STAR Counts 数据...")

query_rna <- GDCquery(
  project           = "TCGA-CHOL",
  data.category     = "Transcriptome Profiling",
  data.type         = "Gene Expression Quantification",
  workflow.type     = "STAR - Counts",
  access            = "open"
)

log_msg(sprintf("查询完成，共 %d 个样本", nrow(query_rna$results[[1]])))

log_msg("Step 2/4: 下载数据（可能需要较长时间）...")
GDCdownload(
  query       = query_rna,
  method      = "api",
  files.per.chunk = 3   # 减小每块文件数，避免下载截断
)

log_msg("Step 3/4: 准备 SummarizedExperiment 对象...")
se_obj <- GDCprepare(
  query = query_rna,
  summarizedExperiment = TRUE
)

log_msg(sprintf("SummarizedExperiment 对象: %d 基因 x %d 样本",
                nrow(se_obj), ncol(se_obj)))

# ---- 保存原始对象 ----
saveRDS(se_obj, "data/raw/TCGA_CHOL/TCGA_CHOL_STAR_counts.rds")
log_msg("原始对象已保存: data/raw/TCGA_CHOL/TCGA_CHOL_STAR_counts.rds")

# ============================================================
# Step 2: 提取 counts 矩阵和基因注释
# ============================================================
log_msg("Step 4/4: 提取表达矩阵和基因注释...")

# counts 矩阵
counts_matrix <- assay(se_obj, "unstranded")
counts_df <- as.data.frame(counts_matrix)
counts_df <- cbind(gene_id = rownames(counts_df), counts_df)

write.csv(counts_df,
          "data/processed/TCGA_CHOL/TCGA_CHOL_counts_matrix.csv",
          row.names = FALSE)
log_msg(sprintf("Counts 矩阵已保存: %d 行 x %d 列 (含 gene_id 列)",
                nrow(counts_df), ncol(counts_df)))

# 基因注释
gene_annot <- as.data.frame(rowRanges(se_obj))
write.csv(gene_annot,
          "data/processed/TCGA_CHOL/TCGA_CHOL_gene_annotation.csv",
          row.names = FALSE)
log_msg(sprintf("基因注释已保存: %d 行 x %d 列",
                nrow(gene_annot), ncol(gene_annot)))

# ============================================================
# Step 3: 下载临床信息
# ============================================================
log_msg("下载临床数据...")

query_clinical <- GDCquery_clinic(
  project    = "TCGA-CHOL",
  type       = "clinical"
)

write.csv(query_clinical,
          "data/clinical/TCGA_CHOL/TCGA_CHOL_clinical.csv",
          row.names = FALSE)
log_msg(sprintf("临床数据已保存: %d 行 x %d 列",
                nrow(query_clinical), ncol(query_clinical)))

# ============================================================
# 最终状态
# ============================================================
log_msg("========================================")
log_msg("TCGA-CHOL 数据下载完成！")
log_msg(sprintf("  data/raw/TCGA_CHOL/TCGA_CHOL_STAR_counts.rds     — %.1f MB",
                file.size("data/raw/TCGA_CHOL/TCGA_CHOL_STAR_counts.rds") / 1e6))
log_msg(sprintf("  data/processed/TCGA_CHOL/TCGA_CHOL_counts_matrix.csv  — %.1f MB",
                file.size("data/processed/TCGA_CHOL/TCGA_CHOL_counts_matrix.csv") / 1e6))
log_msg(sprintf("  data/processed/TCGA_CHOL/TCGA_CHOL_gene_annotation.csv — %.1f MB",
                file.size("data/processed/TCGA_CHOL/TCGA_CHOL_gene_annotation.csv") / 1e6))
log_msg(sprintf("  data/clinical/TCGA_CHOL/TCGA_CHOL_clinical.csv     — %.1f MB",
                file.size("data/clinical/TCGA_CHOL/TCGA_CHOL_clinical.csv") / 1e6))

log_msg("========================================")
log_msg("下一步: 运行 02 脚本，检查数据质量。")
