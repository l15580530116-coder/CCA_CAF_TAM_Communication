# ============================================================
# 分析目的: 下载并预处理 GSE107943 bulk RNA-seq 数据
#          — 作为 TCGA-CHOL 的外部 RNA-seq 验证队列
# 输入/下载来源:
#   GEO accession: GSE107943
#   Supplementary: GSE107943_rawread.txt.gz (raw counts, T/N 命名)
#   Series matrix 含临床信息 (OS, DFS, recurrence, stage 等)
# 输出文件:
#   data/raw/GSE107943/GSE107943_rawread.txt.gz
#   data/raw/GSE107943/GSE107943_series_matrix.txt.gz
#   data/processed/GSE107943/GSE107943_counts_gene_symbol.csv
#   data/processed/GSE107943/GSE107943_log2CPM_gene_symbol.csv
#   data/processed/GSE107943/GSE107943_sample_metadata.csv
#   data/clinical/GSE107943/GSE107943_clinical.csv
#   tables/GSE107943_preprocess_summary.csv
# 主要方法: GEOquery 下载 → raw counts → log2(CPM+1) →
#          gene symbol 聚合 → T/N 分组 → 临床信息提取
# 注: GSE107943_RPKM.txt.gz 的列名为加密 sample ID，
#     无法可靠映射到 T/N 样本名，因此改用从 raw counts
#     计算 CPM 作为归一化表达量。raw counts 采样名为
#     清晰的 T/N 格式 (1T, 2T, 1N...)，可直接分组。
# ============================================================

# ---- 网络设置 ----
options(timeout = 600)

# ---- 加载包 ----
library(GEOquery)
library(data.table)

# ---- 创建目录 ----
for (d in c("data/raw/GSE107943",
            "data/processed/GSE107943",
            "data/clinical/GSE107943",
            "tables")) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- 日志 ----
log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}

# ============================================================
# Step 1: 下载数据
# ============================================================
log_msg("=== Step 1: 下载数据 ===")

rawread_dest <- "data/raw/GSE107943/GSE107943_rawread.txt.gz"
if (!file.exists(rawread_dest)) {
  url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107943/suppl/GSE107943_rawread.txt.gz"
  log_msg(sprintf("下载 rawread: %s", url))
  download.file(url, rawread_dest, mode = "wb")
  log_msg(sprintf("完成: %.1f MB", file.size(rawread_dest) / 1e6))
} else {
  log_msg("rawread 文件已存在")
}

series_mtx_dest <- "data/raw/GSE107943/GSE107943_series_matrix.txt.gz"
if (!file.exists(series_mtx_dest)) {
  log_msg("下载 series matrix...")
  getGEO("GSE107943", GSEMatrix = TRUE, getGPL = FALSE,
         destdir = "data/raw/GSE107943")
} else {
  log_msg("Series matrix 已存在")
}

# ============================================================
# Step 2: 读取 rawread 文件
# ============================================================
log_msg("\n=== Step 2: 读取 rawread 表达矩阵 ===")

# rawread 文件格式:
# 前 7 列: No, Chr, Ensenble, Genesymbol, Start, Stop, CodingLength
# 后续列: 样本表达量 (1T, 2T, 3T, ..., 1N, 2N, ...)

raw <- fread(rawread_dest, data.table = FALSE)
log_msg(sprintf("原始文件: %d 行 x %d 列", nrow(raw), ncol(raw)))

# 提取基因注释列
anno_cols <- c("Chr", "Ensenble", "Genesymbol", "Start", "Stop", "CodingLength")
gene_anno <- raw[, anno_cols, drop = FALSE]
gene_anno$Ensembl <- sub("\\..*", "", gene_anno$Ensenble)

# 提取基因 symbol
gene_symbols <- raw$Genesymbol
n_na_symbol <- sum(is.na(gene_symbols) | gene_symbols == "")
log_msg(sprintf("基因 symbol: %d 个, 缺失: %d", length(gene_symbols), n_na_symbol))

# 提取表达矩阵 (第 8 列起)
expr_col_start <- 8
expr_cols <- colnames(raw)[expr_col_start:ncol(raw)]
log_msg(sprintf("表达列数: %d", length(expr_cols)))
log_msg(sprintf("样本名示例: %s", paste(head(expr_cols, 10), collapse = ", ")))

counts_mat <- as.matrix(raw[, expr_cols, drop = FALSE])
rownames(counts_mat) <- gene_symbols
log_msg(sprintf("Counts 矩阵: %d 基因 x %d 样本", nrow(counts_mat), ncol(counts_mat)))

# ============================================================
# Step 3: 过滤 & 聚合 gene symbol
# ============================================================
log_msg("\n=== Step 3: Gene symbol 过滤与聚合 ===")

# 移除缺失 symbol
keep <- !is.na(gene_symbols) & gene_symbols != ""
counts_mat <- counts_mat[keep, , drop = FALSE]
gs <- gene_symbols[keep]
log_msg(sprintf("有 gene symbol 的基因: %d", length(gs)))

# raw counts: sum by gene symbol
log_msg("聚合 raw counts (sum)...")
counts_symbol <- rowsum(counts_mat, group = gs, reorder = TRUE)
n_genes <- nrow(counts_symbol)
log_msg(sprintf("Gene symbol 数量: %d", n_genes))

# ---- 保存 counts ----
counts_out <- data.frame(gene_symbol = rownames(counts_symbol),
                         counts_symbol, check.names = FALSE)
write.csv(counts_out,
          "data/processed/GSE107943/GSE107943_counts_gene_symbol.csv",
          row.names = FALSE)
log_msg("已保存: GSE107943_counts_gene_symbol.csv")

# ============================================================
# Step 4: 计算 log2(CPM + 1)
# ============================================================
log_msg("\n=== Step 4: 计算 log2(CPM + 1) ===")

lib_sizes <- colSums(counts_symbol)
log_msg(sprintf("文库大小范围: %.0f - %.0f (million reads)",
                min(lib_sizes)/1e6, max(lib_sizes)/1e6))

cpm <- t(t(counts_symbol) / lib_sizes * 1e6)
log2cpm <- log2(cpm + 1)

# ---- 保存 ----
log2cpm_out <- data.frame(gene_symbol = rownames(log2cpm), log2cpm,
                          check.names = FALSE)
write.csv(log2cpm_out,
          "data/processed/GSE107943/GSE107943_log2CPM_gene_symbol.csv",
          row.names = FALSE)
log_msg(sprintf("已保存: GSE107943_log2CPM_gene_symbol.csv (%d x %d)",
                nrow(log2cpm), ncol(log2cpm)))

# ============================================================
# Step 5: 样本分组 T/N
# ============================================================
log_msg("\n=== Step 5: 样本分组 ===")

sample_names <- colnames(counts_symbol)

# T/N 后缀格式识别
ends_T <- grepl("T$", sample_names)
ends_N <- grepl("N$", sample_names)

sample_meta <- data.frame(
  sample_id   = sample_names,
  pair_id     = gsub("[TN]$", "", sample_names),
  sample_type = ifelse(ends_T, "Tumor",
                        ifelse(ends_N, "Normal", "Unknown")),
  stringsAsFactors = FALSE
)

n_tumor  <- sum(sample_meta$sample_type == "Tumor")
n_normal <- sum(sample_meta$sample_type == "Normal")

# 检查配对
paired_numbers <- intersect(
  sample_meta$pair_id[sample_meta$sample_type == "Tumor"],
  sample_meta$pair_id[sample_meta$sample_type == "Normal"]
)
n_paired <- length(paired_numbers)

log_msg(sprintf("Tumor: %d, Normal: %d, Paired: %d 对", n_tumor, n_normal, n_paired))

write.csv(sample_meta,
          "data/processed/GSE107943/GSE107943_sample_metadata.csv",
          row.names = FALSE)
log_msg("已保存: GSE107943_sample_metadata.csv")

# ============================================================
# Step 6: 提取临床信息
# ============================================================
log_msg("\n=== Step 6: 提取临床信息 ===")

# 从 series matrix 直接解析（避免 GEOquery segfault）
# 临床数据格式: !Sample_characteristics_ch1 行，每行一个变量名
#   "tissue: ...", "Sex: ...", "age: ...", "survival(mo): ..." 等
series_lines <- readLines(gzfile(series_mtx_dest))

# 获取样本名
sample_titles <- NULL
title_line <- grep("^!Sample_title\t", series_lines, value = TRUE)
if (length(title_line) > 0) {
  sample_titles <- strsplit(title_line, "\t")[[1]][-1]
  sample_titles <- gsub('"', '', sample_titles)
}
log_msg(sprintf("Series matrix 中的样本数: %d", length(sample_titles)))

# 提取所有 characteristics_ch1 行
char_lines <- grep("^!Sample_characteristics_ch1\t", series_lines, value = TRUE)
log_msg(sprintf("Characteristics 行数: %d", length(char_lines)))

# 解析每行，提取变量名和值
clin_data <- data.frame(sample_name = sample_titles, stringsAsFactors = FALSE)

for (cline in char_lines) {
  vals <- strsplit(cline, "\t")[[1]][-1]
  vals <- gsub('"', '', vals)

  if (length(vals) != length(sample_titles)) next

  # 从第一个值提取变量名 (如 "age: 54" → "age")
  var_name <- sub(": .*", "", vals[1])
  # 清除所有值中的变量名前缀 (如 "age: 54" → "54")
  clean_vals <- sub("^[^:]+: ", "", vals)
  # 清理变量名中的特殊字符
  col_name <- gsub("[^a-zA-Z0-9_.]", "_", var_name)
  col_name <- gsub("_{2,}", "_", col_name)

  clin_data[[col_name]] <- clean_vals
}

log_msg(sprintf("临床表: %d 行 x %d 列", nrow(clin_data), ncol(clin_data)))
log_msg(sprintf("  变量: %s", paste(setdiff(colnames(clin_data), "sample_name"), collapse = ", ")))

# 检查关键变量
key_vars <- c("tissue", "Sex", "age", "survival_mo_", "death", "recurr",
              "dsfree_mo_", "stageajcc")
found_key <- intersect(key_vars, colnames(clin_data))
log_msg(sprintf("  关键变量: %s", paste(found_key, collapse = ", ")))

write.csv(clin_data,
          "data/clinical/GSE107943/GSE107943_clinical.csv",
          row.names = FALSE)
log_msg("已保存: GSE107943_clinical.csv")

# ============================================================
# Step 7: 生成概览表
# ============================================================
log_msg("\n=== Step 7: 生成概览表 ===")

has_os <- "survival_mo_" %in% colnames(clin_data)
has_death <- "death" %in% colnames(clin_data)
has_recurr <- "recurr" %in% colnames(clin_data)
has_dfs <- "dsfree_mo_" %in% colnames(clin_data)

summary_rows <- list(
  c("数据集", "GSE107943"),
  c("癌种", "Intrahepatic cholangiocarcinoma (iCCA)"),
  c("平台", "GPL18573 (Illumina HiSeq)"),
  c("数据类型", "RNA-seq"),
  c("归一化方法", "log2(CPM+1) [从 raw counts 计算]"),
  c("基因数量 (gene symbol)", as.character(n_genes)),
  c("总样本数 (表达)", as.character(ncol(counts_symbol))),
  c("Tumor 样本数", as.character(n_tumor)),
  c("Normal 样本数", as.character(n_normal)),
  c("Paired 样本", sprintf("%d 对", n_paired)),
  c("Counts 矩阵维度", sprintf("%d x %d", nrow(counts_symbol), ncol(counts_symbol))),
  c("log2CPM 矩阵维度", sprintf("%d x %d", nrow(log2cpm), ncol(log2cpm))),
  c("OS 信息 (survival(mo))", if (has_os) "Yes" else "No"),
  c("Death 信息", if (has_death) "Yes" else "No"),
  c("Recurrence 信息", if (has_recurr) "Yes" else "No"),
  c("DFS 信息", if (has_dfs) "Yes" else "No"),
  c("临床样本数", as.character(nrow(clin_data))),
  c("RPKM 文件状态", "已下载但样本ID无法映射，未使用")
)

summary_df <- do.call(rbind, summary_rows)
colnames(summary_df) <- c("指标", "值")

write.csv(summary_df,
          "tables/GSE107943_preprocess_summary.csv",
          row.names = FALSE)
log_msg("已保存: tables/GSE107943_preprocess_summary.csv")

# ============================================================
# 最终报告
# ============================================================
cat(sprintf("\n========================================\n"))
cat(sprintf("GSE107943 下载预处理完成\n"))
cat(sprintf("========================================\n"))
cat(sprintf("归一化方法:          log2(CPM+1) [从 raw counts]\n"))
cat(sprintf("Gene symbol 数量:     %d\n", n_genes))
cat(sprintf("Tumor 样本数:         %d\n", n_tumor))
cat(sprintf("Normal 样本数:        %d\n", n_normal))
cat(sprintf("Paired sample:        %d 对\n", n_paired))
cat(sprintf("Counts 矩阵:          %d x %d\n",
            nrow(counts_symbol), ncol(counts_symbol)))
cat(sprintf("临床 OS 信息:         %s\n", if (has_os) "Yes" else "No"))
cat(sprintf("临床 Death 信息:      %s\n", if (has_death) "Yes" else "No"))
cat(sprintf("临床 Recurr 信息:     %s\n", if (has_recurr) "Yes" else "No"))
cat(sprintf("========================================\n"))
