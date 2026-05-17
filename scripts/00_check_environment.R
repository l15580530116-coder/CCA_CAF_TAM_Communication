# ============================================================
# 分析目的: 检查 R 环境 —— 版本号和核心分析包是否已安装
# 输入文件: 无
# 输出文件: results/R_environment_check.txt
# 主要方法: R base + requireNamespace 检查包可用性
# ============================================================

# 要检查的包列表
required_pkgs <- c(
  "tidyverse",
  "data.table",
  "survival",
  "survminer",
  "TCGAbiolinks",
  "GEOquery",
  "limma",
  "edgeR",
  "DESeq2",
  "clusterProfiler",
  "org.Hs.eg.db",
  "GSVA",
  "ConsensusClusterPlus",
  "NMF",
  "glmnet",
  "randomForestSRC",
  "timeROC",
  "ComplexHeatmap",
  "pheatmap",
  "ggpubr"
)

# 输出文件
output_file <- "results/R_environment_check.txt"

# 开始检查
cat("============================================\n", file = output_file)
cat("R Environment Check\n", file = output_file, append = TRUE)
cat("Date:", as.character(Sys.Date()), "\n", file = output_file, append = TRUE)
cat("============================================\n\n", file = output_file, append = TRUE)

# R 版本
cat("R Version:\n", file = output_file, append = TRUE)
cat(paste0("  ", R.version.string, "\n\n"), file = output_file, append = TRUE)

# 包检查
cat("Package Status:\n", file = output_file, append = TRUE)
cat("--------------------------------------------\n", file = output_file, append = TRUE)

max_name_len <- max(nchar(required_pkgs))

for (pkg in required_pkgs) {
  pkg_padded <- sprintf(paste0("%-", max_name_len + 1, "s"), pkg)
  if (requireNamespace(pkg, quietly = TRUE)) {
    ver <- as.character(packageVersion(pkg))
    status <- paste0("OK (v", ver, ")")
    cat(sprintf("  %s OK (v%s)\n", pkg_padded, ver), file = output_file, append = TRUE)
  } else {
    cat(sprintf("  %s NOT INSTALLED\n", pkg_padded), file = output_file, append = TRUE)
  }
}

cat("\n============================================\n", file = output_file, append = TRUE)

# 汇总
n_ok <- sum(sapply(required_pkgs, function(p) requireNamespace(p, quietly = TRUE)))
n_miss <- length(required_pkgs) - n_ok
cat(sprintf("Summary: %d OK, %d NOT INSTALLED, %d total\n", n_ok, n_miss, length(required_pkgs)),
    file = output_file, append = TRUE)
cat("============================================\n", file = output_file, append = TRUE)

# 同时输出到控制台
cat(readLines(output_file), sep = "\n")

if (n_miss > 0) {
  cat("\nWARNING: Some packages are missing. Install them before proceeding.\n")
} else {
  cat("\nAll required packages are installed. Ready to proceed.\n")
}
