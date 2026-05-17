# ============================================================
# 分析目的: GSVA/ssGSEA 通路评分 → tumor vs normal 比较 →
#          生存分析 + 相关性分析
# 输入文件:
#   data/processed/TCGA_CHOL/TCGA_CHOL_log2TPM_gene_symbol.csv
#   data/processed/TCGA_CHOL/TCGA_CHOL_sample_metadata.csv
#   data/processed/TCGA_CHOL/TCGA_CHOL_tumor_expression_for_survival.csv
#   data/processed/TCGA_CHOL/TCGA_CHOL_matched_survival_metadata.csv
#   data/processed/GSE107943/GSE107943_log2CPM_gene_symbol.csv
#   data/processed/GSE107943/GSE107943_sample_metadata.csv
#   data/clinical/GSE107943/GSE107943_clinical.csv
#   gene_sets/ (6 个基因集文件)
#   tables/enrichment/validated_IM_CAF_TAM_DEGs_{up,down}.csv
# 输出文件:
#   tables/gsva/ (8 个表格)
#   figures/gsva/ (9 张图)
# 主要方法: GSVA::gsva(ssgsea) → Wilcoxon → KM + Cox
# ============================================================

library(GSVA)
library(GSEABase)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(pheatmap)
library(data.table)
library(RColorBrewer)

dir.create("tables/gsva", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/gsva", showWarnings = FALSE, recursive = TRUE)

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

tumor_col  <- "#E41A1C"
normal_col <- "#377EB8"

# ============================================================
# Step 1: 加载数据 & 构建基因集
# ============================================================
log_msg("=== Step 1: 加载数据 & 构建基因集 ===")

# TCGA
tcga_expr <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_log2TPM_gene_symbol.csv",
                      row.names=1, check.names=FALSE)
tcga_meta <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_sample_metadata.csv",
                      stringsAsFactors=FALSE)
rownames(tcga_meta) <- tcga_meta$sample_id

# GSE107943
gse_expr <- read.csv("data/processed/GSE107943/GSE107943_log2CPM_gene_symbol.csv",
                     row.names=1, check.names=FALSE)
# Drop 6N (all NA)
gse_expr <- gse_expr[, colnames(gse_expr) != "6N", drop=FALSE]
gse_meta <- read.csv("data/processed/GSE107943/GSE107943_sample_metadata.csv",
                     stringsAsFactors=FALSE)
rownames(gse_meta) <- gse_meta$sample_id

# Align
ct <- intersect(colnames(tcga_expr), rownames(tcga_meta))
tcga_expr <- as.matrix(tcga_expr[, ct, drop=FALSE])
cg <- intersect(colnames(gse_expr), rownames(gse_meta))
gse_expr <- as.matrix(gse_expr[, cg, drop=FALSE])

# Survival data
tcga_surv <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_matched_survival_metadata.csv",
                      stringsAsFactors=FALSE)
tcga_surv_tumor <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_tumor_expression_for_survival.csv",
                            row.names=1, check.names=FALSE)

gse_clin <- read.csv("data/clinical/GSE107943/GSE107943_clinical.csv",
                     stringsAsFactors=FALSE)

# 基因集
read_gene_set <- function(path) {
  scan(path, what="character", skip=1, quiet=TRUE)
}

im_genes  <- read_gene_set("gene_sets/immunometabolism_genes.txt")
caf_genes <- read_gene_set("gene_sets/CAF_markers.txt")
tam_genes <- read_gene_set("gene_sets/macrophage_TAM_markers.txt")
combined_genes <- read_gene_set("gene_sets/CCA_IM_CAF_TAM_combined_genes.txt")
up_genes_df   <- read.csv("tables/enrichment/validated_IM_CAF_TAM_DEGs_up.csv", stringsAsFactors=FALSE)
down_genes_df <- read.csv("tables/enrichment/validated_IM_CAF_TAM_DEGs_down.csv", stringsAsFactors=FALSE)

gene_sets <- list(
  IM_total       = im_genes,
  CAF            = caf_genes,
  TAM            = tam_genes,
  ECM_up         = up_genes_df$gene_symbol,
  Metabolism_down = down_genes_df$gene_symbol,
  IM_CAF_TAM_combined = combined_genes
)

n_sets <- length(gene_sets)
log_msg(sprintf("Gene sets: %d", n_sets))
for (nm in names(gene_sets)) {
  log_msg(sprintf("  %-25s: %d genes", nm, length(gene_sets[[nm]])))
}

# ============================================================
# Step 2: ssGSEA 评分
# ============================================================
log_msg("\n=== Step 2: ssGSEA 评分 ===")

run_ssgsea <- function(expr_mat, gs_list, label) {
  # Filter genes present in expression matrix
  gs_filt <- lapply(gs_list, function(g) intersect(g, rownames(expr_mat)))
  gs_sizes <- sapply(gs_filt, length)
  log_msg(sprintf("  %s: gene set sizes: %s", label,
    paste(sprintf("%s=%d", names(gs_sizes), gs_sizes), collapse=", ")))

  gs_param <- gsvaParam(expr_mat, gs_filt)
  scores <- gsva(gs_param, verbose=FALSE)
  # scores is a matrix: gene_sets x samples
  return(scores)
}

use_method <- "ssGSEA"

tryCatch({
  tcga_scores <- run_ssgsea(tcga_expr, gene_sets, "TCGA")
  gse_scores  <- run_ssgsea(gse_expr, gene_sets, "GSE107943")
  log_msg("ssGSEA 评分完成")
}, error = function(e) {
  log_msg(sprintf("ssGSEA failed: %s. Using z-score mean fallback.", e$message))
  use_method <<- "z-score mean"

  zscore_mean <- function(expr_mat, gs_list) {
    expr_z <- t(scale(t(expr_mat)))
    expr_z[!is.finite(expr_z)] <- 0
    scores <- matrix(NA, length(gs_list), ncol(expr_mat))
    rownames(scores) <- names(gs_list)
    colnames(scores) <- colnames(expr_mat)
    for (i in seq_along(gs_list)) {
      g <- intersect(gs_list[[i]], rownames(expr_z))
      if (length(g) > 0) scores[i, ] <- colMeans(expr_z[g, , drop=FALSE], na.rm=TRUE)
    }
    return(scores)
  }

  tcga_scores <- zscore_mean(tcga_expr, gene_sets)
  gse_scores  <- zscore_mean(gse_expr, gene_sets)
  log_msg("z-score mean 评分完成")
})

log_msg(sprintf("TCGA scores: %d x %d", nrow(tcga_scores), ncol(tcga_scores)))
log_msg(sprintf("GSE scores:  %d x %d", nrow(gse_scores), ncol(gse_scores)))

# 保存完整评分
write.csv(data.frame(pathway=rownames(tcga_scores), tcga_scores, check.names=FALSE),
          "tables/gsva/TCGA_CHOL_pathway_scores_all_samples.csv", row.names=FALSE)
write.csv(data.frame(pathway=rownames(gse_scores), gse_scores, check.names=FALSE),
          "tables/gsva/GSE107943_pathway_scores_all_samples.csv", row.names=FALSE)
log_msg("已保存 pathway score 表格")

# ============================================================
# Step 3: 构建衍生评分
# ============================================================
log_msg("\n=== Step 3: 构建衍生评分 ===")

# ECM_up score
ecm_idx <- which(rownames(tcga_scores) == "ECM_up")
# Metabolism_down score
met_idx <- which(rownames(tcga_scores) == "Metabolism_down")
# CAF score
caf_idx <- which(rownames(tcga_scores) == "CAF")
# TAM score
tam_idx <- which(rownames(tcga_scores) == "TAM")

# Aggressive microenvironment score
add_aggr_score <- function(scores) {
  aggr <- scale(scores[ecm_idx, ]) + scale(scores[caf_idx, ]) +
          scale(scores[tam_idx, ]) - scale(scores[met_idx, ])
  rbind(scores, CCA_aggressive_microenvironment = aggr[, 1])
}

tcga_scores <- add_aggr_score(tcga_scores)
gse_scores  <- add_aggr_score(gse_scores)
log_msg("已添加 CCA_aggressive_microenvironment 评分")

# ============================================================
# Step 4: Tumor vs Normal 比较
# ============================================================
log_msg("\n=== Step 4: Tumor vs Normal 比较 ===")

compare_tn <- function(scores, meta, label) {
  common <- intersect(colnames(scores), rownames(meta))
  sc <- scores[, common, drop=FALSE]
  mt <- meta[common, , drop=FALSE]

  results <- data.frame(pathway=character(), pvalue=numeric(), FDR=numeric(),
                        mean_Tumor=numeric(), mean_Normal=numeric(),
                        direction=character(), stringsAsFactors=FALSE)

  for (i in seq_len(nrow(sc))) {
    pathway <- rownames(sc)[i]
    vals <- sc[i, ]
    t_vals <- vals[mt$sample_type %in% c("Tumor", "Primary Tumor")]
    n_vals <- vals[mt$sample_type %in% c("Normal", "Solid Tissue Normal")]

    if (length(t_vals) > 0 && length(n_vals) > 0) {
      wt <- wilcox.test(t_vals, n_vals)
      results[i, ] <- c(pathway, wt$p.value, NA,
                        mean(t_vals), mean(n_vals),
                        if (mean(t_vals) > mean(n_vals)) "up_in_tumor" else "down_in_tumor")
    }
  }
  results$pvalue <- as.numeric(results$pvalue)
  results$FDR <- p.adjust(results$pvalue, method="BH")
  results$mean_Tumor <- as.numeric(results$mean_Tumor)
  results$mean_Normal <- as.numeric(results$mean_Normal)

  log_msg(sprintf("%s T vs N: %d pathways, %d FDR<0.05",
    label, nrow(results), sum(results$FDR < 0.05)))
  return(results)
}

tn_tcga <- compare_tn(tcga_scores, tcga_meta, "TCGA")
tn_gse  <- compare_tn(gse_scores, gse_meta, "GSE107943")

# 合并比较
tn_merged <- merge(tn_tcga, tn_gse, by="pathway", suffixes=c("_TCGA", "_GSE"), all=TRUE)
write.csv(tn_merged, "tables/gsva/pathway_score_tumor_normal_comparison.csv", row.names=FALSE)
log_msg("已保存 T vs N 比较表")

# ============================================================
# Step 5: 生存分析
# ============================================================
log_msg("\n=== Step 5: 生存分析 ===")

# Prepare TCGA survival data
tcga_tum_expr <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_tumor_expression_for_survival.csv",
                          row.names=1, check.names=FALSE)
tcga_tum_scores <- tcga_scores[, colnames(tcga_tum_expr), drop=FALSE]
tcga_surv_aligned <- tcga_surv[match(colnames(tcga_tum_scores), tcga_surv$patient_id), ]
# Need to match by sample_id (barcode), not patient_id
# Use the sample metadata to map
tcga_tum_meta <- tcga_meta[colnames(tcga_tum_scores), , drop=FALSE]
tcga_tum_surv <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_matched_survival_metadata.csv",
                          stringsAsFactors=FALSE)

# Build survival data frame with scores
tcga_surv_df <- tcga_tum_surv
rownames(tcga_surv_df) <- tcga_surv_df$patient_id

# Map scores via sample barcode → patient_id
tcga_score_patient <- data.frame(
  sample_id = colnames(tcga_tum_scores),
  patient_id = tcga_tum_meta$patient_id,
  stringsAsFactors = FALSE
)
# Take first sample per patient (if duplicate)
tcga_score_patient <- tcga_score_patient[!duplicated(tcga_score_patient$patient_id), ]
common_pat <- intersect(tcga_score_patient$patient_id, tcga_surv_df$patient_id)

tcga_surv_df <- tcga_surv_df[common_pat, , drop=FALSE]
tcga_score_mat <- tcga_tum_scores[, tcga_score_patient$sample_id[match(common_pat, tcga_score_patient$patient_id)], drop=FALSE]
colnames(tcga_score_mat) <- common_pat

# Add scores to survival
for (i in seq_len(nrow(tcga_score_mat))) {
  nm <- rownames(tcga_score_mat)[i]
  tcga_surv_df[[nm]] <- tcga_score_mat[i, common_pat]
}

log_msg(sprintf("TCGA survival: %d patients", nrow(tcga_surv_df)))

# GSE107943 survival
gse_meta_aligned <- gse_meta[colnames(gse_scores), , drop=FALSE]
gse_tum_scores <- gse_scores[, gse_meta_aligned$sample_type == "Tumor", drop=FALSE]
gse_tum_clin <- gse_clin[gse_clin$tissue == "Tumor (Intrahepatic cholangiocarcinoma)", , drop=FALSE]

# Match by sample_name
common_gse_samp <- intersect(colnames(gse_tum_scores), gse_tum_clin$sample_name)
gse_surv_df <- gse_tum_clin[gse_tum_clin$sample_name %in% common_gse_samp, , drop=FALSE]
rownames(gse_surv_df) <- gse_surv_df$sample_name
gse_surv_df <- gse_surv_df[common_gse_samp, , drop=FALSE]

# Convert survival data
gse_surv_df$OS_time <- as.numeric(gse_surv_df$survival_mo_) * 30.4375  # months to days
gse_surv_df$OS_status <- as.numeric(gse_surv_df$death)
gse_surv_df <- gse_surv_df[!is.na(gse_surv_df$OS_time) & gse_surv_df$OS_time > 0, , drop=FALSE]

# Add scores
common_gse_final <- intersect(common_gse_samp, rownames(gse_surv_df))
gse_surv_df <- gse_surv_df[common_gse_final, , drop=FALSE]
for (i in seq_len(nrow(gse_tum_scores))) {
  nm <- rownames(gse_tum_scores)[i]
  gse_surv_df[[nm]] <- gse_tum_scores[i, common_gse_final]
}

log_msg(sprintf("GSE107943 survival: %d patients", nrow(gse_surv_df)))

# ---- Univariate Cox ----
run_cox <- function(surv_df, score_names, label) {
  results <- data.frame(pathway=character(), HR=numeric(), lower95=numeric(),
                        upper95=numeric(), pvalue=numeric(), stringsAsFactors=FALSE)
  for (i in seq_along(score_names)) {
    nm <- score_names[i]
    if (!nm %in% colnames(surv_df)) next
    cox <- coxph(Surv(OS_time, OS_status) ~ surv_df[[nm]], data=surv_df)
    s <- summary(cox)
    results[i, ] <- c(nm, s$conf.int[1], s$conf.int[3], s$conf.int[4], s$coefficients[5])
  }
  results$HR <- as.numeric(results$HR)
  results$lower95 <- as.numeric(results$lower95)
  results$upper95 <- as.numeric(results$upper95)
  results$pvalue <- as.numeric(results$pvalue)
  results$FDR <- p.adjust(results$pvalue, method="BH")
  results$cohort <- label
  log_msg(sprintf("%s Cox: %d tests, %d FDR<0.05", label, nrow(results), sum(results$FDR < 0.05)))
  return(results)
}

score_names <- rownames(tcga_scores)
cox_tcga <- run_cox(tcga_surv_df, score_names, "TCGA-CHOL")
cox_gse  <- run_cox(gse_surv_df, score_names, "GSE107943")

cox_all <- rbind(cox_tcga, cox_gse)
write.csv(cox_all, "tables/gsva/pathway_score_survival_univariate_cox.csv", row.names=FALSE)
log_msg("已保存 Cox 结果")

# ---- Save tumor scores with survival ----
tcga_surv_out <- tcga_surv_df
write.csv(tcga_surv_out, "tables/gsva/TCGA_CHOL_tumor_scores_with_survival.csv", row.names=FALSE)
gse_surv_out <- gse_surv_df
write.csv(gse_surv_out, "tables/gsva/GSE107943_tumor_scores_with_survival.csv", row.names=FALSE)
log_msg("已保存 tumor scores + survival")

# ============================================================
# Step 6: 相关性分析
# ============================================================
log_msg("\n=== Step 6: 相关性分析 ===")

tcga_tum_sc <- tcga_scores[, colnames(tcga_tum_scores), drop=FALSE]
cor_mat <- cor(t(tcga_tum_sc), method="spearman")
write.csv(data.frame(pathway=rownames(cor_mat), cor_mat, check.names=FALSE),
          "tables/gsva/pathway_score_correlation_matrix.csv", row.names=FALSE)
log_msg("已保存 correlation matrix")

# ============================================================
# Step 7: 可视化
# ============================================================
log_msg("\n=== Step 7: 可视化 ===")

# ---- Fig4A: TCGA boxplot T vs N ----
log_msg("Fig4A: TCGA boxplot")
tcga_tn_plot_data <- data.frame(
  sample = colnames(tcga_scores),
  sample_type = tcga_meta[colnames(tcga_scores), "sample_type"],
  t(as.data.frame(tcga_scores))
)

pdf("figures/gsva/Fig4A_TCGA_score_boxplot_tumor_normal.pdf", width=14, height=10)
par(mfrow=c(3,3), mar=c(5,4,3,1))
for (i in seq_len(nrow(tcga_scores))) {
  nm <- rownames(tcga_scores)[i]
  t_vals <- tcga_scores[i, tcga_meta[colnames(tcga_scores),"sample_type"] == "Primary Tumor"]
  n_vals <- tcga_scores[i, tcga_meta[colnames(tcga_scores),"sample_type"] == "Solid Tissue Normal"]
  boxplot(list(Tumor=t_vals, Normal=n_vals), col=c(tumor_col, normal_col),
          main=nm, ylab="ssGSEA score", las=1)
  pv <- tn_tcga$pvalue[tn_tcga$pathway == nm]
  if (length(pv) > 0) {
    legend("topright", sprintf("p=%.2e", pv), bty="n", cex=0.8)
  }
}
dev.off()
log_msg("  Done")

# ---- Fig4B: GSE107943 boxplot T vs N ----
log_msg("Fig4B: GSE boxplot")
pdf("figures/gsva/Fig4B_GSE107943_score_boxplot_tumor_normal.pdf", width=14, height=10)
par(mfrow=c(3,3), mar=c(5,4,3,1))
for (i in seq_len(nrow(gse_scores))) {
  nm <- rownames(gse_scores)[i]
  t_vals <- gse_scores[i, gse_meta[colnames(gse_scores),"sample_type"] == "Tumor"]
  n_vals <- gse_scores[i, gse_meta[colnames(gse_scores),"sample_type"] == "Normal"]
  boxplot(list(Tumor=t_vals, Normal=n_vals), col=c(tumor_col, normal_col),
          main=nm, ylab="ssGSEA score", las=1)
  pv <- tn_gse$pvalue[tn_gse$pathway == nm]
  if (length(pv) > 0) {
    legend("topright", sprintf("p=%.2e", pv), bty="n", cex=0.8)
  }
}
dev.off()
log_msg("  Done")

# ---- Fig4C: TCGA score heatmap ----
log_msg("Fig4C: TCGA score heatmap")
tcga_z <- t(scale(t(tcga_scores)))
tcga_z[tcga_z > 3] <- 3; tcga_z[tcga_z < -3] <- -3
ann_col <- data.frame(Type=ifelse(tcga_meta[colnames(tcga_z),"sample_type"] == "Primary Tumor", "Tumor", "Normal"),
                      row.names=colnames(tcga_z))
ann_colors <- list(Type=c(Tumor=tumor_col, Normal=normal_col))
pdf("figures/gsva/Fig4C_TCGA_score_heatmap.pdf", width=12, height=6)
pheatmap(tcga_z, annotation_col=ann_col, annotation_colors=ann_colors,
         cluster_rows=TRUE, cluster_cols=TRUE, show_colnames=FALSE,
         main="TCGA-CHOL Pathway Scores (Z-score)", fontsize_row=10)
dev.off()
log_msg("  Done")

# ---- Fig4D: GSE score heatmap ----
log_msg("Fig4D: GSE score heatmap")
gse_z <- t(scale(t(gse_scores)))
gse_z[gse_z > 3] <- 3; gse_z[gse_z < -3] <- -3
ann_col_g <- data.frame(Type=gse_meta[colnames(gse_z),"sample_type"], row.names=colnames(gse_z))
ann_colors_g <- list(Type=c(Tumor=tumor_col, Normal=normal_col))
pdf("figures/gsva/Fig4D_GSE107943_score_heatmap.pdf", width=12, height=6)
pheatmap(gse_z, annotation_col=ann_col_g, annotation_colors=ann_colors_g,
         cluster_rows=TRUE, cluster_cols=TRUE, show_colnames=FALSE,
         main="GSE107943 Pathway Scores (Z-score)", fontsize_row=10)
dev.off()
log_msg("  Done")

# ---- KM plots ----
log_msg("KM plots...")

run_km <- function(surv_df, score_name, title, filename) {
  x <- surv_df[[score_name]]
  cutoff <- median(x, na.rm=TRUE)
  surv_df$group <- ifelse(x > cutoff, "High", "Low")

  fit <- survfit(Surv(OS_time, OS_status) ~ group, data=surv_df)
  # Use custom ggplot for more control
  surv_df_km <- data.frame(
    time = surv_df$OS_time / 30.4375,  # days to months
    status = surv_df$OS_status,
    group = surv_df$group
  )
  fit2 <- survfit(Surv(time, status) ~ group, data=surv_df_km)
  sdiff <- survdiff(Surv(time, status) ~ group, data=surv_df_km)
  pval <- 1 - pchisq(sdiff$chisq, 1)

  cox <- coxph(Surv(OS_time, OS_status) ~ x, data=surv_df)
  hr <- exp(coef(cox))

  pdf(filename, width=7, height=6)
  p <- ggsurvplot(fit2, data=surv_df_km, pval=TRUE,
                  palette=c("High"=tumor_col, "Low"=normal_col),
                  legend.title=score_name, legend.labs=c("High","Low"),
                  xlab="Overall Survival (months)", title=title)
  print(p)
  dev.off()
  cat(sprintf("    %s: p=%.4f, HR=%.2f\n", score_name, pval, hr))
}

# Fig4E: TCGA CAF KM
if ("CAF" %in% colnames(tcga_surv_df)) {
  run_km(tcga_surv_df, "CAF", "TCGA-CHOL: CAF Score",
         "figures/gsva/Fig4E_TCGA_KM_high_low_CAF_score.pdf")
}
# Fig4F: TCGA IM_CAF_TAM combined KM
if ("IM_CAF_TAM_combined" %in% colnames(tcga_surv_df)) {
  run_km(tcga_surv_df, "IM_CAF_TAM_combined",
         "TCGA-CHOL: IM+CAF+TAM Combined Score",
         "figures/gsva/Fig4F_TCGA_KM_high_low_IM_CAF_TAM_score.pdf")
}
# Fig4G: GSE CAF KM
if ("CAF" %in% colnames(gse_surv_df)) {
  run_km(gse_surv_df, "CAF", "GSE107943: CAF Score",
         "figures/gsva/Fig4G_GSE107943_KM_high_low_CAF_score.pdf")
}
# Fig4H: GSE IM_CAF_TAM combined KM
if ("IM_CAF_TAM_combined" %in% colnames(gse_surv_df)) {
  run_km(gse_surv_df, "IM_CAF_TAM_combined",
         "GSE107943: IM+CAF+TAM Combined Score",
         "figures/gsva/Fig4H_GSE107943_KM_high_low_IM_CAF_TAM_score.pdf")
}

# ---- Fig4I: Correlation heatmap ----
log_msg("Fig4I: Correlation heatmap")
pdf("figures/gsva/Fig4I_score_correlation_heatmap.pdf", width=9, height=8)
pheatmap(cor_mat, display_numbers=TRUE, number_format="%.2f",
         color=colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
         breaks=seq(-1, 1, length.out=101),
         main="Spearman Correlation: Pathway Scores (TCGA Tumors)")
dev.off()
log_msg("  Done")

# ============================================================
# Step 8: 汇总表
# ============================================================
log_msg("\n=== Step 8: 汇总表 ===")

# Check direction consistency for aggressive score
tn_tcga_aggr <- tn_tcga[tn_tcga$pathway == "CCA_aggressive_microenvironment", ]
tn_gse_aggr  <- tn_gse[tn_gse$pathway == "CCA_aggressive_microenvironment", ]
aggr_dir_consistent <- if (nrow(tn_tcga_aggr) > 0 && nrow(tn_gse_aggr) > 0) {
  tn_tcga_aggr$direction == tn_gse_aggr$direction
} else NA

summary_rows <- data.frame(
  item=c(
    "评分方法", "总评分通路数",
    "TCGA TvsN 显著数 (FDR<0.05)", "TCGA TvsN 升高", "TCGA TvsN 降低",
    "GSE TvsN 显著数 (FDR<0.05)", "GSE TvsN 升高", "GSE TvsN 降低",
    "TCGA Cox 显著 (FDR<0.05)", "GSE Cox 显著 (FDR<0.05)",
    "Aggressive score 方向一致"
  ),
  value=c(
    use_method, as.character(nrow(tcga_scores)),
    as.character(sum(tn_tcga$FDR < 0.05)),
    paste(tn_tcga$pathway[tn_tcga$direction=="up_in_tumor" & tn_tcga$FDR<0.05], collapse="; "),
    paste(tn_tcga$pathway[tn_tcga$direction=="down_in_tumor" & tn_tcga$FDR<0.05], collapse="; "),
    as.character(sum(tn_gse$FDR < 0.05)),
    paste(tn_gse$pathway[tn_gse$direction=="up_in_tumor" & tn_gse$FDR<0.05], collapse="; "),
    paste(tn_gse$pathway[tn_gse$direction=="down_in_tumor" & tn_gse$FDR<0.05], collapse="; "),
    paste(cox_tcga$pathway[cox_tcga$FDR < 0.05], collapse="; "),
    paste(cox_gse$pathway[cox_gse$FDR < 0.05], collapse="; "),
    as.character(aggr_dir_consistent)
  ),
  stringsAsFactors=FALSE
)

write.csv(summary_rows, "tables/gsva/GSVA_analysis_summary.csv", row.names=FALSE)

# ============================================================
# 最终报告
# ============================================================
cat(sprintf("\n========================================\n"))
cat(sprintf("GSVA/ssGSEA 评分与生存分析完成\n"))
cat(sprintf("========================================\n"))
cat(sprintf("评分方法: %s\n", use_method))
cat(sprintf("TCGA TvsN 显著: %d (up=%d, down=%d)\n",
  sum(tn_tcga$FDR < 0.05),
  sum(tn_tcga$direction=="up_in_tumor" & tn_tcga$FDR < 0.05),
  sum(tn_tcga$direction=="down_in_tumor" & tn_tcga$FDR < 0.05)))
cat(sprintf("GSE TvsN 显著:  %d (up=%d, down=%d)\n",
  sum(tn_gse$FDR < 0.05),
  sum(tn_gse$direction=="up_in_tumor" & tn_gse$FDR < 0.05),
  sum(tn_gse$direction=="down_in_tumor" & tn_gse$FDR < 0.05)))
cat(sprintf("TCGA Cox 显著: %s\n",
  paste(cox_tcga$pathway[cox_tcga$FDR < 0.05], collapse=", ")))
cat(sprintf("GSE Cox 显著:  %s\n",
  paste(cox_gse$pathway[cox_gse$FDR < 0.05], collapse=", ")))
cat(sprintf("Aggressive score 方向一致: %s\n", aggr_dir_consistent))
cat(sprintf("========================================\n"))
