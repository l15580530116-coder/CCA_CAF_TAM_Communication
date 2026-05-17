# ============================================================
# 分析目的: 免疫浸润 + 免疫检查点 + aggressive score 关联分析
# 输入文件:
#   data/processed/TCGA_CHOL/TCGA_CHOL_log2TPM_gene_symbol.csv
#   data/processed/TCGA_CHOL/TCGA_CHOL_tumor_expression_for_survival.csv
#   data/processed/TCGA_CHOL/TCGA_CHOL_sample_metadata.csv
#   data/processed/TCGA_CHOL/TCGA_CHOL_matched_survival_metadata.csv
#   data/processed/GSE107943/GSE107943_log2CPM_gene_symbol.csv
#   data/processed/GSE107943/GSE107943_sample_metadata.csv
#   data/clinical/GSE107943/GSE107943_clinical.csv
#   tables/gsva/TCGA_CHOL_tumor_scores_with_survival.csv
#   tables/gsva/GSE107943_tumor_scores_with_survival.csv
# 输出文件:
#   tables/immune/ (11 个表格)
#   figures/immune/ (10 张图)
# 主要方法: marker-based ssGSEA immune scores + Spearman correlation +
#          survival + checkpoint analysis
# 注: immunedeconv 可用但使用内置 marker ssGSEA 以获得一致的全谱系覆盖
# ============================================================

library(GSVA)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(pheatmap)
library(RColorBrewer)
library(data.table)

dir.create("tables/immune", showWarnings=FALSE, recursive=TRUE)
dir.create("figures/immune", showWarnings=FALSE, recursive=TRUE)

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

use_immunedeconv <- requireNamespace("immunedeconv", quietly=TRUE)
log_msg(sprintf("immunedeconv available: %s", use_immunedeconv))

# ============================================================
# Step 1: 加载数据
# ============================================================
log_msg("=== Step 1: 加载数据 ===")

tcga_expr <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_log2TPM_gene_symbol.csv",
                      row.names=1, check.names=FALSE)
tcga_tum <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_tumor_expression_for_survival.csv",
                     row.names=1, check.names=FALSE)
tcga_meta <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_sample_metadata.csv", stringsAsFactors=FALSE)
rownames(tcga_meta) <- tcga_meta$sample_id
tcga_surv <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_matched_survival_metadata.csv", stringsAsFactors=FALSE)

gse_expr <- read.csv("data/processed/GSE107943/GSE107943_log2CPM_gene_symbol.csv",
                     row.names=1, check.names=FALSE)
gse_expr <- gse_expr[, colnames(gse_expr)!="6N", drop=FALSE]
gse_meta <- read.csv("data/processed/GSE107943/GSE107943_sample_metadata.csv", stringsAsFactors=FALSE)
rownames(gse_meta) <- gse_meta$sample_id
gse_clin <- read.csv("data/clinical/GSE107943/GSE107943_clinical.csv", stringsAsFactors=FALSE)

# Align
ct <- intersect(colnames(tcga_expr), rownames(tcga_meta))
tcga_expr <- as.matrix(tcga_expr[, ct, drop=FALSE])
cg <- intersect(colnames(gse_expr), rownames(gse_meta))
gse_expr <- as.matrix(gse_expr[, cg, drop=FALSE])

# Aggressive scores
tcga_gsva <- read.csv("tables/gsva/TCGA_CHOL_tumor_scores_with_survival.csv", stringsAsFactors=FALSE)
gse_gsva  <- read.csv("tables/gsva/GSE107943_tumor_scores_with_survival.csv", stringsAsFactors=FALSE)

# ============================================================
# Step 2: 定义 marker gene sets
# ============================================================
log_msg("\n=== Step 2: 定义 immune marker gene sets ===")

immune_markers <- list(
  Macrophage_TAM = c("CD68","CD163","MRC1","CSF1R","LYZ","C1QA","C1QB","C1QC","APOE","TREM2","SPP1","MARCO","MSR1"),
  M1_macrophage  = c("IL1B","TNF","CXCL9","CXCL10","CXCL11","CD86","NOS2","IRF5"),
  M2_macrophage  = c("CD163","MRC1","MSR1","IL10","TGFB1","CCL18","MAF","VSIG4"),
  CAF            = c("ACTA2","FAP","PDGFRB","PDGFRA","COL1A1","COL1A2","COL3A1","COL6A1","DCN","LUM","POSTN","MMP11","CXCL12"),
  CD8_T_cell     = c("CD8A","CD8B","GZMA","GZMB","PRF1","NKG7","GNLY"),
  CD4_T_cell     = c("CD4","IL7R","CCR7","LTB","TCF7"),
  Treg           = c("FOXP3","IL2RA","CTLA4","TIGIT","IKZF2","CCR8"),
  NK_cell        = c("NKG7","GNLY","KLRD1","KLRK1","FCGR3A","PRF1"),
  B_cell         = c("MS4A1","CD79A","CD79B","CD19","CD22","BANK1"),
  Dendritic_cell = c("ITGAX","CD1C","CLEC9A","LILRA4","FCER1A","CST3"),
  Neutrophil     = c("S100A8","S100A9","FCGR3B","CXCR2","CSF3R","MPO"),
  Endothelial_cell = c("PECAM1","VWF","KDR","ENG","CDH5","EMCN"),
  Exhausted_T_cell = c("PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","TOX","CXCL13")
)

# Filter: keep only cell types with >=3 markers in expression data
immune_filtered <- list()
for (nm in names(immune_markers)) {
  g <- intersect(immune_markers[[nm]], rownames(tcga_expr))
  g <- intersect(g, rownames(gse_expr))
  if (length(g) >= 3) {
    immune_filtered[[nm]] <- g
    log_msg(sprintf("  %-20s: %d markers", nm, length(g)))
  } else {
    log_msg(sprintf("  %-20s: SKIPPED (only %d markers)", nm, length(g)))
  }
}

# Checkpoint genes
checkpoint_genes <- c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT",
                      "CD80","CD86","CD276","VSIR","IDO1","ICOS","TNFRSF9","TNFRSF4")
checkpoint_avail <- intersect(checkpoint_genes, intersect(rownames(tcga_expr), rownames(gse_expr)))
log_msg(sprintf("Checkpoint genes available: %d/%d", length(checkpoint_avail), length(checkpoint_genes)))

# ============================================================
# Step 3: ssGSEA immune scores
# ============================================================
log_msg("\n=== Step 3: ssGSEA immune scores ===")

run_immune_ssgsea <- function(expr_mat, marker_list, label) {
  gs_param <- gsvaParam(expr_mat, marker_list)
  scores <- gsva(gs_param, verbose=FALSE)
  log_msg(sprintf("  %s: %d cell types x %d samples", label, nrow(scores), ncol(scores)))
  return(scores)
}

tcga_immune <- run_immune_ssgsea(tcga_expr, immune_filtered, "TCGA")
gse_immune  <- run_immune_ssgsea(gse_expr, immune_filtered, "GSE107943")

# Save
write.csv(data.frame(cell_type=rownames(tcga_immune), tcga_immune, check.names=FALSE),
          "tables/immune/TCGA_immune_scores.csv", row.names=FALSE)
write.csv(data.frame(cell_type=rownames(gse_immune), gse_immune, check.names=FALSE),
          "tables/immune/GSE107943_immune_scores.csv", row.names=FALSE)
log_msg("Immune scores saved")

# ============================================================
# Step 4: Tumor vs Normal comparison
# ============================================================
log_msg("\n=== Step 4: Tumor vs Normal ===")

compare_immune_tn <- function(scores, meta, label) {
  common <- intersect(colnames(scores), rownames(meta))
  sc <- scores[, common, drop=FALSE]
  mt <- meta[common, , drop=FALSE]
  res <- data.frame(stringsAsFactors=FALSE)
  for (i in seq_len(nrow(sc))) {
    ct <- rownames(sc)[i]
    t_vals <- sc[i, mt$sample_type %in% c("Tumor","Primary Tumor")]
    n_vals <- sc[i, mt$sample_type %in% c("Normal","Solid Tissue Normal")]
    if (length(t_vals)>0 && length(n_vals)>0) {
      wt <- wilcox.test(t_vals, n_vals)
      res <- rbind(res, data.frame(cell_type=ct, mean_Tumor=mean(t_vals),
        mean_Normal=mean(n_vals), log2FC=mean(t_vals)-mean(n_vals),
        pvalue=wt$p.value, stringsAsFactors=FALSE))
    }
  }
  res$FDR <- p.adjust(res$pvalue, method="BH")
  res$direction <- ifelse(res$log2FC>0, "up_in_tumor", "down_in_tumor")
  log_msg(sprintf("  %s: %d/%d FDR<0.05", label, sum(res$FDR<0.05), nrow(res)))
  return(res)
}

imm_tn_tcga <- compare_immune_tn(tcga_immune, tcga_meta, "TCGA")
imm_tn_gse  <- compare_immune_tn(gse_immune, gse_meta, "GSE107943")

write.csv(imm_tn_tcga, "tables/immune/TCGA_immune_score_tumor_normal_comparison.csv", row.names=FALSE)
write.csv(imm_tn_gse,  "tables/immune/GSE107943_immune_score_tumor_normal_comparison.csv", row.names=FALSE)

# ============================================================
# Step 5: Aggressive score vs immune correlation (tumor only)
# ============================================================
log_msg("\n=== Step 5: Aggressive score vs immune correlation ===")

tumor_samples_tcga <- colnames(tcga_tum)
tumor_samples_gse <- colnames(gse_expr)[gse_meta[colnames(gse_expr),"sample_type"]=="Tumor"]

tcga_imm_tum <- tcga_immune[, intersect(tumor_samples_tcga, colnames(tcga_immune)), drop=FALSE]
gse_imm_tum  <- gse_immune[, intersect(tumor_samples_gse, colnames(gse_immune)), drop=FALSE]

# Get aggressive scores for tumor samples
tcga_gsva_rn <- tcga_gsva
if ("patient_id" %in% colnames(tcga_gsva_rn)) rownames(tcga_gsva_rn) <- tcga_gsva_rn$patient_id

# Aggressive score from ssGSEA pathway scores
tcga_path <- read.csv("tables/gsva/TCGA_CHOL_pathway_scores_all_samples.csv", row.names=1, check.names=FALSE)
gse_path  <- read.csv("tables/gsva/GSE107943_pathway_scores_all_samples.csv", row.names=1, check.names=FALSE)

score_vars <- c("CAF","ECM_up","TAM","Metabolism_down","IM_total","IM_CAF_TAM_combined","CCA_aggressive_microenvironment")

# Compute aggressive score if not in CSV
compute_aggr <- function(path_mat) {
  path_mat <- as.matrix(path_mat)
  # Ensure numeric
  mode(path_mat) <- "numeric"
  ecm <- path_mat["ECM_up", ]; caf <- path_mat["CAF", ]
  tam <- path_mat["TAM", ]; met <- path_mat["Metabolism_down", ]
  aggr <- scale(as.numeric(ecm)) + scale(as.numeric(caf)) +
          scale(as.numeric(tam)) - scale(as.numeric(met))
  return(aggr[,1])
}

# Get aggr score for tumor samples
tcga_aggr <- compute_aggr(tcga_path[, tumor_samples_tcga, drop=FALSE])
names(tcga_aggr) <- tumor_samples_tcga
gse_aggr <- compute_aggr(gse_path[, tumor_samples_gse, drop=FALSE])
names(gse_aggr) <- tumor_samples_gse

# Add the computed aggr to score matrices
tcga_score_vec <- compute_aggr(tcga_path)
tcga_all_scores <- rbind(tcga_path[intersect(score_vars, rownames(tcga_path)), , drop=FALSE],
                         CCA_aggressive_microenvironment = tcga_score_vec)
tcga_all_scores <- tcga_all_scores[, tumor_samples_tcga, drop=FALSE]

gse_score_vec <- compute_aggr(gse_path)
gse_all_scores <- rbind(gse_path[intersect(score_vars, rownames(gse_path)), , drop=FALSE],
                        CCA_aggressive_microenvironment = gse_score_vec)
gse_all_scores <- gse_all_scores[, tumor_samples_gse, drop=FALSE]

# Also get individual pathway scores
score_vars <- c("CAF","ECM_up","TAM","Metabolism_down","IM_total","IM_CAF_TAM_combined","CCA_aggressive_microenvironment")
# (score matrices now built above with aggr score included)

# Spearman correlation
run_spearman <- function(immune_mat, score_vec, label) {
  common_s <- intersect(colnames(immune_mat), names(score_vec))
  if (length(common_s) < 5) return(NULL)
  res <- data.frame(stringsAsFactors=FALSE)
  for (i in seq_len(nrow(immune_mat))) {
    ct <- rownames(immune_mat)[i]
    ct_vals <- as.numeric(immune_mat[i, common_s])
    sv_vals <- as.numeric(score_vec[common_s])
    valid <- is.finite(ct_vals) & is.finite(sv_vals)
    if (sum(valid) >= 5) {
      sp <- cor.test(ct_vals[valid], sv_vals[valid], method="spearman")
      res <- rbind(res, data.frame(cell_type=ct, rho=sp$estimate, pvalue=sp$p.value,
                                   stringsAsFactors=FALSE))
    }
  }
  if (nrow(res) > 0) res$FDR <- p.adjust(res$pvalue, method="BH")
  return(res)
}

tcga_aggr_cor <- run_spearman(tcga_imm_tum, tcga_aggr, "TCGA_aggr")
gse_aggr_cor  <- run_spearman(gse_imm_tum, gse_aggr, "GSE_aggr")

# Full correlation matrix
full_cor <- function(score_mat, immune_mat, label) {
  common <- intersect(colnames(score_mat), colnames(immune_mat))
  sc <- score_mat[, common, drop=FALSE]
  im <- immune_mat[, common, drop=FALSE]
  cor_mat <- cor(t(im), t(sc), method="spearman")
  write.csv(data.frame(cell_type=rownames(cor_mat), cor_mat, check.names=FALSE),
            paste0("tables/immune/", label, "_score_immune_correlation.csv"), row.names=FALSE)
  return(cor_mat)
}

tcga_cor_mat <- full_cor(tcga_all_scores, tcga_imm_tum, "TCGA")
gse_cor_mat  <- full_cor(gse_all_scores, gse_imm_tum, "GSE107943")

# Also single-score results
write.csv(tcga_aggr_cor, "tables/immune/TCGA_aggr_vs_immune_correlation.csv", row.names=FALSE)
write.csv(gse_aggr_cor,  "tables/immune/GSE107943_aggr_vs_immune_correlation.csv", row.names=FALSE)

# ============================================================
# Step 6: Checkpoint correlation
# ============================================================
log_msg("\n=== Step 6: Checkpoint correlation ===")

checkpoint_cor <- function(expr_mat, score_vec, cp_genes, label) {
  common <- intersect(colnames(expr_mat), names(score_vec))
  res <- data.frame(stringsAsFactors=FALSE)
  for (cg in cp_genes) {
    if (!cg %in% rownames(expr_mat)) next
    cp_vals <- as.numeric(expr_mat[cg, common])
    sv_vals <- as.numeric(score_vec[common])
    valid <- is.finite(cp_vals) & is.finite(sv_vals)
    if (sum(valid) >= 5) {
      sp <- cor.test(cp_vals[valid], sv_vals[valid], method="spearman")
      res <- rbind(res, data.frame(checkpoint=cg, rho=sp$estimate, pvalue=sp$p.value,
                                   stringsAsFactors=FALSE))
    }
  }
  if (nrow(res) > 0) res$FDR <- p.adjust(res$pvalue, method="BH")
  write.csv(res, paste0("tables/immune/", label, "_checkpoint_expression_correlation.csv"), row.names=FALSE)
  return(res)
}

tcga_cp_cor <- checkpoint_cor(tcga_expr, tcga_aggr, checkpoint_avail, "TCGA")
gse_cp_cor  <- checkpoint_cor(gse_expr, gse_aggr, checkpoint_avail, "GSE107943")

# ============================================================
# Step 7: Survival exploration
# ============================================================
log_msg("\n=== Step 7: Survival exploration ===")

immune_surv_cox <- function(immune_mat, surv_df, meta_data, label) {
  # Map samples to patient survival
  common_s <- intersect(colnames(immune_mat), meta_data$sample_id)
  im <- immune_mat[, common_s, drop=FALSE]
  # Get patient IDs
  pat_map <- setNames(meta_data$patient_id, meta_data$sample_id)
  pat_ids <- pat_map[common_s]
  common_p <- intersect(pat_ids, surv_df$patient_id)
  im_pat <- im[, common_s[pat_ids %in% common_p], drop=FALSE]
  colnames(im_pat) <- pat_ids[pat_ids %in% common_p]

  res <- data.frame(stringsAsFactors=FALSE)
  for (i in seq_len(nrow(im_pat))) {
    ct <- rownames(im_pat)[i]
    vals <- as.numeric(im_pat[i, ])
    names(vals) <- colnames(im_pat)
    surv_sub <- surv_df[surv_df$patient_id %in% names(vals), ]
    surv_vals <- vals[surv_sub$patient_id]
    valid <- is.finite(surv_vals) & !is.na(surv_sub$OS_time) & surv_sub$OS_time > 0
    if (sum(valid) >= 10) {
      cox <- tryCatch(
        coxph(Surv(OS_time, OS_status) ~ surv_vals, data=surv_sub[valid,]),
        error=function(e) NULL)
      if (!is.null(cox)) {
        s <- summary(cox)
        res <- rbind(res, data.frame(cell_type=ct, HR=s$conf.int[1],
          lower95=s$conf.int[3], upper95=s$conf.int[4],
          pvalue=s$coefficients[5], stringsAsFactors=FALSE))
      }
    }
  }
  if (nrow(res) > 0) res$FDR <- p.adjust(res$pvalue, method="BH")
  write.csv(res, paste0("tables/immune/", label, "_immune_survival_univariate_cox.csv"), row.names=FALSE)
  log_msg(sprintf("  %s: %d cell types, %d FDR<0.05", label, nrow(res), sum(res$FDR<0.05)))
  return(res)
}

tcga_meta_full <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_sample_metadata.csv", stringsAsFactors=FALSE)
cox_tcga <- immune_surv_cox(tcga_immune, tcga_surv, tcga_meta_full, "TCGA")

# GSE survival
gse_tum_immune <- gse_immune[, tumor_samples_gse, drop=FALSE]
gse_tum_clin <- gse_clin[gse_clin$tissue == "Tumor (Intrahepatic cholangiocarcinoma)", ]
gse_surv_immune <- gse_tum_clin
gse_surv_immune$OS_time <- as.numeric(gse_surv_immune$survival_mo_) * 30.4375
gse_surv_immune$OS_status <- as.numeric(gse_surv_immune$death)
gse_surv_immune <- gse_surv_immune[!is.na(gse_surv_immune$OS_time) & gse_surv_immune$OS_time > 0, ]

res_gse <- data.frame(stringsAsFactors=FALSE)
for (i in seq_len(nrow(gse_tum_immune))) {
  ct <- rownames(gse_tum_immune)[i]
  vals <- as.numeric(gse_tum_immune[i, ])
  common_s <- intersect(names(vals), gse_surv_immune$sample_name)
  if (length(common_s) < 10) next
  surv_sub <- gse_surv_immune[match(common_s, gse_surv_immune$sample_name), ]
  sv_vals <- vals[common_s]
  valid <- is.finite(sv_vals) & !is.na(surv_sub$OS_time) & surv_sub$OS_time > 0
  if (sum(valid) < 10) next
  cox <- tryCatch(coxph(Surv(OS_time, OS_status) ~ sv_vals, data=surv_sub[valid,]), error=function(e)NULL)
  if (!is.null(cox)) {
    s <- summary(cox)
    res_gse <- rbind(res_gse, data.frame(cell_type=ct, HR=s$conf.int[1],
      lower95=s$conf.int[3], upper95=s$conf.int[4], pvalue=s$coefficients[5], stringsAsFactors=FALSE))
  }
}
if (nrow(res_gse) > 0) res_gse$FDR <- p.adjust(res_gse$pvalue, method="BH")
write.csv(res_gse, "tables/immune/GSE107943_immune_survival_univariate_cox.csv", row.names=FALSE)

# ============================================================
# Step 8: 可视化
# ============================================================
log_msg("\n=== Step 8: 可视化 ===")

tumor_col <- "#E41A1C"; normal_col <- "#377EB8"

# ---- Fig6A/B: Immune boxplot T vs N ----
plot_immune_box <- function(scores, meta, filename, title) {
  common <- intersect(colnames(scores), rownames(meta))
  sc <- scores[, common, drop=FALSE]
  mt <- meta[common, ]

  pdf(filename, width=14, height=10)
  par(mfrow=c(3,4), mar=c(5,4,3,1))
  for (i in seq_len(min(nrow(sc), 12))) {
    nm <- rownames(sc)[i]
    t_vals <- sc[i, mt$sample_type %in% c("Tumor","Primary Tumor")]
    n_vals <- sc[i, mt$sample_type %in% c("Normal","Solid Tissue Normal")]
    boxplot(list(Tumor=t_vals, Normal=n_vals), col=c(tumor_col, normal_col),
            main=nm, ylab="ssGSEA", las=1)
  }
  dev.off()
}

log_msg("Fig6A/B: Immune boxplots")
plot_immune_box(tcga_immune, tcga_meta, "figures/immune/Fig6A_TCGA_immune_score_boxplot_tumor_normal.pdf", "TCGA")
plot_immune_box(gse_immune, gse_meta, "figures/immune/Fig6B_GSE107943_immune_score_boxplot_tumor_normal.pdf", "GSE")

# ---- Fig6C/D: Aggressive vs immune correlation heatmap ----
log_msg("Fig6C/D: Correlation heatmaps")

plot_cor_heatmap <- function(cor_mat, filename, title) {
  pdf(filename, width=10, height=6)
  pheatmap(cor_mat, display_numbers=TRUE, number_format="%.2f",
           color=colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
           breaks=seq(-1,1,length.out=101), main=title, fontsize=10)
  dev.off()
}

plot_cor_heatmap(tcga_cor_mat, "figures/immune/Fig6C_TCGA_aggressive_immune_correlation_heatmap.pdf",
                 "TCGA-CHOL: Pathway Scores vs Immune Cells")
plot_cor_heatmap(gse_cor_mat, "figures/immune/Fig6D_GSE107943_aggressive_immune_correlation_heatmap.pdf",
                 "GSE107943: Pathway Scores vs Immune Cells")

# ---- Fig6E/F: Checkpoint correlation heatmap ----
log_msg("Fig6E/F: Checkpoint heatmaps")

cp_cor_mat_tcga <- matrix(NA, 1, length(checkpoint_avail))
colnames(cp_cor_mat_tcga) <- checkpoint_avail
rownames(cp_cor_mat_tcga) <- "Aggressive"
for (i in seq_len(nrow(tcga_cp_cor))) {
  cp_cor_mat_tcga[1, tcga_cp_cor$checkpoint[i]] <- tcga_cp_cor$rho[i]
}

cp_cor_mat_gse <- matrix(NA, 1, length(checkpoint_avail))
colnames(cp_cor_mat_gse) <- checkpoint_avail
rownames(cp_cor_mat_gse) <- "Aggressive"
for (i in seq_len(nrow(gse_cp_cor))) {
  cp_cor_mat_gse[1, gse_cp_cor$checkpoint[i]] <- gse_cp_cor$rho[i]
}

# Replace NA with 0 for heatmap
cp_cor_mat_tcga[is.na(cp_cor_mat_tcga)] <- 0
cp_cor_mat_gse[is.na(cp_cor_mat_gse)] <- 0

pdf("figures/immune/Fig6E_TCGA_checkpoint_correlation_heatmap.pdf", width=10, height=3)
pheatmap(cp_cor_mat_tcga, display_numbers=TRUE, number_format="%.2f",
         color=colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
         breaks=seq(-1,1,length.out=101), main="TCGA: Aggressive Score vs Checkpoint Genes",
         cluster_rows=FALSE, cluster_cols=FALSE, fontsize=10)
dev.off()

pdf("figures/immune/Fig6F_GSE107943_checkpoint_correlation_heatmap.pdf", width=10, height=3)
pheatmap(cp_cor_mat_gse, display_numbers=TRUE, number_format="%.2f",
         color=colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
         breaks=seq(-1,1,length.out=101), main="GSE107943: Aggressive Score vs Checkpoint Genes",
         cluster_rows=FALSE, cluster_cols=FALSE, fontsize=10)
dev.off()

# ---- Fig6G/H: Aggressive vs Macrophage scatter ----
log_msg("Fig6G/H: Scatter plots")

plot_scatter <- function(score_vec, immune_mat, cell_type, filename, title) {
  common <- intersect(names(score_vec), colnames(immune_mat))
  if (length(common) < 5) return()
  sv <- as.numeric(score_vec[common])
  iv <- as.numeric(immune_mat[cell_type, common])
  df <- data.frame(aggressive=sv, immune=iv)
  sp <- cor.test(sv, iv, method="spearman")

  pdf(filename, width=7, height=6)
  p <- ggplot(df, aes(x=aggressive, y=immune)) + geom_point(size=2, alpha=0.7) +
    geom_smooth(method="lm", se=TRUE, color="red") +
    labs(x="Aggressive Microenvironment Score", y=paste(cell_type, "Score"),
         title=title, subtitle=sprintf("Spearman rho=%.3f, p=%.4f", sp$estimate, sp$p.value)) +
    theme_pubr(base_size=12)
  print(p); dev.off()
}

# Macrophage — wrapped in tryCatch
tryCatch({
  plot_scatter(tcga_aggr, tcga_imm_tum, "Macrophage_TAM",
    "figures/immune/Fig6G_aggressive_vs_macrophage_scatter_TCGA.pdf",
    "TCGA: Aggressive vs Macrophage_TAM")
}, error=function(e) { pdf("figures/immune/Fig6G_aggressive_vs_macrophage_scatter_TCGA.pdf"); plot.new(); text(0.5,0.5,paste("Error:",e$message)); dev.off() })

tryCatch({
  plot_scatter(gse_aggr, gse_imm_tum, "Macrophage_TAM",
    "figures/immune/Fig6H_aggressive_vs_macrophage_scatter_GSE107943.pdf",
    "GSE107943: Aggressive vs Macrophage_TAM")
}, error=function(e) { pdf("figures/immune/Fig6H_aggressive_vs_macrophage_scatter_GSE107943.pdf"); plot.new(); text(0.5,0.5,paste("Error:",e$message)); dev.off() })

# ---- Fig6I/J: Aggressive vs key checkpoint scatter ----
log_msg("Fig6I/J: Checkpoint scatter")
# Use PDCD1 or CTLA4 as representative
key_cp <- intersect(c("PDCD1","CTLA4","HAVCR2","CD274"), rownames(tcga_expr))[1]
tryCatch({
  if (!is.na(key_cp)) {
    common_t <- intersect(colnames(tcga_expr), names(tcga_aggr))
    df_t <- data.frame(aggressive=as.numeric(tcga_aggr[common_t]),
                       checkpoint=as.numeric(tcga_expr[key_cp, common_t]))
    df_t <- df_t[is.finite(df_t$aggressive) & is.finite(df_t$checkpoint),]
    if (nrow(df_t) >= 5) {
      sp_t <- cor.test(df_t$aggressive, df_t$checkpoint, method="spearman")
      pdf("figures/immune/Fig6I_aggressive_vs_checkpoint_scatter_TCGA.pdf", width=7, height=6)
      p <- ggplot(df_t, aes(x=aggressive, y=checkpoint)) + geom_point(size=2, alpha=0.7) +
        geom_smooth(method="lm", se=TRUE, color="red") +
        labs(x="Aggressive Score", y=paste(key_cp,"Expression"),
             title=sprintf("TCGA: Aggressive vs %s",key_cp),
             subtitle=sprintf("rho=%.3f, p=%.4f",sp_t$estimate,sp_t$p.value)) +
        theme_pubr(base_size=12)
      print(p); dev.off()
    }
    common_g <- intersect(colnames(gse_expr), names(gse_aggr))
    df_g <- data.frame(aggressive=as.numeric(gse_aggr[common_g]),
                       checkpoint=as.numeric(gse_expr[key_cp, common_g]))
    df_g <- df_g[is.finite(df_g$aggressive) & is.finite(df_g$checkpoint),]
    if (nrow(df_g) >= 5) {
      sp_g <- cor.test(df_g$aggressive, df_g$checkpoint, method="spearman")
      pdf("figures/immune/Fig6J_aggressive_vs_checkpoint_scatter_GSE107943.pdf", width=7, height=6)
      p <- ggplot(df_g, aes(x=aggressive, y=checkpoint)) + geom_point(size=2, alpha=0.7) +
        geom_smooth(method="lm", se=TRUE, color="red") +
        labs(x="Aggressive Score", y=paste(key_cp,"Expression"),
             title=sprintf("GSE107943: Aggressive vs %s",key_cp),
             subtitle=sprintf("rho=%.3f, p=%.4f",sp_g$estimate,sp_g$p.value)) +
        theme_pubr(base_size=12)
      print(p); dev.off()
    }
  }
}, error=function(e) {
  pdf("figures/immune/Fig6I_aggressive_vs_checkpoint_scatter_TCGA.pdf"); plot.new(); text(0.5,0.5,e$message); dev.off()
  pdf("figures/immune/Fig6J_aggressive_vs_checkpoint_scatter_GSE107943.pdf"); plot.new(); text(0.5,0.5,e$message); dev.off()
})

# ============================================================
# Step 9: 汇总表
# ============================================================
log_msg("\n=== Step 9: 汇总 ===")

# Find consistent signals
if (!is.null(tcga_aggr_cor) && !is.null(gse_aggr_cor)) {
  common_ct <- intersect(tcga_aggr_cor$cell_type, gse_aggr_cor$cell_type)
  consistent <- data.frame(stringsAsFactors=FALSE)
  for (ct in common_ct) {
    r_t <- tcga_aggr_cor$rho[tcga_aggr_cor$cell_type==ct]
    r_g <- gse_aggr_cor$rho[gse_aggr_cor$cell_type==ct]
    if (length(r_t)>0 && length(r_g)>0 && sign(r_t)==sign(r_g)) {
      consistent <- rbind(consistent, data.frame(cell_type=ct, TCGA_rho=r_t, GSE_rho=r_g))
    }
  }
  n_consistent <- nrow(consistent)
} else {
  n_consistent <- 0
}

summary_rows <- data.frame(
  item = c("免疫浸润方法","Cell types scored","TCGA TvsN FDR<0.05","GSE TvsN FDR<0.05",
           "TCGA aggr+immune significant","GSE aggr+immune significant",
           "Consistent immune signals","TCGA checkpoint significant",
           "GSE checkpoint significant","TCGA immune Cox FDR<0.05",
           "GSE immune Cox FDR<0.05","immunedeconv used"),
  value = c("Marker-based ssGSEA (GSVA)", as.character(length(immune_filtered)),
    paste(imm_tn_tcga$cell_type[imm_tn_tcga$FDR<0.05], collapse="; "),
    paste(imm_tn_gse$cell_type[imm_tn_gse$FDR<0.05], collapse="; "),
    if (!is.null(tcga_aggr_cor)) paste(tcga_aggr_cor$cell_type[tcga_aggr_cor$FDR<0.05], collapse="; ") else "N/A",
    if (!is.null(gse_aggr_cor)) paste(gse_aggr_cor$cell_type[gse_aggr_cor$FDR<0.05], collapse="; ") else "N/A",
    as.character(n_consistent),
    if (!is.null(tcga_cp_cor)) paste(tcga_cp_cor$checkpoint[tcga_cp_cor$FDR<0.05], collapse="; ") else "N/A",
    if (!is.null(gse_cp_cor)) paste(gse_cp_cor$checkpoint[gse_cp_cor$FDR<0.05], collapse="; ") else "N/A",
    if (!is.null(cox_tcga)) paste(cox_tcga$cell_type[cox_tcga$FDR<0.05], collapse="; ") else "N/A",
    if (nrow(res_gse)>0) paste(res_gse$cell_type[res_gse$FDR<0.05], collapse="; ") else "N/A",
    if (use_immunedeconv) "Available but not used (marker ssGSEA for consistency)" else "Not available"
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_rows, "tables/immune/immune_analysis_summary.csv", row.names=FALSE)

# ============================================================
# 最终报告
# ============================================================
cat(sprintf("\n========================================\n"))
cat(sprintf("免疫浸润分析完成\n"))
cat(sprintf("========================================\n"))
cat(sprintf("方法: Marker-based ssGSEA\n"))
cat(sprintf("Cell types: %d\n", length(immune_filtered)))
cat(sprintf("TCGA TvsN FDR<0.05: %d\n", sum(imm_tn_tcga$FDR<0.05)))
cat(sprintf("GSE  TvsN FDR<0.05: %d\n", sum(imm_tn_gse$FDR<0.05)))
if (!is.null(tcga_aggr_cor)) cat(sprintf("TCGA aggr+immune sig: %d\n", sum(tcga_aggr_cor$FDR<0.05)))
if (!is.null(gse_aggr_cor))  cat(sprintf("GSE  aggr+immune sig: %d\n", sum(gse_aggr_cor$FDR<0.05)))
cat(sprintf("Consistent immune signals: %d\n", n_consistent))
if (!is.null(tcga_cp_cor)) cat(sprintf("TCGA checkpoint sig: %d\n", sum(tcga_cp_cor$FDR<0.05)))
if (!is.null(gse_cp_cor))  cat(sprintf("GSE  checkpoint sig: %d\n", sum(gse_cp_cor$FDR<0.05)))
cat(sprintf("========================================\n"))
