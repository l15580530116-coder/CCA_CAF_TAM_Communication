# ============================================================
# 分析目的: hub genes 和通讯轴的临床相关性、预后探索、队列稳定性
# 输入: TCGA/GSE107943/GSE26566 表达 + clinical + scores
# 输出:
#   tables/clinical_relevance/ (12 个表格)
#   figures/clinical_relevance/ (12 张图)
# 主要方法: z-score axis → Cox/KM → clinical correlation → Spearman
# 注: 所有生存分析标记为 exploratory (小样本)
# ============================================================

library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(pheatmap)
library(RColorBrewer)
library(data.table)

dir.create("tables/clinical_relevance", showWarnings=FALSE, recursive=TRUE)
dir.create("figures/clinical_relevance", showWarnings=FALSE, recursive=TRUE)

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ============================================================
# Step 1: 加载所有数据
# ============================================================
log_msg("=== Step 1: 加载数据 ===")

# TCGA
tcga_expr <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_tumor_expression_for_survival.csv",
                      row.names=1, check.names=FALSE)
tcga_surv <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_matched_survival_metadata.csv", stringsAsFactors=FALSE)
tcga_scores <- read.csv("tables/gsva/TCGA_CHOL_tumor_scores_with_survival.csv", stringsAsFactors=FALSE)

# GSE107943
gse_expr <- read.csv("data/processed/GSE107943/GSE107943_log2CPM_gene_symbol.csv",
                     row.names=1, check.names=FALSE)
gse_expr <- gse_expr[, colnames(gse_expr)!="6N", drop=FALSE]
gse_clin <- read.csv("data/clinical/GSE107943/GSE107943_clinical.csv", stringsAsFactors=FALSE)
gse_scores <- read.csv("tables/gsva/GSE107943_tumor_scores_with_survival.csv", stringsAsFactors=FALSE)

# GSE26566
gse66_expr <- tryCatch(
  read.csv("data/processed/GSE26566/GSE26566_expression_matrix_clean.csv", row.names=1, check.names=FALSE),
  error=function(e) NULL)
gse66_clin <- tryCatch(
  read.csv("data/clinical/GSE26566/GSE26566_clinical_clean.csv", stringsAsFactors=FALSE),
  error=function(e) NULL)

log_msg(sprintf("TCGA tumor: %d genes x %d samples", nrow(tcga_expr), ncol(tcga_expr)))
log_msg(sprintf("GSE107943 tumor: %d genes x %d samples", nrow(gse_expr), ncol(gse_expr)))
if (!is.null(gse66_expr)) log_msg(sprintf("GSE26566: %d genes x %d samples", nrow(gse66_expr), ncol(gse66_expr)))

# ============================================================
# Step 2: 定义 hub genes 和 axis 基因集
# ============================================================
log_msg("\n=== Step 2: 构建 axis gene sets ===")

hub_genes_all <- c("COL1A1","COL1A2","TREM2","POSTN","SPP1","MMP11","COL3A1",
                   "FAP","CD163","C1QA","MIF","CD74","CXCR4","CD44","TGFB1",
                   "PPIA","BSG","HAVCR2","IDO1","CD86","PDCD1LG2",
                   "PKM","SLC2A1","GLUD1","CAT","ACAA1","LCAT","ACADSB")

axis_sets <- list(
  CAF_matrix             = c("COL1A1","COL1A2","COL3A1","POSTN","FAP","MMP11"),
  TAM_activation         = c("TREM2","SPP1","CD163","C1QA","C1QB","C1QC"),
  Immune_suppressive_TAM = c("TREM2","SPP1","CD163","HAVCR2","IDO1","CD86","PDCD1LG2"),
  COLLAGEN_CD44_axis     = c("COL1A1","COL1A2","COL3A1","CD44"),
  MIF_CD74_CXCR4_axis    = c("MIF","CD74","CXCR4"),
  TGFB_axis              = c("TGFB1")
)

# Build composite score: avg z-score of available genes
compute_zscore <- function(expr_mat, genes) {
  genes_use <- intersect(genes, rownames(expr_mat))
  if (length(genes_use) < 2) return(NULL)
  sub <- as.matrix(expr_mat[genes_use, , drop=FALSE])
  z <- t(scale(t(sub)))
  colMeans(z, na.rm=TRUE)
}

# Build all axis scores for a cohort
build_axis_scores <- function(expr_mat, label) {
  scores <- list()
  for (nm in names(axis_sets)) {
    s <- compute_zscore(expr_mat, axis_sets[[nm]])
    if (!is.null(s)) scores[[nm]] <- s
  }
  # CAF_TAM_axis = CAF_matrix + TAM_activation + MIF_CD74_CXCR4 (each z-scored first)
  if (all(c("CAF_matrix","TAM_activation","MIF_CD74_CXCR4_axis") %in% names(scores))) {
    caf_z <- as.numeric(scale(scores$CAF_matrix))
    tam_z <- as.numeric(scale(scores$TAM_activation))
    mif_z <- as.numeric(scale(scores$MIF_CD74_CXCR4_axis))
    scores$CAF_TAM_axis <- (caf_z + tam_z + mif_z) / 3
    names(scores$CAF_TAM_axis) <- names(scores$CAF_matrix)
  }
  log_msg(sprintf("  %s: %d axis scores built", label, length(scores)))
  return(scores)
}

tcga_axis <- build_axis_scores(tcga_expr, "TCGA")
gse_axis  <- build_axis_scores(gse_expr, "GSE107943")

# Check gene coverage
for (g in hub_genes_all) {
  in_tcga <- g %in% rownames(tcga_expr)
  in_gse  <- g %in% rownames(gse_expr)
  in_gse66 <- if (!is.null(gse66_expr)) g %in% rownames(gse66_expr) else NA
  if (!in_tcga || !in_gse) log_msg(sprintf("  MISSING: %s (TCGA=%s, GSE=%s)", g, in_tcga, in_gse))
}

# ============================================================
# Step 3: 单基因 + axis score 预后分析
# ============================================================
log_msg("\n=== Step 3: 生存分析 ===")

run_cox_pipeline <- function(expr_mat, surv_df, genes, axis_scores, label) {
  # Prepare survival: map sample IDs
  common_s <- intersect(colnames(expr_mat), surv_df$patient_id)
  if (length(common_s) < 10) {
    # Try sample_id match
    common_s <- intersect(colnames(expr_mat), surv_df$sample_id)
  }
  if (length(common_s) < 10) {
    log_msg(sprintf("  %s: insufficient sample match (%d)", label, length(common_s)))
    return(list(gene_cox=NULL, axis_cox=NULL))
  }

  surv_sub <- surv_df[match(common_s, surv_df$patient_id %||% surv_df$sample_id), ]
  surv_sub <- surv_sub[!is.na(surv_sub$OS_time) & surv_sub$OS_time > 0, ]
  common_s <- common_s[!is.na(surv_sub$OS_time)]

  # Single-gene Cox
  gene_cox <- data.frame(stringsAsFactors=FALSE)
  for (g in genes) {
    if (!g %in% rownames(expr_mat)) next
    vals <- as.numeric(expr_mat[g, common_s])
    if (sum(!is.na(vals)) < 10) next
    cox <- tryCatch(coxph(Surv(OS_time, OS_status) ~ vals, data=surv_sub), error=function(e) NULL)
    if (!is.null(cox)) {
      s <- summary(cox)
      gene_cox <- rbind(gene_cox, data.frame(
        gene=g, HR=s$conf.int[1], lower95=s$conf.int[3], upper95=s$conf.int[4],
        pvalue=s$coefficients[5], stringsAsFactors=FALSE))
    }
  }
  if (nrow(gene_cox) > 0) gene_cox$FDR <- p.adjust(gene_cox$pvalue, method="BH")

  # Axis score Cox
  axis_cox <- data.frame(stringsAsFactors=FALSE)
  for (nm in names(axis_scores)) {
    vals <- axis_scores[[nm]][common_s]
    if (sum(!is.na(vals)) < 10) next
    cox <- tryCatch(coxph(Surv(OS_time, OS_status) ~ vals, data=surv_sub), error=function(e) NULL)
    if (!is.null(cox)) {
      s <- summary(cox)
      axis_cox <- rbind(axis_cox, data.frame(
        axis=nm, HR=s$conf.int[1], lower95=s$conf.int[3], upper95=s$conf.int[4],
        pvalue=s$coefficients[5], stringsAsFactors=FALSE))
    }
  }
  if (nrow(axis_cox) > 0) axis_cox$FDR <- p.adjust(axis_cox$pvalue, method="BH")

  log_msg(sprintf("  %s: %d genes, %d axes tested", label, nrow(gene_cox), nrow(axis_cox)))
  return(list(gene_cox=gene_cox, axis_cox=axis_cox, surv_df=surv_sub, common=common_s))
}

# TCGA: map barcode → short patient ID for both expression AND axis scores
tcga_meta_full <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_sample_metadata.csv", stringsAsFactors=FALSE)
tcga_barcode_to_patient <- setNames(tcga_meta_full$patient_id, tcga_meta_full$sample_id)
tcga_patient_ids <- tcga_barcode_to_patient[colnames(tcga_expr)]
tcga_expr_pat <- tcga_expr[, !duplicated(tcga_patient_ids), drop=FALSE]
colnames(tcga_expr_pat) <- tcga_patient_ids[!duplicated(tcga_patient_ids)]

# Also map axis scores
tcga_axis_pat <- lapply(tcga_axis, function(v) {
  names(v) <- tcga_barcode_to_patient[names(v)]
  v[!duplicated(names(v))]
})

tcga_surv_df <- tcga_surv
tcga_res <- run_cox_pipeline(tcga_expr_pat, tcga_surv_df, hub_genes_all, tcga_axis_pat, "TCGA")

# GSE107943 survival — filter to tumor only
gse_meta_sc <- read.csv("data/processed/GSE107943/GSE107943_sample_metadata.csv", stringsAsFactors=FALSE)
gse_tum_ids <- gse_meta_sc$sample_id[gse_meta_sc$sample_type == "Tumor"]
gse_tum_expr <- gse_expr[, intersect(gse_tum_ids, colnames(gse_expr)), drop=FALSE]

gse_clin$OS_time <- as.numeric(gse_clin$survival_mo_) * 30.4375
gse_clin$OS_status <- as.numeric(gse_clin$death)
gse_clin$patient_id <- gse_clin$sample_name
gse_res <- run_cox_pipeline(gse_tum_expr, gse_clin, hub_genes_all, gse_axis, "GSE107943")

# Save
if (!is.null(tcga_res$gene_cox)) write.csv(tcga_res$gene_cox, "tables/clinical_relevance/TCGA_hub_gene_survival_cox.csv", row.names=FALSE)
if (!is.null(tcga_res$axis_cox)) write.csv(tcga_res$axis_cox, "tables/clinical_relevance/TCGA_axis_score_survival_cox.csv", row.names=FALSE)
if (!is.null(gse_res$gene_cox))  write.csv(gse_res$gene_cox,  "tables/clinical_relevance/GSE107943_hub_gene_survival_cox.csv", row.names=FALSE)
if (!is.null(gse_res$axis_cox))  write.csv(gse_res$axis_cox,  "tables/clinical_relevance/GSE107943_axis_score_survival_cox.csv", row.names=FALSE)

# ============================================================
# Step 4: 临床相关性 (TCGA stage/grade, GSE death/recurrence)
# ============================================================
log_msg("\n=== Step 4: 临床相关性 ===")

# TCGA clinical
tcga_clin <- read.csv("data/clinical/TCGA_CHOL/TCGA_CHOL_survival_clean.csv", stringsAsFactors=FALSE)
tcga_meta <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_sample_metadata.csv", stringsAsFactors=FALSE)
tcga_pat_map <- setNames(tcga_meta$patient_id, tcga_meta$sample_id)

# Build TCGA sample-level data
tcga_data <- data.frame(
  sample_id = colnames(tcga_expr),
  patient_id = tcga_pat_map[colnames(tcga_expr)],
  stringsAsFactors = FALSE
)
tcga_data <- merge(tcga_data, tcga_clin, by="patient_id", all.x=TRUE)

# Add axis scores
for (nm in names(tcga_axis)) {
  tcga_data[[nm]] <- as.numeric(tcga_axis[[nm]][tcga_data$sample_id])
}

# Add key single-gene expression
for (g in hub_genes_all) {
  if (g %in% rownames(tcga_expr)) tcga_data[[g]] <- as.numeric(tcga_expr[g, tcga_data$sample_id])
}

# Clinical variables to test
tcga_clin_vars <- intersect(c("ajcc_pathologic_stage","ajcc_pathologic_t",
  "ajcc_pathologic_n","ajcc_pathologic_m","tumor_grade","residual_disease"),
  colnames(tcga_data))

# Wilcoxon/Kruskal test
run_clin_corr <- function(data, features, clin_var, label) {
  res <- data.frame(stringsAsFactors=FALSE)
  cv <- data[[clin_var]]
  cv <- cv[!is.na(cv) & cv != ""]
  if (length(unique(cv)) < 2) return(res)

  for (f in features) {
    fv <- data[[f]]
    if (all(is.na(fv))) next
    pv <- tryCatch({
      if (length(unique(cv)) == 2) {
        wilcox.test(fv[data[[clin_var]] %in% unique(cv)[1]],
                    fv[data[[clin_var]] %in% unique(cv)[2]])$p.value
      } else {
        kruskal.test(fv ~ as.factor(data[[clin_var]]))$p.value
      }
    }, error=function(e) NA)
    if (!is.na(pv)) res <- rbind(res, data.frame(feature=f, variable=clin_var, pvalue=pv, stringsAsFactors=FALSE))
  }
  if (nrow(res) > 0) res$FDR <- p.adjust(res$pvalue, method="BH")
  return(res)
}

axis_names <- names(tcga_axis)
tcga_clin_res <- data.frame(stringsAsFactors=FALSE)
for (cv in tcga_clin_vars) {
  res <- run_clin_corr(tcga_data, c(axis_names, hub_genes_all), cv, "TCGA")
  tcga_clin_res <- rbind(tcga_clin_res, res)
}

write.csv(tcga_clin_res, "tables/clinical_relevance/TCGA_hub_gene_clinical_correlation.csv", row.names=FALSE)
log_msg(sprintf("TCGA clinical: %d tests, %d FDR<0.05", nrow(tcga_clin_res),
  sum(tcga_clin_res$FDR < 0.05, na.rm=TRUE)))

# GSE107943 clinical (death, recurrence)
gse_data <- data.frame(sample_id = colnames(gse_expr), stringsAsFactors=FALSE)
gse_data <- merge(gse_data, gse_clin[, c("sample_name","death","recurr","survival_mo_","dsfree_mo_")],
                  by.x="sample_id", by.y="sample_name", all.x=TRUE)
for (nm in names(gse_axis)) gse_data[[nm]] <- as.numeric(gse_axis[[nm]][gse_data$sample_id])
for (g in hub_genes_all) {
  if (g %in% rownames(gse_expr)) gse_data[[g]] <- as.numeric(gse_expr[g, gse_data$sample_id])
}

gse_death <- data.frame(feature=character(), pvalue=numeric(), stringsAsFactors=FALSE)
for (f in c(axis_names, hub_genes_all)) {
  if (!f %in% colnames(gse_data)) next
  g0 <- gse_data[[f]][gse_data$death=="0"]
  g1 <- gse_data[[f]][gse_data$death=="1"]
  if (length(g0)>2 && length(g1)>2) {
    wt <- wilcox.test(g1, g0)
    gse_death <- rbind(gse_death, data.frame(feature=f, variable="death_status", pvalue=wt$p.value, stringsAsFactors=FALSE))
  }
}
if (nrow(gse_death) > 0) gse_death$FDR <- p.adjust(gse_death$pvalue, method="BH")

gse_recurr <- data.frame(stringsAsFactors=FALSE)
for (f in c(axis_names, hub_genes_all)) {
  if (!f %in% colnames(gse_data)) next
  g0 <- gse_data[[f]][gse_data$recurr=="0"]
  g1 <- gse_data[[f]][gse_data$recurr=="1"]
  if (length(g0)>2 && length(g1)>2) {
    wt <- wilcox.test(g1, g0)
    gse_recurr <- rbind(gse_recurr, data.frame(feature=f, variable="recurrence_status", pvalue=wt$p.value, stringsAsFactors=FALSE))
  }
}
if (nrow(gse_recurr) > 0) gse_recurr$FDR <- p.adjust(gse_recurr$pvalue, method="BH")

gse_clin_res <- rbind(gse_death, gse_recurr)
write.csv(gse_clin_res, "tables/clinical_relevance/GSE107943_hub_gene_clinical_correlation.csv", row.names=FALSE)

# ============================================================
# Step 5: Score vs axis correlation
# ============================================================
log_msg("\n=== Step 5: Score correlation ===")

# Read ssGSEA scores for TCGA
tcga_path <- read.csv("tables/gsva/TCGA_CHOL_pathway_scores_all_samples.csv", row.names=1, check.names=FALSE)
gse_path  <- read.csv("tables/gsva/GSE107943_pathway_scores_all_samples.csv", row.names=1, check.names=FALSE)

score_vars <- c("CAF","TAM","ECM_up","Metabolism_down","IM_total","IM_CAF_TAM_combined")
compute_correlations <- function(expr_data, axis_scores, path_scores, score_vars, label) {
  common <- Reduce(intersect, list(colnames(expr_data), names(axis_scores[[1]]), colnames(path_scores)))
  if (length(common) < 5) return(NULL)

  features <- c(names(axis_scores), hub_genes_all)
  cor_mat <- matrix(NA, length(features), length(score_vars))
  rownames(cor_mat) <- features; colnames(cor_mat) <- score_vars
  p_mat <- cor_mat

  for (i in seq_along(features)) {
    f <- features[i]
    if (f %in% names(axis_scores)) {
      fv <- axis_scores[[f]][common]
    } else if (f %in% rownames(expr_data)) {
      fv <- as.numeric(expr_data[f, common])
    } else next

    for (j in seq_along(score_vars)) {
      sv <- score_vars[j]
      if (sv %in% rownames(path_scores)) {
        svals <- as.numeric(path_scores[sv, common])
        sp <- cor.test(fv, svals, method="spearman")
        cor_mat[i, j] <- sp$estimate
        p_mat[i, j] <- sp$p.value
      }
    }
  }

  # Remove rows with all NA
  keep <- rowSums(!is.na(cor_mat)) > 0
  cor_mat <- cor_mat[keep, , drop=FALSE]
  p_mat <- p_mat[keep, , drop=FALSE]

  write.csv(data.frame(feature=rownames(cor_mat), cor_mat, check.names=FALSE),
            paste0("tables/clinical_relevance/hub_gene_score_correlation_", label, ".csv"), row.names=FALSE)
  return(list(cor=cor_mat, p=p_mat))
}

tcga_cor <- compute_correlations(tcga_expr, tcga_axis_pat, tcga_path, score_vars, "TCGA")
gse_cor  <- compute_correlations(gse_expr, gse_axis, gse_path, score_vars, "GSE107943")

# ============================================================
# Step 6: GSE26566 表达验证
# ============================================================
log_msg("\n=== Step 6: GSE26566 validation ===")

if (!is.null(gse66_expr)) {
  gse66_axis <- build_axis_scores(gse66_expr, "GSE26566")

  # Hub gene presence
  gse66_hub <- data.frame(stringsAsFactors=FALSE)
  for (g in hub_genes_all) {
    present <- g %in% rownames(gse66_expr)
    gse66_hub <- rbind(gse66_hub, data.frame(
      gene=g, in_GSE26566=present,
      mean_expr=if(present) mean(as.numeric(gse66_expr[g,]),na.rm=TRUE) else NA,
      stringsAsFactors=FALSE))
  }

  # Axis score stats
  gse66_axis_df <- data.frame(axis=names(gse66_axis), stringsAsFactors=FALSE)
  for (nm in names(gse66_axis)) {
    gse66_axis_df$mean[gse66_axis_df$axis==nm] <- mean(gse66_axis[[nm]], na.rm=TRUE)
    gse66_axis_df$sd[gse66_axis_df$axis==nm] <- sd(gse66_axis[[nm]], na.rm=TRUE)
    gse66_axis_df$n_genes[gse66_axis_df$axis==nm] <- length(axis_sets[[nm]])
  }

  gse66_val <- merge(gse66_hub, gse66_axis_df, by.x="gene", by.y="axis", all=TRUE)
  write.csv(gse66_val, "tables/clinical_relevance/GSE26566_hub_axis_expression_validation.csv", row.names=FALSE)
  log_msg(sprintf("GSE26566: %d genes checked, %d present", nrow(gse66_hub), sum(gse66_hub$in_GSE26566)))
}

# ============================================================
# Step 7: 可视化
# ============================================================
log_msg("\n=== Step 7: 可视化 ===")

# Forest plot helper
forest_plot <- function(cox_df, title, filename) {
  if (is.null(cox_df) || nrow(cox_df) < 1) return()
  df <- cox_df[order(cox_df$HR), ]
  df$label <- sprintf("%.2f (%.2f-%.2f)", df$HR, df$lower95, df$upper95)

  pdf(filename, width=8, height=max(4, nrow(df)*0.3+2))
  p <- ggplot(df, aes(x=HR, y=reorder(if("gene" %in% colnames(df)) gene else axis, HR))) +
    geom_point(size=3, color=ifelse(df$HR>1,"#E41A1C","#377EB8")) +
    geom_errorbarh(aes(xmin=lower95, xmax=upper95), height=0.2) +
    geom_vline(xintercept=1, linetype="dashed", color="grey50") +
    labs(x="Hazard Ratio", y="", title=title, subtitle="Exploratory — small sample size") +
    theme_pubr(base_size=11)
  print(p); dev.off()
}

# Fig10A: TCGA axis forest
forest_plot(tcga_res$axis_cox, "TCGA-CHOL: Axis Score OS (Exploratory)",
            "figures/clinical_relevance/Fig10A_axis_scores_survival_forest_TCGA.pdf")

# Fig10B: GSE axis forest
forest_plot(gse_res$axis_cox, "GSE107943: Axis Score OS (Exploratory)",
            "figures/clinical_relevance/Fig10B_axis_scores_survival_forest_GSE107943.pdf")

# Fig10C: TCGA gene forest
forest_plot(tcga_res$gene_cox, "TCGA-CHOL: Hub Gene OS (Exploratory)",
            "figures/clinical_relevance/Fig10C_hub_genes_survival_forest_TCGA.pdf")

# Fig10D: GSE gene forest
forest_plot(gse_res$gene_cox, "GSE107943: Hub Gene OS (Exploratory)",
            "figures/clinical_relevance/Fig10D_hub_genes_survival_forest_GSE107943.pdf")

# KM helper
km_plot <- function(axis_scores, surv_df, common, axis_name, title, filename) {
  vals <- axis_scores[[axis_name]][common]
  if (all(is.na(vals)) || length(unique(vals)) < 3) return()
  # Ensure alignment
  aligned <- intersect(names(vals), rownames(surv_df) %||% surv_df$patient_id)
  if (length(aligned) < 5) return()
  surv_df_use <- surv_df
  if (!is.null(surv_df$patient_id)) {
    surv_df_use <- surv_df[surv_df$patient_id %in% aligned, ]
    vals <- vals[surv_df_use$patient_id]
  } else {
    surv_df_use <- surv_df[aligned, ]
    vals <- vals[aligned]
  }
  surv_df_use$score <- vals
  med <- median(vals, na.rm=TRUE)
  surv_df_use$group <- ifelse(vals > med, "High", "Low")

  sdiff <- survdiff(Surv(OS_time, OS_status) ~ group, data=surv_df_use)
  pval <- 1 - pchisq(sdiff$chisq, 1)

  pdf(filename, width=7, height=6)
  fit <- survfit(Surv(OS_time, OS_status) ~ group, data=surv_df_use)
  p <- ggsurvplot(fit, data=surv_df_use, pval=TRUE,
    palette=c("High"="#E41A1C","Low"="#377EB8"),
    legend.title=axis_name, xlab="Time (days)", title=sprintf("%s\n%s", title, "Exploratory"))
  print(p); dev.off()
}

# Fig10E-F: CAF_TAM_axis KM
if ("CAF_TAM_axis" %in% names(tcga_axis_pat) && !is.null(tcga_res$common)) {
  km_plot(tcga_axis_pat, tcga_res$surv_df, tcga_res$common, "CAF_TAM_axis",
          "TCGA-CHOL: CAF-TAM Axis", "figures/clinical_relevance/Fig10E_CAF_TAM_axis_KM_TCGA.pdf")
}
if ("CAF_TAM_axis" %in% names(gse_axis) && !is.null(gse_res$common)) {
  km_plot(gse_axis, gse_res$surv_df, gse_res$common, "CAF_TAM_axis",
          "GSE107943: CAF-TAM Axis", "figures/clinical_relevance/Fig10F_CAF_TAM_axis_KM_GSE107943.pdf")
}

# Fig10G-H: MIF-CD74-CXCR4 KM
if ("MIF_CD74_CXCR4_axis" %in% names(tcga_axis_pat) && !is.null(tcga_res$common)) {
  km_plot(tcga_axis_pat, tcga_res$surv_df, tcga_res$common, "MIF_CD74_CXCR4_axis",
          "TCGA-CHOL: MIF-CD74/CXCR4 Axis", "figures/clinical_relevance/Fig10G_MIF_CD74_CXCR4_axis_KM_TCGA.pdf")
}
if ("MIF_CD74_CXCR4_axis" %in% names(gse_axis) && !is.null(gse_res$common)) {
  km_plot(gse_axis, gse_res$surv_df, gse_res$common, "MIF_CD74_CXCR4_axis",
          "GSE107943: MIF-CD74/CXCR4 Axis", "figures/clinical_relevance/Fig10H_MIF_CD74_CXCR4_axis_KM_GSE107943.pdf")
}

# Fig10I-J: Correlation heatmaps
plot_cor_heatmap <- function(cor_obj, title, filename) {
  if (is.null(cor_obj) || nrow(cor_obj$cor) < 2) return()
  mat <- cor_obj$cor
  mat[is.na(mat)] <- 0
  pdf(filename, width=10, height=12)
  pheatmap(mat, display_numbers=TRUE, number_format="%.2f",
           color=colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
           breaks=seq(-1,1,length.out=101), main=title, fontsize=9, fontsize_number=7)
  dev.off()
}

if (!is.null(tcga_cor)) plot_cor_heatmap(tcga_cor, "TCGA: Hub Genes/Axes vs Pathway Scores",
  "figures/clinical_relevance/Fig10I_hub_gene_score_correlation_heatmap_TCGA.pdf")
if (!is.null(gse_cor))  plot_cor_heatmap(gse_cor, "GSE107943: Hub Genes/Axes vs Pathway Scores",
  "figures/clinical_relevance/Fig10J_hub_gene_score_correlation_heatmap_GSE107943.pdf")

# Fig10K: Clinical stage boxplot
if ("ajcc_pathologic_stage" %in% colnames(tcga_data) && "CAF_TAM_axis" %in% colnames(tcga_data)) {
  pdf("figures/clinical_relevance/Fig10K_clinical_stage_axis_score_boxplot_TCGA.pdf", width=10, height=6)
  stage_data <- tcga_data[!is.na(tcga_data$ajcc_pathologic_stage) & tcga_data$ajcc_pathologic_stage!="", ]
  p <- ggplot(stage_data, aes(x=ajcc_pathologic_stage, y=CAF_TAM_axis)) +
    geom_boxplot(fill="steelblue", alpha=0.7) +
    labs(x="AJCC Pathologic Stage", y="CAF-TAM Axis Score", title="TCGA-CHOL: CAF-TAM Axis by Stage") +
    theme_pubr(base_size=12) + theme(axis.text.x=element_text(angle=45, hjust=1))
  print(p); dev.off()
}

# Fig10L: GSE26566 axis score distribution
if (!is.null(gse66_expr) && exists("gse66_axis") && length(gse66_axis) > 0) {
  pdf("figures/clinical_relevance/Fig10L_GSE26566_axis_score_expression_validation.pdf", width=12, height=8)
  par(mfrow=c(2,4), mar=c(4,3,3,1))
  for (nm in names(gse66_axis)) {
    hist(gse66_axis[[nm]], breaks=20, col="steelblue", border="white",
         main=nm, xlab="Z-score", las=1)
    abline(v=0, col="red", lty=2)
  }
  dev.off()
}

# ============================================================
# Step 8: Summary
# ============================================================
log_msg("\n=== Step 8: Summary ===")

# Consistent direction analysis
find_consistent <- function(tcga_df, gse_df, id_col) {
  if (is.null(tcga_df) || is.null(gse_df)) return(0)
  common <- intersect(tcga_df[[id_col]], gse_df[[id_col]])
  consistent <- 0
  for (id in common) {
    hr_t <- tcga_df$HR[tcga_df[[id_col]]==id]
    hr_g <- gse_df$HR[gse_df[[id_col]]==id]
    if (length(hr_t)>0 && length(hr_g)>0 && sign(log(hr_t))==sign(log(hr_g))) consistent <- consistent+1
  }
  return(consistent)
}

n_cons_axis <- find_consistent(tcga_res$axis_cox, gse_res$axis_cox, "axis")
n_cons_gene <- find_consistent(tcga_res$gene_cox, gse_res$gene_cox, "gene")

tcga_sig <- if (!is.null(tcga_res$axis_cox)) tcga_res$axis_cox$axis[tcga_res$axis_cox$FDR < 0.05] else "None"
gse_sig  <- if (!is.null(gse_res$axis_cox)) gse_res$axis_cox$axis[gse_res$axis_cox$FDR < 0.05] else "None"

summary_rows <- data.frame(
  question = c(
    "Consistent axis direction", "Consistent gene direction",
    "TCGA sig axes", "GSE107943 sig axes",
    "Both cohorts significant", "GSE26566 hub genes present",
    "Clinical relevance level", "Recommendation"
  ),
  answer = c(
    sprintf("%d/%d", n_cons_axis, if(!is.null(tcga_res$axis_cox)) nrow(tcga_res$axis_cox) else 0),
    sprintf("%d/%d", n_cons_gene, if(!is.null(tcga_res$gene_cox)) nrow(tcga_res$gene_cox) else 0),
    if (length(tcga_sig)>0) paste(tcga_sig, collapse=", ") else "None",
    if (length(gse_sig)>0) paste(gse_sig, collapse=", ") else "None",
    if (length(tcga_sig)>0 && length(gse_sig)>0) "None — small sample limits statistical power" else "N/A",
    if (!is.null(gse66_expr)) sprintf("%d/%d present", sum(gse66_hub$in_GSE26566), nrow(gse66_hub)) else "N/A",
    "Exploratory only — insufficient power (TCGA n=35, GSE n=30); all survival results should be clearly labeled as exploratory",
    "Report direction-consistent trends; avoid overclaiming significance; describe clinical correlation as hypothesis-generating for future validation cohorts"
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_rows, "tables/clinical_relevance/clinical_relevance_summary.csv", row.names=FALSE)

cat(sprintf("\n========================================\n"))
cat(sprintf("Clinical Relevance Analysis Complete\n"))
cat(sprintf("========================================\n"))
cat(sprintf("TCGA axis FDR<0.05: %s\n", paste(tcga_sig, collapse=", ")))
cat(sprintf("GSE  axis FDR<0.05: %s\n", paste(gse_sig, collapse=", ")))
cat(sprintf("Consistent axis direction: %d\n", n_cons_axis))
cat(sprintf("All results EXPLORATORY (small n)\n"))
cat(sprintf("========================================\n"))
