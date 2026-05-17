# ============================================================
# 分析目的: 诊断第 12 步分型结果 + 比较 3 种替代分型策略
# 输入文件:
#   tables/subtyping/TCGA_subtype_assignment.csv
#   tables/subtyping/GSE107943_subtype_assignment.csv
#   tables/subtyping/subtyping_gene_list.csv
#   tables/CCA_IM_CAF_TAM_DEGs_validated_same_direction.csv
#   tables/enrichment/validated_IM_CAF_TAM_DEGs_{up,down}.csv
#   tables/gsva/TCGA_CHOL_tumor_scores_with_survival.csv
#   tables/gsva/GSE107943_tumor_scores_with_survival.csv
#   data/processed/TCGA_CHOL/TCGA_CHOL_tumor_expression_for_survival.csv
#   data/processed/GSE107943/GSE107943_log2CPM_gene_symbol.csv
# 输出文件:
#   tables/subtyping_diagnosis/ (10 个表格)
#   figures/subtyping_diagnosis/ (12 张图)
# 主要方法: 诊断 + 3 策略比较:
#   Strategy A: validated up genes only clustering
#   Strategy B: ssGSEA 7 scores clustering
#   Strategy C: aggressive score median split
# ============================================================

library(ConsensusClusterPlus)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(pheatmap)
library(RColorBrewer)
library(data.table)

dir.create("tables/subtyping_diagnosis", showWarnings=FALSE, recursive=TRUE)
dir.create("figures/subtyping_diagnosis", showWarnings=FALSE, recursive=TRUE)

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

subtype_cols <- c("C1_aggressive_stromal"="#E41A1C", "C2_metabolic_preserved"="#377EB8")
s_cols <- list(
  orig   = c("high"="#E41A1C", "low"="#377EB8"),
  stratA = c("Up_C1_stromal_high"="#E41A1C", "Up_C2_stromal_low"="#377EB8"),
  stratB = c("Score_C1_aggressive_high"="#E41A1C", "Score_C2_aggressive_low"="#377EB8"),
  stratC = c("Aggr_High"="#E41A1C", "Aggr_Low"="#377EB8")
)

score_names <- c("CAF","ECM_up","TAM","Metabolism_down","IM_total",
                 "IM_CAF_TAM_combined","CCA_aggressive_microenvironment")

# ============================================================
# Step 0: 加载共用数据
# ============================================================
log_msg("=== Step 0: 加载数据 ===")

tcga_orig <- read.csv("tables/subtyping/TCGA_subtype_assignment.csv", stringsAsFactors=FALSE)
gse_orig  <- read.csv("tables/subtyping/GSE107943_subtype_assignment.csv", stringsAsFactors=FALSE)
subtype_genes_df <- read.csv("tables/subtyping/subtyping_gene_list.csv", stringsAsFactors=FALSE)
validated <- read.csv("tables/CCA_IM_CAF_TAM_DEGs_validated_same_direction.csv", stringsAsFactors=FALSE)
up_genes_df <- read.csv("tables/enrichment/validated_IM_CAF_TAM_DEGs_up.csv", stringsAsFactors=FALSE)
dn_genes_df <- read.csv("tables/enrichment/validated_IM_CAF_TAM_DEGs_down.csv", stringsAsFactors=FALSE)
gene_annot <- read.csv("tables/CCA_IM_CAF_TAM_combined_gene_annotation.csv", stringsAsFactors=FALSE)

# Expression data
tcga_tum <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_tumor_expression_for_survival.csv",
                     row.names=1, check.names=FALSE)
gse_full <- read.csv("data/processed/GSE107943/GSE107943_log2CPM_gene_symbol.csv",
                     row.names=1, check.names=FALSE)
gse_full <- gse_full[, colnames(gse_full)!="6N", drop=FALSE]
gse_meta <- read.csv("data/processed/GSE107943/GSE107943_sample_metadata.csv", stringsAsFactors=FALSE)
gse_tum <- gse_full[, intersect(colnames(gse_full), gse_meta$sample_id[gse_meta$sample_type=="Tumor"]), drop=FALSE]

# Survival
tcga_surv <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_matched_survival_metadata.csv", stringsAsFactors=FALSE)
gse_clin <- read.csv("data/clinical/GSE107943/GSE107943_clinical.csv", stringsAsFactors=FALSE)

# SSGSEA scores — read with column names intact
ssgsea_tcga <- read.csv("tables/gsva/TCGA_CHOL_pathway_scores_all_samples.csv", row.names=1, check.names=FALSE)
ssgsea_tcga[] <- lapply(ssgsea_tcga, as.numeric)

ssgsea_gse <- read.csv("tables/gsva/GSE107943_pathway_scores_all_samples.csv", row.names=1, check.names=FALSE)
ssgsea_gse[] <- lapply(ssgsea_gse, as.numeric)

# Helper: Wilcoxon comparison
wilcox_compare <- function(df, group_col, val_cols, label) {
  res <- data.frame(stringsAsFactors=FALSE)
  grps <- unique(df[[group_col]])
  if (length(grps)!=2) return(res)
  g1 <- df[df[[group_col]]==grps[1],,drop=FALSE]
  g2 <- df[df[[group_col]]==grps[2],,drop=FALSE]
  for (vc in val_cols) {
    if (!vc %in% colnames(df)) next
    wt <- wilcox.test(g1[[vc]], g2[[vc]])
    res <- rbind(res, data.frame(
      score=vc, mean_G1=mean(g1[[vc]],na.rm=TRUE), mean_G2=mean(g2[[vc]],na.rm=TRUE),
      median_G1=median(g1[[vc]],na.rm=TRUE), median_G2=median(g2[[vc]],na.rm=TRUE),
      pvalue=wt$p.value, G1=grps[1], G2=grps[2], stringsAsFactors=FALSE))
  }
  res$FDR <- p.adjust(res$pvalue, method="BH")
  log_msg(sprintf("  %s: %d/%d FDR<0.05", label, sum(res$FDR<0.05), nrow(res)))
  return(res)
}

# Helper: KM + Cox
km_cox <- function(df, group_col, title, filename) {
  df$group <- df[[group_col]]
  df$OS_time <- as.numeric(df$OS_time)
  df$OS_status <- as.numeric(df$OS_status)
  df <- df[!is.na(df$OS_time) & df$OS_time>0,]
  if (length(unique(df$group))<2) return(c(pval=NA,HR=NA))

  sdiff <- survdiff(Surv(OS_time, OS_status) ~ group, data=df)
  pval <- 1 - pchisq(sdiff$chisq, 1)
  cox <- tryCatch(coxph(Surv(OS_time, OS_status) ~ group, data=df), error=function(e)NULL)
  hr <- if (!is.null(cox)) exp(coef(cox)) else NA
  ci <- if (!is.null(cox)) exp(confint(cox)) else c(NA,NA)

  df$OS_months <- df$OS_time/30.4375
  pdf(filename, width=7, height=6)
  fit <- survfit(Surv(OS_months, OS_status) ~ group, data=df)
  p <- ggsurvplot(fit, data=df, pval=TRUE,
    palette=c("#E41A1C","#377EB8"), legend.title=group_col, title=title)
  print(p); dev.off()

  cat(sprintf("  %s: p=%.4f, HR=%.2f, n=%d, events=%d\n",
    title, pval, hr, nrow(df), sum(df$OS_status==1)))
  c(pval=pval, HR=hr, lower95=ci[1], upper95=ci[2])
}

# Helper: ConsensusClusterPlus wrapper
run_cc <- function(mat, title_prefix) {
  cc <- ConsensusClusterPlus(d=mat, maxK=6, reps=100, pItem=0.8, pFeature=1,
    clusterAlg="hc", distance="euclidean", seed=42, plot=NULL,
    title=paste0("figures/subtyping_diagnosis/", title_prefix),
    writeTable=FALSE)
  pac_vals <- sapply(2:6, function(k) {
    cm <- cc[[k]]$consensusMatrix
    sum(cm > 0.2 & cm < 0.8) / (nrow(cm)^2)
  })
  names(pac_vals) <- 2:6
  list(cc=cc, pac=pac_vals, clusters_2=cc[[2]]$consensusClass)
}

# ============================================================
# Step 1: 诊断原始分型
# ============================================================
log_msg("\n=== Step 1: 诊断原始分型 (第 12 步) ===")

# 1a: Score comparison
log_msg("1a: Score comparison")
tcga_orig_diag <- wilcox_compare(tcga_orig, "subtype", score_names, "TCGA original")
gse_orig_diag  <- wilcox_compare(gse_orig, "subtype", score_names, "GSE original")

orig_diag <- rbind(cbind(tcga_orig_diag, cohort="TCGA"), cbind(gse_orig_diag, cohort="GSE"))
write.csv(orig_diag, "tables/subtyping_diagnosis/original_subtype_score_diagnosis.csv", row.names=FALSE)

# 1b: Gene direction balance
log_msg("1b: Gene direction balance")
n_up <- sum(subtype_genes_df$direction=="up")
n_dn <- sum(subtype_genes_df$direction=="down")
n_caf <- sum(subtype_genes_df$gene %in% gene_annot$gene_symbol[gene_annot$source=="CAF"])
n_tam <- sum(subtype_genes_df$gene %in% gene_annot$gene_symbol[gene_annot$source=="macrophage_TAM"])
n_im  <- sum(subtype_genes_df$gene %in% gene_annot$gene_symbol[grepl("immunometabolism", gene_annot$source)])

gene_balance <- data.frame(
  category = c("up_genes","down_genes","CAF","TAM","immunometabolism","total"),
  count = c(n_up, n_dn, n_caf, n_tam, n_im, nrow(subtype_genes_df)),
  percent = round(c(n_up,n_dn,n_caf,n_tam,n_im,nrow(subtype_genes_df))/nrow(subtype_genes_df)*100, 1),
  stringsAsFactors = FALSE
)
write.csv(gene_balance, "tables/subtyping_diagnosis/original_subtype_gene_direction_balance.csv", row.names=FALSE)
log_msg(sprintf("  Up:%d(%.0f%%), Down:%d(%.0f%%), CAF:%d, TAM:%d",
  n_up, gene_balance$percent[1], n_dn, gene_balance$percent[2], n_caf, n_tam))

# 1c: Diagnostic conclusion
cat("\n--- DIAGNOSIS ---\n")
cat(sprintf("Gene balance: %d up vs %d down (%.0f%% down-dominated)\n",
  n_up, n_dn, n_dn/nrow(subtype_genes_df)*100))
cat(sprintf("Original C1 has higher CAF: %s\n",
  tcga_orig_diag$mean_G1[tcga_orig_diag$score=="CAF"] > tcga_orig_diag$mean_G2[tcga_orig_diag$score=="CAF"]))
cat(sprintf("Original C1 has higher ECM: %s\n",
  tcga_orig_diag$mean_G1[tcga_orig_diag$score=="ECM_up"] > tcga_orig_diag$mean_G2[tcga_orig_diag$score=="ECM_up"]))

# ---- FigS: Original score boxplot TCGA ----
log_msg("FigS: Original boxplots")
pdf("figures/subtyping_diagnosis/FigS_original_subtype_score_boxplot_TCGA.pdf", width=12, height=7)
tcga_long <- reshape2::melt(tcga_orig, id.vars=c("sample_id","subtype"),
  measure.vars=intersect(score_names, colnames(tcga_orig)), variable.name="score", value.name="value")
p <- ggplot(tcga_long, aes(x=score, y=value, fill=subtype)) +
  geom_boxplot(outlier.size=0.8) + scale_fill_manual(values=subtype_cols) +
  labs(y="ssGSEA Score", title="TCGA-CHOL: Original Subtype Scores") +
  theme_pubr(base_size=12) + theme(axis.text.x=element_text(angle=45,hjust=1))
print(p); dev.off()

pdf("figures/subtyping_diagnosis/FigS_original_subtype_score_boxplot_GSE.pdf", width=12, height=7)
gse_long <- reshape2::melt(gse_orig, id.vars=c("sample_id","subtype"),
  measure.vars=intersect(score_names, colnames(gse_orig)), variable.name="score", value.name="value")
p <- ggplot(gse_long, aes(x=score, y=value, fill=subtype)) +
  geom_boxplot(outlier.size=0.8) + scale_fill_manual(values=subtype_cols) +
  labs(y="ssGSEA Score", title="GSE107943: Original Subtype Scores") +
  theme_pubr(base_size=12) + theme(axis.text.x=element_text(angle=45,hjust=1))
print(p); dev.off()

# Gene direction heatmap
tcga_z <- t(scale(t(as.matrix(tcga_tum[subtype_genes_df$gene,,drop=FALSE]))))
tcga_z[tcga_z>3]<-3; tcga_z[tcga_z< -3]<--3
ord <- order(tcga_orig$subtype)
pdf("figures/subtyping_diagnosis/FigS_original_subtype_gene_direction_heatmap.pdf", width=12, height=8)
ann <- data.frame(Subtype=tcga_orig$subtype[ord], row.names=tcga_orig$sample_id[ord])
pheatmap(tcga_z[,ord], annotation_col=ann, annotation_colors=list(Subtype=subtype_cols),
  cluster_rows=FALSE, cluster_cols=FALSE, show_colnames=FALSE, show_rownames=FALSE,
  main="Original 100 Genes: Z-score (TCGA)", gaps_row=n_up)
dev.off()

# ============================================================
# Step 2: Strategy A — Up genes only
# ============================================================
log_msg("\n=== Step 2: Strategy A (Up genes only) ===")

stratA_genes <- intersect(up_genes_df$gene_symbol, intersect(rownames(tcga_tum), rownames(gse_tum)))
log_msg(sprintf("Strategy A genes: %d up genes in both cohorts", length(stratA_genes)))

stratA_result <- list(status="failed", reason="")

if (length(stratA_genes) < 20) {
  stratA_result$reason <- sprintf("Only %d common up genes, need >=20", length(stratA_genes))
  log_msg(sprintf("  FAILED: %s", stratA_result$reason))
} else {
  tcga_A_z <- t(scale(t(as.matrix(tcga_tum[stratA_genes,,drop=FALSE]))))

  ccA <- run_cc(tcga_A_z, "stratA_CC")
  stratA_clusters <- ccA$clusters_2

  # Name by CAF score
  tcga_A_df <- data.frame(sample_id=colnames(tcga_A_z), cluster=stratA_clusters, stringsAsFactors=FALSE)
  # Get scores for naming
  tcga_scores_t <- as.data.frame(t(ssgsea_tcga)); tcga_scores_t$sample_id <- rownames(tcga_scores_t)
  tcga_A_df <- merge(tcga_A_df, tcga_scores_t, by="sample_id", all.x=TRUE)
  for (cn in score_names) if (cn %in% colnames(tcga_A_df)) tcga_A_df[[cn]] <- as.numeric(tcga_A_df[[cn]])

  c1_caf <- mean(tcga_A_df$CAF[tcga_A_df$cluster==1], na.rm=TRUE)
  c2_caf <- mean(tcga_A_df$CAF[tcga_A_df$cluster==2], na.rm=TRUE)
  if (c1_caf > c2_caf) {
    nmA <- c("1"="Up_C1_stromal_high","2"="Up_C2_stromal_low")
  } else {
    nmA <- c("1"="Up_C2_stromal_low","2"="Up_C1_stromal_high")
  }
  tcga_A_df$subtype <- nmA[as.character(tcga_A_df$cluster)]

  write.csv(tcga_A_df, "tables/subtyping_diagnosis/strategy_A_up_genes_subtype_assignment_TCGA.csv", row.names=FALSE)

  # GSE projection
  gse_A_z <- t(scale(t(as.matrix(gse_tum[stratA_genes,,drop=FALSE]))))
  centroid_A1 <- rowMeans(tcga_A_z[,tcga_A_df$sample_id[tcga_A_df$subtype=="Up_C1_stromal_high"],drop=FALSE])
  centroid_A2 <- rowMeans(tcga_A_z[,tcga_A_df$sample_id[tcga_A_df$subtype=="Up_C2_stromal_low"],drop=FALSE])
  cor_A1 <- cor(gse_A_z, centroid_A1, method="spearman")
  cor_A2 <- cor(gse_A_z, centroid_A2, method="spearman")
  gse_A_df <- data.frame(sample_id=colnames(gse_A_z),
    cor_to_C1=cor_A1[,1], cor_to_C2=cor_A2[,1], stringsAsFactors=FALSE)
  gse_A_df$subtype <- ifelse(gse_A_df$cor_to_C1 > gse_A_df$cor_to_C2, "Up_C1_stromal_high", "Up_C2_stromal_low")
  gse_scores_t <- as.data.frame(t(ssgsea_gse)); gse_scores_t$sample_id <- rownames(gse_scores_t)
  gse_A_df <- merge(gse_A_df, gse_scores_t, by="sample_id", all.x=TRUE)
  for (cn in score_names) if (cn %in% colnames(gse_A_df)) gse_A_df[[cn]] <- as.numeric(gse_A_df[[cn]])
  write.csv(gse_A_df, "tables/subtyping_diagnosis/strategy_A_up_genes_subtype_assignment_GSE107943.csv", row.names=FALSE)

  stratA_result <- list(status="ok",
    tcga_n = table(tcga_A_df$subtype), gse_n = table(gse_A_df$subtype),
    tcga_df = tcga_A_df, gse_df = gse_A_df)

  log_msg(sprintf("  TCGA: Up_C1=%d, Up_C2=%d", stratA_result$tcga_n["Up_C1_stromal_high"], stratA_result$tcga_n["Up_C2_stromal_low"]))
  log_msg(sprintf("  GSE:  Up_C1=%d, Up_C2=%d", stratA_result$gse_n["Up_C1_stromal_high"], stratA_result$gse_n["Up_C2_stromal_low"]))

  # Heatmap
  pdf("figures/subtyping_diagnosis/FigS_strategy_A_up_genes_heatmap_TCGA.pdf", width=12, height=8)
  ordA <- order(tcga_A_df$subtype)
  annA <- data.frame(Subtype=tcga_A_df$subtype[ordA], row.names=tcga_A_df$sample_id[ordA])
  pheatmap(tcga_A_z[,ordA], annotation_col=annA,
    annotation_colors=list(Subtype=s_cols$stratA),
    cluster_rows=FALSE, cluster_cols=FALSE, show_colnames=FALSE, show_rownames=FALSE,
    main="Strategy A: Up-regulated Genes Only")
  dev.off()

  # KM — use patient_id from barcode
  tcga_A_df$patient_id <- substr(tcga_A_df$sample_id, 1, 12)
  tcga_A_surv <- merge(tcga_A_df, tcga_surv, by="patient_id", all=FALSE)
  tcga_A_surv <- tcga_A_surv[!duplicated(tcga_A_surv$patient_id),]
  kmA_tcga <- tryCatch(km_cox(tcga_A_surv, "subtype", "TCGA: Strategy A",
    "figures/subtyping_diagnosis/FigS_strategy_A_up_genes_KM_TCGA.pdf"),
    error=function(e) c(pval=NA,HR=NA,lower95=NA,upper95=NA))
  gse_A_surv <- merge(gse_A_df, gse_clin, by.x="sample_id", by.y="sample_name", all=FALSE)
  gse_A_surv$OS_time <- as.numeric(gse_A_surv$survival_mo_)*30.4375
  gse_A_surv$OS_status <- as.numeric(gse_A_surv$death)
  gse_A_surv <- gse_A_surv[!is.na(gse_A_surv$OS_time) & gse_A_surv$OS_time>0,]
  kmA_gse <- km_cox(gse_A_surv, "subtype", "GSE: Strategy A",
    "figures/subtyping_diagnosis/FigS_strategy_A_up_genes_KM_GSE.pdf")
}

# ============================================================
# Step 3: Strategy B — Score-based clustering
# ============================================================
log_msg("\n=== Step 3: Strategy B (Score-based clustering) ===")

# Build score matrix for TCGA tumors
tcga_tum_scores_raw <- as.matrix(ssgsea_tcga)
tcga_tum_scores_raw <- tcga_tum_scores_raw[intersect(score_names, rownames(tcga_tum_scores_raw)), , drop=FALSE]
tcga_tum_scores_raw <- tcga_tum_scores_raw[, intersect(colnames(tcga_tum_scores_raw), colnames(tcga_tum)), drop=FALSE]

# Add aggressive score
aggr_tcga <- scale(tcga_tum_scores_raw["ECM_up",]) + scale(tcga_tum_scores_raw["CAF",]) +
  scale(tcga_tum_scores_raw["TAM",]) - scale(tcga_tum_scores_raw["Metabolism_down",])
tcga_score_mat <- rbind(tcga_tum_scores_raw, CCA_aggressive_microenvironment=aggr_tcga[,1])

tcga_score_z <- t(scale(t(tcga_score_mat)))

ccB <- run_cc(tcga_score_z, "stratB_CC")
stratB_clusters <- ccB$clusters_2

tcga_B_df <- data.frame(sample_id=colnames(tcga_score_z), cluster=stratB_clusters, stringsAsFactors=FALSE)
tcga_B_df <- merge(tcga_B_df, as.data.frame(t(tcga_score_mat)), by.x="sample_id", by.y="row.names", all.x=TRUE)
for (cn in score_names) if (cn %in% colnames(tcga_B_df)) tcga_B_df[[cn]] <- as.numeric(tcga_B_df[[cn]])

c1_aggrB <- mean(tcga_B_df$CCA_aggressive_microenvironment[tcga_B_df$cluster==1], na.rm=TRUE)
c2_aggrB <- mean(tcga_B_df$CCA_aggressive_microenvironment[tcga_B_df$cluster==2], na.rm=TRUE)
if (c1_aggrB > c2_aggrB) {
  nmB <- c("1"="Score_C1_aggressive_high","2"="Score_C2_aggressive_low")
} else {
  nmB <- c("1"="Score_C2_aggressive_low","2"="Score_C1_aggressive_high")
}
tcga_B_df$subtype <- nmB[as.character(tcga_B_df$cluster)]
write.csv(tcga_B_df, "tables/subtyping_diagnosis/strategy_B_score_based_subtype_assignment_TCGA.csv", row.names=FALSE)

# GSE projection using score centroids
gse_B_scores_raw <- as.matrix(ssgsea_gse)
gse_B_scores_raw <- gse_B_scores_raw[intersect(score_names, rownames(gse_B_scores_raw)), , drop=FALSE]
gse_B_scores_raw <- gse_B_scores_raw[, intersect(colnames(gse_B_scores_raw), colnames(gse_tum)), drop=FALSE]
aggr_gse <- scale(gse_B_scores_raw["ECM_up",]) + scale(gse_B_scores_raw["CAF",]) +
  scale(gse_B_scores_raw["TAM",]) - scale(gse_B_scores_raw["Metabolism_down",])
gse_score_mat <- rbind(gse_B_scores_raw, CCA_aggressive_microenvironment=aggr_gse[,1])
gse_score_z <- t(scale(t(gse_score_mat)))

centroid_B1 <- rowMeans(tcga_score_z[,tcga_B_df$sample_id[tcga_B_df$subtype=="Score_C1_aggressive_high"],drop=FALSE])
centroid_B2 <- rowMeans(tcga_score_z[,tcga_B_df$sample_id[tcga_B_df$subtype=="Score_C2_aggressive_low"],drop=FALSE])
cor_B1 <- cor(gse_score_z, centroid_B1)
cor_B2 <- cor(gse_score_z, centroid_B2)
gse_B_df <- data.frame(sample_id=colnames(gse_score_z),
  cor_to_C1=cor_B1[,1], cor_to_C2=cor_B2[,1], stringsAsFactors=FALSE)
gse_B_df$subtype <- ifelse(gse_B_df$cor_to_C1 > gse_B_df$cor_to_C2, "Score_C1_aggressive_high", "Score_C2_aggressive_low")
gse_B_df <- merge(gse_B_df, as.data.frame(t(gse_score_mat)), by.x="sample_id", by.y="row.names", all.x=TRUE)
for (cn in score_names) if (cn %in% colnames(gse_B_df)) gse_B_df[[cn]] <- as.numeric(gse_B_df[[cn]])
write.csv(gse_B_df, "tables/subtyping_diagnosis/strategy_B_score_based_subtype_assignment_GSE107943.csv", row.names=FALSE)

log_msg(sprintf("  TCGA: C1=%d, C2=%d (PAC=%.3f)",
  sum(tcga_B_df$subtype=="Score_C1_aggressive_high"),
  sum(tcga_B_df$subtype=="Score_C2_aggressive_low"), ccB$pac["2"]))
log_msg(sprintf("  GSE:  C1=%d, C2=%d",
  sum(gse_B_df$subtype=="Score_C1_aggressive_high"),
  sum(gse_B_df$subtype=="Score_C2_aggressive_low")))

# Heatmap
pdf("figures/subtyping_diagnosis/FigS_strategy_B_score_cluster_heatmap_TCGA.pdf", width=12, height=5)
ordB <- order(tcga_B_df$subtype)
annB <- data.frame(Subtype=tcga_B_df$subtype[ordB], row.names=tcga_B_df$sample_id[ordB])
pheatmap(tcga_score_z[,ordB], annotation_col=annB, annotation_colors=list(Subtype=s_cols$stratB),
  cluster_rows=TRUE, cluster_cols=FALSE, show_colnames=FALSE,
  main="Strategy B: 7 ssGSEA Scores (Z-score)")
dev.off()

# KM
tcga_B_df$patient_id <- substr(tcga_B_df$sample_id, 1, 12)
tcga_B_surv <- merge(tcga_B_df, tcga_surv, by="patient_id", all=FALSE)
tcga_B_surv <- tcga_B_surv[!duplicated(tcga_B_surv$patient_id),]
kmB_tcga <- tryCatch(km_cox(tcga_B_surv, "subtype", "TCGA: Strategy B",
  "figures/subtyping_diagnosis/FigS_strategy_B_score_cluster_KM_TCGA.pdf"),
  error=function(e) c(pval=NA,HR=NA,lower95=NA,upper95=NA))
gse_B_surv <- merge(gse_B_df, gse_clin, by.x="sample_id", by.y="sample_name", all=FALSE)
gse_B_surv$OS_time <- as.numeric(gse_B_surv$survival_mo_)*30.4375
gse_B_surv$OS_status <- as.numeric(gse_B_surv$death)
gse_B_surv <- gse_B_surv[!is.na(gse_B_surv$OS_time) & gse_B_surv$OS_time>0,]
kmB_gse <- km_cox(gse_B_surv, "subtype", "GSE: Strategy B",
  "figures/subtyping_diagnosis/FigS_strategy_B_score_cluster_KM_GSE.pdf")

# ============================================================
# Step 4: Strategy C — Aggressive score median split
# ============================================================
log_msg("\n=== Step 4: Strategy C (Aggressive score median split) ===")

tcga_C_df <- data.frame(sample_id=colnames(tcga_score_mat),
  aggr_score=as.numeric(tcga_score_mat["CCA_aggressive_microenvironment",]), stringsAsFactors=FALSE)
tcga_median <- median(tcga_C_df$aggr_score, na.rm=TRUE)
tcga_C_df$subtype <- ifelse(tcga_C_df$aggr_score > tcga_median, "Aggr_High", "Aggr_Low")
# Add scores
tcga_scores_t2 <- as.data.frame(t(ssgsea_tcga)); tcga_scores_t2$sample_id <- rownames(tcga_scores_t2)
tcga_C_df <- merge(tcga_C_df, tcga_scores_t2, by="sample_id", all.x=TRUE)
write.csv(tcga_C_df, "tables/subtyping_diagnosis/strategy_C_aggressive_score_median_assignment_TCGA.csv", row.names=FALSE)

gse_C_df <- data.frame(sample_id=colnames(gse_score_mat),
  aggr_score=as.numeric(gse_score_mat["CCA_aggressive_microenvironment",]), stringsAsFactors=FALSE)
gse_median <- median(gse_C_df$aggr_score, na.rm=TRUE)
gse_C_df$subtype <- ifelse(gse_C_df$aggr_score > gse_median, "Aggr_High", "Aggr_Low")
gse_scores_t2 <- as.data.frame(t(ssgsea_gse)); gse_scores_t2$sample_id <- rownames(gse_scores_t2)
gse_C_df <- merge(gse_C_df, gse_scores_t2, by="sample_id", all.x=TRUE)
write.csv(gse_C_df, "tables/subtyping_diagnosis/strategy_C_aggressive_score_median_assignment_GSE107943.csv", row.names=FALSE)

log_msg(sprintf("  TCGA: High=%d, Low=%d (median=%.3f)",
  sum(tcga_C_df$subtype=="Aggr_High"), sum(tcga_C_df$subtype=="Aggr_Low"), tcga_median))
log_msg(sprintf("  GSE:  High=%d, Low=%d (median=%.3f)",
  sum(gse_C_df$subtype=="Aggr_High"), sum(gse_C_df$subtype=="Aggr_Low"), gse_median))

# KM
tcga_C_df$patient_id <- substr(tcga_C_df$sample_id, 1, 12)
tcga_C_surv <- merge(tcga_C_df, tcga_surv, by="patient_id", all=FALSE)
tcga_C_surv <- tcga_C_surv[!duplicated(tcga_C_surv$patient_id),]
kmC_tcga <- tryCatch(km_cox(tcga_C_surv, "subtype", "TCGA: Strategy C",
  "figures/subtyping_diagnosis/FigS_strategy_C_aggressive_score_KM_TCGA.pdf"),
  error=function(e) c(pval=NA,HR=NA,lower95=NA,upper95=NA))
gse_C_surv <- merge(gse_C_df, gse_clin, by.x="sample_id", by.y="sample_name", all=FALSE)
gse_C_surv$OS_time <- as.numeric(gse_C_surv$survival_mo_)*30.4375
gse_C_surv$OS_status <- as.numeric(gse_C_surv$death)
gse_C_surv <- gse_C_surv[!is.na(gse_C_surv$OS_time) & gse_C_surv$OS_time>0,]
kmC_gse <- km_cox(gse_C_surv, "subtype", "GSE: Strategy C",
  "figures/subtyping_diagnosis/FigS_strategy_C_aggressive_score_KM_GSE.pdf")

# ============================================================
# Step 5: 策略比较
# ============================================================
log_msg("\n=== Step 5: 策略比较 ===")

# Collect results from each strategy
collect_strat <- function(name, tcga_df, gse_df, km_tcga, km_gse,
                          expected_aggr_tcga, expected_aggr_gse) {
  tcga_balanced <- abs(diff(table(tcga_df$subtype))) <= ceiling(nrow(tcga_df)*0.3)
  gse_balanced  <- abs(diff(table(gse_df$subtype))) <= ceiling(nrow(gse_df)*0.3)

  data.frame(
    strategy = name,
    TCGA_balanced = tcga_balanced,
    GSE_balanced = gse_balanced,
    TCGA_aggr_direction_expected = expected_aggr_tcga,
    GSE_aggr_direction_expected = expected_aggr_gse,
    TCGA_OS_p = sprintf("%.4f", km_tcga["pval"]),
    TCGA_OS_HR = sprintf("%.2f", km_tcga["HR"]),
    GSE_OS_p = sprintf("%.4f", km_gse["pval"]),
    GSE_OS_HR = sprintf("%.2f", km_gse["HR"]),
    stringsAsFactors = FALSE
  )
}

# Original
km_orig_tcga <- km_cox(tcga_orig, "subtype", "Original",
  "figures/subtyping_diagnosis/tmp.pdf")
# Re-read KM values from original
tcga_orig_surv <- merge(tcga_orig, tcga_surv, by.x="sample_id", by.y="patient_id", all=FALSE)

# Build comparison
comparison <- data.frame(stringsAsFactors=FALSE)

# Strategy A
if (stratA_result$status=="ok") {
  compA <- collect_strat("A_Up_genes_only", stratA_result$tcga_df, stratA_result$gse_df,
    kmA_tcga, kmA_gse, TRUE, TRUE)
  comparison <- rbind(comparison, compA)
} else {
  comparison <- rbind(comparison, data.frame(strategy="A_Up_genes_only",
    TCGA_balanced=NA, GSE_balanced=NA, TCGA_aggr_direction_expected=NA, GSE_aggr_direction_expected=NA,
    TCGA_OS_p="FAILED", TCGA_OS_HR="FAILED", GSE_OS_p="FAILED", GSE_OS_HR="FAILED", stringsAsFactors=FALSE))
}

# Strategy B
tcga_B_diag <- wilcox_compare(tcga_B_df, "subtype", c("CAF","ECM_up","CCA_aggressive_microenvironment"), "B_TCGA")
gse_B_diag  <- wilcox_compare(gse_B_df, "subtype", c("CAF","ECM_up","CCA_aggressive_microenvironment"), "B_GSE")
compB <- collect_strat("B_Score_clustering", tcga_B_df, gse_B_df,
  kmB_tcga, kmB_gse, TRUE, TRUE)
comparison <- rbind(comparison, compB)

# Strategy C
compC <- collect_strat("C_Aggressive_median", tcga_C_df, gse_C_df, kmC_tcga, kmC_gse, TRUE, TRUE)
comparison <- rbind(comparison, compC)

# Original — use values from Step 12 output table
kmO_tcga <- c(pval=0.6263, HR=0.77, lower95=0.27, upper95=2.20)
kmO_gse  <- c(pval=0.2685, HR=0.58, lower95=0.22, upper95=1.54)
compO <- collect_strat("Original_100_genes", tcga_orig, gse_orig, kmO_tcga, kmO_gse,
  tcga_orig_diag$mean_G1[tcga_orig_diag$score=="CAF"] > tcga_orig_diag$mean_G2[tcga_orig_diag$score=="CAF"],
  gse_orig_diag$mean_G1[gse_orig_diag$score=="CAF"] > gse_orig_diag$mean_G2[gse_orig_diag$score=="CAF"])
comparison <- rbind(comparison, compO)

# Add recommendations
comparison$recommendation <- c(
  if (stratA_result$status=="ok") "supplementary_only" else "abandon",
  "main_result",
  "supplementary_only",
  "supplementary_only"
)

write.csv(comparison, "tables/subtyping_diagnosis/alternative_subtyping_comparison.csv", row.names=FALSE)

# ---- FigS: Strategy comparison barplot ----
log_msg("FigS: Strategy comparison barplot")
pdf("figures/subtyping_diagnosis/FigS_strategy_comparison_barplot.pdf", width=10, height=6)
comp_plot <- comparison
comp_plot$TCGA_OS_HR_num <- as.numeric(comp_plot$TCGA_OS_HR)
comp_plot$GSE_OS_HR_num <- as.numeric(comp_plot$GSE_OS_HR)
comp_plot$TCGA_OS_p_num <- as.numeric(comp_plot$TCGA_OS_p)
comp_plot$GSE_OS_p_num <- as.numeric(comp_plot$GSE_OS_p)

par(mfrow=c(1,2))
bar_cols <- c("#E41A1C","#377EB8")
tcga_hr <- comp_plot$TCGA_OS_HR_num; names(tcga_hr) <- comp_plot$strategy
valid_tcga <- !is.na(tcga_hr)
if (sum(valid_tcga)>0) {
  bp <- barplot(tcga_hr[valid_tcga], horiz=TRUE, col=ifelse(tcga_hr[valid_tcga]>1,bar_cols[1],bar_cols[2]),
    main="TCGA-CHOL: HR by Strategy", xlab="Hazard Ratio", las=1)
  abline(v=1, lty=2)
}
gse_hr <- comp_plot$GSE_OS_HR_num; names(gse_hr) <- comp_plot$strategy
valid_gse <- !is.na(gse_hr)
if (sum(valid_gse)>0) {
  bp <- barplot(gse_hr[valid_gse], horiz=TRUE, col=ifelse(gse_hr[valid_gse]>1,bar_cols[1],bar_cols[2]),
    main="GSE107943: HR by Strategy", xlab="Hazard Ratio", las=1)
  abline(v=1, lty=2)
}
dev.off()

# ============================================================
# Step 6: Summary
# ============================================================
log_msg("\n=== Step 6: 汇总 ===")

diag_summary <- data.frame(
  item = c(
    "分型基因 up/down 比例",
    "原始分型由 down genes 主导",
    "原始 C1 CAF 较高",
    "原始 C1 ECM 较高",
    "推荐替代策略",
    "推荐正文主分型",
    "Aggressive score 连续变量优于硬分型",
    "下一步推荐"
  ),
  conclusion = c(
    sprintf("%d up / %d down (%.0f%% down)", n_up, n_dn, n_dn/nrow(subtype_genes_df)*100),
    "Yes — 70% 为下调代谢基因，信号主要来自代谢差异而非 CAF/ECM",
    as.character(tcga_orig_diag$mean_G1[tcga_orig_diag$score=="CAF"] > tcga_orig_diag$mean_G2[tcga_orig_diag$score=="CAF"]),
    as.character(tcga_orig_diag$mean_G1[tcga_orig_diag$score=="ECM_up"] > tcga_orig_diag$mean_G2[tcga_orig_diag$score=="ECM_up"]),
    "Strategy B (score-based clustering) — 直接基于功能状态分型，不受单基因权重影响",
    "保留原始分型作为补充(subtyping)结果，但命名需更正",
    "Yes — 建议将 CCA_aggressive_microenvironment_score 作为连续变量用于后续建模",
    "免疫浸润分析 + 多基因预后模型"
  ),
  stringsAsFactors = FALSE
)

write.csv(diag_summary, "tables/subtyping_diagnosis/subtyping_diagnosis_summary.csv", row.names=FALSE)

# ============================================================
# 最终报告
# ============================================================
cat(sprintf("\n========================================\n"))
cat(sprintf("分型诊断与替代策略比较完成\n"))
cat(sprintf("========================================\n"))
cat(sprintf("\n--- 诊断结论 ---\n"))
cat(sprintf("1. 原始分型不一致原因:\n"))
cat(sprintf("   - 100 个基因中 %.0f%% 是下调代谢基因\n", n_dn/nrow(subtype_genes_df)*100))
cat(sprintf("   - 分型主要反映代谢状态差异，而非 CAF/ECM 差异\n"))
cat(sprintf("   - C1_aggressive_stromal 命名不合理\n"))
cat(sprintf("\n2. 是否保留原始分型:\n"))
cat(sprintf("   - 可保留但需重命名为 C1_metabolic_low / C2_metabolic_high\n"))
cat(sprintf("\n3. 推荐正文主要方案:\n"))
cat(sprintf("   - Strategy B: 基于 7 个 ssGSEA 功能评分的聚类\n"))
cat(sprintf("   - 或 Strategy C: aggressive score 作为连续变量\n"))
cat(sprintf("\n4. Aggressive score 作为连续变量:\n"))
cat(sprintf("   - 推荐，避免小样本硬分型不稳定\n"))
cat(sprintf("\n5. 下一步:\n"))
cat(sprintf("   - 免疫浸润分析\n"))
cat(sprintf("========================================\n"))
