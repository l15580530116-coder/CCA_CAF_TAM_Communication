# ============================================================
# 分析目的: 基于 validated IM/CAF/TAM DEGs 和 ssGSEA scores
#          进行胆管癌分子亚型分型 (TCGA discovery → GSE projection)
# 输入文件:
#   tables/CCA_IM_CAF_TAM_DEGs_validated_same_direction.csv
#   data/processed/TCGA_CHOL/TCGA_CHOL_tumor_expression_for_survival.csv
#   data/processed/TCGA_CHOL/TCGA_CHOL_matched_survival_metadata.csv
#   data/processed/GSE107943/GSE107943_log2CPM_gene_symbol.csv
#   data/processed/GSE107943/GSE107943_sample_metadata.csv
#   tables/gsva/TCGA_CHOL_pathway_scores_all_samples.csv
#   tables/gsva/GSE107943_pathway_scores_all_samples.csv
# 输出文件:
#   tables/subtyping/ (8 个表格)
#   figures/subtyping/ (9 张图)
# 主要方法: ConsensusClusterPlus (k=2-6) → nearest centroid
#          projection (GSE107943) → subtype comparison + KM survival
# ============================================================

library(ConsensusClusterPlus)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(pheatmap)
library(RColorBrewer)
library(data.table)

dir.create("tables/subtyping", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/subtyping", showWarnings = FALSE, recursive = TRUE)

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

subtype_cols <- c("C1_aggressive_stromal" = "#E41A1C",
                  "C2_metabolic_preserved" = "#377EB8")

score_names <- c("CAF", "ECM_up", "TAM", "Metabolism_down", "IM_total",
                 "IM_CAF_TAM_combined", "CCA_aggressive_microenvironment")

# ============================================================
# Step 1: 准备分型基因
# ============================================================
log_msg("=== Step 1: 准备分型基因 ===")

validated <- read.csv("tables/CCA_IM_CAF_TAM_DEGs_validated_same_direction.csv", stringsAsFactors=FALSE)
tcga_tum <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_tumor_expression_for_survival.csv",
                     row.names=1, check.names=FALSE)
gse_full <- read.csv("data/processed/GSE107943/GSE107943_log2CPM_gene_symbol.csv",
                     row.names=1, check.names=FALSE)
gse_full <- gse_full[, colnames(gse_full) != "6N", drop=FALSE]
gse_meta <- read.csv("data/processed/GSE107943/GSE107943_sample_metadata.csv", stringsAsFactors=FALSE)
gse_tum <- gse_full[, intersect(colnames(gse_full), gse_meta$sample_id[gse_meta$sample_type == "Tumor"]), drop=FALSE]

# 共同基因
common_genes <- intersect(validated$gene_symbol, intersect(rownames(tcga_tum), rownames(gse_tum)))
log_msg(sprintf("Common validated genes: %d", length(common_genes)))

validated_common <- validated[validated$gene_symbol %in% common_genes, ]
up_genes <- validated_common$gene_symbol[validated_common$direction_TCGA == "up"]
dn_genes <- validated_common$gene_symbol[validated_common$direction_TCGA == "down"]
log_msg(sprintf("  Up: %d, Down: %d", length(up_genes), length(dn_genes)))

# Top variable genes, balanced up/down
tcga_var <- apply(tcga_tum[common_genes, , drop=FALSE], 1, var, na.rm=TRUE)
up_var <- sort(tcga_var[up_genes], decreasing=TRUE)
dn_var <- sort(tcga_var[dn_genes], decreasing=TRUE)

n_up_keep <- min(30, length(up_var))
n_dn_keep <- min(70, length(dn_var))
top_up <- names(up_var)[1:n_up_keep]
top_dn <- names(dn_var)[1:n_dn_keep]
subtype_genes <- c(top_up, top_dn)

log_msg(sprintf("Subtyping genes: %d (up=%d, down=%d)", length(subtype_genes), n_up_keep, n_dn_keep))

write.csv(data.frame(gene=subtype_genes,
                     direction=ifelse(subtype_genes %in% top_up, "up", "down"),
                     stringsAsFactors=FALSE),
          "tables/subtyping/subtyping_gene_list.csv", row.names=FALSE)

# ============================================================
# Step 2: TCGA ConsensusClusterPlus
# ============================================================
log_msg("\n=== Step 2: TCGA ConsensusClusterPlus ===")

tcga_z <- as.matrix(tcga_tum[subtype_genes, , drop=FALSE])
tcga_z <- t(scale(t(tcga_z)))

# 使用 ConsensusClusterPlus
cc_results <- ConsensusClusterPlus(
  d = tcga_z,
  maxK = 6,
  reps = 100,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "hc",
  distance = "euclidean",
  seed = 42,
  plot = NULL,
  title = "figures/subtyping/TCGA_CC",
  writeTable = FALSE
)

log_msg("ConsensusClusterPlus 完成")
for (k in 2:6) {
  if (!is.null(cc_results[[k]])) {
    icl <- cc_results[[k]][["consensusMatrix"]]
    log_msg(sprintf("  k=%d: PAC=%.3f", k,
      if (!is.null(attr(cc_results[[k]],"PAC"))) attr(cc_results[[k]],"PAC") else NA))
  }
}

# 评估 k 的选择：使用 PAC (proportion of ambiguous clustering) 或 delta area
icl <- sapply(2:6, function(k) {
  if (!is.null(cc_results[[k]])) {
    cm <- cc_results[[k]]$consensusMatrix
    sum(cm > 0.2 & cm < 0.8) / (nrow(cm)^2)  # simplified PAC
  } else NA
})
names(icl) <- 2:6
log_msg(sprintf("PAC by k: %s", paste(sprintf("k%d=%.3f", 2:6, icl), collapse=", ")))

# 选择 k=2
k_opt <- 2
tcga_clusters <- cc_results[[k_opt]]$consensusClass
log_msg(sprintf("Selected k=%d", k_opt))

# ---- Fig5A: CDF plot ----
log_msg("Fig5A: Consensus CDF")
pdf("figures/subtyping/Fig5A_TCGA_consensus_CDF.pdf", width=7, height=5)
# Plot CDFs for each k
cols <- brewer.pal(6, "Set1")[-6]
plot(ecdf(cc_results[[2]]$consensusMatrix[lower.tri(cc_results[[2]]$consensusMatrix)]),
     col=cols[1], lwd=2, main="Consensus CDF", xlab="Consensus Index", ylab="CDF",
     xlim=c(0,1), do.points=FALSE, verticals=TRUE)
for (k in 3:6) {
  if (!is.null(cc_results[[k]])) {
    plot(ecdf(cc_results[[k]]$consensusMatrix[lower.tri(cc_results[[k]]$consensusMatrix)]),
         col=cols[k-1], lwd=2, do.points=FALSE, verticals=TRUE, add=TRUE)
  }
}
legend("bottomright", legend=paste0("k=", 2:6), col=cols[1:5], lwd=2, bty="n")
dev.off()
log_msg("  Done")

# ---- Fig5B: Consensus matrix k=2 ----
log_msg("Fig5B: Consensus matrix k=2")
cm2 <- cc_results[[2]]$consensusMatrix
ord <- order(tcga_clusters)
cm2_ord <- cm2[ord, ord]

pdf("figures/subtyping/Fig5B_TCGA_consensus_matrix_k2.pdf", width=8, height=7)
# Simplified consensus matrix plot using base heatmap
heatmap(cm2_ord, Colv=NA, Rowv=NA, scale="none",
        col=colorRampPalette(c("white", "steelblue", "darkblue"))(50),
        main="TCGA-CHOL: Consensus Matrix (k=2)", labRow=FALSE, labCol=FALSE,
        ColSideColors=ifelse(tcga_clusters[ord]==1, "pink", "lightblue"),
        RowSideColors=ifelse(tcga_clusters[ord]==1, "pink", "lightblue"))
legend("topright", legend=c("Cluster 1", "Cluster 2"), fill=c("pink", "lightblue"), bty="n")
dev.off()
log_msg("  Done")

# ---- Assign subtype names based on pathway scores ----
ssgsea_tcga_raw <- read.csv("tables/gsva/TCGA_CHOL_pathway_scores_all_samples.csv", check.names=FALSE)
rownames(ssgsea_tcga_raw) <- ssgsea_tcga_raw$pathway
ssgsea_tcga_raw$pathway <- NULL
# Convert all columns to numeric
ssgsea_tcga <- as.data.frame(sapply(ssgsea_tcga_raw, function(x) as.numeric(as.character(x))))
rownames(ssgsea_tcga) <- rownames(ssgsea_tcga_raw)
ssgsea_tcga_t <- as.data.frame(t(ssgsea_tcga))
ssgsea_tcga_t$sample_id <- rownames(ssgsea_tcga_t)

tcga_subtype <- data.frame(
  sample_id = colnames(tcga_z),
  cluster = tcga_clusters,
  stringsAsFactors = FALSE
)
tcga_subtype <- merge(tcga_subtype, ssgsea_tcga_t, by="sample_id", all.x=TRUE)

# Use CAF + ECM_up to determine which cluster is "aggressive stromal"
c1_caf <- mean(as.numeric(tcga_subtype$CAF[tcga_subtype$cluster == 1]), na.rm=TRUE)
c2_caf <- mean(as.numeric(tcga_subtype$CAF[tcga_subtype$cluster == 2]), na.rm=TRUE)
c1_ecm <- mean(as.numeric(tcga_subtype$ECM_up[tcga_subtype$cluster == 1]), na.rm=TRUE)
c2_ecm <- mean(as.numeric(tcga_subtype$ECM_up[tcga_subtype$cluster == 2]), na.rm=TRUE)
log_msg(sprintf("Cluster 1: CAF=%.3f, ECM=%.3f", c1_caf, c1_ecm))
log_msg(sprintf("Cluster 2: CAF=%.3f, ECM=%.3f", c2_caf, c2_ecm))

if (c1_caf > c2_caf) {
  cluster_name <- c("1"="C1_aggressive_stromal", "2"="C2_metabolic_preserved")
} else {
  cluster_name <- c("1"="C2_metabolic_preserved", "2"="C1_aggressive_stromal")
}
tcga_subtype$subtype <- cluster_name[as.character(tcga_subtype$cluster)]

# Ensure all score columns are numeric
for (cn in score_names) {
  if (cn %in% colnames(tcga_subtype)) tcga_subtype[[cn]] <- as.numeric(tcga_subtype[[cn]])
}

c1_n <- sum(tcga_subtype$subtype == "C1_aggressive_stromal")
c2_n <- sum(tcga_subtype$subtype == "C2_metabolic_preserved")
log_msg(sprintf("TCGA: C1=%d, C2=%d", c1_n, c2_n))

write.csv(tcga_subtype, "tables/subtyping/TCGA_subtype_assignment.csv", row.names=FALSE)
log_msg("已保存 TCGA subtype")

# ============================================================
# Step 3: GSE107943 Nearest Centroid Projection
# ============================================================
log_msg("\n=== Step 3: GSE107943 nearest centroid projection ===")

# centroids
tcga_c1_mat <- tcga_z[, tcga_subtype$sample_id[tcga_subtype$subtype == "C1_aggressive_stromal"], drop=FALSE]
tcga_c2_mat <- tcga_z[, tcga_subtype$sample_id[tcga_subtype$subtype == "C2_metabolic_preserved"], drop=FALSE]
centroid_c1 <- rowMeans(tcga_c1_mat)
centroid_c2 <- rowMeans(tcga_c2_mat)

# GSE z-score
gse_z <- as.matrix(gse_tum[subtype_genes, , drop=FALSE])
gse_z <- t(scale(t(gse_z)))

cor_c1 <- cor(gse_z, centroid_c1, method="spearman")
cor_c2 <- cor(gse_z, centroid_c2, method="spearman")

gse_subtype <- data.frame(
  sample_id = colnames(gse_z),
  cor_to_C1 = cor_c1[,1],
  cor_to_C2 = cor_c2[,1],
  stringsAsFactors = FALSE
)
gse_subtype$subtype <- ifelse(gse_subtype$cor_to_C1 > gse_subtype$cor_to_C2,
                               "C1_aggressive_stromal", "C2_metabolic_preserved")

# Add scores
ssgsea_gse_raw <- read.csv("tables/gsva/GSE107943_pathway_scores_all_samples.csv", check.names=FALSE)
rownames(ssgsea_gse_raw) <- ssgsea_gse_raw$pathway
ssgsea_gse_raw$pathway <- NULL
ssgsea_gse <- as.data.frame(sapply(ssgsea_gse_raw, function(x) as.numeric(as.character(x))))
rownames(ssgsea_gse) <- rownames(ssgsea_gse_raw)
ssgsea_gse_t <- as.data.frame(t(ssgsea_gse))
ssgsea_gse_t$sample_id <- rownames(ssgsea_gse_t)
gse_subtype <- merge(gse_subtype, ssgsea_gse_t, by="sample_id", all.x=TRUE)
for (cn in score_names) {
  if (cn %in% colnames(gse_subtype)) gse_subtype[[cn]] <- as.numeric(gse_subtype[[cn]])
}

gse_c1_n <- sum(gse_subtype$subtype == "C1_aggressive_stromal")
gse_c2_n <- sum(gse_subtype$subtype == "C2_metabolic_preserved")
log_msg(sprintf("GSE107943: C1=%d, C2=%d", gse_c1_n, gse_c2_n))

write.csv(gse_subtype, "tables/subtyping/GSE107943_subtype_assignment.csv", row.names=FALSE)
log_msg("已保存 GSE subtype")

# ============================================================
# Step 4: 亚型评分比较
# ============================================================
log_msg("\n=== Step 4: 亚型评分比较 ===")

compare_subtypes <- function(df, label) {
  results <- data.frame(stringsAsFactors=FALSE)
  c1 <- df[df$subtype == "C1_aggressive_stromal", , drop=FALSE]
  c2 <- df[df$subtype == "C2_metabolic_preserved", , drop=FALSE]
  for (sc in score_names) {
    if (!sc %in% colnames(df)) next
    wt <- wilcox.test(c1[[sc]], c2[[sc]])
    results <- rbind(results, data.frame(
      score=sc, mean_C1=mean(c1[[sc]],na.rm=TRUE), mean_C2=mean(c2[[sc]],na.rm=TRUE),
      pvalue=wt$p.value, stringsAsFactors=FALSE))
  }
  results$FDR <- p.adjust(results$pvalue, method="BH")
  log_msg(sprintf("%s: %d scores, %d FDR<0.05", label, nrow(results), sum(results$FDR<0.05)))
  return(results)
}

tcga_comp_sc <- compare_subtypes(tcga_subtype, "TCGA")
gse_comp_sc  <- compare_subtypes(gse_subtype, "GSE107943")

write.csv(tcga_comp_sc, "tables/subtyping/TCGA_subtype_score_comparison.csv", row.names=FALSE)
write.csv(gse_comp_sc,  "tables/subtyping/GSE107943_subtype_score_comparison.csv", row.names=FALSE)

# ============================================================
# Step 5: 生存分析
# ============================================================
log_msg("\n=== Step 5: 生存分析 ===")

# TCGA
tcga_surv <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_matched_survival_metadata.csv", stringsAsFactors=FALSE)
tcga_meta <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_sample_metadata.csv", stringsAsFactors=FALSE)
tcga_meta_map <- setNames(tcga_meta$patient_id, tcga_meta$sample_id)
tcga_subtype$patient_id <- tcga_meta_map[tcga_subtype$sample_id]
tcga_surv_sub <- merge(tcga_subtype, tcga_surv, by="patient_id", all=FALSE)
tcga_surv_sub <- tcga_surv_sub[!duplicated(tcga_surv_sub$patient_id), ]

# GSE
gse_clin <- read.csv("data/clinical/GSE107943/GSE107943_clinical.csv", stringsAsFactors=FALSE)
gse_surv_sub <- merge(gse_subtype, gse_clin, by.x="sample_id", by.y="sample_name", all=FALSE)
gse_surv_sub$OS_time <- as.numeric(gse_surv_sub$survival_mo_) * 30.4375
gse_surv_sub$OS_status <- as.numeric(gse_surv_sub$death)
gse_surv_sub <- gse_surv_sub[!is.na(gse_surv_sub$OS_time) & gse_surv_sub$OS_time > 0, ]

km_subtype <- function(df, title, filename) {
  sdiff <- survdiff(Surv(OS_time, OS_status) ~ subtype, data=df)
  pval <- 1 - pchisq(sdiff$chisq, 1)
  cox <- coxph(Surv(OS_time, OS_status) ~ subtype, data=df)
  hr <- exp(coef(cox))
  ci <- exp(confint(cox))
  df$OS_months <- df$OS_time / 30.4375

  cat(sprintf("  %s: p=%.4f, HR=%.2f (%.2f-%.2f), C1=%d, C2=%d, events=%d\n",
    title, pval, hr, ci[1], ci[2],
    sum(df$subtype=="C1_aggressive_stromal"),
    sum(df$subtype=="C2_metabolic_preserved"),
    sum(df$OS_status==1)))

  pdf(filename, width=7, height=6)
  fit2 <- survfit(Surv(OS_months, OS_status) ~ subtype, data=df)
  p <- ggsurvplot(fit2, data=df, pval=TRUE, palette=subtype_cols,
    legend.title="Subtype", xlab="Overall Survival (months)", title=title)
  print(p); dev.off()
  c(pval=pval, HR=hr, lower95=ci[1], upper95=ci[2])
}

km_tcga <- km_subtype(tcga_surv_sub, "TCGA-CHOL", "figures/subtyping/Fig5E_TCGA_subtype_KM_OS.pdf")
km_gse  <- km_subtype(gse_surv_sub, "GSE107943", "figures/subtyping/Fig5H_GSE107943_subtype_KM_OS.pdf")

# Survival summary tables
make_surv_sum <- function(df) {
  data.frame(
    subtype = names(table(df$subtype)),
    n = as.integer(table(df$subtype)),
    events = as.integer(tapply(df$OS_status, df$subtype, sum)),
    median_OS_mo = round(as.numeric(tapply(df$OS_time/30.4375, df$subtype, median)), 1),
    stringsAsFactors = FALSE
  )
}
write.csv(make_surv_sum(tcga_surv_sub), "tables/subtyping/TCGA_subtype_survival_summary.csv", row.names=FALSE)
write.csv(make_surv_sum(gse_surv_sub),  "tables/subtyping/GSE107943_subtype_survival_summary.csv", row.names=FALSE)

# ============================================================
# Step 6: 可视化
# ============================================================
log_msg("\n=== Step 6: 可视化 ===")

# ---- Fig5C: TCGA subtype heatmap ----
log_msg("Fig5C: TCGA heatmap")
tcga_heat_z <- tcga_z
tcga_heat_z[tcga_heat_z > 3] <- 3; tcga_heat_z[tcga_heat_z < -3] <- -3
ord_tcga <- order(tcga_subtype$subtype)
tcga_heat_z <- tcga_heat_z[, tcga_subtype$sample_id[ord_tcga], drop=FALSE]

ann_tcga <- data.frame(
  Subtype = tcga_subtype$subtype[ord_tcga],
  row.names = tcga_subtype$sample_id[ord_tcga]
)

pdf("figures/subtyping/Fig5C_TCGA_subtype_heatmap.pdf", width=12, height=8)
pheatmap(tcga_heat_z, annotation_col=ann_tcga,
         annotation_colors=list(Subtype=subtype_cols),
         cluster_rows=FALSE, cluster_cols=FALSE, show_colnames=FALSE,
         show_rownames=(nrow(tcga_heat_z)<=100),
         main="TCGA-CHOL: Subtyping Genes (Z-score)", fontsize_row=6)
dev.off()
log_msg("  Done")

# ---- Fig5D: TCGA score boxplot ----
log_msg("Fig5D: TCGA boxplot")
tcga_score_cols <- intersect(score_names, colnames(tcga_subtype))
tcga_long <- reshape2::melt(tcga_subtype, id.vars=c("sample_id","subtype"),
                            measure.vars=tcga_score_cols, variable.name="score", value.name="value")
pdf("figures/subtyping/Fig5D_TCGA_subtype_score_boxplot.pdf", width=12, height=7)
p <- ggplot(tcga_long, aes(x=score, y=value, fill=subtype)) +
  geom_boxplot(outlier.size=0.8) + scale_fill_manual(values=subtype_cols) +
  labs(y="ssGSEA Score", title="TCGA-CHOL: Subtype Score Comparison") +
  theme_pubr(base_size=12) + theme(axis.text.x=element_text(angle=45, hjust=1))
# Add p-values
tcga_pvals <- tcga_comp_sc[, c("score","FDR")]
print(p); dev.off()
log_msg("  Done")

# ---- Fig5F: GSE heatmap ----
log_msg("Fig5F: GSE heatmap")
gse_heat_z <- gse_z
gse_heat_z[gse_heat_z > 3] <- 3; gse_heat_z[gse_heat_z < -3] <- -3
ord_gse <- order(gse_subtype$subtype)
gse_heat_z <- gse_heat_z[, gse_subtype$sample_id[ord_gse], drop=FALSE]

ann_gse <- data.frame(
  Subtype = gse_subtype$subtype[ord_gse],
  row.names = gse_subtype$sample_id[ord_gse]
)

pdf("figures/subtyping/Fig5F_GSE107943_subtype_heatmap.pdf", width=12, height=8)
pheatmap(gse_heat_z, annotation_col=ann_gse,
         annotation_colors=list(Subtype=subtype_cols),
         cluster_rows=FALSE, cluster_cols=FALSE, show_colnames=FALSE,
         show_rownames=(nrow(gse_heat_z)<=100),
         main="GSE107943: Subtyping Genes (Z-score)", fontsize_row=6)
dev.off()
log_msg("  Done")

# ---- Fig5G: GSE boxplot ----
log_msg("Fig5G: GSE boxplot")
gse_score_cols <- intersect(score_names, colnames(gse_subtype))
gse_long <- reshape2::melt(gse_subtype, id.vars=c("sample_id","subtype"),
                           measure.vars=gse_score_cols, variable.name="score", value.name="value")
pdf("figures/subtyping/Fig5G_GSE107943_subtype_score_boxplot.pdf", width=12, height=7)
p <- ggplot(gse_long, aes(x=score, y=value, fill=subtype)) +
  geom_boxplot(outlier.size=0.8) + scale_fill_manual(values=subtype_cols) +
  labs(y="ssGSEA Score", title="GSE107943: Subtype Score Comparison") +
  theme_pubr(base_size=12) + theme(axis.text.x=element_text(angle=45, hjust=1))
print(p); dev.off()
log_msg("  Done")

# ---- Fig5I: Consistency barplot ----
log_msg("Fig5I: Signature consistency")
tcga_means <- aggregate(tcga_long$value, by=list(score=tcga_long$score, subtype=tcga_long$subtype), mean)
names(tcga_means)[3] <- "mean_score"
gse_means <- aggregate(gse_long$value, by=list(score=gse_long$score, subtype=gse_long$subtype), mean)
names(gse_means)[3] <- "mean_score"

p1 <- ggplot(tcga_means, aes(x=score, y=mean_score, fill=subtype)) +
  geom_bar(stat="identity", position="dodge") + scale_fill_manual(values=subtype_cols) +
  labs(title="TCGA-CHOL", y="Mean ssGSEA") + theme_pubr(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
p2 <- ggplot(gse_means, aes(x=score, y=mean_score, fill=subtype)) +
  geom_bar(stat="identity", position="dodge") + scale_fill_manual(values=subtype_cols) +
  labs(title="GSE107943", y="Mean ssGSEA") + theme_pubr(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

pdf("figures/subtyping/Fig5I_TCGA_GSE_subtype_signature_consistency.pdf", width=11, height=6)
print(ggarrange(p1, p2, common.legend=TRUE, legend="bottom"))
dev.off()
log_msg("  Done")

# ============================================================
# Step 7: 汇总表
# ============================================================
log_msg("\n=== Step 7: 汇总表 ===")

summary_rows <- data.frame(
  item = c(
    "分型方法", "分型基因数",
    "TCGA C1 数", "TCGA C2 数",
    "GSE C1 数", "GSE C2 数",
    "TCGA OS p", "TCGA OS HR", "TCGA OS 95CI",
    "GSE OS p", "GSE OS HR", "GSE OS 95CI",
    "C1 higher CAF (TCGA)", "C1 higher CAF (GSE)",
    "C1 higher aggressive (TCGA)", "C1 higher aggressive (GSE)"
  ),
  value = c(
    "ConsensusClusterPlus k=2 + nearest centroid projection",
    length(subtype_genes),
    c1_n, c2_n, gse_c1_n, gse_c2_n,
    sprintf("%.4f", km_tcga["pval"]), sprintf("%.2f", km_tcga["HR"]),
    sprintf("%.2f-%.2f", km_tcga["lower95"], km_tcga["upper95"]),
    sprintf("%.4f", km_gse["pval"]), sprintf("%.2f", km_gse["HR"]),
    sprintf("%.2f-%.2f", km_gse["lower95"], km_gse["upper95"]),
    as.character(tcga_comp_sc$mean_C1[tcga_comp_sc$score=="CAF"] > tcga_comp_sc$mean_C2[tcga_comp_sc$score=="CAF"]),
    as.character(gse_comp_sc$mean_C1[gse_comp_sc$score=="CAF"] > gse_comp_sc$mean_C2[gse_comp_sc$score=="CAF"]),
    as.character(tcga_comp_sc$mean_C1[tcga_comp_sc$score=="CCA_aggressive_microenvironment"] > tcga_comp_sc$mean_C2[tcga_comp_sc$score=="CCA_aggressive_microenvironment"]),
    as.character(gse_comp_sc$mean_C1[gse_comp_sc$score=="CCA_aggressive_microenvironment"] > gse_comp_sc$mean_C2[gse_comp_sc$score=="CCA_aggressive_microenvironment"])
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_rows, "tables/subtyping/subtyping_analysis_summary.csv", row.names=FALSE)

# Generate DFS if available
log_msg("DFS check...")
if ("dsfree_mo_" %in% colnames(gse_surv_sub)) {
  gse_dfs <- gse_surv_sub[!is.na(gse_surv_sub$dsfree_mo_) & gse_surv_sub$dsfree_mo_ != "NA", ]
  gse_dfs$DFS_time <- as.numeric(gse_dfs$dsfree_mo_) * 30.4375
  gse_dfs$DFS_status <- as.numeric(gse_dfs$recurr)
  gse_dfs <- gse_dfs[!is.na(gse_dfs$DFS_time) & gse_dfs$DFS_time > 0, ]
  if (nrow(gse_dfs) > 10) {
    log_msg(sprintf("GSE DFS: %d patients available", nrow(gse_dfs)))
    km_subtype(gse_dfs, "GSE107943 DFS", "figures/subtyping/Fig5H2_GSE107943_subtype_KM_DFS.pdf")
    cat(sprintf("GSE DFS available: %d patients\n", nrow(gse_dfs)))
  }
} else {
  cat("GSE DFS: not available in clinical data\n")
}

# ============================================================
# 最终报告
# ============================================================
cat(sprintf("\n========================================\n"))
cat(sprintf("分子分型分析完成\n"))
cat(sprintf("========================================\n"))
cat(sprintf("分型基因数: %d\n", length(subtype_genes)))
cat(sprintf("TCGA:  C1=%d, C2=%d\n", c1_n, c2_n))
cat(sprintf("GSE:   C1=%d, C2=%d\n", gse_c1_n, gse_c2_n))
cat(sprintf("TCGA C1 CAF↑/Aggressive↑: %s\n",
  tcga_comp_sc$mean_C1[tcga_comp_sc$score=="CAF"] > tcga_comp_sc$mean_C2[tcga_comp_sc$score=="CAF"]))
cat(sprintf("GSE  C1 CAF↑/Aggressive↑: %s\n",
  gse_comp_sc$mean_C1[gse_comp_sc$score=="CAF"] > gse_comp_sc$mean_C2[gse_comp_sc$score=="CAF"]))
cat(sprintf("TCGA OS: p=%.4f, HR=%.2f\n", km_tcga["pval"], km_tcga["HR"]))
cat(sprintf("GSE  OS: p=%.4f, HR=%.2f\n", km_gse["pval"], km_gse["HR"]))
cat(sprintf("========================================\n"))
