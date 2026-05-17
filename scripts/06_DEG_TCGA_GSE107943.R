# ============================================================
# 分析目的: TCGA-CHOL 与 GSE107943 差异表达分析 + 跨队列验证
#          — 聚焦免疫代谢/CAF/巨噬细胞基因集
# 输入文件:
#   data/processed/TCGA_CHOL/TCGA_CHOL_counts_gene_symbol.csv
#   data/processed/TCGA_CHOL/TCGA_CHOL_log2TPM_gene_symbol.csv
#   data/processed/TCGA_CHOL/TCGA_CHOL_sample_metadata.csv
#   data/processed/GSE107943/GSE107943_counts_gene_symbol.csv
#   data/processed/GSE107943/GSE107943_log2CPM_gene_symbol.csv
#   data/processed/GSE107943/GSE107943_sample_metadata.csv
#   gene_sets/CCA_IM_CAF_TAM_combined_genes.txt
#   tables/CCA_IM_CAF_TAM_combined_gene_annotation.csv
# 输出文件:
#   tables/DEG_TCGA_CHOL_all.csv / significant.csv
#   tables/DEG_GSE107943_paired_all.csv / significant.csv
#   tables/DEG_TCGA_GSE107943_overlap_same_direction.csv
#   tables/CCA_IM_CAF_TAM_DEGs_TCGA.csv
#   tables/CCA_IM_CAF_TAM_DEGs_GSE107943.csv
#   tables/CCA_IM_CAF_TAM_DEGs_validated_same_direction.csv
#   tables/DEG_analysis_summary.csv
#   figures/DEG/Fig2A-G (7 张图)
# 主要方法:
#   TCGA: DESeq2 (Wald test, ~ sample_type)
#   GSE107943: edgeR paired design (~ pair_id + sample_type)
#   交叉验证: intersection + log2FC direction consistency
#   IM/CAF/TAM 基因集筛选
# ============================================================

library(DESeq2)
library(edgeR)
library(limma)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(ComplexHeatmap)
library(data.table)
library(RColorBrewer)
library(gridExtra)

dir.create("figures/DEG", showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ---- 通用配色 ----
tumor_col  <- "#E41A1C"
normal_col <- "#377EB8"
fc_cut  <- 1      # abs(log2FC) > 1
fdr_cut <- 0.05

# ============================================================
# Step 1: 加载数据
# ============================================================
log_msg("=== Step 1: 加载数据 ===")

# TCGA
tcga_counts <- fread("data/processed/TCGA_CHOL/TCGA_CHOL_counts_gene_symbol.csv", data.table = FALSE)
rownames(tcga_counts) <- tcga_counts$gene_symbol; tcga_counts$gene_symbol <- NULL

tcga_log2tpm <- fread("data/processed/TCGA_CHOL/TCGA_CHOL_log2TPM_gene_symbol.csv", data.table = FALSE)
rownames(tcga_log2tpm) <- tcga_log2tpm$gene_symbol; tcga_log2tpm$gene_symbol <- NULL

tcga_meta <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_sample_metadata.csv", stringsAsFactors = FALSE)
rownames(tcga_meta) <- tcga_meta$sample_id

# 对齐
common_tcga <- intersect(colnames(tcga_counts), rownames(tcga_meta))
tcga_counts <- tcga_counts[, common_tcga, drop = FALSE]
tcga_log2tpm <- tcga_log2tpm[, common_tcga, drop = FALSE]
tcga_meta <- tcga_meta[common_tcga, , drop = FALSE]

log_msg(sprintf("TCGA counts: %d genes x %d samples", nrow(tcga_counts), ncol(tcga_counts)))

# GSE107943
gse_counts <- fread("data/processed/GSE107943/GSE107943_counts_gene_symbol.csv", data.table = FALSE)
rownames(gse_counts) <- gse_counts$gene_symbol; gse_counts$gene_symbol <- NULL

gse_log2cpm <- fread("data/processed/GSE107943/GSE107943_log2CPM_gene_symbol.csv", data.table = FALSE)
rownames(gse_log2cpm) <- gse_log2cpm$gene_symbol; gse_log2cpm$gene_symbol <- NULL

gse_meta <- read.csv("data/processed/GSE107943/GSE107943_sample_metadata.csv", stringsAsFactors = FALSE)
rownames(gse_meta) <- gse_meta$sample_id

common_gse <- intersect(colnames(gse_counts), rownames(gse_meta))
gse_counts <- gse_counts[, common_gse, drop = FALSE]
gse_log2cpm <- gse_log2cpm[, common_gse, drop = FALSE]
gse_meta <- gse_meta[common_gse, , drop = FALSE]

log_msg(sprintf("GSE107943 counts: %d genes x %d samples", nrow(gse_counts), ncol(gse_counts)))

# 基因集
im_genes <- scan("gene_sets/CCA_IM_CAF_TAM_combined_genes.txt", what = "character",
                 skip = 1, quiet = TRUE)
gene_annot <- read.csv("tables/CCA_IM_CAF_TAM_combined_gene_annotation.csv", stringsAsFactors = FALSE)

# ============================================================
# Step 2: TCGA-CHOL DESeq2 DEG
# ============================================================
log_msg("\n=== Step 2: TCGA-CHOL DESeq2 DEG ===")

tcga_meta$sample_type <- factor(tcga_meta$sample_type,
                                levels = c("Primary Tumor", "Solid Tissue Normal"))

# 过滤低表达
tcga_counts_int <- round(as.matrix(tcga_counts))
storage.mode(tcga_counts_int) <- "integer"
keep <- rowSums(tcga_counts_int >= 10) >= ceiling(0.2 * ncol(tcga_counts_int))
tcga_counts_filt <- tcga_counts_int[keep, , drop = FALSE]
log_msg(sprintf("过滤: %d → %d 基因", nrow(tcga_counts_int), nrow(tcga_counts_filt)))

# DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = tcga_counts_filt,
  colData   = tcga_meta,
  design    = ~ sample_type
)
dds$sample_type <- relevel(dds$sample_type, ref = "Solid Tissue Normal")
dds <- DESeq(dds)

tcga_res <- results(dds, contrast = c("sample_type", "Primary Tumor", "Solid Tissue Normal"),
                    alpha = fdr_cut)
tcga_res <- tcga_res[order(tcga_res$pvalue), , drop = FALSE]
tcga_df <- as.data.frame(tcga_res)
tcga_df$gene_symbol <- rownames(tcga_df)
tcga_df <- tcga_df[, c("gene_symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
tcga_df$significant <- !is.na(tcga_df$padj) & tcga_df$padj < fdr_cut & abs(tcga_df$log2FoldChange) > fc_cut
tcga_df$direction <- ifelse(tcga_df$log2FoldChange > 0, "up", "down")

n_tcga_up   <- sum(tcga_df$significant & tcga_df$direction == "up")
n_tcga_down <- sum(tcga_df$significant & tcga_df$direction == "down")
n_tcga_sig  <- n_tcga_up + n_tcga_down
log_msg(sprintf("TCGA DEG: %d (up=%d, down=%d)", n_tcga_sig, n_tcga_up, n_tcga_down))

# 保存
write.csv(tcga_df, "tables/DEG_TCGA_CHOL_all.csv", row.names = FALSE)
write.csv(tcga_df[tcga_df$significant, , drop = FALSE],
          "tables/DEG_TCGA_CHOL_significant.csv", row.names = FALSE)
log_msg("已保存 TCGA DEG 表格")

# ============================================================
# Step 3: GSE107943 edgeR paired DEG
# ============================================================
log_msg("\n=== Step 3: GSE107943 edgeR paired DEG ===")

# 只取 paired samples
paired_ids <- names(which(table(gse_meta$pair_id) == 2))
gse_paired_meta <- gse_meta[gse_meta$pair_id %in% paired_ids, , drop = FALSE]
log_msg(sprintf("Paired 样本: %d (来自 %d 对)", nrow(gse_paired_meta), length(paired_ids)))

gse_paired_counts <- gse_counts[, rownames(gse_paired_meta), drop = FALSE]
gse_paired_counts <- round(as.matrix(gse_paired_counts))
# 移除含 NA 的行
na_rows <- apply(gse_paired_counts, 1, function(x) any(is.na(x)))
if (sum(na_rows) > 0) {
  log_msg(sprintf("移除 %d 个含 NA 基因: %s",
    sum(na_rows), paste(rownames(gse_paired_counts)[na_rows], collapse=", ")))
  gse_paired_counts <- gse_paired_counts[!na_rows, , drop = FALSE]
}
storage.mode(gse_paired_counts) <- "integer"

log_msg(sprintf("Paired counts: %d genes x %d samples", nrow(gse_paired_counts), ncol(gse_paired_counts)))

# edgeR 过滤
gse_paired_meta$sample_type <- factor(gse_paired_meta$sample_type,
                                      levels = c("Tumor", "Normal"))
gse_paired_meta$pair_id <- factor(gse_paired_meta$pair_id)

y <- DGEList(counts = gse_paired_counts, samples = gse_paired_meta,
             group = gse_paired_meta$sample_type)
keep <- filterByExpr(y, group = gse_paired_meta$sample_type)
y <- y[keep, ]
log_msg(sprintf("过滤: %d → %d 基因", nrow(gse_paired_counts), nrow(y$counts)))

# 归一化
y <- calcNormFactors(y)

# paired design
design <- model.matrix(~ pair_id + sample_type, data = gse_paired_meta)
# 找到 sample_type 对应的列
sample_coef <- grep("sample_type", colnames(design), value = TRUE)
log_msg(sprintf("Design matrix sample_type column: %s", sample_coef))

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = sample_coef)

# 修正方向：如果 coef 是 sample_typeNormal，logFC 符号需反转
gse_df <- as.data.frame(topTags(qlf, n = Inf))
gse_df$gene_symbol <- rownames(gse_df)
if (grepl("Normal", sample_coef)) {
  log_msg("翻转 GSE log2FC 符号 (coef = Normal, 转变为 Tumor vs Normal)")
  gse_df$logFC <- -gse_df$logFC
}
gse_df <- gse_df[, c("gene_symbol", "logCPM", "logFC", "PValue", "FDR")]
colnames(gse_df) <- c("gene_symbol", "baseMean", "log2FoldChange", "pvalue", "padj")
gse_df$significant <- !is.na(gse_df$padj) & gse_df$padj < fdr_cut & abs(gse_df$log2FoldChange) > fc_cut
gse_df$direction <- ifelse(gse_df$log2FoldChange > 0, "up", "down")

n_gse_up   <- sum(gse_df$significant & gse_df$direction == "up")
n_gse_down <- sum(gse_df$significant & gse_df$direction == "down")
n_gse_sig  <- n_gse_up + n_gse_down
log_msg(sprintf("GSE107943 DEG: %d (up=%d, down=%d)", n_gse_sig, n_gse_up, n_gse_down))

# 保存
write.csv(gse_df, "tables/DEG_GSE107943_paired_all.csv", row.names = FALSE)
write.csv(gse_df[gse_df$significant, , drop = FALSE],
          "tables/DEG_GSE107943_paired_significant.csv", row.names = FALSE)
log_msg("已保存 GSE107943 DEG 表格")

# ============================================================
# Step 4: 跨队列交叉验证
# ============================================================
log_msg("\n=== Step 4: 跨队列交叉验证 ===")

tcga_sig <- tcga_df[tcga_df$significant, c("gene_symbol", "log2FoldChange", "padj", "direction")]
gse_sig  <- gse_df[gse_df$significant, c("gene_symbol", "log2FoldChange", "padj", "direction")]

# 交集
common_genes <- intersect(tcga_sig$gene_symbol, gse_sig$gene_symbol)
log_msg(sprintf("共同显著 DEG: %d", length(common_genes)))

tcga_sig_c <- tcga_sig[tcga_sig$gene_symbol %in% common_genes, ]
gse_sig_c  <- gse_sig[gse_sig$gene_symbol %in% common_genes, ]

# 合并
overlap_df <- merge(tcga_sig_c, gse_sig_c, by = "gene_symbol", suffixes = c("_TCGA", "_GSE"))
overlap_df$same_direction <- overlap_df$direction_TCGA == overlap_df$direction_GSE
n_same_dir <- sum(overlap_df$same_direction)
log_msg(sprintf("方向一致: %d / %d", n_same_dir, nrow(overlap_df)))

overlap_df <- overlap_df[order(overlap_df$padj_TCGA), ]
write.csv(overlap_df, "tables/DEG_TCGA_GSE107943_overlap_same_direction.csv", row.names = FALSE)
log_msg("已保存 overlap 表格")

# ============================================================
# Step 5: 聚焦 IM/CAF/TAM 基因集
# ============================================================
log_msg("\n=== Step 5: 聚焦 IM/CAF/TAM 基因集 ===")

# TCGA DEG ∩ IM/CAF/TAM
tcga_im <- tcga_df[tcga_df$gene_symbol %in% im_genes & tcga_df$significant, , drop = FALSE]
log_msg(sprintf("TCGA IM/CAF/TAM DEGs: %d", nrow(tcga_im)))
write.csv(tcga_im, "tables/CCA_IM_CAF_TAM_DEGs_TCGA.csv", row.names = FALSE)

# GSE DEG ∩ IM/CAF/TAM
gse_im <- gse_df[gse_df$gene_symbol %in% im_genes & gse_df$significant, , drop = FALSE]
log_msg(sprintf("GSE107943 IM/CAF/TAM DEGs: %d", nrow(gse_im)))
write.csv(gse_im, "tables/CCA_IM_CAF_TAM_DEGs_GSE107943.csv", row.names = FALSE)

# validated IM/CAF/TAM
overlap_same <- overlap_df[overlap_df$same_direction, , drop = FALSE]
overlap_im <- overlap_same[overlap_same$gene_symbol %in% im_genes, , drop = FALSE]
overlap_im <- merge(overlap_im, gene_annot, by = "gene_symbol", all.x = TRUE)
overlap_im <- overlap_im[order(overlap_im$padj_TCGA), ]
log_msg(sprintf("Validated IM/CAF/TAM DEGs (same direction): %d", nrow(overlap_im)))
write.csv(overlap_im, "tables/CCA_IM_CAF_TAM_DEGs_validated_same_direction.csv", row.names = FALSE)

# ============================================================
# Step 6: 可视化
# ============================================================
log_msg("\n=== Step 6: 可视化 ===")

# ---- Fig2A: TCGA PCA ----
log_msg("Fig2A: TCGA PCA")
tcga_sds <- apply(tcga_log2tpm, 1, sd, na.rm = TRUE)
tcga_sds <- tcga_sds[is.finite(tcga_sds) & tcga_sds > 0]
top_var_tcga <- names(sort(tcga_sds, decreasing = TRUE))[1:min(5000, length(tcga_sds))]
tcga_log2tpm_pca <- tcga_log2tpm[top_var_tcga, , drop = FALSE]
tcga_log2tpm_pca <- tcga_log2tpm_pca[complete.cases(tcga_log2tpm_pca), , drop = FALSE]
log_msg(sprintf("TCGA PCA 基因: %d", nrow(tcga_log2tpm_pca)))
tcga_pca <- prcomp(t(tcga_log2tpm_pca), scale. = TRUE)
pca_df <- data.frame(PC1 = tcga_pca$x[,1], PC2 = tcga_pca$x[,2],
                     sample_type = tcga_meta$sample_type)
pca_var <- round(100 * summary(tcga_pca)$importance[2, 1:2], 1)

pdf("figures/DEG/Fig2A_TCGA_PCA_tumor_normal.pdf", width = 7, height = 6)
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = sample_type)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = c("Primary Tumor" = tumor_col,
                                "Solid Tissue Normal" = normal_col)) +
  labs(x = sprintf("PC1 (%s%%)", pca_var[1]), y = sprintf("PC2 (%s%%)", pca_var[2]),
       title = "TCGA-CHOL PCA (log2 TPM+1)", color = "") +
  theme_pubr(base_size = 14) + theme(legend.position = "top")
print(p)
dev.off()
log_msg("  Done")

# ---- Fig2B: TCGA Volcano ----
log_msg("Fig2B: TCGA Volcano")
tcga_df$neg_log10_p <- -log10(tcga_df$pvalue)
tcga_df$label <- ifelse(tcga_df$significant &
                          rank(tcga_df$neg_log10_p, ties.method = "random") > (nrow(tcga_df) - 10),
                        tcga_df$gene_symbol, "")

pdf("figures/DEG/Fig2B_TCGA_volcano.pdf", width = 8, height = 7)
p <- ggplot(tcga_df, aes(x = log2FoldChange, y = neg_log10_p)) +
  geom_point(aes(color = significant), size = 0.8, alpha = 0.6) +
  scale_color_manual(values = c("TRUE" = tumor_col, "FALSE" = "grey70")) +
  geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", color = "grey50") +
  ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = 20) +
  labs(x = "log2 Fold Change (Tumor vs Normal)", y = "-log10(p-value)",
       title = "TCGA-CHOL DEG Volcano", color = "Significant") +
  theme_pubr(base_size = 14)
print(p)
dev.off()
log_msg("  Done")

# ---- Fig2C: TCGA Top IM/CAF/TAM DEG Heatmap ----
log_msg("Fig2C: TCGA IM/CAF/TAM DEG Heatmap")
tcga_im_genes <- tcga_im$gene_symbol
if (length(tcga_im_genes) > 50) {
  tcga_im_genes <- tcga_im$gene_symbol[1:min(50, nrow(tcga_im))]
}

mat_tcga_heat <- as.matrix(tcga_log2tpm[tcga_im_genes, , drop = FALSE])
mat_tcga_z <- t(scale(t(mat_tcga_heat)))

# 注释
ha <- HeatmapAnnotation(
  Type = tcga_meta$sample_type,
  col = list(Type = c("Primary Tumor" = tumor_col, "Solid Tissue Normal" = normal_col)),
  annotation_legend_param = list(Type = list(title = "Sample Type"))
)

pdf("figures/DEG/Fig2C_TCGA_top_IM_CAF_TAM_DEG_heatmap.pdf", width = 12, height = 10)
ht <- Heatmap(mat_tcga_z, name = "Z-score",
              top_annotation = ha,
              show_row_names = (nrow(mat_tcga_z) <= 50),
              show_column_names = FALSE,
              cluster_rows = TRUE, cluster_columns = TRUE,
              row_names_gp = gpar(fontsize = 7),
              column_title = sprintf("TCGA-CHOL: Top %d IM/CAF/TAM DEGs", nrow(mat_tcga_z)))
draw(ht, annotation_legend_side = "bottom")
dev.off()
log_msg("  Done")

# ---- Fig2D: GSE107943 PCA ----
log_msg("Fig2D: GSE107943 PCA")
# 取 top 5000 最高方差的基因做 PCA
gse_sds <- apply(gse_log2cpm, 1, sd, na.rm = TRUE)
gse_sds <- gse_sds[is.finite(gse_sds) & gse_sds > 0]
top_var_gse <- names(sort(gse_sds, decreasing = TRUE))[1:min(5000, length(gse_sds))]
gse_log2cpm_pca <- gse_log2cpm[top_var_gse, , drop = FALSE]
gse_log2cpm_pca <- gse_log2cpm_pca[complete.cases(gse_log2cpm_pca), , drop = FALSE]
log_msg(sprintf("GSE PCA 基因: %d", nrow(gse_log2cpm_pca)))
gse_pca <- prcomp(t(gse_log2cpm_pca), scale. = TRUE)
pca_gse <- data.frame(PC1 = gse_pca$x[,1], PC2 = gse_pca$x[,2],
                      sample_type = gse_meta$sample_type)
gse_pca_var <- round(100 * summary(gse_pca)$importance[2, 1:2], 1)

pdf("figures/DEG/Fig2D_GSE107943_PCA_tumor_normal.pdf", width = 7, height = 6)
p <- ggplot(pca_gse, aes(x = PC1, y = PC2, color = sample_type)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = c("Tumor" = tumor_col, "Normal" = normal_col)) +
  labs(x = sprintf("PC1 (%s%%)", gse_pca_var[1]), y = sprintf("PC2 (%s%%)", gse_pca_var[2]),
       title = "GSE107943 PCA (log2 CPM+1)", color = "") +
  theme_pubr(base_size = 14) + theme(legend.position = "top")
print(p)
dev.off()
log_msg("  Done")

# ---- Fig2E: GSE107943 Paired Volcano ----
log_msg("Fig2E: GSE107943 Volcano")
gse_df$neg_log10_p <- -log10(gse_df$pvalue)
gse_df$label <- ifelse(gse_df$significant &
                         rank(gse_df$neg_log10_p, ties.method = "random") > (nrow(gse_df) - 10),
                       gse_df$gene_symbol, "")

pdf("figures/DEG/Fig2E_GSE107943_paired_volcano.pdf", width = 8, height = 7)
p <- ggplot(gse_df, aes(x = log2FoldChange, y = neg_log10_p)) +
  geom_point(aes(color = significant), size = 0.8, alpha = 0.6) +
  scale_color_manual(values = c("TRUE" = tumor_col, "FALSE" = "grey70")) +
  geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", color = "grey50") +
  ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = 20) +
  labs(x = "log2 Fold Change (Tumor vs Normal)", y = "-log10(p-value)",
       title = "GSE107943 Paired DEG Volcano", color = "Significant") +
  theme_pubr(base_size = 14)
print(p)
dev.off()
log_msg("  Done")

# ---- Fig2F: Overlap UpSet / Venn ----
log_msg("Fig2F: Overlap diagram")
tcga_sig_set <- tcga_df$gene_symbol[tcga_df$significant]
gse_sig_set  <- gse_df$gene_symbol[gse_df$significant]
tcga_up_set  <- tcga_df$gene_symbol[tcga_df$significant & tcga_df$direction == "up"]
tcga_dn_set  <- tcga_df$gene_symbol[tcga_df$significant & tcga_df$direction == "down"]
gse_up_set   <- gse_df$gene_symbol[gse_df$significant & gse_df$direction == "up"]
gse_dn_set   <- gse_df$gene_symbol[gse_df$significant & gse_df$direction == "down"]

pdf("figures/DEG/Fig2F_TCGA_GSE107943_overlap_upset_or_venn.pdf", width = 10, height = 6)
par(mfrow = c(1, 2))
# Left: all DEG overlap
all_sets <- list(
  TCGA_UP = tcga_up_set,
  TCGA_DOWN = tcga_dn_set,
  GSE_UP = gse_up_set,
  GSE_DOWN = gse_dn_set
)

# Venn: TCGA vs GSE DEG overlap
par(mar = c(3, 3, 4, 2))
tcga_only <- length(setdiff(tcga_sig_set, gse_sig_set))
gse_only  <- length(setdiff(gse_sig_set, tcga_sig_set))
both_deg  <- length(intersect(tcga_sig_set, gse_sig_set))

# Simple Venn via draw.pairwise.venn
if (requireNamespace("VennDiagram", quietly = TRUE)) {
  VennDiagram::draw.pairwise.venn(
    area1 = length(tcga_sig_set), area2 = length(gse_sig_set),
    cross.area = both_deg,
    category = c("TCGA-CHOL", "GSE107943"),
    cat.cex = 1.2, cex = 1.2,
    fill = c(tumor_col, "darkorange"), alpha = c(0.5, 0.5)
  )
} else {
  # fallback: text only
  plot.new()
  text(0.5, 0.8, "DEG OVERLAP", cex = 1.5, font = 2)
  text(0.5, 0.6, sprintf("TCGA only: %d", tcga_only), cex = 1.1)
  text(0.5, 0.5, sprintf("GSE only: %d", gse_only), cex = 1.1)
  text(0.5, 0.4, sprintf("Both: %d", both_deg), cex = 1.1, font = 2)
}

# Right: validated IM DEG
valid_im_genes <- overlap_im$gene_symbol
valid_im_up <- overlap_im$gene_symbol[overlap_im$direction_TCGA == "up"]
valid_im_dn <- overlap_im$gene_symbol[overlap_im$direction_TCGA == "down"]

# Text summary
plot.new()
text(0.5, 0.8, "CROSS-COHORT VALIDATION", cex = 1.5, font = 2)
text(0.5, 0.65, sprintf("Common DEGs: %d", length(common_genes)), cex = 1.2)
text(0.5, 0.55, sprintf("Same direction: %d", n_same_dir), cex = 1.2)
text(0.5, 0.45, sprintf("Validated IM/CAF/TAM: %d", nrow(overlap_im)), cex = 1.2)
text(0.5, 0.35, sprintf("  Up-regulated: %d", sum(overlap_im$direction_TCGA == "up")), cex = 1)
text(0.5, 0.25, sprintf("  Down-regulated: %d", sum(overlap_im$direction_TCGA == "down")), cex = 1)
dev.off()
log_msg("  Done")

# ---- Fig2G: Validated IM/CAF/TAM DEGs Heatmap (both cohorts) ----
log_msg("Fig2G: Validated IM/CAF/TAM heatmap")
valid_genes <- overlap_im$gene_symbol

if (length(valid_genes) > 0) {
  # TCGA side
  valid_tcga <- intersect(valid_genes, rownames(tcga_log2tpm))
  mat_tcga_v <- as.matrix(tcga_log2tpm[valid_tcga, , drop = FALSE])
  mat_tcga_z <- t(scale(t(mat_tcga_v)))

  # GSE side
  valid_gse <- intersect(valid_genes, rownames(gse_log2cpm))
  mat_gse_v <- as.matrix(gse_log2cpm[valid_gse, , drop = FALSE])
  mat_gse_z <- t(scale(t(mat_gse_v)))

  # 合并两个矩阵（只取共同基因）
  common_valid <- intersect(valid_tcga, valid_gse)
  if (length(common_valid) > 0) {
    mat_combined <- cbind(mat_tcga_z[common_valid, , drop = FALSE],
                          mat_gse_z[common_valid, , drop = FALSE])

    # 注释
    cohort_ann <- c(rep("TCGA-CHOL", ncol(mat_tcga_z)),
                    rep("GSE107943", ncol(mat_gse_z)))
    type_ann <- c(as.character(tcga_meta$sample_type),
                  as.character(gse_meta$sample_type))

    ha_top <- HeatmapAnnotation(
      Cohort = cohort_ann,
      Type = type_ann,
      col = list(
        Cohort = c("TCGA-CHOL" = "#4DAF4A", "GSE107943" = "#984EA3"),
        Type = c("Primary Tumor" = tumor_col, "Solid Tissue Normal" = normal_col,
                 "Tumor" = tumor_col, "Normal" = normal_col)
      ),
      annotation_legend_param = list(Cohort = list(title = "Cohort"),
                                     Type = list(title = "Sample Type"))
    )

    n_show <- min(length(common_valid), 60)
    if (length(common_valid) > n_show) {
      # 选择 FDR 最小的 n_show 个
      top_valid <- overlap_im$gene_symbol[order(overlap_im$padj_TCGA)][1:n_show]
      top_valid <- intersect(top_valid, common_valid)
      mat_combined <- mat_combined[top_valid, , drop = FALSE]
    }

    pdf("figures/DEG/Fig2G_validated_IM_CAF_TAM_DEGs_heatmap.pdf", width = 16, height = min(14, nrow(mat_combined) * 0.18 + 3))
    ht <- Heatmap(mat_combined, name = "Z-score",
                  top_annotation = ha_top,
                  show_row_names = (nrow(mat_combined) <= 60),
                  show_column_names = FALSE,
                  cluster_rows = TRUE, cluster_columns = FALSE,
                  column_split = cohort_ann,
                  row_names_gp = gpar(fontsize = 7),
                  column_title = "Validated IM/CAF/TAM DEGs (Same Direction)")
    draw(ht, annotation_legend_side = "bottom")
    dev.off()
    log_msg(sprintf("  Heatmap: %d genes", nrow(mat_combined)))
  } else {
    pdf("figures/DEG/Fig2G_validated_IM_CAF_TAM_DEGs_heatmap.pdf", width = 8, height = 4)
    plot.new(); text(0.5, 0.5, "No common validated genes between cohorts", cex = 1.2)
    dev.off()
  }
} else {
  pdf("figures/DEG/Fig2G_validated_IM_CAF_TAM_DEGs_heatmap.pdf", width = 8, height = 4)
  plot.new(); text(0.5, 0.5, "No validated IM/CAF/TAM DEGs found", cex = 1.2)
  dev.off()
}
log_msg("  Done")

# ============================================================
# Step 7: 汇总表
# ============================================================
log_msg("\n=== Step 7: 生成汇总表 ===")

summary_rows <- data.frame(
  item = c(
    "TCGA 原始基因数", "TCGA 过滤后基因数", "TCGA 显著 DEG 数",
    "TCGA 上调数", "TCGA 下调数",
    "GSE107943 原始基因数", "GSE107943 过滤后基因数",
    "GSE107943 显著 DEG 数", "GSE107943 上调数", "GSE107943 下调数",
    "共同 DEG 数", "方向一致 validated DEG 数",
    "IM/CAF/TAM TCGA DEG", "IM/CAF/TAM GSE DEG",
    "IM/CAF/TAM validated same-direction DEG",
    "FDR cutoff", "log2FC cutoff",
    "TCGA 方法", "GSE107943 方法"
  ),
  value = c(
    nrow(tcga_counts_int), nrow(tcga_counts_filt), n_tcga_sig,
    n_tcga_up, n_tcga_down,
    nrow(gse_paired_counts), nrow(y$counts), n_gse_sig,
    n_gse_up, n_gse_down,
    length(common_genes), n_same_dir,
    nrow(tcga_im), nrow(gse_im), nrow(overlap_im),
    fdr_cut, fc_cut,
    "DESeq2 Wald test", "edgeR glmQLF paired (~ pair_id + sample_type)"
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_rows, "tables/DEG_analysis_summary.csv", row.names = FALSE)
log_msg("已保存: tables/DEG_analysis_summary.csv")

# ---- 打印前 20 validated IM/CAF/TAM DEGs ----
cat(sprintf("\n=== Top 20 Validated IM/CAF/TAM DEGs ===\n"))
if (nrow(overlap_im) > 0) {
  top20 <- head(overlap_im, 20)
  print(data.frame(
    gene = top20$gene_symbol,
    TCGA_log2FC = round(top20$log2FoldChange_TCGA, 3),
    GSE_log2FC  = round(top20$log2FoldChange_GSE, 3),
    direction   = top20$direction_TCGA,
    source      = top20$source
  ), row.names = FALSE)
} else {
  cat("  No validated IM/CAF/TAM DEGs found.\n")
}

# ============================================================
# 最终报告
# ============================================================
cat(sprintf("\n========================================\n"))
cat(sprintf("DEG 分析和交叉验证完成\n"))
cat(sprintf("========================================\n"))
cat(sprintf("TCGA DEG: %d (up=%d, down=%d)\n", n_tcga_sig, n_tcga_up, n_tcga_down))
cat(sprintf("GSE107943 DEG: %d (up=%d, down=%d)\n", n_gse_sig, n_gse_up, n_gse_down))
cat(sprintf("Common DEGs: %d\n", length(common_genes)))
cat(sprintf("Same-direction validated: %d\n", n_same_dir))
cat(sprintf("IM/CAF/TAM validated: %d\n", nrow(overlap_im)))
cat(sprintf("Figures: figures/DEG/Fig2A-G.pdf\n"))
cat(sprintf("========================================\n"))
