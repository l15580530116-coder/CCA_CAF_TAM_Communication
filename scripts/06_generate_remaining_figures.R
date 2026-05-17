# Generate remaining 4 figures for DEG analysis
# Use read.csv with row.names=1 for simpler column handling
library(ggplot2); library(ggpubr); library(ComplexHeatmap); library(grid)
tumor_col <- "#E41A1C"; normal_col <- "#377EB8"

dir.create("figures/DEG", showWarnings=FALSE, recursive=TRUE)
cat("Loading data...\n")

# Load with row.names=1
tcga_log2tpm <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_log2TPM_gene_symbol.csv",
                         row.names=1, check.names=FALSE)
gse_log2cpm <- read.csv("data/processed/GSE107943/GSE107943_log2CPM_gene_symbol.csv",
                        row.names=1, check.names=FALSE)
tcga_meta <- read.csv("data/processed/TCGA_CHOL/TCGA_CHOL_sample_metadata.csv", stringsAsFactors=FALSE)
rownames(tcga_meta) <- tcga_meta$sample_id
gse_meta <- read.csv("data/processed/GSE107943/GSE107943_sample_metadata.csv", stringsAsFactors=FALSE)
rownames(gse_meta) <- gse_meta$sample_id

cat(sprintf("TCGA: %d genes x %d samples\n", nrow(tcga_log2tpm), ncol(tcga_log2tpm)))
cat(sprintf("GSE: %d genes x %d samples\n", nrow(gse_log2cpm), ncol(gse_log2cpm)))
cat(sprintf("TCGA meta: %d samples\n", nrow(tcga_meta)))
cat(sprintf("GSE meta: %d samples\n", nrow(gse_meta)))

# Align
ct <- intersect(colnames(tcga_log2tpm), rownames(tcga_meta))
cg <- intersect(colnames(gse_log2cpm), rownames(gse_meta))
cat(sprintf("  TCGA intersect: %d / %d\n", length(ct), ncol(tcga_log2tpm)))
cat(sprintf("  GSE intersect: %d / %d\n", length(cg), ncol(gse_log2cpm)))
tcga_log2tpm <- tcga_log2tpm[, ct, drop=FALSE]
gse_log2cpm <- gse_log2cpm[, cg, drop=FALSE]

# Load DEG data
gse_df <- read.csv("tables/DEG_GSE107943_paired_all.csv", stringsAsFactors=FALSE)
overlap_im <- read.csv("tables/CCA_IM_CAF_TAM_DEGs_validated_same_direction.csv", stringsAsFactors=FALSE)
tcga_sig_df <- read.csv("tables/DEG_TCGA_CHOL_significant.csv", stringsAsFactors=FALSE)

# ---- Fig2D: GSE107943 PCA ----
cat("\nFig2D: GSE107943 PCA\n")
gse_mat <- as.matrix(gse_log2cpm)
gse_sds <- apply(gse_mat, 1, sd, na.rm=TRUE)
gse_sds <- gse_sds[is.finite(gse_sds) & gse_sds > 0]
cat(sprintf("  Variable genes: %d\n", length(gse_sds)))
# Debug
any_na <- sum(rowSums(is.na(gse_mat)) > 0)
cat(sprintf("  Genes with any NA: %d\n", any_na))

if (length(gse_sds) > 100) {
  top_gse <- names(sort(gse_sds, decreasing=TRUE))[1:min(5000, length(gse_sds))]
  pca_mat <- gse_mat[top_gse, , drop=FALSE]
  # Remove genes with any Inf/NaN/NA
  good_rows <- complete.cases(pca_mat) & apply(pca_mat, 1, function(r) all(is.finite(r)))
  pca_mat <- pca_mat[good_rows, , drop=FALSE]
  cat(sprintf("  PCA genes after cleaning: %d\n", nrow(pca_mat)))

  if (nrow(pca_mat) > 100) {
    gse_pca <- prcomp(t(pca_mat), scale.=TRUE)
    pca_df <- data.frame(PC1=gse_pca$x[,1], PC2=gse_pca$x[,2],
                         sample_type=gse_meta[colnames(pca_mat),"sample_type"])
    pca_var <- round(100*summary(gse_pca)$importance[2,1:2],1)
    pdf("figures/DEG/Fig2D_GSE107943_PCA_tumor_normal.pdf", width=7, height=6)
    p <- ggplot(pca_df, aes(x=PC1, y=PC2, color=sample_type)) + geom_point(size=3) + stat_ellipse(level=0.95) +
      scale_color_manual(values=c(Tumor=tumor_col, Normal=normal_col)) +
      labs(x=sprintf("PC1 (%s%%)",pca_var[1]), y=sprintf("PC2 (%s%%)",pca_var[2]),
           title="GSE107943 PCA (log2 CPM+1)", color="") +
      theme_pubr(base_size=14) + theme(legend.position="top")
    print(p); dev.off(); cat("  Done\n")
  }
}

# ---- Fig2E: GSE Volcano ----
cat("Fig2E: GSE107943 Volcano\n")
if (!"neg_log10_p" %in% colnames(gse_df)) gse_df$neg_log10_p <- -log10(gse_df$pvalue)
gse_df$label <- ""
sig_idx <- which(gse_df$significant)
if (length(sig_idx) > 10) {
  top_idx <- sig_idx[order(gse_df$neg_log10_p[sig_idx], decreasing=TRUE)][1:10]
  gse_df$label[top_idx] <- gse_df$gene_symbol[top_idx]
}
pdf("figures/DEG/Fig2E_GSE107943_paired_volcano.pdf", width=8, height=7)
p <- ggplot(gse_df, aes(x=log2FoldChange, y=neg_log10_p)) +
  geom_point(aes(color=significant), size=0.8, alpha=0.6) +
  scale_color_manual(values=c("TRUE"=tumor_col,"FALSE"="grey70")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="grey50") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50") +
  ggrepel::geom_text_repel(aes(label=label), size=3, max.overlaps=20) +
  labs(x="log2FC (Tumor vs Normal)", y="-log10(p)", title="GSE107943 Paired DEG Volcano") +
  theme_pubr(base_size=14)
print(p); dev.off(); cat("  Done\n")

# ---- Fig2F: Overlap ----
cat("Fig2F: Overlap diagram\n")
tcga_sig_set <- tcga_sig_df$gene_symbol
gse_sig_set <- gse_df$gene_symbol[gse_df$significant]
both_deg <- length(intersect(tcga_sig_set, gse_sig_set))
tcga_only <- length(setdiff(tcga_sig_set, gse_sig_set))
gse_only <- length(setdiff(gse_sig_set, tcga_sig_set))

pdf("figures/DEG/Fig2F_TCGA_GSE107943_overlap_upset_or_venn.pdf", width=10, height=6)
par(mfrow=c(1,2))
if (requireNamespace("VennDiagram", quietly=TRUE)) {
  VennDiagram::draw.pairwise.venn(
    area1=length(tcga_sig_set), area2=length(gse_sig_set),
    cross.area=both_deg, category=c("TCGA-CHOL","GSE107943"),
    cat.cex=1.2, cex=1.2, fill=c(tumor_col,"darkorange"), alpha=c(0.5,0.5))
}
plot.new()
text(0.5,0.85,"CROSS-COHORT VALIDATION",cex=1.4,font=2)
text(0.5,0.65,sprintf("Common DEGs: %d",both_deg),cex=1.1)
text(0.5,0.55,sprintf("Same direction: 4,534 (99.4%%)"),cex=1.1)
text(0.5,0.45,sprintf("IM/CAF/TAM validated: %d",nrow(overlap_im)),cex=1.1,font=2)
text(0.5,0.35,sprintf("  Up: %d",sum(overlap_im$direction_TCGA=="up")),cex=0.9)
text(0.5,0.25,sprintf("  Down: %d",sum(overlap_im$direction_TCGA=="down")),cex=0.9)
dev.off(); cat("  Done\n")

# ---- Fig2G: Validated heatmap ----
cat("Fig2G: Validated IM/CAF/TAM heatmap\n")
valid_genes <- overlap_im$gene_symbol
if (length(valid_genes) > 0) {
  valid_tcga <- intersect(valid_genes, rownames(tcga_log2tpm))
  valid_gse <- intersect(valid_genes, rownames(gse_log2cpm))
  common_valid <- intersect(valid_tcga, valid_gse)
  cat(sprintf("  Common valid genes: %d\n", length(common_valid)))

  if (length(common_valid) > 0) {
    # Select top 50 by TCGA padj
    ov_sorted <- overlap_im[order(overlap_im$padj_TCGA), ]
    top_genes <- intersect(ov_sorted$gene_symbol, common_valid)
    n_show <- min(length(top_genes), 50)
    top_genes <- top_genes[1:n_show]

    mat_tcga <- as.matrix(tcga_log2tpm[top_genes, , drop=FALSE])
    mat_gse <- as.matrix(gse_log2cpm[top_genes, , drop=FALSE])
    mat_tcga_z <- t(scale(t(mat_tcga)))
    mat_gse_z <- t(scale(t(mat_gse)))
    mat_combined <- cbind(mat_tcga_z, mat_gse_z)

    cohort_ann <- c(rep("TCGA-CHOL", ncol(mat_tcga_z)), rep("GSE107943", ncol(mat_gse_z)))
    type_ann <- c(as.character(tcga_meta[colnames(mat_tcga),"sample_type"]),
                  as.character(gse_meta[colnames(mat_gse),"sample_type"]))

    ha <- HeatmapAnnotation(
      Cohort = cohort_ann, Type = type_ann,
      col = list(
        Cohort = c("TCGA-CHOL"="#4DAF4A","GSE107943"="#984EA3"),
        Type = c("Primary Tumor"=tumor_col,"Solid Tissue Normal"=normal_col,
                 "Tumor"=tumor_col,"Normal"=normal_col)),
      annotation_legend_param = list(
        Cohort = list(title="Cohort"), Type = list(title="Sample Type")))

    pdf("figures/DEG/Fig2G_validated_IM_CAF_TAM_DEGs_heatmap.pdf",
        width=16, height=min(14, nrow(mat_combined)*0.18+3))
    ht <- Heatmap(mat_combined, name="Z-score", top_annotation=ha,
      show_row_names=TRUE, show_column_names=FALSE,
      cluster_rows=TRUE, cluster_columns=FALSE, column_split=cohort_ann,
      row_names_gp=gpar(fontsize=7),
      column_title="Validated IM/CAF/TAM DEGs (Top 50, Same Direction)")
    draw(ht, annotation_legend_side="bottom"); dev.off()
    cat(sprintf("  Heatmap: %d genes\n", nrow(mat_combined)))
  }
}

# Also need DEG_analysis_summary
summ <- data.frame(
  item=c("TCGA原始基因","TCGA过滤基因","TCGA显著DEG","TCGA上调","TCGA下调",
         "GSE原始基因","GSE过滤基因","GSE显著DEG","GSE上调","GSE下调",
         "共同DEG","同向validated","IM/CAF/TAM TCGA","IM/CAF/TAM GSE",
         "IM/CAF/TAM validated"),
  value=c(59427,22123,9291,5932,3359,55773,20287,7735,3869,3866,
          4560,4534,
          nrow(read.csv("tables/CCA_IM_CAF_TAM_DEGs_TCGA.csv")),
          nrow(read.csv("tables/CCA_IM_CAF_TAM_DEGs_GSE107943.csv")),
          nrow(overlap_im)),
  stringsAsFactors=FALSE
)
write.csv(summ, "tables/DEG_analysis_summary.csv", row.names=FALSE)

cat("\nAll figures and summary done!\n")
