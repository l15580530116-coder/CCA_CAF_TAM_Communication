# Generate GSE107943 PCA figure
library(ggplot2); library(ggpubr)
tumor_col <- "#E41A1C"; normal_col <- "#377EB8"

gse_log2cpm <- read.csv("data/processed/GSE107943/GSE107943_log2CPM_gene_symbol.csv",
                        row.names=1, check.names=FALSE)
gse_meta <- read.csv("data/processed/GSE107943/GSE107943_sample_metadata.csv", stringsAsFactors=FALSE)
rownames(gse_meta) <- gse_meta$sample_id

# Drop 6N (all NA)
gse_log2cpm <- gse_log2cpm[, colnames(gse_log2cpm) != "6N", drop=FALSE]
cg <- intersect(colnames(gse_log2cpm), rownames(gse_meta))
gse_log2cpm <- gse_log2cpm[, cg, drop=FALSE]

gse_mat <- as.matrix(gse_log2cpm)
gse_sds <- apply(gse_mat, 1, sd, na.rm=TRUE)
gse_sds <- gse_sds[is.finite(gse_sds) & gse_sds > 0]
top_gse <- names(sort(gse_sds, decreasing=TRUE))[1:5000]
pca_mat <- gse_mat[top_gse, , drop=FALSE]
pca_mat <- pca_mat[complete.cases(pca_mat), , drop=FALSE]
cat(sprintf("PCA genes: %d\n", nrow(pca_mat)))

gse_pca <- prcomp(t(pca_mat), scale.=TRUE)
pca_df <- data.frame(PC1=gse_pca$x[,1], PC2=gse_pca$x[,2],
                     sample_type=gse_meta[colnames(pca_mat),"sample_type"])
pca_var <- round(100*summary(gse_pca)$importance[2,1:2],1)

pdf("figures/DEG/Fig2D_GSE107943_PCA_tumor_normal.pdf", width=7, height=6)
p <- ggplot(pca_df, aes(x=PC1, y=PC2, color=sample_type)) + geom_point(size=3) +
  stat_ellipse(level=0.95) +
  scale_color_manual(values=c(Tumor=tumor_col, Normal=normal_col)) +
  labs(x=sprintf("PC1 (%s%%)",pca_var[1]), y=sprintf("PC2 (%s%%)",pca_var[2]),
       title="GSE107943 PCA (log2 CPM+1)", color="") +
  theme_pubr(base_size=14) + theme(legend.position="top")
print(p); dev.off()
cat(sprintf("Fig2D Done: PC1=%s%%, PC2=%s%%\n", pca_var[1], pca_var[2]))
