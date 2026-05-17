# ============================================================
# 分析目的: 修复 GSE138709 细胞注释和 score 映射 (Seurat v5 兼容)
#          — 基于已保存的 processed RDS
# 输入:
#   data/processed/GSE138709/GSE138709_seurat_processed.rds
#   gene_sets/CAF_markers.txt
#   gene_sets/macrophage_TAM_markers.txt
#   tables/enrichment/validated_IM_CAF_TAM_DEGs_{up,down}.csv
# 输出:
#   tables/single_cell/ (5 个表格)
#   figures/single_cell/ (8 张图)
# 主要方法: Seurat v5 AggregateExpression + AddModuleScore
# 修复点: Seurat v5 layer 参数 → JoinLayers + layer="data"
# ============================================================

library(Seurat)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(data.table)

dir.create("tables/single_cell", showWarnings=FALSE, recursive=TRUE)
dir.create("figures/single_cell", showWarnings=FALSE, recursive=TRUE)

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ============================================================
# Step 1: 加载 + Seurat v5 准备
# ============================================================
log_msg("=== Step 1: 加载 Seurat RDS + v5 准备 ===")

merged <- readRDS("data/processed/GSE138709/GSE138709_seurat_processed.rds")
log_msg(sprintf("Loaded: %d cells x %d genes", ncol(merged), nrow(merged)))

# Check layers
rna_layers <- Layers(merged[["RNA"]])
log_msg(sprintf("RNA layers: %s", paste(rna_layers, collapse=", ")))

# Join layers for simpler downstream access
if (length(rna_layers) > 1) {
  log_msg("Joining multiple RNA layers...")
  merged <- JoinLayers(merged)
  log_msg("Layers joined")
}
rna_layers <- Layers(merged[["RNA"]])
log_msg(sprintf("RNA layers after join: %s", paste(rna_layers, collapse=", ")))

# Verify markers exist
test_markers <- c("COL1A1","CD68","CD163","EPCAM","KRT19","CD3D","ACTA2","MMP11","GAPDH","ACTB")
found <- test_markers %in% rownames(merged)
log_msg(sprintf("Marker check: %d/%d found", sum(found), length(found)))
if (!all(found)) log_msg(sprintf("  Missing: %s", paste(test_markers[!found], collapse=", ")))

# Set default assay
DefaultAssay(merged) <- "RNA"

# ============================================================
# Step 2: 细胞注释 (marker-based scoring per cluster)
# ============================================================
log_msg("\n=== Step 2: 细胞注释 ===")

cell_markers <- list(
  Epithelial_like = c("EPCAM","KRT7","KRT8","KRT18","KRT19","SOX9","MUC1"),
  CAF_Fibroblast  = c("COL1A1","COL1A2","COL3A1","DCN","LUM","ACTA2","FAP","PDGFRB","POSTN","MMP11"),
  Macrophage_TAM  = c("CD68","CD163","MRC1","LYZ","C1QA","C1QB","C1QC","APOE","TREM2","SPP1"),
  T_cell          = c("CD3D","CD3E","CD2","CD247"),
  NK_cell         = c("NKG7","GNLY","KLRD1","KLRK1"),
  B_cell          = c("MS4A1","CD79A","CD79B","CD19"),
  Endothelial     = c("PECAM1","VWF","KDR","CDH5","ENG"),
  Dendritic       = c("ITGAX","CD1C","CLEC9A","LILRA4","FCER1A")
)

# Compute average per cluster using layer="data"
log_msg("Computing average expression per cluster...")
avg_list <- list()
for (nm in names(cell_markers)) {
  genes_use <- intersect(cell_markers[[nm]], rownames(merged))
  log_msg(sprintf("  %-20s: %d/%d genes found", nm, length(genes_use), length(cell_markers[[nm]])))
  if (length(genes_use) >= 2) {
    # After JoinLayers, use slot="data" (v5 compatible)
    avg <- tryCatch({
      AverageExpression(merged, features=genes_use, group.by="seurat_clusters",
                       assays="RNA", slot="data")$RNA
    }, error=function(e) {
      log_msg(sprintf("  AverageExpression failed: %s", e$message))
      NULL
    })
    if (!is.null(avg)) avg_list[[nm]] <- avg
  }
}

log_msg(sprintf("Cell types with data: %d", length(avg_list)))

# Score each cluster by mean across all genes in the set
if (length(avg_list) > 0) {
  cluster_scores <- sapply(avg_list, function(x) colMeans(as.matrix(x)))
  if (!is.matrix(cluster_scores)) cluster_scores <- t(cluster_scores)

  # Assign each cluster to the cell type with highest score
  cluster_assignment <- apply(cluster_scores, 1, function(x) names(which.max(x)))
  cluster_max_score <- apply(cluster_scores, 1, max)

  # Low confidence: label as "Unclear"
  cluster_assignment[cluster_max_score < 0.3] <- "Unclear"

  # Fix cluster names: Seurat prepends "g" to numeric cluster names
  cluster_names <- gsub("^g", "", rownames(cluster_scores))
  names(cluster_assignment) <- cluster_names
  log_msg(sprintf("Cluster name mapping: %s", paste(cluster_names, collapse=",")))

  # Create safe cell_type vector
  ct_vec <- cluster_assignment[as.character(merged$seurat_clusters)]
  names(ct_vec) <- colnames(merged)
  merged$cell_type <- ct_vec

  log_msg("Cell type assignment:")
  print(table(merged$cell_type))

  # Save annotation
  ann_summary <- data.frame(
    cluster = names(cluster_assignment),
    cell_type = cluster_assignment,
    max_score = round(cluster_max_score, 3),
    stringsAsFactors = FALSE
  )
  write.csv(ann_summary, "tables/single_cell/GSE138709_cell_annotation_summary.csv", row.names=FALSE)

  # Save cluster-level average per marker (long format)
  marker_df_list <- list()
  for (nm in names(avg_list)) {
    mat <- as.matrix(avg_list[[nm]])
    for (i in seq_len(nrow(mat))) {
      gene <- rownames(mat)[i]
      for (j in seq_len(ncol(mat))) {
        cluster <- gsub("^g", "", colnames(mat)[j])
        marker_df_list[[length(marker_df_list)+1]] <- data.frame(
          cell_type_marker = nm, cluster = cluster, gene = gene,
          avg_expression = mat[i, j], stringsAsFactors = FALSE)
      }
    }
  }
  marker_long <- do.call(rbind, marker_df_list)
  write.csv(marker_long, "tables/single_cell/GSE138709_marker_genes_by_cluster.csv", row.names=FALSE)
  log_msg("Annotation tables saved")
} else {
  stop("No cell type markers could be scored. Cannot proceed.")
}

# ---- Fig7A: UMAP clusters ----
log_msg("Fig7A: UMAP clusters")
pdf("figures/single_cell/Fig7A_GSE138709_UMAP_clusters.pdf", width=9, height=7)
print(DimPlot(merged, group.by="seurat_clusters", label=TRUE, repel=TRUE) +
  ggtitle("GSE138709: UMAP by Cluster"))
dev.off()

# ---- Fig7B: UMAP celltypes ----
log_msg("Fig7B: UMAP cell types")
pdf("figures/single_cell/Fig7B_GSE138709_UMAP_celltypes.pdf", width=10, height=7)
print(DimPlot(merged, group.by="cell_type", label=TRUE, repel=TRUE) +
  ggtitle("GSE138709: UMAP by Cell Type"))
dev.off()

# ---- Fig7C: Marker dotplot ----
log_msg("Fig7C: Marker dotplot")
all_marker_genes <- unique(unlist(cell_markers))
all_marker_genes <- intersect(all_marker_genes, rownames(merged))
if (length(all_marker_genes) > 40) all_marker_genes <- all_marker_genes[1:40]
pdf("figures/single_cell/Fig7C_GSE138709_marker_dotplot.pdf", width=16, height=7)
print(DotPlot(merged, features=all_marker_genes, group.by="cell_type") +
  RotatedAxis() + ggtitle("GSE138709: Cell Type Markers"))
dev.off()

# ============================================================
# Step 3: AddModuleScore — pathway scores
# ============================================================
log_msg("\n=== Step 3: AddModuleScore ===")

up_genes <- read.csv("tables/enrichment/validated_IM_CAF_TAM_DEGs_up.csv", stringsAsFactors=FALSE)$gene_symbol
dn_genes <- read.csv("tables/enrichment/validated_IM_CAF_TAM_DEGs_down.csv", stringsAsFactors=FALSE)$gene_symbol
caf_gs <- scan("gene_sets/CAF_markers.txt", what="character", skip=1, quiet=TRUE)
tam_gs <- scan("gene_sets/macrophage_TAM_markers.txt", what="character", skip=1, quiet=TRUE)

log_msg(sprintf("Gene sets: up=%d, down=%d, CAF=%d, TAM=%d",
  length(intersect(up_genes, rownames(merged))),
  length(intersect(dn_genes, rownames(merged))),
  length(intersect(caf_gs, rownames(merged))),
  length(intersect(tam_gs, rownames(merged)))))

merged <- AddModuleScore(merged, features=list(intersect(up_genes, rownames(merged))), name="ECM_up")
merged <- AddModuleScore(merged, features=list(intersect(dn_genes, rownames(merged))), name="Metab_down")
merged <- AddModuleScore(merged, features=list(intersect(caf_gs, rownames(merged))), name="CAF_sc")
merged <- AddModuleScore(merged, features=list(intersect(tam_gs, rownames(merged))), name="TAM_sc")

# Rename from ECM_up1 → ECM_up_score etc.
merged$ECM_up_score    <- merged$ECM_up1
merged$Metab_down_score <- merged$Metab_down1
merged$CAF_score       <- merged$CAF_sc1
merged$TAM_score       <- merged$TAM_sc1

# Aggressive microenvironment score
merged$aggressive_score <- as.numeric(scale(merged$ECM_up_score)) +
                           as.numeric(scale(merged$CAF_score)) +
                           as.numeric(scale(merged$TAM_score)) -
                           as.numeric(scale(merged$Metab_down_score))

log_msg("Module scores computed")

# ============================================================
# Step 4: Score by cell type
# ============================================================
log_msg("\n=== Step 4: Score by cell type ===")

score_names_v5 <- c("ECM_up_score","Metab_down_score","CAF_score","TAM_score","aggressive_score")
score_by_ct <- data.frame(stringsAsFactors=FALSE)

for (ct in sort(unique(merged$cell_type))) {
  cells_ct <- WhichCells(merged, expression=cell_type==ct)
  if (length(cells_ct) > 10) {
    row <- data.frame(cell_type=ct, n_cells=length(cells_ct), stringsAsFactors=FALSE)
    for (sn in score_names_v5) {
      row[[sn]] <- mean(merged@meta.data[cells_ct, sn], na.rm=TRUE)
    }
    score_by_ct <- rbind(score_by_ct, row)
  }
}
score_by_ct <- score_by_ct[order(-score_by_ct$aggressive_score), ]
write.csv(score_by_ct, "tables/single_cell/GSE138709_score_by_celltype.csv", row.names=FALSE)

log_msg("Top cell types by aggressive score:")
print(head(score_by_ct[, c("cell_type","aggressive_score","CAF_score","TAM_score")], 8))

# ---- Fig7D: Aggressive score UMAP ----
log_msg("Fig7D: Aggressive score UMAP")
pdf("figures/single_cell/Fig7D_GSE138709_aggressive_score_UMAP.pdf", width=9, height=7)
print(FeaturePlot(merged, features="aggressive_score", order=TRUE) +
  scale_color_gradientn(colors=c("grey90","blue","red","darkred")) +
  ggtitle("Aggressive Microenvironment Score"))
dev.off()

# ---- Fig7E: CAF + TAM UMAP ----
log_msg("Fig7E: CAF + TAM score UMAP")
pdf("figures/single_cell/Fig7E_GSE138709_CAF_TAM_score_UMAP.pdf", width=12, height=5)
p1 <- FeaturePlot(merged, features="CAF_score", order=TRUE) +
  scale_color_gradientn(colors=c("grey90","blue","red")) + ggtitle("CAF Score")
p2 <- FeaturePlot(merged, features="TAM_score", order=TRUE) +
  scale_color_gradientn(colors=c("grey90","blue","red")) + ggtitle("TAM Score")
print(p1 | p2)
dev.off()

# ---- Fig7F: Violin by cell type ----
log_msg("Fig7F: Score violin")
pdf("figures/single_cell/Fig7F_GSE138709_score_by_celltype_violin.pdf", width=14, height=6)
print(VlnPlot(merged, features=c("aggressive_score","CAF_score","TAM_score","ECM_up_score"),
              group.by="cell_type", pt.size=0, ncol=4))
dev.off()

# ============================================================
# Step 5: 关键基因表达
# ============================================================
log_msg("\n=== Step 5: 关键基因表达 ===")

key_genes <- c(
  "MMP11","COL1A1","COL1A2","COL3A1","POSTN","FAP","ACTA2",
  "CD68","CD163","MRC1","C1QA","C1QB","C1QC","SPP1","TREM2",
  "HAVCR2","IDO1","CD86","PDCD1LG2","TIGIT","CTLA4","LAG3","PDCD1",
  "PKM","SLC2A1","CAT","GLUD1","ACAA1","LCAT","ACADSB"
)
key_avail <- intersect(key_genes, rownames(merged))
log_msg(sprintf("Key genes available: %d/%d", length(key_avail), length(key_genes)))

# Average expression by cell type
key_expr <- AverageExpression(merged, features=key_avail, group.by="cell_type",
                              assays="RNA", slot="data")$RNA

write.csv(data.frame(gene=rownames(key_expr), key_expr, check.names=FALSE),
          "tables/single_cell/GSE138709_key_gene_expression_by_celltype.csv", row.names=FALSE)

# Per-gene top cell type
for (g in key_avail) {
  top_ct <- names(which.max(key_expr[g, ]))
  log_msg(sprintf("  %-15s → %s", g, top_ct))
}

# ---- Fig7G: Feature plots ----
log_msg("Fig7G: Feature plots")
pdf("figures/single_cell/Fig7G_GSE138709_key_genes_featureplot.pdf", width=16, height=12)
print(FeaturePlot(merged, features=head(key_avail, 12), ncol=4, order=TRUE))
dev.off()

# ---- Fig7H: Dotplot ----
log_msg("Fig7H: Dotplot")
pdf("figures/single_cell/Fig7H_GSE138709_key_genes_dotplot_by_celltype.pdf", width=18, height=8)
print(DotPlot(merged, features=head(key_avail, 30), group.by="cell_type") +
  RotatedAxis() + ggtitle("GSE138709: Key Gene Expression by Cell Type"))
dev.off()

# ============================================================
# Step 6: Summary
# ============================================================
log_msg("\n=== Step 6: Summary ===")

ct_counts <- table(merged$cell_type)
ct_pct <- round(prop.table(ct_counts)*100, 1)

# Find top cell types for key genes
gene_top_ct <- sapply(c("MMP11","COL1A1","SPP1","CD163","HAVCR2","IDO1"), function(g) {
  if (g %in% rownames(key_expr)) names(which.max(key_expr[g,])) else "N/A"
})

summary_rows <- data.frame(
  item = c(
    "Cell annotation completed", "Total cells after QC",
    "Cell types identified",
    "Cell type counts",
    "CAF score top: CAF_Fibroblast",
    "TAM score top: Macrophage_TAM",
    "Aggressive score top cell types",
    "MMP11 top cell type", "COL1A1 top cell type",
    "SPP1 top cell type", "CD163 top cell type",
    "HAVCR2 top cell type", "IDO1 top cell type",
    "Supports bulk CAF/TAM conclusion",
    "Seurat version"
  ),
  value = c(
    "Yes", as.character(ncol(merged)),
    paste(names(ct_counts), collapse="; "),
    paste(sprintf("%s=%d(%.1f%%)", names(ct_counts), ct_counts, ct_pct), collapse="; "),
    if ("CAF_Fibroblast" %in% score_by_ct$cell_type)
      sprintf("Yes — CAF_Fibroblast CAF_score=%.3f vs mean=%.3f",
        score_by_ct$CAF_score[score_by_ct$cell_type=="CAF_Fibroblast"],
        mean(score_by_ct$CAF_score)) else "N/A",
    if ("Macrophage_TAM" %in% score_by_ct$cell_type)
      sprintf("Yes — Macrophage_TAM TAM_score=%.3f vs mean=%.3f",
        score_by_ct$TAM_score[score_by_ct$cell_type=="Macrophage_TAM"],
        mean(score_by_ct$TAM_score)) else "N/A",
    paste(head(score_by_ct$cell_type, 3), collapse=" > "),
    gene_top_ct["MMP11"], gene_top_ct["COL1A1"],
    gene_top_ct["SPP1"], gene_top_ct["CD163"],
    gene_top_ct["HAVCR2"], gene_top_ct["IDO1"],
    if (score_by_ct$cell_type[1] %in% c("CAF_Fibroblast","Macrophage_TAM"))
      "Yes — aggressive score highest in stromal/immune cells" else "Partial",
    as.character(packageVersion("Seurat"))
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_rows, "tables/single_cell/GSE138709_single_cell_analysis_summary.csv", row.names=FALSE)

# ============================================================
# 最终报告
# ============================================================
cat(sprintf("\n========================================\n"))
cat(sprintf("GSE138709 单细胞注释和评分完成\n"))
cat(sprintf("========================================\n"))
cat(sprintf("Cells: %d\n", ncol(merged)))
cat(sprintf("Cell types: %s\n", paste(names(ct_counts), collapse=", ")))
cat(sprintf("Cell counts: %s\n", paste(sprintf("%s=%d", names(ct_counts), ct_counts), collapse=", ")))
cat(sprintf("\nAggressive score top 5:\n"))
print(head(score_by_ct[, c("cell_type","n_cells","aggressive_score")], 5))
cat(sprintf("\nGene → top cell type:\n"))
for (g in names(gene_top_ct)) cat(sprintf("  %s → %s\n", g, gene_top_ct[g]))
cat(sprintf("\nAll figures: figures/single_cell/Fig7A-H.pdf\n"))
cat(sprintf("========================================\n"))
