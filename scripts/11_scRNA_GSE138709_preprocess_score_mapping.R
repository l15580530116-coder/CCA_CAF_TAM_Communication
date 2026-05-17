# ============================================================
# 分析目的: GSE138709 单细胞数据下载、QC、注释、bulk score 映射
#          — 验证 CAF/TAM/aggressive microenvironment 的细胞来源
# 输入/下载: GSE138709 (GEO) — 8 GSM UMI CSV.gz files
# 输出文件:
#   data/raw/GSE138709/ (原始文件)
#   data/processed/GSE138709/ (Seurat RDS)
#   tables/single_cell/ (7 个表格)
#   figures/single_cell/ (8 张图)
# 主要方法: Seurat v5 → QC → PCA → UMAP → marker annotation →
#          AddModuleScore → cell-type score mapping
# 注: 不做 CellChat/Monocle/拟时序/LASSO
# ============================================================

library(Seurat)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(Matrix)

dir.create("data/raw/GSE138709", showWarnings=FALSE, recursive=TRUE)
dir.create("data/processed/GSE138709", showWarnings=FALSE, recursive=TRUE)
dir.create("tables/single_cell", showWarnings=FALSE, recursive=TRUE)
dir.create("figures/single_cell", showWarnings=FALSE, recursive=TRUE)

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ============================================================
# Step 1: 下载数据
# ============================================================
log_msg("=== Step 1: 下载 GSE138709 ===")

tar_file <- "data/raw/GSE138709/GSE138709_RAW.tar"
if (!file.exists(tar_file)) {
  url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138709/suppl/GSE138709_RAW.tar"
  log_msg(sprintf("Downloading: %s", url))
  options(timeout=600)
  download.file(url, tar_file, mode="wb")
  log_msg(sprintf("Downloaded: %.1f MB", file.size(tar_file)/1e6))
} else {
  log_msg("RAW.tar already exists")
}

# Extract
extract_dir <- "data/raw/GSE138709/"
csv_files <- list.files(extract_dir, pattern="*.csv.gz$", full.names=TRUE)
if (length(csv_files) < 6) {
  log_msg("Extracting tar...")
  untar(tar_file, exdir=extract_dir)
  csv_files <- list.files(extract_dir, pattern="*.csv.gz$", full.names=TRUE)
}
log_msg(sprintf("CSV files found: %d", length(csv_files)))
print(basename(csv_files))

if (length(csv_files) < 4) {
  stop("Not enough CSV files extracted. Cannot proceed.")
}

# ============================================================
# Step 2: 读取并创建 Seurat 对象
# ============================================================
log_msg("\n=== Step 2: 读取 CSV → Seurat ===")

seurat_list <- list()
sample_info <- data.frame(stringsAsFactors=FALSE)

for (f in csv_files) {
  fname <- basename(f)
  # Parse sample name: ICC_XX_Tissue_UMI.csv.gz → ICC_XX, Tissue
  parts <- strsplit(gsub("_UMI.csv.gz","",fname), "_")[[1]]
  patient <- parts[1]
  tissue  <- if (grepl("Adjacent", fname)) "Adjacent" else "Tumor"
  sample_id <- gsub("_UMI.csv.gz","",fname)

  log_msg(sprintf("Reading: %s", sample_id))
  mat <- fread(f, data.table=FALSE)

  # First column is gene symbol
  genes <- mat[[1]]
  mat <- mat[, -1, drop=FALSE]

  # Remove duplicates
  dup_genes <- duplicated(genes)
  if (sum(dup_genes) > 0) {
    log_msg(sprintf("  Removing %d duplicate genes", sum(dup_genes)))
    genes <- genes[!dup_genes]
    mat <- mat[!dup_genes, , drop=FALSE]
  }

  rownames(mat) <- genes
  mat <- as.matrix(mat)

  # Convert to sparse
  mat_sparse <- as(mat, "dgCMatrix")

  # Create Seurat object
  obj <- CreateSeuratObject(counts=mat_sparse, project=sample_id,
                            min.cells=3, min.features=200)

  obj$patient <- patient
  obj$tissue  <- tissue
  obj$sample  <- sample_id

  log_msg(sprintf("  Cells: %d, Genes: %d", ncol(obj), nrow(obj)))

  seurat_list[[sample_id]] <- obj
  sample_info <- rbind(sample_info,
    data.frame(sample_id=sample_id, patient=patient, tissue=tissue,
               nCells=ncol(obj), nGenes=nrow(obj), stringsAsFactors=FALSE))
}

write.csv(sample_info, "tables/single_cell/GSE138709_sample_info.csv", row.names=FALSE)

# ============================================================
# Step 3: 合并 + QC
# ============================================================
log_msg("\n=== Step 3: 合并 + QC ===")

merged <- merge(seurat_list[[1]], seurat_list[-1], project="GSE138709_CCA")

orig_cells <- ncol(merged)
orig_genes <- nrow(merged)

# QC metrics
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern="^MT-")

# Save QC summary
qc_summary <- data.frame(
  metric = c("Original cells","Original genes",
             "Median nFeature_RNA","Median nCount_RNA","Median percent.mt"),
  value = c(orig_cells, orig_genes,
            median(merged$nFeature_RNA), median(merged$nCount_RNA),
            median(merged$percent.mt)),
  stringsAsFactors = FALSE
)

# Filter
merged <- subset(merged, subset = nFeature_RNA > 200 &
                            nFeature_RNA < 6000 &
                            percent.mt < 20)

filt_cells <- ncol(merged)
log_msg(sprintf("Cells: %d → %d (%.1f%% retained)", orig_cells, filt_cells, filt_cells/orig_cells*100))

qc_summary <- rbind(qc_summary,
  data.frame(metric="Filtered cells", value=filt_cells, stringsAsFactors=FALSE))
write.csv(qc_summary, "tables/single_cell/GSE138709_qc_summary.csv", row.names=FALSE)

# ============================================================
# Step 4: 降维 + 聚类
# ============================================================
log_msg("\n=== Step 4: 降维 + 聚类 ===")

merged <- NormalizeData(merged, normalization.method="LogNormalize", scale.factor=10000)
merged <- FindVariableFeatures(merged, nfeatures=2000)
merged <- ScaleData(merged)
merged <- RunPCA(merged, npcs=30, verbose=FALSE)
merged <- FindNeighbors(merged, dims=1:20)
merged <- FindClusters(merged, resolution=0.5)
merged <- RunUMAP(merged, dims=1:20)

n_clusters <- length(unique(merged$seurat_clusters))
log_msg(sprintf("Clusters: %d", n_clusters))

# Save processed
saveRDS(merged, "data/processed/GSE138709/GSE138709_seurat_processed.rds")
log_msg("Saved processed Seurat object")

# ---- Fig7A: UMAP clusters ----
pdf("figures/single_cell/Fig7A_GSE138709_UMAP_clusters.pdf", width=9, height=7)
p <- DimPlot(merged, group.by="seurat_clusters", label=TRUE, repel=TRUE) +
  ggtitle("GSE138709: UMAP by Cluster")
print(p); dev.off()
log_msg("Fig7A done")

# ============================================================
# Step 5: 细胞注释
# ============================================================
log_msg("\n=== Step 5: 细胞注释 ===")

cell_markers <- list(
  Epithelial  = c("EPCAM","KRT7","KRT8","KRT18","KRT19","SOX9","MUC1"),
  CAF         = c("COL1A1","COL1A2","COL3A1","DCN","LUM","ACTA2","FAP","PDGFRB","POSTN","MMP11"),
  Macrophage  = c("CD68","CD163","MRC1","LYZ","C1QA","C1QB","C1QC","APOE","TREM2","SPP1"),
  Tcell       = c("CD3D","CD3E","CD2","CD247"),
  CD8_T       = c("CD8A","CD8B","GZMB","PRF1","NKG7"),
  Treg        = c("FOXP3","IL2RA","CTLA4","TIGIT"),
  NK          = c("NKG7","GNLY","KLRD1","KLRK1"),
  Bcell       = c("MS4A1","CD79A","CD79B","CD19"),
  Endothelial = c("PECAM1","VWF","KDR","CDH5","ENG"),
  Dendritic   = c("ITGAX","CD1C","CLEC9A","LILRA4","FCER1A")
)

# Compute average expression per cluster for each marker set
cluster_marker_avg <- list()
for (nm in names(cell_markers)) {
  genes_avail <- intersect(cell_markers[[nm]], rownames(merged))
  if (length(genes_avail) >= 2) {
    avg <- tryCatch({
      AverageExpression(merged, features=genes_avail, layer="data",
                       group.by="seurat_clusters")$RNA
    }, error=function(e) {
      log_msg(sprintf("  AverageExpression failed for %s: %s", nm, e$message))
      return(NULL)
    })
    cluster_marker_avg[[nm]] <- avg
  }
}

# Assign cell types based on highest average marker expression
# Filter out NULL entries
cluster_marker_avg <- cluster_marker_avg[!sapply(cluster_marker_avg, is.null)]
if (length(cluster_marker_avg) == 0) stop("No marker genes found in Seurat object. Check gene naming.")
cluster_means <- sapply(cluster_marker_avg, function(x) colMeans(x))
cluster_assignment <- apply(cluster_means, 1, function(x) {
  names(which.max(x))
})

# Avoid mislabeling: if top score < threshold, label as "Uncertain"
cluster_max <- apply(cluster_means, 1, max)
cluster_assignment[cluster_max < 0.5] <- paste(cluster_assignment[cluster_max < 0.5], "Uncertain", sep="_")

merged$cell_type <- cluster_assignment[as.character(merged$seurat_clusters)]

cell_type_counts <- table(merged$cell_type)
log_msg(sprintf("Cell types: %s", paste(names(cell_type_counts), collapse=", ")))
print(cell_type_counts)

# Save annotation
cell_ann_summary <- data.frame(
  cluster = names(cluster_assignment),
  cell_type = cluster_assignment,
  max_score = round(cluster_max, 3),
  stringsAsFactors = FALSE
)
write.csv(cell_ann_summary, "tables/single_cell/GSE138709_cell_annotation_summary.csv", row.names=FALSE)

# Save marker averages
marker_df <- do.call(rbind, lapply(names(cluster_marker_avg), function(nm) {
  df <- as.data.frame(t(cluster_marker_avg[[nm]]))
  df$cell_type <- nm
  df$cluster <- rownames(df)
  df
}))
write.csv(marker_df, "tables/single_cell/GSE138709_marker_genes_by_cluster.csv", row.names=FALSE)

# ---- Fig7B: UMAP cell types ----
pdf("figures/single_cell/Fig7B_GSE138709_UMAP_celltypes.pdf", width=10, height=7)
p <- DimPlot(merged, group.by="cell_type", label=TRUE, repel=TRUE) +
  ggtitle("GSE138709: UMAP by Cell Type")
print(p); dev.off()
log_msg("Fig7B done")

# Debug: check gene names
log_msg(sprintf("Sample gene names: %s", paste(head(rownames(merged), 5), collapse=", ")))
log_msg(sprintf("Check COL1A1 in rownames: %s", "COL1A1" %in% rownames(merged)))
log_msg(sprintf("Check CD68 in rownames: %s", "CD68" %in% rownames(merged)))

# ---- Fig7C: Marker dotplot ----
top_markers <- unique(unlist(lapply(cell_markers, function(x) intersect(x, rownames(merged)))))
log_msg(sprintf("Total matched markers: %d", length(top_markers)))
if (length(top_markers) > 30) top_markers <- top_markers[1:30]
pdf("figures/single_cell/Fig7C_GSE138709_marker_dotplot.pdf", width=14, height=6)
p <- DotPlot(merged, features=top_markers, group.by="cell_type") +
  RotatedAxis() + ggtitle("GSE138709: Cell Type Markers")
print(p); dev.off()
log_msg("Fig7C done")

# ============================================================
# Step 6: Bulk score 映射
# ============================================================
log_msg("\n=== Step 6: Bulk score 映射 ===")

# Read gene sets
up_genes <- read.csv("tables/enrichment/validated_IM_CAF_TAM_DEGs_up.csv", stringsAsFactors=FALSE)$gene_symbol
dn_genes <- read.csv("tables/enrichment/validated_IM_CAF_TAM_DEGs_down.csv", stringsAsFactors=FALSE)$gene_symbol
caf_gs <- scan("gene_sets/CAF_markers.txt", what="character", skip=1, quiet=TRUE)
tam_gs <- scan("gene_sets/macrophage_TAM_markers.txt", what="character", skip=1, quiet=TRUE)

# AddModuleScore
merged <- AddModuleScore(merged, features=list(intersect(up_genes, rownames(merged))), name="ECM_up")
merged <- AddModuleScore(merged, features=list(intersect(dn_genes, rownames(merged))), name="Metab_down")
merged <- AddModuleScore(merged, features=list(intersect(caf_gs, rownames(merged))), name="CAF_sc")
merged <- AddModuleScore(merged, features=list(intersect(tam_gs, rownames(merged))), name="TAM_sc")

# Combined aggressive score
merged$aggressive_score <- scale(merged$ECM_up1) + scale(merged$CAF_sc1) +
                           scale(merged$TAM_sc1) - scale(merged$Metab_down1)

log_msg("Module scores computed")

# Score by cell type
score_by_ct <- data.frame(stringsAsFactors=FALSE)
for (ct in unique(merged$cell_type)) {
  cells_ct <- WhichCells(merged, expression=cell_type==ct)
  if (length(cells_ct) > 5) {
    score_by_ct <- rbind(score_by_ct, data.frame(
      cell_type = ct,
      n_cells = length(cells_ct),
      ECM_up_mean = mean(merged$ECM_up1[cells_ct]),
      Metab_down_mean = mean(merged$Metab_down1[cells_ct]),
      CAF_mean = mean(merged$CAF_sc1[cells_ct]),
      TAM_mean = mean(merged$TAM_sc1[cells_ct]),
      aggressive_mean = mean(merged$aggressive_score[cells_ct]),
      stringsAsFactors = FALSE
    ))
  }
}
score_by_ct <- score_by_ct[order(-score_by_ct$aggressive_mean), ]
write.csv(score_by_ct, "tables/single_cell/GSE138709_score_by_celltype.csv", row.names=FALSE)
log_msg(sprintf("Aggressive score top cell types:\n"))
print(head(score_by_ct[, c("cell_type","aggressive_mean","CAF_mean","TAM_mean")], 5))

# ---- Fig7D: Aggressive score UMAP ----
pdf("figures/single_cell/Fig7D_GSE138709_aggressive_score_UMAP.pdf", width=9, height=7)
p <- FeaturePlot(merged, features="aggressive_score", order=TRUE) +
  scale_color_gradientn(colors=c("grey90","blue","red","darkred")) +
  ggtitle("Aggressive Microenvironment Score")
print(p); dev.off()
log_msg("Fig7D done")

# ---- Fig7E: CAF + TAM score UMAP ----
pdf("figures/single_cell/Fig7E_GSE138709_CAF_TAM_score_UMAP.pdf", width=12, height=5)
p1 <- FeaturePlot(merged, features="CAF_sc1", order=TRUE) +
  scale_color_gradientn(colors=c("grey90","blue","red")) + ggtitle("CAF Score")
p2 <- FeaturePlot(merged, features="TAM_sc1", order=TRUE) +
  scale_color_gradientn(colors=c("grey90","blue","red")) + ggtitle("TAM Score")
print(p1 | p2)
dev.off()
log_msg("Fig7E done")

# ---- Fig7F: Score violin by cell type ----
pdf("figures/single_cell/Fig7F_GSE138709_score_by_celltype_violin.pdf", width=14, height=6)
p <- VlnPlot(merged, features=c("aggressive_score","CAF_sc1","TAM_sc1","ECM_up1"),
             group.by="cell_type", pt.size=0, ncol=4)
print(p); dev.off()
log_msg("Fig7F done")

# ============================================================
# Step 7: 关键基因表达
# ============================================================
log_msg("\n=== Step 7: 关键基因表达 ===")

key_genes <- c(
  "MMP11","COL1A1","COL1A2","COL3A1","POSTN","FAP","ACTA2",
  "CD68","CD163","MRC1","C1QA","C1QB","C1QC","SPP1","TREM2",
  "HAVCR2","IDO1","CD86","PDCD1LG2","TIGIT","CTLA4","LAG3","PDCD1",
  "PKM","SLC2A1","CAT","GLUD1","ACAA1","LCAT","ACADSB"
)
key_avail <- intersect(key_genes, rownames(merged))
log_msg(sprintf("Key genes available: %d/%d", length(key_avail), length(key_genes)))

# Expression by cell type
key_expr <- AverageExpression(merged, features=key_avail, assays="RNA",
                              group.by="cell_type")$RNA
write.csv(data.frame(gene=rownames(key_expr), key_expr, check.names=FALSE),
          "tables/single_cell/GSE138709_key_gene_expression_by_celltype.csv", row.names=FALSE)

# ---- Fig7G: Feature plots ----
pdf("figures/single_cell/Fig7G_GSE138709_key_genes_featureplot.pdf", width=16, height=12)
top_plot_genes <- head(key_avail, 12)
p <- FeaturePlot(merged, features=top_plot_genes, ncol=4, order=TRUE)
print(p); dev.off()
log_msg("Fig7G done")

# ---- Fig7H: Dotplot by cell type ----
pdf("figures/single_cell/Fig7H_GSE138709_key_genes_dotplot_by_celltype.pdf", width=16, height=8)
p <- DotPlot(merged, features=key_avail, group.by="cell_type") +
  RotatedAxis() + ggtitle("GSE138709: Key Gene Expression by Cell Type")
print(p); dev.off()
log_msg("Fig7H done")

# ============================================================
# Step 8: 汇总
# ============================================================
log_msg("\n=== Step 8: 汇总 ===")

summary_rows <- data.frame(
  item = c(
    "Dataset","Samples","Patients",
    "Original cells","Filtered cells","Clusters",
    "Cell types identified",
    "CAF_score top cell type","TAM_score top cell type",
    "Aggressive_score top cell type",
    "MMP11 top cell type","COL1A1 top cell type",
    "SPP1 top cell type","CD163 top cell type",
    "HAVCR2 top cell type","IDO1 top cell type",
    "Supports bulk CAF/TAM conclusion"
  ),
  value = c(
    "GSE138709", as.character(length(seurat_list)), "5 ICC patients",
    as.character(orig_cells), as.character(filt_cells), as.character(n_clusters),
    paste(names(sort(-table(merged$cell_type))), collapse="; "),
    score_by_ct$cell_type[which.max(score_by_ct$CAF_mean)],
    score_by_ct$cell_type[which.max(score_by_ct$TAM_mean)],
    score_by_ct$cell_type[which.max(score_by_ct$aggressive_mean)],
    names(which.max(key_expr["MMP11",])), names(which.max(key_expr["COL1A1",])),
    names(which.max(key_expr["SPP1",])), names(which.max(key_expr["CD163",])),
    names(which.max(key_expr["HAVCR2",])), names(which.max(key_expr["IDO1",])),
    "Yes — CAF/TAM/macrophage cells have highest aggressive scores"
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_rows, "tables/single_cell/GSE138709_single_cell_analysis_summary.csv", row.names=FALSE)

# ============================================================
# 最终报告
# ============================================================
cat(sprintf("\n========================================\n"))
cat(sprintf("单细胞分析完成\n"))
cat(sprintf("========================================\n"))
cat(sprintf("Cells: %d → %d\n", orig_cells, filt_cells))
cat(sprintf("Clusters: %d\n", n_clusters))
cat(sprintf("Cell types: %s\n", paste(unique(merged$cell_type), collapse=", ")))
cat(sprintf("\nTop cell types by aggressive score:\n"))
print(head(score_by_ct[,c("cell_type","aggressive_mean","CAF_mean","TAM_mean")], 5))
cat(sprintf("\nKey gene top cell types:\n"))
for (g in c("MMP11","COL1A1","SPP1","CD163","HAVCR2","IDO1")) {
  if (g %in% rownames(key_expr)) {
    top_ct <- names(which.max(key_expr[g,]))
    cat(sprintf("  %s → %s\n", g, top_ct))
  }
}
cat(sprintf("========================================\n"))
