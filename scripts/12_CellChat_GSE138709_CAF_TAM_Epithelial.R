# ============================================================
# 分析目的: GSE138709 CellChat 细胞通讯分析
#          — 重点: CAF ↔ TAM ↔ Epithelial_like 通讯网络
# 输入:
#   data/processed/GSE138709/GSE138709_seurat_processed.rds (已注释)
# 输出:
#   tables/cellchat/ (7 个表格)
#   figures/cellchat/ (8 张图)
# 主要方法: CellChat v2 (CellChatDB.human)
# ============================================================

library(CellChat)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(data.table)

dir.create("tables/cellchat", showWarnings=FALSE, recursive=TRUE)
dir.create("figures/cellchat", showWarnings=FALSE, recursive=TRUE)

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ============================================================
# Step 1: 加载 + 准备数据
# ============================================================
log_msg("=== Step 1: 加载 Seurat 对象 ===")

merged <- readRDS("data/processed/GSE138709/GSE138709_seurat_processed.rds")

# Check cell_type — if missing, annotate in-place
if (!"cell_type" %in% colnames(merged@meta.data)) {
  log_msg("cell_type not in meta.data — annotating...")

  # Join layers if needed
  if (length(Layers(merged[["RNA"]])) > 3) merged <- JoinLayers(merged)

  # Marker-based annotation
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

  avg_list <- list()
  for (nm in names(cell_markers)) {
    genes_use <- intersect(cell_markers[[nm]], rownames(merged))
    if (length(genes_use) >= 2) {
      avg <- tryCatch(
        AverageExpression(merged, features=genes_use, group.by="seurat_clusters",
                         assays="RNA", slot="data")$RNA,
        error=function(e) NULL)
      if (!is.null(avg)) avg_list[[nm]] <- avg
    }
  }

  cluster_scores <- sapply(avg_list, function(x) colMeans(as.matrix(x)))
  if (!is.matrix(cluster_scores)) cluster_scores <- t(cluster_scores)
  cluster_assignment <- apply(cluster_scores, 1, function(x) names(which.max(x)))
  cluster_max_score <- apply(cluster_scores, 1, max)
  cluster_assignment[cluster_max_score < 0.3] <- "Unclear"
  cluster_names <- gsub("^g", "", rownames(cluster_scores))
  names(cluster_assignment) <- cluster_names
  ct_vec <- cluster_assignment[as.character(merged$seurat_clusters)]
  names(ct_vec) <- colnames(merged)
  merged$cell_type <- ct_vec
  log_msg(sprintf("Annotated: %s", paste(unique(merged$cell_type), collapse=", ")))

  # Save annotated RDS
  saveRDS(merged, "data/processed/GSE138709/GSE138709_seurat_annotated.rds")
  log_msg("Annotated RDS saved")
}
log_msg(sprintf("Metadata columns: %s", paste(colnames(merged@meta.data), collapse=", ")))
log_msg(sprintf("Cell types: %s", paste(unique(merged$cell_type), collapse=", ")))

# Filter to tumor samples
if ("tissue" %in% colnames(merged@meta.data)) {
  log_msg(sprintf("Tissue types: %s", paste(unique(merged$tissue), collapse=", ")))
  merged_tumor <- subset(merged, subset = tissue == "Tumor")
  log_msg(sprintf("Tumor-only cells: %d", ncol(merged_tumor)))
} else {
  log_msg("No 'tissue' column — checking sample names...")
  # Check if sample contains "Tumor"
  if ("sample" %in% colnames(merged@meta.data)) {
    is_tumor <- grepl("Tumor", merged$sample)
    merged_tumor <- subset(merged, cells = colnames(merged)[is_tumor])
    log_msg(sprintf("Using %d cells with 'Tumor' in sample name", ncol(merged_tumor)))
  } else {
    merged_tumor <- merged
    log_msg("Using all cells (no tissue/sample filter available)")
  }
}

# Cell type counts for CellChat
ct_counts <- table(merged_tumor$cell_type)
log_msg("Cell type counts:")
print(ct_counts)

# Keep cell types with >= 30 cells
keep_ct <- names(ct_counts)[ct_counts >= 30]
log_msg(sprintf("Cell types for CellChat (>=30 cells): %s", paste(keep_ct, collapse=", ")))
merged_tumor <- subset(merged_tumor, subset = cell_type %in% keep_ct)

ct_counts_final <- table(merged_tumor$cell_type)
write.csv(data.frame(cell_type=names(ct_counts_final), n_cells=as.integer(ct_counts_final),
                     stringsAsFactors=FALSE),
          "tables/cellchat/GSE138709_cellchat_celltype_counts.csv", row.names=FALSE)

# Ensure RNA assay layers are joined
if (length(Layers(merged_tumor[["RNA"]])) > 3) {
  log_msg("Joining RNA layers...")
  merged_tumor <- JoinLayers(merged_tumor)
}

# ============================================================
# Step 2: 创建 CellChat 对象
# ============================================================
log_msg("\n=== Step 2: 创建 CellChat 对象 ===")

data_input <- GetAssayData(merged_tumor, assay="RNA", layer="data")
meta_input <- merged_tumor@meta.data

cellchat <- createCellChat(object=data_input, meta=meta_input, group.by="cell_type")
log_msg(sprintf("CellChat object created: %d genes x %d cells", nrow(cellchat@data), ncol(cellchat@data)))

# Set cell identities
cellchat <- setIdent(cellchat, ident.use="cell_type")
log_msg(sprintf("Cell groups: %s", paste(levels(cellchat@idents), collapse=", ")))

# ============================================================
# Step 3: 运行 CellChat 通路
# ============================================================
log_msg("\n=== Step 3: CellChat pipeline ===")

# Use human database
cellchat@DB <- CellChatDB.human
log_msg(sprintf("CellChatDB.human loaded: %d interactions",
  nrow(cellchat@DB$interaction)))

# Subset to secreted signaling + ECM-receptor + cell-cell contact
cellchat <- subsetData(cellchat)
log_msg("subsetData done")

# Identify over-expressed genes and interactions
cellchat <- identifyOverExpressedGenes(cellchat, do.fast=FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)
log_msg("Over-expressed interactions identified")

# Compute communication probabilities
cellchat <- computeCommunProb(cellchat, type="triMean")
log_msg("CommunProb computed")

# Filter: remove weak interactions
cellchat <- filterCommunication(cellchat, min.cells=10)

# Pathway-level
cellchat <- computeCommunProbPathway(cellchat)
log_msg("CommunProbPathway computed")

# Aggregate networks
cellchat <- aggregateNet(cellchat)
log_msg("AggregateNet computed")

# ============================================================
# Step 4: 保存交互和通路表格
# ============================================================
log_msg("\n=== Step 4: 保存表格 ===")

# Overall interactions
net_count <- cellchat@net$count
net_weight <- cellchat@net$weight
overall_df <- data.frame(
  source = rep(rownames(net_count), ncol(net_count)),
  target = rep(colnames(net_count), each=nrow(net_count)),
  interaction_count = as.vector(net_count),
  interaction_weight = as.vector(net_weight),
  stringsAsFactors = FALSE
)
overall_df <- overall_df[overall_df$interaction_count > 0, ]
write.csv(overall_df, "tables/cellchat/GSE138709_cellchat_overall_interactions.csv", row.names=FALSE)

# Pathway summary
pathway_df <- cellchat@netP$pathways
if (length(pathway_df) > 0) {
  pathway_summary <- data.frame(pathway=names(pathway_df), stringsAsFactors=FALSE)
  write.csv(pathway_summary, "tables/cellchat/GSE138709_cellchat_pathway_summary.csv", row.names=FALSE)
}

# Top outgoing/incoming signals
# Use centrality if available
tryCatch({
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name="netP")
  outgoing <- cellchat@netP$centrality$outdegree
  incoming <- cellchat@netP$centrality$indegree
  if (!is.null(outgoing)) {
    out_df <- as.data.frame(outgoing); out_df$cell_type <- rownames(out_df)
    write.csv(out_df, "tables/cellchat/GSE138709_top_outgoing_signals_by_celltype.csv", row.names=FALSE)
    in_df <- as.data.frame(incoming); in_df$cell_type <- rownames(in_df)
    write.csv(in_df, "tables/cellchat/GSE138709_top_incoming_signals_by_celltype.csv", row.names=FALSE)
  }
}, error=function(e) log_msg(sprintf("Centrality failed: %s", e$message)))

# ============================================================
# Step 5: Extract CAF-TAM-Epithelial LR pairs
# ============================================================
log_msg("\n=== Step 5: CAF-TAM-Epithelial LR pairs ===")

focus_types <- c("CAF_Fibroblast", "Macrophage_TAM", "Epithelial_like")
focus_types <- intersect(focus_types, levels(cellchat@idents))
log_msg(sprintf("Focus cell types: %s", paste(focus_types, collapse=", ")))

# Extract all LR pairs between focus cell types
lr_df_list <- list()
for (src in focus_types) {
  for (tgt in focus_types) {
    if (src == tgt) next  # skip autocrine for now
    tryCatch({
      lr_pairs <- subsetCommunication(cellchat, sources.use=src, targets.use=tgt)
      if (nrow(lr_pairs) > 0) {
        lr_pairs$source_cell <- src
        lr_pairs$target_cell <- tgt
        lr_df_list[[paste(src, tgt, sep="->")]] <- lr_pairs
      }
    }, error=function(e) NULL)
  }
}

if (length(lr_df_list) > 0) {
  all_lr <- do.call(rbind, lr_df_list)
  all_lr <- all_lr[order(-all_lr$prob), ]
  write.csv(all_lr, "tables/cellchat/GSE138709_CAF_TAM_Epithelial_LR_pairs.csv", row.names=FALSE)
  log_msg(sprintf("Total LR pairs between focus types: %d", nrow(all_lr)))

  # Top pairs
  cat("\nTop 10 LR pairs:\n")
  print(head(all_lr[, c("source_cell","target_cell","ligand","receptor","pathway_name","prob")], 10))
}

# ============================================================
# Step 6: 可视化
# ============================================================
log_msg("\n=== Step 6: 可视化 ===")

# ---- Fig8A: Interaction number network ----
log_msg("Fig8A: Interaction number")
pdf("figures/cellchat/Fig8A_CellChat_overall_interaction_number.pdf", width=8, height=7)
netVisual_circle(cellchat@net$count, vertex.weight=as.numeric(table(cellchat@idents)),
                 weight.scale=TRUE, label.edge=FALSE,
                 title.name="Interaction counts")
dev.off()

# ---- Fig8B: Interaction weight network ----
log_msg("Fig8B: Interaction weight")
pdf("figures/cellchat/Fig8B_CellChat_overall_interaction_weight.pdf", width=8, height=7)
netVisual_circle(cellchat@net$weight, vertex.weight=as.numeric(table(cellchat@idents)),
                 weight.scale=TRUE, label.edge=FALSE,
                 title.name="Interaction weights/strength")
dev.off()

# ---- Fig8C: CAF/TAM/Epithelial sub-network ----
log_msg("Fig8C: CAF-TAM-Epithelial network")
if (length(focus_types) >= 2) {
  pdf("figures/cellchat/Fig8C_CellChat_CAF_TAM_Epithelial_network.pdf", width=10, height=8)
  tryCatch({
    netVisual_circle(cellchat@net$weight,
                     sources.use=focus_types, targets.use=focus_types,
                     vertex.weight=as.numeric(table(cellchat@idents)[focus_types]),
                     weight.scale=TRUE,
                     title.name="CAF ↔ TAM ↔ Epithelial")
  }, error=function(e) {
    plot.new(); text(0.5, 0.5, paste("Network viz failed:", e$message), cex=0.8)
  })
  dev.off()
}

# ---- Fig8D: Outgoing/incoming heatmap ----
log_msg("Fig8D: Outgoing/incoming heatmap")
pdf("figures/cellchat/Fig8D_CellChat_outgoing_incoming_heatmap.pdf", width=12, height=8)
tryCatch({
  netAnalysis_signalingRole_heatmap(cellchat, pattern="outgoing", width=8, height=12)
}, error=function(e) {
  plot.new(); text(0.5, 0.5, paste("Heatmap failed:", e$message), cex=0.8)
})
dev.off()

# ---- Fig8E: Top pathways bubble ----
log_msg("Fig8E: Bubble plot")
pdf("figures/cellchat/Fig8E_CellChat_top_pathways_bubble.pdf", width=12, height=8)
tryCatch({
  pathways_show <- head(cellchat@netP$pathways, 15)
  netAnalysis_dot(cellchat, pattern="outgoing", signaling=pathways_show) +
    ggtitle("Top Signaling Pathways — All Cell Types")
}, error=function(e) {
  plot.new(); text(0.5, 0.5, paste("Bubble failed:", e$message), cex=0.8)
})
dev.off()

# ---- Fig8F: CAF → TAM bubble ----
log_msg("Fig8F: CAF → TAM bubble")
if ("CAF_Fibroblast" %in% focus_types && "Macrophage_TAM" %in% focus_types) {
  pdf("figures/cellchat/Fig8F_CellChat_CAF_to_TAM_bubble.pdf", width=10, height=7)
  tryCatch({
    netVisual_bubble(cellchat, sources.use="CAF_Fibroblast",
                     targets.use="Macrophage_TAM", remove.isolate=FALSE) +
      ggtitle("CAF_Fibroblast → Macrophage_TAM")
  }, error=function(e) {
    plot.new(); text(0.5, 0.5, paste("CAF→TAM bubble failed:", e$message), cex=0.8)
  })
  dev.off()
}

# ---- Fig8G: TAM → CAF bubble ----
log_msg("Fig8G: TAM → CAF bubble")
if ("Macrophage_TAM" %in% focus_types && "CAF_Fibroblast" %in% focus_types) {
  pdf("figures/cellchat/Fig8G_CellChat_TAM_to_CAF_bubble.pdf", width=10, height=7)
  tryCatch({
    netVisual_bubble(cellchat, sources.use="Macrophage_TAM",
                     targets.use="CAF_Fibroblast", remove.isolate=FALSE) +
      ggtitle("Macrophage_TAM → CAF_Fibroblast")
  }, error=function(e) {
    plot.new(); text(0.5, 0.5, paste("TAM→CAF bubble failed:", e$message), cex=0.8)
  })
  dev.off()
}

# ---- Fig8H: Key pathway chord/heatmap for CAF-TAM-Epithelial ----
log_msg("Fig8H: Key pathways")
pdf("figures/cellchat/Fig8H_CellChat_CAF_TAM_Epithelial_key_pathways.pdf", width=14, height=10)
tryCatch({
  key_pathways <- intersect(c("COLLAGEN","FN1","THBS","TGFb","SPP1","MIF",
                              "CCL","CXCL","GALECTIN","COMPLEMENT","IL1","TNF"),
                            cellchat@netP$pathways)
  if (length(key_pathways) > 0) {
    netAnalysis_dot(cellchat, pattern="outgoing", signaling=key_pathways,
                    sources.use=focus_types, targets.use=focus_types) +
      ggtitle("Key Pathways: CAF ↔ TAM ↔ Epithelial")
  } else {
    plot.new(); text(0.5, 0.5, "No key pathways found among focus cell types", cex=1)
  }
}, error=function(e) {
  plot.new(); text(0.5, 0.5, paste("Key pathways failed:", e$message), cex=0.8)
})
dev.off()

# ============================================================
# Step 7: Summary
# ============================================================
log_msg("\n=== Step 7: Summary ===")

# Find main outgoing/incoming cell types
total_out <- rowSums(cellchat@net$count)
total_in <- colSums(cellchat@net$count)
top_out <- names(sort(total_out, decreasing=TRUE))
top_in <- names(sort(total_in, decreasing=TRUE))

# CAF ↔ TAM connections
caf_to_tam <- tryCatch({
  lr <- subsetCommunication(cellchat, sources.use="CAF_Fibroblast", targets.use="Macrophage_TAM")
  nrow(lr)
}, error=function(e) 0)

tam_to_caf <- tryCatch({
  lr <- subsetCommunication(cellchat, sources.use="Macrophage_TAM", targets.use="CAF_Fibroblast")
  nrow(lr)
}, error=function(e) 0)

summary_rows <- data.frame(
  item = c(
    "CellChat 是否成功运行", "使用 tumor-only 分析",
    "纳入细胞数", "纳入细胞类型数",
    "主要 outgoing source", "主要 incoming receiver",
    "CAF→TAM LR pairs", "TAM→CAF LR pairs",
    "Total pathways identified", "Focus cell types"
  ),
  value = c(
    "Yes", if ("tissue" %in% colnames(merged@meta.data)) "Yes" else "All cells (no tissue filter)",
    as.character(ncol(merged_tumor)),
    as.character(length(levels(cellchat@idents))),
    paste(top_out[1:2], collapse=", "),
    paste(top_in[1:2], collapse=", "),
    as.character(caf_to_tam), as.character(tam_to_caf),
    as.character(length(cellchat@netP$pathways)),
    paste(focus_types, collapse=", ")
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_rows, "tables/cellchat/GSE138709_cellchat_analysis_summary.csv", row.names=FALSE)

cat(sprintf("\n========================================\n"))
cat(sprintf("CellChat 分析完成\n"))
cat(sprintf("========================================\n"))
cat(sprintf("Cells: %d\n", ncol(merged_tumor)))
cat(sprintf("Cell types: %s\n", paste(levels(cellchat@idents), collapse=", ")))
cat(sprintf("Pathways: %d\n", length(cellchat@netP$pathways)))
cat(sprintf("Top outgoing: %s\n", paste(top_out[1:3], collapse=", ")))
cat(sprintf("Top incoming: %s\n", paste(top_in[1:3], collapse=", ")))
cat(sprintf("CAF→TAM LR pairs: %d\n", caf_to_tam))
cat(sprintf("TAM→CAF LR pairs: %d\n", tam_to_caf))
cat(sprintf("========================================\n"))
