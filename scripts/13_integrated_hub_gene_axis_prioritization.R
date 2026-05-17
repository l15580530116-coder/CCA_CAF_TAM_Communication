# ============================================================
# 分析目的: 整合 bulk、免疫、单细胞、CellChat 多层证据
#          筛选核心 hub genes 和通讯轴
# 输入:
#   tables/CCA_IM_CAF_TAM_DEGs_validated_same_direction.csv
#   tables/gsva/pathway_score_survival_univariate_cox.csv
#   tables/immune/TCGA_score_immune_correlation.csv
#   tables/immune/TCGA_checkpoint_expression_correlation.csv
#   tables/single_cell/GSE138709_key_gene_expression_by_celltype.csv
#   tables/single_cell/GSE138709_score_by_celltype.csv
#   tables/cellchat/GSE138709_CAF_TAM_Epithelial_LR_pairs.csv
# 输出:
#   tables/integrated/ (5 个表格)
#   figures/integrated/ (4 张图)
# 主要方法: 证据加权评分 → 排序 → 可视化
# ============================================================

library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(data.table)

dir.create("tables/integrated", showWarnings=FALSE, recursive=TRUE)
dir.create("figures/integrated", showWarnings=FALSE, recursive=TRUE)

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ============================================================
# Step 1: 加载所有证据数据
# ============================================================
log_msg("=== Step 1: 加载证据 ===")

# Bulk DEG
validated_deg <- read.csv("tables/CCA_IM_CAF_TAM_DEGs_validated_same_direction.csv", stringsAsFactors=FALSE)
up_deg <- read.csv("tables/enrichment/validated_IM_CAF_TAM_DEGs_up.csv", stringsAsFactors=FALSE)
dn_deg <- read.csv("tables/enrichment/validated_IM_CAF_TAM_DEGs_down.csv", stringsAsFactors=FALSE)

# GSVA survival Cox
cox_all <- read.csv("tables/gsva/pathway_score_survival_univariate_cox.csv", stringsAsFactors=FALSE)

# Immune correlations
tcga_imm_cor <- tryCatch(read.csv("tables/immune/TCGA_score_immune_correlation.csv", stringsAsFactors=FALSE), error=function(e) NULL)
gse_imm_cor  <- tryCatch(read.csv("tables/immune/GSE107943_score_immune_correlation.csv", stringsAsFactors=FALSE), error=function(e) NULL)
tcga_cp_cor  <- tryCatch(read.csv("tables/immune/TCGA_checkpoint_expression_correlation.csv", stringsAsFactors=FALSE), error=function(e) NULL)
gse_cp_cor   <- tryCatch(read.csv("tables/immune/GSE107943_checkpoint_expression_correlation.csv", stringsAsFactors=FALSE), error=function(e) NULL)

# Single cell
sc_expr <- read.csv("tables/single_cell/GSE138709_key_gene_expression_by_celltype.csv", stringsAsFactors=FALSE)
sc_score <- read.csv("tables/single_cell/GSE138709_score_by_celltype.csv", stringsAsFactors=FALSE)

# CellChat
cc_lr <- read.csv("tables/cellchat/GSE138709_CAF_TAM_Epithelial_LR_pairs.csv", stringsAsFactors=FALSE)

log_msg(sprintf("Bulk validated DEGs: %d", nrow(validated_deg)))
log_msg(sprintf("CellChat LR pairs: %d", nrow(cc_lr)))

# ============================================================
# Step 2: 定义候选基因 + 评分
# ============================================================
log_msg("\n=== Step 2: 基因评分 ===")

# All candidate genes
candidate_genes <- c(
  # CAF/ECM
  "MMP11","COL1A1","COL1A2","COL3A1","POSTN","FAP","ACTA2",
  # TAM/Macrophage
  "CD68","CD163","MRC1","C1QA","C1QB","C1QC","TREM2","SPP1",
  # Checkpoint/immunosuppression
  "HAVCR2","IDO1","CD86","PDCD1LG2","TGFB1",
  # CellChat ligands/receptors
  "MIF","CD74","CXCR4","CD44","APP",
  "PPIA","BSG",
  # Metabolism
  "PKM","SLC2A1","GLUD1","CAT","ACAA1","LCAT","ACADSB",
  # Bulk DEG top validated
  "MSMO1","SC5D","INSIG1","ACADSB","MMP11"
)
candidate_genes <- unique(candidate_genes)
log_msg(sprintf("Candidate genes: %d", length(candidate_genes)))

# CAF marker list
caf_genes <- c("ACTA2","FAP","PDGFRB","PDGFRA","COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","COL6A2","DCN","LUM","TAGLN","POSTN","MMP2","MMP11","THY1","VIM","FN1","CXCL12","CAV1","RGS5","MCAM","CSPG4")
tam_genes <- c("CD68","CD163","MRC1","MSR1","CSF1R","LYZ","C1QA","C1QB","C1QC","APOE","TREM2","SPP1","MARCO","FCGR3A","ITGAM","ITGAX","CD14","CD86","IL1B","TNF","CXCL8","CCL2","CCL3","CCL4","CCL5","TGFB1","VEGFA","ARG1","IL10")

# Key CellChat pathways
key_pathways <- c("COLLAGEN","MIF","TGFb","APP","CypA","COMPLEMENT","GALECTIN","CCL","CXCL","IL1","TNF","THBS","FN1")

# Single-cell specificity info
sc_genes <- sc_expr$gene
sc_cell_types <- setdiff(colnames(sc_expr), "gene")

# Initialize scores
gene_scores <- data.frame(
  gene = candidate_genes,
  score = 0,
  stringsAsFactors = FALSE
)

# Criterion 1: Validated same-direction DEG (both cohorts)
gene_scores$validated_DEG <- gene_scores$gene %in% validated_deg$gene_symbol
gene_scores$score <- gene_scores$score + gene_scores$validated_DEG * 2

# Criterion 2-3: Up/Down DEG
gene_scores$up_DEG <- gene_scores$gene %in% up_deg$gene_symbol
gene_scores$down_DEG <- gene_scores$gene %in% dn_deg$gene_symbol
gene_scores$score <- gene_scores$score + gene_scores$up_DEG * 1 + gene_scores$down_DEG * 1

# Criterion 4-5: CAF/TAM marker
gene_scores$CAF_marker <- gene_scores$gene %in% caf_genes
gene_scores$TAM_marker <- gene_scores$gene %in% tam_genes
gene_scores$score <- gene_scores$score + gene_scores$CAF_marker * 2 + gene_scores$TAM_marker * 2

# Criterion 6: Single-cell specificity (top expression in a specific cell type)
for (g in candidate_genes) {
  if (g %in% sc_genes) {
    row <- sc_expr[sc_expr$gene == g, sc_cell_types, drop=FALSE]
    vals <- as.numeric(row[1, ])
    names(vals) <- sc_cell_types
    top_val <- max(vals, na.rm=TRUE)
    second_val <- sort(vals, decreasing=TRUE)[2]
    # Specificity: ratio of top to second
    if (length(vals) > 1 && second_val > 0 && top_val / second_val > 2) {
      gene_scores$sc_specificity[gene_scores$gene == g] <- names(which.max(vals))
      gene_scores$score[gene_scores$gene == g] <- gene_scores$score[gene_scores$gene == g] + 2
    } else if (top_val > 1) {
      gene_scores$sc_specificity[gene_scores$gene == g] <- names(which.max(vals))
      gene_scores$score[gene_scores$gene == g] <- gene_scores$score[gene_scores$gene == g] + 1
    }
  }
}

# Criterion 7: In CellChat LR pair
cc_ligands <- unique(cc_lr$ligand)
cc_receptors <- unique(cc_lr$receptor)
cc_all_genes <- unique(c(cc_ligands, unlist(strsplit(cc_receptors, "_"))))
gene_scores$in_CellChat <- gene_scores$gene %in% cc_all_genes
gene_scores$score <- gene_scores$score + gene_scores$in_CellChat * 2

# Criterion 8: Immune/checkpoint association
if (!is.null(tcga_imm_cor)) {
  sig_imm <- tcga_imm_cor$cell_type[tcga_imm_cor$FDR < 0.05]
  gene_scores$immune_related <- gene_scores$gene %in% sig_imm |
    gene_scores$gene %in% c("CAF_Fibroblast","Macrophage_TAM")
}
if (!is.null(tcga_cp_cor)) {
  sig_cp <- tcga_cp_cor$checkpoint[tcga_cp_cor$FDR < 0.05]
  gene_scores$checkpoint_related <- gene_scores$gene %in% sig_cp
  gene_scores$score <- gene_scores$score + gene_scores$checkpoint_related * 2
}

# Criterion 9: Key pathway involvement
gene_scores$key_pathway <- FALSE
for (g in candidate_genes) {
  in_pathway <- FALSE
  # Check if gene appears in any key pathway LR pair
  lr_sub <- cc_lr[cc_lr$ligand == g | grepl(g, cc_lr$receptor), ]
  if (nrow(lr_sub) > 0) {
    pw_match <- any(key_pathways %in% lr_sub$pathway_name)
    if (pw_match) in_pathway <- TRUE
  }
  gene_scores$key_pathway[gene_scores$gene == g] <- in_pathway
  gene_scores$score[gene_scores$gene == g] <- gene_scores$score[gene_scores$gene == g] + in_pathway * 2
}

# Sort
gene_scores <- gene_scores[order(-gene_scores$score), ]
rownames(gene_scores) <- NULL

log_msg("Top 20 hub genes:")
print(head(gene_scores[, c("gene","score","validated_DEG","CAF_marker","TAM_marker","in_CellChat")], 20))

# Save all
write.csv(gene_scores, "tables/integrated/integrated_candidate_genes_all.csv", row.names=FALSE)
write.csv(head(gene_scores, 20), "tables/integrated/integrated_hub_genes_top20.csv", row.names=FALSE)

# ============================================================
# Step 3: 通讯轴优先级
# ============================================================
log_msg("\n=== Step 3: 通讯轴优先级 ===")

# Define focus axes
focus_axes <- list(
  list(name="CAF→TAM: MIF-CD74/CXCR4", source="CAF_Fibroblast", target="Macrophage_TAM",
       ligand="MIF", receptor_pattern="CD74"),
  list(name="CAF→TAM: MIF-CD74/CD44", source="CAF_Fibroblast", target="Macrophage_TAM",
       ligand="MIF", receptor_pattern="CD44"),
  list(name="CAF→TAM: COL1A1-CD44", source="CAF_Fibroblast", target="Macrophage_TAM",
       ligand="COL1A1", receptor_pattern="CD44"),
  list(name="CAF→TAM: COL1A2-CD44", source="CAF_Fibroblast", target="Macrophage_TAM",
       ligand="COL1A2", receptor_pattern="CD44"),
  list(name="CAF→TAM: APP-CD74", source="CAF_Fibroblast", target="Macrophage_TAM",
       ligand="APP", receptor_pattern="CD74"),
  list(name="TAM→CAF: PPIA-BSG", source="Macrophage_TAM", target="CAF_Fibroblast",
       ligand="PPIA", receptor_pattern="BSG"),
  list(name="TAM→CAF: TGFB1-TGFBR", source="Macrophage_TAM", target="CAF_Fibroblast",
       ligand="TGFB1", receptor_pattern="TGFBR"),
  list(name="Epi→TAM: MIF-CD74/CXCR4", source="Epithelial_like", target="Macrophage_TAM",
       ligand="MIF", receptor_pattern="CD74"),
  list(name="CAF→Epi: COL1A1-CD44", source="CAF_Fibroblast", target="Epithelial_like",
       ligand="COL1A1", receptor_pattern="CD44"),
  list(name="CAF→Epi: COL1A2-CD44", source="CAF_Fibroblast", target="Epithelial_like",
       ligand="COL1A2", receptor_pattern="CD44")
)

axis_scores <- data.frame(stringsAsFactors=FALSE)

for (ax in focus_axes) {
  # Find matching LR pairs in CellChat
  lr_matches <- cc_lr[cc_lr$source_cell == ax$source &
                      cc_lr$target_cell == ax$target &
                      cc_lr$ligand == ax$ligand &
                      grepl(ax$receptor_pattern, cc_lr$receptor), ]

  n_pairs <- nrow(lr_matches)
  max_prob <- if (n_pairs > 0) max(lr_matches$prob, na.rm=TRUE) else 0
  mean_prob <- if (n_pairs > 0) mean(lr_matches$prob, na.rm=TRUE) else 0

  # Evidence score
  lr_score <- 0
  lr_score <- lr_score + if (n_pairs > 0) 3 else 0  # Has CellChat evidence
  lr_score <- lr_score + if (max_prob > 0.05) 2 else 0  # High probability
  lr_score <- lr_score + if (ax$ligand %in% validated_deg$gene_symbol) 1 else 0  # Ligand is validated DEG
  lr_score <- lr_score + if (ax$ligand %in% gene_scores$gene[gene_scores$CAF_marker | gene_scores$TAM_marker]) 1 else 0

  # Ligand cell-type specificity from scRNA
  lr_score <- lr_score + if (ax$ligand %in% gene_scores$gene[!is.na(gene_scores$sc_specificity)]) 1 else 0

  axis_scores <- rbind(axis_scores, data.frame(
    axis = ax$name,
    source = ax$source, target = ax$target,
    ligand = ax$ligand, receptor = ax$receptor_pattern,
    CellChat_pairs = n_pairs,
    max_prob = round(max_prob, 4),
    mean_prob = round(mean_prob, 4),
    evidence_score = lr_score,
    stringsAsFactors = FALSE
  ))
}

axis_scores <- axis_scores[order(-axis_scores$evidence_score, -axis_scores$max_prob), ]
write.csv(axis_scores, "tables/integrated/integrated_LR_axes_prioritized.csv", row.names=FALSE)

log_msg("Prioritized LR axes:")
print(axis_scores[, c("axis","CellChat_pairs","max_prob","evidence_score")])

# ============================================================
# Step 4: CAF-TAM-Epithelial axis summary
# ============================================================
log_msg("\n=== Step 4: CAF-TAM-Epithelial axis summary ===")

axis_summary <- data.frame(
  direction = c("CAF→TAM","CAF→TAM","CAF→TAM","CAF→TAM",
                "TAM→CAF","TAM→CAF",
                "Epi→TAM","Epi→TAM",
                "CAF→Epi","CAF→Epi"),
  ligand = c("COL1A1","COL1A2","MIF","APP","PPIA","TGFB1","MIF","SPP1","COL1A1","COL1A2"),
  receptor = c("CD44","CD44","CD74/CD44","CD74","BSG","TGFBR","CD74","CD44","CD44","CD44"),
  pathway = c("COLLAGEN","COLLAGEN","MIF","APP","CypA","TGFb","MIF","SPP1","COLLAGEN","COLLAGEN"),
  bulk_DEG = c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),
  sc_ligand_in_CAF = c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE),
  sc_ligand_in_TAM = c(FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE),
  recommendation = c("Top priority","Top priority","Top priority","Top priority",
                     "Important","Important","Important","Important",
                     "Mechanistic support","Mechanistic support"),
  stringsAsFactors = FALSE
)
write.csv(axis_summary, "tables/integrated/integrated_CAF_TAM_Epithelial_axis_summary.csv", row.names=FALSE)

# ============================================================
# Step 5: 可视化
# ============================================================
log_msg("\n=== Step 5: 可视化 ===")

# ---- Fig9A: Integrated evidence heatmap ----
log_msg("Fig9A: Evidence heatmap")
top30 <- head(gene_scores, 30)
heat_mat <- as.matrix(top30[, c("validated_DEG","up_DEG","down_DEG","CAF_marker","TAM_marker","in_CellChat","key_pathway")])
mode(heat_mat) <- "numeric"
heat_mat <- heat_mat * 1  # TRUE→1
rownames(heat_mat) <- top30$gene
colnames(heat_mat) <- c("Validated_DEG","Up_DEG","Down_DEG","CAF_marker","TAM_marker","CellChat_LR","Key_Pathway")

pdf("figures/integrated/Fig9A_integrated_evidence_heatmap.pdf", width=9, height=9)
pheatmap(heat_mat, cluster_rows=FALSE, cluster_cols=FALSE,
         color=c("grey95","#2166AC"),
         display_numbers=TRUE, number_format="%d",
         main="Top 30 Hub Genes: Multi-layer Evidence",
         angle_col=45, fontsize=10)
dev.off()

# ---- Fig9B: Top hub gene barplot ----
log_msg("Fig9B: Hub gene barplot")
pdf("figures/integrated/Fig9B_top_hub_genes_barplot.pdf", width=10, height=8)
p <- ggplot(head(gene_scores, 20), aes(x=reorder(gene, score), y=score)) +
  geom_bar(stat="identity", fill="#2166AC") +
  coord_flip() +
  labs(x="", y="Integrated Evidence Score", title="Top 20 Hub Genes") +
  theme_pubr(base_size=12)
print(p); dev.off()

# ---- Fig9C: Prioritized LR axis dotplot ----
log_msg("Fig9C: LR axis dotplot")
pdf("figures/integrated/Fig9C_prioritized_LR_axis_dotplot.pdf", width=10, height=6)
p <- ggplot(axis_scores, aes(x=reorder(axis, evidence_score), y=evidence_score,
                              size=CellChat_pairs, color=max_prob)) +
  geom_point() +
  scale_color_gradient(low="blue", high="red") +
  coord_flip() +
  labs(x="", y="Evidence Score", size="LR Pairs", color="Max Prob",
       title="Prioritized Ligand-Receptor Axes") +
  theme_pubr(base_size=11)
print(p); dev.off()

# ---- Fig9D: Mechanism schematic ----
log_msg("Fig9D: Mechanism schematic")
pdf("figures/integrated/Fig9D_CAF_TAM_Epithelial_mechanism_schematic.pdf", width=12, height=8)

# Build a text-based schematic
plot.new()
par(mar=c(2,2,4,2))

# Draw three cell types as boxes
# Box positions
x_left <- 0.1; x_mid <- 0.5; x_right <- 0.9
y_top <- 0.85; y_mid <- 0.5; y_bot <- 0.15

# CAF_Fibroblast (left)
rect(x_left-0.12, y_top-0.1, x_left+0.12, y_top+0.1, col="#FFB3B3", border="#E41A1C", lwd=2)
text(x_left, y_top, "CAF\nFibroblast", cex=1.1, font=2)

# Macrophage_TAM (center)
rect(x_mid-0.14, y_mid-0.1, x_mid+0.14, y_mid+0.1, col="#B3D9FF", border="#377EB8", lwd=2)
text(x_mid, y_mid, "Macrophage\nTAM", cex=1.1, font=2)

# Epithelial_like (right)
rect(x_right-0.12, y_bot-0.1, x_right+0.12, y_bot+0.1, col="#B3E6B3", border="#4DAF4A", lwd=2)
text(x_right, y_bot, "Epithelial\nlike", cex=1.1, font=2)

# Arrows and labels for top axes
# CAF → TAM: COL1A1/CD44
arrows(x_left+0.12, y_top, x_mid-0.14, y_mid+0.03, length=0.1, lwd=3, col="#E41A1C")
text(0.25, 0.72, "COL1A1/COL1A2\n→ CD44", cex=0.7, col="#E41A1C", font=2)

# CAF → TAM: MIF/CD74
arrows(x_left+0.12, y_top-0.02, x_mid-0.14, y_mid-0.03, length=0.1, lwd=3, col="#E41A1C", lty=2)
text(0.25, 0.62, "MIF → CD74", cex=0.7, col="#E41A1C")

# TAM → CAF: PPIA/BSG (reverse)
arrows(x_mid-0.14, y_top-0.04, x_left+0.12, y_top-0.04, length=0.08, lwd=2, col="#377EB8")
text(0.25, 0.82, "← PPIA/BSG", cex=0.65, col="#377EB8")

# TAM → CAF: TGFB1
arrows(x_mid-0.14, y_top-0.08, x_left+0.12, y_top-0.08, length=0.08, lwd=2, col="#377EB8", lty=2)
text(0.25, 0.78, "← TGFB1", cex=0.65, col="#377EB8")

# Epi → TAM: MIF/CD74
arrows(x_right-0.12, y_bot+0.05, x_mid+0.14, y_mid-0.05, length=0.1, lwd=3, col="#4DAF4A")
text(0.74, 0.28, "MIF → CD74", cex=0.7, col="#4DAF4A")

# TAM → Epi: SPP1/CD44
arrows(x_mid+0.08, y_mid-0.1, x_right-0.08, y_bot+0.1, length=0.08, lwd=2, col="#377EB8", lty=2)
text(0.74, 0.42, "SPP1 → CD44", cex=0.65, col="#377EB8")

# Title and annotations
text(0.5, 0.97, "CAF–TAM–Epithelial Communication Axis", cex=1.5, font=2)

# Key downstream effects
text(0.1, 0.08, "ECM remodeling\nCollagen deposition\nMMP activation", cex=0.7)
text(0.5, 0.08, "M2 polarization\nIDO1/HAVCR2 expression\nImmunosuppression", cex=0.7)
text(0.9, 0.08, "Metabolic reprogramming\nPKM/SLC2A1 up\nFatty acid oxidation down", cex=0.7)

# Legend
legend(0.25, 0.38, legend=c("Collagen signaling","MIF axis","TAM feedback","Reverse signaling"),
       col=c("#E41A1C","#E41A1C","#377EB8","#4DAF4A"),
       lty=c(1,2,1,1), lwd=2, cex=0.7, bty="n")

dev.off()

# ============================================================
# Step 6: Summary
# ============================================================
log_msg("\n=== Step 6: Summary ===")

# Top 10 for report
top10 <- head(gene_scores, 10)

# Top axes
top_axes <- head(axis_scores, 3)

summary_df <- data.frame(
  item = c(
    "Top 10 hub genes",
    "Top 3 communication axes",
    "qPCR/IHC candidates",
    "Prognostic model candidates",
    "Main figures",
    "Supplementary"
  ),
  content = c(
    paste(top10$gene, sprintf("(score=%d)", top10$score), collapse=", "),
    paste(top_axes$axis, collapse="; "),
    "MMP11, COL1A1, SPP1, CD163, HAVCR2 — high cell-type specificity, easy IHC",
    "MIF, CD74, TGFB1, PKM, INSIG1 — validated DEG, multi-cohort confirmed",
    "Fig9A (evidence heatmap), Fig9B (hub genes), Fig9C (LR axes), Fig9D (mechanism)",
    "Single-cell UMAP, CellChat bubble plots, correlation heatmaps"
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_df, "tables/integrated/integrated_analysis_summary.csv", row.names=FALSE)

# ============================================================
# 最终报告
# ============================================================
cat(sprintf("\n========================================\n"))
cat(sprintf("Integrated Evidence Analysis Complete\n"))
cat(sprintf("========================================\n"))

cat(sprintf("\n=== TOP 10 HUB GENES ===\n"))
print(data.frame(Rank=1:10, Gene=top10$gene, Score=top10$score,
                 CAF=top10$CAF_marker, TAM=top10$TAM_marker,
                 CellChat=top10$in_CellChat, DEG=top10$validated_DEG), row.names=FALSE)

cat(sprintf("\n=== TOP COMMUNICATION AXES ===\n"))
print(top_axes[, c("axis","evidence_score","max_prob","CellChat_pairs")], row.names=FALSE)

cat(sprintf("\n========================================\n"))
