# ============================================================
# 分析目的: 对 TCGA-CHOL × GSE107943 共同验证的 IM/CAF/TAM DEGs
#          做功能富集分析 (GO BP / KEGG)，分上调和下调
# 输入文件:
#   tables/CCA_IM_CAF_TAM_DEGs_validated_same_direction.csv
#   tables/CCA_IM_CAF_TAM_combined_gene_annotation.csv
# 输出文件:
#   tables/enrichment/validated_IM_CAF_TAM_DEGs_up.csv
#   tables/enrichment/validated_IM_CAF_TAM_DEGs_down.csv
#   tables/enrichment/GO_BP_{up,down}_enrichment.csv
#   tables/enrichment/KEGG_{up,down}_enrichment.csv
#   tables/enrichment/enrichment_summary.csv
#   figures/enrichment/Fig3A-F (6 张图)
# 主要方法: clusterProfiler enrichGO / enrichKEGG, org.Hs.eg.db
# 注: ReactomePA 未安装，Reactome 富集跳过
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(ggpubr)
library(data.table)

dir.create("tables/enrichment", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/enrichment", showWarnings = FALSE, recursive = TRUE)

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ---- 配色 ----
up_col   <- "#E41A1C"
down_col <- "#377EB8"

# ============================================================
# Step 1: 读取数据并分组
# ============================================================
log_msg("=== Step 1: 读取数据并分组 ===")

validated <- read.csv("tables/CCA_IM_CAF_TAM_DEGs_validated_same_direction.csv",
                      stringsAsFactors = FALSE)
log_msg(sprintf("Validated IM/CAF/TAM DEGs: %d", nrow(validated)))

# 分组
up_genes <- validated[validated$direction_TCGA == "up", , drop = FALSE]
down_genes <- validated[validated$direction_TCGA == "down", , drop = FALSE]

log_msg(sprintf("  上调: %d", nrow(up_genes)))
log_msg(sprintf("  下调: %d", nrow(down_genes)))

# 保存
write.csv(up_genes, "tables/enrichment/validated_IM_CAF_TAM_DEGs_up.csv", row.names = FALSE)
write.csv(down_genes, "tables/enrichment/validated_IM_CAF_TAM_DEGs_down.csv", row.names = FALSE)
log_msg("已保存 up/down 基因列表")

# 来源分布
log_msg(sprintf("来源: CAF=%d, macrophage_TAM=%d, immunometabolism=%d, 多来源=%d",
  sum(validated$source == "CAF"),
  sum(validated$source == "macrophage_TAM"),
  sum(validated$source == "immunometabolism"),
  sum(grepl(";", validated$source))))

# ============================================================
# Step 2: Entrez ID 映射
# ============================================================
log_msg("\n=== Step 2: Entrez ID 映射 ===")

map_to_entrez <- function(genes) {
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}

up_entrez   <- map_to_entrez(up_genes$gene_symbol)
down_entrez <- map_to_entrez(down_genes$gene_symbol)

log_msg(sprintf("上调基因 Entrez 映射: %d / %d", nrow(up_entrez), nrow(up_genes)))
log_msg(sprintf("下调基因 Entrez 映射: %d / %d", nrow(down_entrez), nrow(down_genes)))

# 背景: 所有 validated DEGs 中能映射到 Entrez 的基因
all_entrez <- map_to_entrez(validated$gene_symbol)
log_msg(sprintf("背景基因 (universe): %d", nrow(all_entrez)))

# ============================================================
# Step 3: GO BP 富集
# ============================================================
log_msg("\n=== Step 3: GO BP 富集 ===")

run_go <- function(entrez_ids, label) {
  ego <- enrichGO(
    gene          = entrez_ids$ENTREZID,
    universe      = all_entrez$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  if (!is.null(ego) && nrow(ego) > 0) {
    log_msg(sprintf("GO BP %s: %d 显著条目", label, nrow(ego)))
  } else {
    log_msg(sprintf("GO BP %s: 无显著条目", label))
  }
  return(ego)
}

go_up   <- run_go(up_entrez, "up")
go_down <- run_go(down_entrez, "down")

if (!is.null(go_up) && nrow(go_up) > 0) {
  write.csv(as.data.frame(go_up), "tables/enrichment/GO_BP_up_enrichment.csv", row.names = FALSE)
  log_msg("已保存 GO_BP_up")
}
if (!is.null(go_down) && nrow(go_down) > 0) {
  write.csv(as.data.frame(go_down), "tables/enrichment/GO_BP_down_enrichment.csv", row.names = FALSE)
  log_msg("已保存 GO_BP_down")
}

# ============================================================
# Step 4: KEGG 富集
# ============================================================
log_msg("\n=== Step 4: KEGG 富集 ===")

run_kegg <- function(entrez_ids, label) {
  ekegg <- enrichKEGG(
    gene          = entrez_ids$ENTREZID,
    universe      = all_entrez$ENTREZID,
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  if (!is.null(ekegg) && nrow(ekegg) > 0) {
    ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    log_msg(sprintf("KEGG %s: %d 显著条目", label, nrow(ekegg)))
  } else {
    log_msg(sprintf("KEGG %s: 无显著条目", label))
  }
  return(ekegg)
}

kegg_up   <- run_kegg(up_entrez, "up")
kegg_down <- run_kegg(down_entrez, "down")

if (!is.null(kegg_up) && nrow(kegg_up) > 0) {
  write.csv(as.data.frame(kegg_up), "tables/enrichment/KEGG_up_enrichment.csv", row.names = FALSE)
  log_msg("已保存 KEGG_up")
}
if (!is.null(kegg_down) && nrow(kegg_down) > 0) {
  write.csv(as.data.frame(kegg_down), "tables/enrichment/KEGG_down_enrichment.csv", row.names = FALSE)
  log_msg("已保存 KEGG_down")
}

# ============================================================
# Step 5: 重点关注通路标记
# ============================================================
log_msg("\n=== Step 5: 关注通路标记 ===")

focus_pathways <- c(
  "extracellular matrix", "collagen", "focal adhesion",
  "PI3K-Akt", "cytokine", "NF-kappa", "glycolysis", "gluconeogenesis",
  "HIF-1", "HIF1",
  "fatty acid", "bile acid", "PPAR", "xenobiotic", "cytochrome P450",
  "oxidative phosphorylation", "glutathione", "cholesterol",
  "complement and coagulation", "ECM"
)

focus_terms <- function(obj) {
  # Return just the focus logical vector
  if (is.null(obj) || nrow(obj) == 0) return(logical(0))
  focus <- rep(FALSE, nrow(obj))
  for (kw in focus_pathways) {
    focus <- focus | grepl(kw, obj@result$Description, ignore.case = TRUE)
  }
  return(focus)
}

go_focus_up   <- focus_terms(go_up)
go_focus_down <- focus_terms(go_down)
kegg_focus_up   <- focus_terms(kegg_up)
kegg_focus_down <- focus_terms(kegg_down)

# ============================================================
# Step 6: 汇总表
# ============================================================
log_msg("\n=== Step 6: 生成汇总表 ===")

n_go_up   <- if (!is.null(go_up) && !is.null(nrow(go_up))) nrow(go_up) else 0
n_go_down <- if (!is.null(go_down) && !is.null(nrow(go_down))) nrow(go_down) else 0
n_kegg_up   <- if (!is.null(kegg_up) && !is.null(nrow(kegg_up))) nrow(kegg_up) else 0
n_kegg_down <- if (!is.null(kegg_down) && !is.null(nrow(kegg_down))) nrow(kegg_down) else 0

summary_rows <- data.frame(
  item = c(
    "Validated IM/CAF/TAM DEGs total", "上调基因数", "下调基因数",
    "GO BP up 显著条目", "GO BP down 显著条目",
    "KEGG up 显著条目", "KEGG down 显著条目",
    "ReactomePA 状态",
    "universe 基因数",
    "来源: immunometabolism", "来源: CAF", "来源: macrophage_TAM", "来源: multiple"
  ),
  value = c(
    nrow(validated), nrow(up_genes), nrow(down_genes),
    n_go_up, n_go_down, n_kegg_up, n_kegg_down,
    "未安装，已跳过",
    nrow(all_entrez),
    sum(validated$source == "immunometabolism"),
    sum(validated$source == "CAF"),
    sum(validated$source == "macrophage_TAM"),
    sum(grepl(";", validated$source))
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_rows, "tables/enrichment/enrichment_summary.csv", row.names = FALSE)
log_msg("已保存 enrichment_summary.csv")

# ============================================================
# Step 7: 可视化
# ============================================================
log_msg("\n=== Step 7: 可视化 ===")

# Helper: safe dotplot
safe_dotplot <- function(obj, title, n = 15) {
  if (!is.null(obj) && nrow(obj) > 0) {
    n_show <- min(n, nrow(obj))
    dotplot(obj, showCategory = n_show, title = title)
  } else {
    ggplot() + annotate("text", x=0, y=0, label="No significant enrichment") +
      ggtitle(title) + theme_void()
  }
}

# ---- Fig3A: GO BP up dotplot ----
log_msg("Fig3A: GO BP up dotplot")
pdf("figures/enrichment/Fig3A_GO_BP_up_dotplot.pdf", width=10, height=8)
print(safe_dotplot(go_up, "GO BP: Up-regulated in Tumor", 15))
dev.off()

# ---- Fig3B: GO BP down dotplot ----
log_msg("Fig3B: GO BP down dotplot")
pdf("figures/enrichment/Fig3B_GO_BP_down_dotplot.pdf", width=10, height=8)
print(safe_dotplot(go_down, "GO BP: Down-regulated in Tumor", 15))
dev.off()

# ---- Fig3C: KEGG up dotplot ----
log_msg("Fig3C: KEGG up dotplot")
pdf("figures/enrichment/Fig3C_KEGG_up_dotplot.pdf", width=10, height=7)
print(safe_dotplot(kegg_up, "KEGG: Up-regulated in Tumor", 15))
dev.off()

# ---- Fig3D: KEGG down dotplot ----
log_msg("Fig3D: KEGG down dotplot")
pdf("figures/enrichment/Fig3D_KEGG_down_dotplot.pdf", width=10, height=7)
print(safe_dotplot(kegg_down, "KEGG: Down-regulated in Tumor", 15))
dev.off()

# ---- Fig3E: Up/Down pathway barplot ----
log_msg("Fig3E: Up/Down pathway barplot")

build_bar_data <- function(go_obj, kegg_obj, direction) {
  rows <- list()
  if (!is.null(go_obj) && nrow(go_obj) > 0) {
    go_top <- go_obj[order(go_obj$p.adjust), ]
    if (nrow(go_top) > 8) go_top <- go_top[1:8, ]
    for (i in seq_len(nrow(go_top))) {
      rows[[length(rows)+1]] <- c(go_top$Description[i], "GO_BP", go_top$Count[i], -log10(go_top$p.adjust[i]), direction)
    }
  }
  if (!is.null(kegg_obj) && nrow(kegg_obj) > 0) {
    kegg_top <- kegg_obj[order(kegg_obj$p.adjust), ]
    if (nrow(kegg_top) > 8) kegg_top <- kegg_top[1:8, ]
    for (i in seq_len(nrow(kegg_top))) {
      rows[[length(rows)+1]] <- c(kegg_top$Description[i], "KEGG", kegg_top$Count[i], -log10(kegg_top$p.adjust[i]), direction)
    }
  }
  if (length(rows) == 0) return(NULL)
  df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  colnames(df) <- c("Pathway", "Database", "Count", "negLog10P", "Direction")
  df$Count <- as.integer(df$Count)
  df$negLog10P <- as.numeric(df$negLog10P)
  df
}

bar_up   <- build_bar_data(go_up, kegg_up, "Up")
bar_down <- build_bar_data(go_down, kegg_down, "Down")

if (!is.null(bar_up) || !is.null(bar_down)) {
  bar_all <- rbind(if (!is.null(bar_up)) bar_up else NULL,
                   if (!is.null(bar_down)) bar_down else NULL)

  if (nrow(bar_all) > 0) {
    # Truncate long pathway names
    bar_all$Pathway <- sapply(bar_all$Pathway, function(x) {
      if (nchar(x) > 50) paste0(substr(x, 1, 47), "...") else x
    })
    bar_all$Pathway <- factor(bar_all$Pathway, levels = rev(unique(bar_all$Pathway)))

    pdf("figures/enrichment/Fig3E_up_down_pathway_barplot.pdf", width=12, height=8)
    p <- ggplot(bar_all, aes(x = Pathway, y = negLog10P, fill = Direction)) +
      geom_bar(stat = "identity", position = "dodge") +
      coord_flip() +
      scale_fill_manual(values = c("Up" = up_col, "Down" = down_col)) +
      labs(x = "", y = "-log10(adjusted p-value)",
           title = "Enriched Pathways: Up vs Down in Tumor") +
      facet_wrap(~ Database, scales = "free_y") +
      theme_pubr(base_size = 12)
    print(p); dev.off()
    log_msg("Fig3E Done")
  }
}

# ---- Fig3F: Source distribution ----
log_msg("Fig3F: Source distribution")

# Parse source
src_counts <- c(
  immunometabolism = sum(validated$in_immunometabolism == "TRUE" | validated$source == "immunometabolism"),
  CAF = sum(validated$source == "CAF" | validated$in_CAF == "TRUE"),
  macrophage_TAM = sum(validated$source == "macrophage_TAM" | validated$in_macrophage_TAM == "TRUE")
)

# Handle multiple source — use grepl
src_counts <- c(
  immunometabolism = sum(grepl("immunometabolism", validated$source)),
  CAF = sum(grepl("CAF", validated$source)),
  macrophage_TAM = sum(grepl("macrophage_TAM", validated$source))
)
# Multiple occurs when more than one source
src_counts["multiple_source"] <- sum(grepl(";", validated$source))

src_df <- data.frame(source = names(src_counts), count = as.integer(src_counts),
                     stringsAsFactors = FALSE)
src_df$label <- sprintf("%s\n(n=%d)", src_df$source, src_df$count)

pdf("figures/enrichment/Fig3F_validated_DEG_source_distribution.pdf", width=7, height=6)
p <- ggplot(src_df, aes(x = "", y = count, fill = source)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = sprintf("Validated DEGs Source Distribution (n=%d)", nrow(validated)),
       fill = "Source") +
  theme_void() + theme(legend.position = "right")
print(p); dev.off()
log_msg("Fig3F Done")

# ============================================================
# 最终报告
# ============================================================
cat(sprintf("\n========================================\n"))
cat(sprintf("功能富集分析完成\n"))
cat(sprintf("========================================\n"))
cat(sprintf("Validated DEGs total:    %d\n", nrow(validated)))
cat(sprintf("  上调: %d\n", nrow(up_genes)))
cat(sprintf("  下调: %d\n", nrow(down_genes)))
cat(sprintf("GO BP up:   %d\n", n_go_up))
cat(sprintf("GO BP down: %d\n", n_go_down))
cat(sprintf("KEGG up:    %d\n", n_kegg_up))
cat(sprintf("KEGG down:  %d\n", n_kegg_down))
cat(sprintf("ReactomePA: 未安装，已跳过\n"))
cat(sprintf("========================================\n"))

# Top terms
print_top <- function(obj, label, n = 10) {
  cat(sprintf("\n--- %s Top %d ---\n", label, n))
  if (!is.null(obj) && nrow(obj) > 0) {
    top <- head(obj[order(obj$p.adjust), ], n)
    print(data.frame(
      Description = top$Description,
      Count = top$Count,
      p.adjust = format(top$p.adjust, digits = 3, scientific = TRUE)
    ), row.names = FALSE)
  } else {
    cat("  No significant terms\n")
  }
}

print_top(go_up, "GO BP UP")
print_top(go_down, "GO BP DOWN")
print_top(kegg_up, "KEGG UP")
print_top(kegg_down, "KEGG DOWN")

cat(sprintf("\nFigures: figures/enrichment/Fig3A-F.pdf\n"))
cat(sprintf("========================================\n"))
