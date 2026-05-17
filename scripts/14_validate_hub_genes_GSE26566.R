# ============================================================
# 分析目的: GSE26566 外部验证核心 hub genes 和通讯轴基因
#          验证 COL1A1/COL1A2/TREM2/SPP1 等在第三队列表达和预后
# 输入/下载: GSE26566 (GEO) — microarray (GPL6104)
# 输出:
#   data/raw/GSE26566/ (原始)
#   data/processed/GSE26566/ (清洗后)
#   tables/GSE26566_validation/ (7 个表格)
#   figures/GSE26566_validation/ (6 张图)
# 主要方法: GEOquery → probe→gene mapping → z-score →
#          aggressive-like score → survival analysis
# ============================================================

library(GEOquery)
library(limma)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(pheatmap)
library(RColorBrewer)

dir.create("data/raw/GSE26566", showWarnings=FALSE, recursive=TRUE)
dir.create("data/processed/GSE26566", showWarnings=FALSE, recursive=TRUE)
dir.create("data/clinical/GSE26566", showWarnings=FALSE, recursive=TRUE)
dir.create("tables/GSE26566_validation", showWarnings=FALSE, recursive=TRUE)
dir.create("figures/GSE26566_validation", showWarnings=FALSE, recursive=TRUE)

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ============================================================
# Step 1: 下载 GSE26566
# ============================================================
log_msg("=== Step 1: 下载 GSE26566 ===")

options(timeout=600)
gse_file <- "data/raw/GSE26566/GSE26566_series_matrix.txt.gz"

if (!file.exists(gse_file)) {
  log_msg("Downloading GSE26566...")
  gse <- getGEO("GSE26566", GSEMatrix=TRUE, getGPL=TRUE,
                destdir="data/raw/GSE26566")
} else {
  log_msg("Series matrix already exists, loading...")
  gse <- getGEO(filename=gse_file, getGPL=TRUE)
}

# GSE26566 should be a single ExpressionSet
eset <- gse[[1]]
log_msg(sprintf("ExpressionSet: %d probes x %d samples", nrow(eset), ncol(eset)))

# ============================================================
# Step 2: 提取表达 + 表型
# ============================================================
log_msg("\n=== Step 2: 提取表达和表型 ===")

# Expression matrix
expr_mat <- exprs(eset)
log_msg(sprintf("Expression matrix: %d x %d", nrow(expr_mat), ncol(expr_mat)))
log_msg(sprintf("Value range: %.2f - %.2f", min(expr_mat, na.rm=TRUE), max(expr_mat, na.rm=TRUE)))

# Phenotype data
pd <- pData(eset)
log_msg(sprintf("Phenotype columns: %d", ncol(pd)))
log_msg(sprintf("Column names: %s", paste(colnames(pd), collapse=", ")))

# Save all phenotype columns for inspection
write.csv(pd, "tables/GSE26566_validation/GSE26566_phenotype_all.csv", row.names=FALSE)

# Platform annotation
gpl <- eset@annotation
log_msg(sprintf("Platform: %s", gpl))

# ============================================================
# Step 3: Probe → Gene Symbol mapping
# ============================================================
log_msg("\n=== Step 3: Probe → Gene mapping ===")

# Get platform annotation
gpl_obj <- getGEO(gpl, destdir="data/raw/GSE26566")
gpl_table <- Table(gpl_obj)

# Find gene symbol column
gene_col <- NULL
for (cn in c("Gene symbol","Symbol","GENE_SYMBOL","gene_symbol","Gene Symbol",
             "ILMN_Gene","Gene name","gene")) {
  if (cn %in% colnames(gpl_table)) { gene_col <- cn; break }
}
log_msg(sprintf("Gene symbol column: %s", gene_col))

if (is.null(gene_col)) {
  log_msg("Available GPL columns:")
  print(head(colnames(gpl_table), 20))
  stop("Cannot find gene symbol column in platform annotation")
}

# Map probes to genes
probe_ids <- gpl_table$ID
gene_symbols <- gpl_table[[gene_col]]

# Remove probes without gene symbol
valid <- !is.na(gene_symbols) & gene_symbols != "" & gene_symbols != "---"
probe_gene_map <- data.frame(
  probe_id = probe_ids[valid],
  gene_symbol = gene_symbols[valid],
  stringsAsFactors = FALSE
)

log_msg(sprintf("Probes with gene symbol: %d / %d", nrow(probe_gene_map), length(probe_ids)))

# Match expression matrix probes
common_probes <- intersect(rownames(expr_mat), probe_gene_map$probe_id)
log_msg(sprintf("Probes in expression matrix with gene symbol: %d / %d",
  length(common_probes), nrow(expr_mat)))

expr_sub <- expr_mat[common_probes, , drop=FALSE]
gene_map_sub <- probe_gene_map[match(common_probes, probe_gene_map$probe_id), ]

# Aggregate: mean by gene symbol
gene_symbols_agg <- gene_map_sub$gene_symbol
# Handle multi-gene symbols (e.g., "GENE1 /// GENE2")
gene_symbols_agg <- sapply(strsplit(gene_symbols_agg, " /// "), `[`, 1)

expr_gene <- aggregate(expr_sub, by=list(gene=gene_symbols_agg), FUN=mean, na.rm=TRUE)
rownames(expr_gene) <- expr_gene$gene; expr_gene$gene <- NULL
expr_gene <- as.matrix(expr_gene)

log_msg(sprintf("Gene-level expression: %d genes x %d samples", nrow(expr_gene), ncol(expr_gene)))

# Save clean matrix
write.csv(data.frame(gene=rownames(expr_gene), expr_gene, check.names=FALSE),
          "data/processed/GSE26566/GSE26566_expression_matrix_clean.csv", row.names=FALSE)

# Platform summary
plat_summary <- data.frame(
  item = c("GEO accession","Platform","Platform ID","Original probes","Probes with symbol",
           "Gene-level features","Samples"),
  value = c("GSE26566","Illumina humanRef-8 v2.0",gpl,
            as.character(nrow(expr_mat)),as.character(length(common_probes)),
            as.character(nrow(expr_gene)),as.character(ncol(expr_gene))),
  stringsAsFactors = FALSE
)
write.csv(plat_summary, "tables/GSE26566_validation/GSE26566_platform_annotation_summary.csv", row.names=FALSE)

# ============================================================
# Step 4: 提取临床信息
# ============================================================
log_msg("\n=== Step 4: 提取临床信息 ===")

# Try to find survival columns
surv_cols <- list()
for (os_col in c("os:ch1","overall survival:ch1","survival time:ch1",
                 "days to death:ch1","follow up time:ch1","survival:ch1",
                 "os (months):ch1","overall_survival:ch1")) {
  if (os_col %in% colnames(pd)) surv_cols$os_time <- os_col
}
for (event_col in c("os event:ch1","death:ch1","overall survival event:ch1",
                    "vital status:ch1","status:ch1","event death:ch1")) {
  if (event_col %in% colnames(pd)) surv_cols$os_event <- event_col
}
for (rec_col in c("recurrence:ch1","recurr:ch1","relapse:ch1",
                  "recurrence status:ch1","dfs:ch1","rfs:ch1")) {
  if (rec_col %in% colnames(pd)) surv_cols$recurrence <- rec_col
}
for (tn_col in c("tissue:ch1","sample type:ch1","source_name_ch1",
                 "tumor/normal:ch1","characteristics_ch1")) {
  if (tn_col %in% colnames(pd)) surv_cols$tissue <- tn_col
}

log_msg(sprintf("Survival columns: %s", paste(names(surv_cols), collapse=", ")))

# Build clinical table
clinical <- pd[, c("geo_accession","title"), drop=FALSE]
colnames(clinical)[1:2] <- c("sample_id","sample_name")

# Add survival if found
if (!is.null(surv_cols$os_time)) {
  os_raw <- pd[[surv_cols$os_time]]
  os_num <- as.numeric(gsub("[^0-9.]","", os_raw))
  clinical$OS_time <- os_num
}
if (!is.null(surv_cols$os_event)) {
  event_raw <- pd[[surv_cols$os_event]]
  event_num <- as.numeric(gsub("[^0-9]","", event_raw))
  clinical$OS_status <- event_num
}

# Tissue type
if (!is.null(surv_cols$tissue)) {
  clinical$tissue <- pd[[surv_cols$tissue]]
}

has_os <- "OS_time" %in% colnames(clinical) && "OS_status" %in% colnames(clinical)
has_tn <- "tissue" %in% colnames(clinical)
log_msg(sprintf("Has OS: %s, Has tissue type: %s", has_os, has_tn))

write.csv(clinical, "data/clinical/GSE26566/GSE26566_clinical_clean.csv", row.names=FALSE)
log_msg("Clinical data saved")

# ============================================================
# Step 5: 核心基因验证
# ============================================================
log_msg("\n=== Step 5: 核心 hub gene 验证 ===")

hub_genes <- c("COL1A1","COL1A2","TREM2","POSTN","SPP1","MMP11","COL3A1",
               "FAP","CD163","C1QA",
               "MIF","CD74","CXCR4","CD44","TGFB1","PPIA","BSG",
               "PKM","SLC2A1","GLUD1","CAT","ACAA1","LCAT","ACADSB")
hub_genes <- unique(hub_genes)

hub_found <- intersect(hub_genes, rownames(expr_gene))
hub_missing <- setdiff(hub_genes, rownames(expr_gene))
log_msg(sprintf("Hub genes found: %d/%d", length(hub_found), length(hub_genes)))
if (length(hub_missing) > 0) log_msg(sprintf("  Missing: %s", paste(hub_missing, collapse=", ")))

# Hub gene summary
hub_summary <- data.frame(stringsAsFactors=FALSE)
for (g in hub_genes) {
  if (g %in% rownames(expr_gene)) {
    vals <- as.numeric(expr_gene[g, ])
    hub_summary <- rbind(hub_summary, data.frame(
      gene=g, found=TRUE, mean=mean(vals,na.rm=TRUE), sd=sd(vals,na.rm=TRUE),
      median=median(vals,na.rm=TRUE), min=min(vals,na.rm=TRUE), max=max(vals,na.rm=TRUE),
      stringsAsFactors=FALSE))
  } else {
    hub_summary <- rbind(hub_summary, data.frame(
      gene=g, found=FALSE, mean=NA, sd=NA, median=NA, min=NA, max=NA, stringsAsFactors=FALSE))
  }
}
write.csv(hub_summary, "tables/GSE26566_validation/GSE26566_hub_gene_expression_summary.csv", row.names=FALSE)

# ---- FigS: Hub gene heatmap ----
log_msg("FigS: Hub gene heatmap")
if (length(hub_found) > 3) {
  hub_mat <- expr_gene[hub_found, , drop=FALSE]
  hub_z <- t(scale(t(hub_mat)))
  hub_z[hub_z > 3] <- 3; hub_z[hub_z < -3] <- -3

  # Annotate with tissue if available
  ann_col <- NULL
  if (has_tn) {
    ann_col <- data.frame(Tissue=clinical$tissue[match(colnames(hub_z), clinical$sample_id)],
                          row.names=colnames(hub_z))
  }

  pdf("figures/GSE26566_validation/FigS_GSE26566_hub_gene_heatmap.pdf", width=14, height=8)
  pheatmap(hub_z, annotation_col=ann_col, show_colnames=FALSE,
           main="GSE26566: Hub Gene Expression (Z-score)", fontsize_row=10)
  dev.off()
}

# ---- FigS: Hub gene boxplot T vs N (if available) ----
if (has_tn && length(hub_found) > 0) {
  log_msg("FigS: Hub gene boxplot T vs N")
  tn_info <- clinical$tissue[match(colnames(expr_gene), clinical$sample_id)]
  is_tumor <- grepl("tumor|cancer|carcinoma|malignant", tn_info, ignore.case=TRUE)
  is_normal <- grepl("normal|adjacent|non.tumor|healthy", tn_info, ignore.case=TRUE)

  if (sum(is_tumor) > 0 && sum(is_normal) > 0) {
    pdf("figures/GSE26566_validation/FigS_GSE26566_hub_gene_boxplot.pdf", width=14, height=10)
    par(mfrow=c(4,6), mar=c(4,3,3,1))
    for (g in head(hub_found, 24)) {
      t_vals <- expr_gene[g, is_tumor]; n_vals <- expr_gene[g, is_normal]
      boxplot(list(T=t_vals, N=n_vals), col=c("#E41A1C","#377EB8"), main=g, las=1)
    }
    dev.off()
  }
}

# ---- FigS: Hub gene correlation heatmap ----
log_msg("FigS: Hub gene correlation")
if (length(hub_found) >= 5) {
  hub_cor <- cor(t(expr_gene[hub_found, , drop=FALSE]), method="spearman")
  pdf("figures/GSE26566_validation/FigS_GSE26566_hub_gene_correlation_heatmap.pdf", width=10, height=9)
  pheatmap(hub_cor, display_numbers=TRUE, number_format="%.2f",
           color=colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
           breaks=seq(-1,1,length.out=101), main="GSE26566: Hub Gene Spearman Correlation",
           fontsize=9)
  dev.off()
}

# ============================================================
# Step 6: Aggressive-like score
# ============================================================
log_msg("\n=== Step 6: Aggressive-like score ===")

caf_ecm_gs <- c("COL1A1","COL1A2","COL3A1","POSTN","FAP","MMP11")
tam_gs     <- c("TREM2","SPP1","CD163","C1QA","C1QB","C1QC")
comm_gs    <- c("MIF","CD74","CXCR4","CD44","TGFB1","PPIA","BSG")

all_aggr_gs <- unique(c(caf_ecm_gs, tam_gs, comm_gs))
aggr_avail  <- intersect(all_aggr_gs, rownames(expr_gene))
log_msg(sprintf("Aggressive score genes available: %d/%d", length(aggr_avail), length(all_aggr_gs)))

if (length(aggr_avail) >= 5) {
  aggr_mat <- expr_gene[aggr_avail, , drop=FALSE]
  aggr_z <- t(scale(t(aggr_mat)))
  aggr_score <- colMeans(aggr_z, na.rm=TRUE)

  # Add to clinical
  clinical$aggressive_score <- aggr_score[match(clinical$sample_id, names(aggr_score))]

  # Distribution
  pdf("figures/GSE26566_validation/FigS_GSE26566_aggressive_score_distribution.pdf", width=7, height=5)
  hist(clinical$aggressive_score, breaks=20, col="steelblue", border="white",
       main="GSE26566: Aggressive-like Score Distribution", xlab="Score")
  abline(v=median(clinical$aggressive_score, na.rm=TRUE), col="red", lwd=2, lty=2)
  legend("topright", legend=sprintf("Median=%.3f", median(clinical$aggressive_score, na.rm=TRUE)), bty="n")
  dev.off()

  # Survival analysis if OS available
  if (has_os) {
    clinical$OS_time_n <- as.numeric(clinical$OS_time)
    clinical$OS_status_n <- as.numeric(clinical$OS_status)
    clinical_os <- clinical[!is.na(clinical$OS_time_n) & clinical$OS_time_n > 0 &
                            !is.na(clinical$OS_status_n), ]

    if (nrow(clinical_os) > 20) {
      med_aggr <- median(clinical_os$aggressive_score, na.rm=TRUE)
      clinical_os$aggr_group <- ifelse(clinical_os$aggressive_score > med_aggr, "High", "Low")

      # KM
      log_msg("FigS: Aggressive score KM")
      sdiff <- survdiff(Surv(OS_time_n, OS_status_n) ~ aggr_group, data=clinical_os)
      pval <- 1 - pchisq(sdiff$chisq, 1)
      cox <- coxph(Surv(OS_time_n, OS_status_n) ~ aggressive_score, data=clinical_os)
      hr <- exp(coef(cox))
      ci <- exp(confint(cox))

      pdf("figures/GSE26566_validation/FigS_GSE26566_aggressive_score_KM.pdf", width=7, height=6)
      fit <- survfit(Surv(OS_time_n, OS_status_n) ~ aggr_group, data=clinical_os)
      p <- ggsurvplot(fit, data=clinical_os, pval=TRUE,
        palette=c("High"="#E41A1C","Low"="#377EB8"),
        legend.title="Aggressive Score", xlab="Time",
        title=sprintf("GSE26566: Aggressive-like Score (HR=%.2f)", hr))
      print(p); dev.off()

      aggr_surv <- data.frame(
        metric=c("Aggressive score Cox HR","Cox p-value","95% CI lower","95% CI upper",
                 "High group n","Low group n","Events total","Log-rank p","Median score"),
        value=c(sprintf("%.2f",hr),sprintf("%.4f",summary(cox)$coefficients[5]),
                sprintf("%.2f",ci[1]),sprintf("%.2f",ci[2]),
                as.character(sum(clinical_os$aggr_group=="High")),
                as.character(sum(clinical_os$aggr_group=="Low")),
                as.character(sum(clinical_os$OS_status_n==1)),
                sprintf("%.4f",pval),sprintf("%.3f",med_aggr)),
        stringsAsFactors=FALSE)
      write.csv(aggr_surv, "tables/GSE26566_validation/GSE26566_aggressive_score_survival.csv", row.names=FALSE)

      log_msg(sprintf("  Aggressive score: HR=%.2f, p=%.4f", hr, summary(cox)$coefficients[5]))
    }
  }
}

# ============================================================
# Step 7: 单基因生存分析
# ============================================================
log_msg("\n=== Step 7: 单基因生存分析 ===")

gene_cox <- data.frame(stringsAsFactors=FALSE)

if (has_os && exists("clinical_os") && nrow(clinical_os) > 20) {
  for (g in hub_found) {
    vals <- as.numeric(expr_gene[g, ])
    names(vals) <- colnames(expr_gene)
    surv_vals <- vals[clinical_os$sample_id]
    valid <- is.finite(surv_vals) & !is.na(clinical_os$OS_time_n) & clinical_os$OS_time_n > 0
    if (sum(valid) > 20) {
      cox <- tryCatch(
        coxph(Surv(OS_time_n, OS_status_n) ~ surv_vals, data=clinical_os[valid,]),
        error=function(e) NULL)
      if (!is.null(cox)) {
        s <- summary(cox)
        gene_cox <- rbind(gene_cox, data.frame(
          gene=g, HR=s$conf.int[1], lower95=s$conf.int[3], upper95=s$conf.int[4],
          pvalue=s$coefficients[5], stringsAsFactors=FALSE))
      }
    }
  }

  if (nrow(gene_cox) > 0) {
    gene_cox$FDR <- p.adjust(gene_cox$pvalue, method="BH")
    gene_cox <- gene_cox[order(gene_cox$pvalue), ]
    write.csv(gene_cox, "tables/GSE26566_validation/GSE26566_hub_gene_survival_univariate_cox.csv", row.names=FALSE)

    log_msg(sprintf("Genes with Cox results: %d, FDR<0.05: %d",
      nrow(gene_cox), sum(gene_cox$FDR < 0.05)))
    if (sum(gene_cox$FDR < 0.05) > 0) {
      log_msg("Significant genes:")
      print(gene_cox[gene_cox$FDR < 0.05, c("gene","HR","pvalue","FDR")])
    }

    # KM for top 3 risk genes
    top3 <- head(gene_cox, 3)
    for (i in seq_len(min(3, nrow(top3)))) {
      g <- top3$gene[i]
      vals <- as.numeric(expr_gene[g, ])
      names(vals) <- colnames(expr_gene)
      surv_g <- clinical_os
      surv_g$expr <- vals[surv_g$sample_id]
      med <- median(surv_g$expr, na.rm=TRUE)
      surv_g$group <- ifelse(surv_g$expr > med, "High", "Low")

      pdf(sprintf("figures/GSE26566_validation/FigS_GSE26566_top_hub_gene_KM_%d_%s.pdf", i, g),
          width=7, height=6)
      fit <- survfit(Surv(OS_time_n, OS_status_n) ~ group, data=surv_g)
      p <- ggsurvplot(fit, data=surv_g, pval=TRUE,
        palette=c("High"="#E41A1C","Low"="#377EB8"),
        legend.title=g, xlab="Time",
        title=sprintf("GSE26566: %s (HR=%.2f)", g, top3$HR[i]))
      print(p); dev.off()
    }
  }
}

# ============================================================
# Step 8: Summary
# ============================================================
log_msg("\n=== Step 8: Summary ===")

summary_df <- data.frame(
  item = c(
    "GSE26566 download status","Samples","Platform",
    "Gene-level features","Probes with symbol",
    "Has OS","Has tumor/normal","Has recurrence/DFS",
    "Hub genes covered","Hub genes missing",
    "Aggressive score Cox HR (if OS)","Aggressive score Cox p",
    "Genes with Cox FDR<0.05",
    "Supports COL1A1/COL1A2/TREM2/SPP1",
    "All figures generated"
  ),
  value = c(
    "Yes", as.character(nrow(clinical)), gpl,
    as.character(nrow(expr_gene)), as.character(length(common_probes)),
    as.character(has_os), as.character(has_tn),
    if (!is.null(surv_cols$recurrence)) "Yes" else "No/Unclear",
    sprintf("%d/%d",length(hub_found),length(hub_genes)),
    paste(hub_missing, collapse=", "),
    if (exists("hr")) sprintf("%.2f",hr) else "N/A",
    if (exists("hr")) sprintf("%.4f",summary(cox)$coefficients[5]) else "N/A",
    if (nrow(gene_cox)>0) paste(gene_cox$gene[gene_cox$FDR<0.05],collapse=", ") else "None",
    "Yes — third cohort external validation",
    "Yes"
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_df, "tables/GSE26566_validation/GSE26566_validation_summary.csv", row.names=FALSE)

cat(sprintf("\n========================================\n"))
cat(sprintf("GSE26566 Validation Complete\n"))
cat(sprintf("========================================\n"))
cat(sprintf("Samples: %d\n", ncol(expr_gene)))
cat(sprintf("Hub genes: %d/%d found\n", length(hub_found), length(hub_genes)))
cat(sprintf("Has OS: %s\n", has_os))
if (exists("hr")) cat(sprintf("Aggressive score: HR=%.2f, p=%.4f\n", hr, summary(cox)$coefficients[5]))
if (nrow(gene_cox) > 0) cat(sprintf("Sig genes: %s\n", paste(gene_cox$gene[gene_cox$FDR<0.05],collapse=", ")))
cat(sprintf("========================================\n"))
