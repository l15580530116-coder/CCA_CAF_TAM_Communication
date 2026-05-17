# ============================================================
# 分析目的: 诊断 GSE138709 gene identifier 和矩阵方向问题
# 输入:
#   data/raw/GSE138709/ (8 CSV.gz files)
#   data/processed/GSE138709/GSE138709_seurat_processed.rds
# 输出:
#   tables/single_cell_diagnosis/ (6 个诊断表格)
#   figures/single_cell_diagnosis/ (诊断图)
# 主要方法: 逐文件检查 → marker 搜索 → Seurat feature 扫描
# ============================================================

library(Seurat)
library(data.table)

dir.create("tables/single_cell_diagnosis", showWarnings=FALSE, recursive=TRUE)
dir.create("figures/single_cell_diagnosis", showWarnings=FALSE, recursive=TRUE)

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# Marker genes to search
epithelial  <- c("EPCAM","KRT7","KRT8","KRT18","KRT19","SOX9","MUC1")
caf_markers <- c("COL1A1","COL1A2","COL3A1","DCN","LUM","ACTA2","FAP","PDGFRB","POSTN","MMP11")
macro_genes <- c("CD68","CD163","MRC1","LYZ","C1QA","C1QB","C1QC","APOE","TREM2","SPP1")
tcell_genes <- c("CD3D","CD3E","CD8A","CD8B","FOXP3","CTLA4")
checkpoint  <- c("HAVCR2","IDO1","CD86","PDCD1LG2","TIGIT","LAG3","PDCD1")
metabolism  <- c("PKM","SLC2A1","CAT","GLUD1","ACAA1","LCAT","ACADSB")
housekeeper <- c("GAPDH","ACTB","B2M","RPLP0","HPRT1")

all_markers <- unique(c(epithelial, caf_markers, macro_genes, tcell_genes, checkpoint, metabolism, housekeeper))

# ============================================================
# Step 1: 检查原始文件结构
# ============================================================
log_msg("=== Step 1: 检查原始文件结构 ===")

csv_files <- list.files("data/raw/GSE138709", pattern="*.csv.gz$", full.names=TRUE)
log_msg(sprintf("CSV files: %d", length(csv_files)))

raw_structure <- data.frame(stringsAsFactors=FALSE)

for (f in csv_files) {
  fname <- basename(f)
  log_msg(sprintf("\n--- %s ---", fname))

  # Read using R's gzfile
  con <- gzfile(f, "r")
  lines_100 <- readLines(con, n=101)
  close(con)

  header <- strsplit(lines_100[1], ",")[[1]]
  ncols <- length(header)

  # Gene names from column 1 (skip header)
  gene_col <- sapply(lines_100[-1], function(l) strsplit(l, ",")[[1]][1], USE.NAMES=FALSE)
  gene_col <- as.character(gene_col)

  # Total rows: count lines in gz file
  con2 <- gzfile(f, "r")
  nrows <- length(readLines(con2)) - 1  # minus header
  close(con2)
  nrows <- nrows  # already computed above via readLines

  log_msg(sprintf("  Rows (excl header): %d, Cols: %d", nrows, ncols))
  log_msg(sprintf("  Header[1]: '%s'", header[1]))
  log_msg(sprintf("  Header[2]: '%s'", header[2]))
  log_msg(sprintf("  First genes: %s", paste(head(gene_col, 10), collapse=", ")))
  log_msg(sprintf("  Last genes: %s", paste(tail(gene_col, 10), collapse=", ")))

  # Determine orientation
  if (header[1] == "" || header[1] == '""') {
    orientation <- "genes_x_cells (column 1 empty = gene names)"
  } else if (nrows > ncols * 5) {
    orientation <- "genes_x_cells (many more rows than columns)"
  } else if (ncols > nrows * 5) {
    orientation <- "cells_x_genes (many more columns than rows)"
  } else {
    orientation <- "unclear"
  }

  # Search markers in gene column
  marker_hits <- all_markers[all_markers %in% gene_col]
  marker_hits_ci <- all_markers[toupper(all_markers) %in% toupper(gene_col)]

  raw_structure <- rbind(raw_structure, data.frame(
    file = fname, rows = nrows, cols = ncols,
    header1 = substr(header[1], 1, 50), header2 = substr(header[2], 1, 50),
    first_genes = paste(head(gene_col, 5), collapse="; "),
    orientation = orientation,
    exact_markers_found = length(marker_hits),
    case_insensitive_found = length(marker_hits_ci),
    found_markers = paste(marker_hits, collapse="; "),
    stringsAsFactors = FALSE
  ))
}

write.csv(raw_structure, "tables/single_cell_diagnosis/GSE138709_raw_file_structure.csv", row.names=FALSE)
log_msg("\nSaved raw file structure table")

# ============================================================
# Step 2: 在 raw files 中深度搜索 marker
# ============================================================
log_msg("\n=== Step 2: 在 raw files 中深度搜索 marker ===")

marker_search_raw <- data.frame(stringsAsFactors=FALSE)

# Use the first tumor CSV for thorough search
tumor_csv <- csv_files[grepl("Tumor", basename(csv_files))][1]
log_msg(sprintf("Using: %s", basename(tumor_csv)))

# Read all gene names (column 1) from this file
con_t <- gzfile(tumor_csv, "r")
all_lines <- readLines(con_t)
close(con_t)
all_genes_raw <- sapply(all_lines[-1], function(l) strsplit(l, ",")[[1]][1], USE.NAMES=FALSE)
all_genes_raw <- as.character(all_genes_raw)
log_msg(sprintf("Total unique genes in raw file: %d", length(unique(all_genes_raw))))

# Show all gene name formats
log_msg("Gene name examples (first 50):")
print(head(all_genes_raw, 50))

# Check for special patterns
has_ensg <- sum(grepl("^ENSG", all_genes_raw))
has_dot_version <- sum(grepl("\\.\\d+$", all_genes_raw))
has_pipe <- sum(grepl("\\|", all_genes_raw))
has_dash <- sum(grepl("-", all_genes_raw))
has_underscore <- sum(grepl("_", all_genes_raw))
is_allcaps <- sum(grepl("^[A-Z][A-Z0-9]+$", all_genes_raw))

log_msg(sprintf("ENSG prefix: %d", has_ensg))
log_msg(sprintf("Dot version (.1): %d", has_dot_version))
log_msg(sprintf("Pipe separator: %d", has_pipe))
log_msg(sprintf("Contains dash: %d", has_dash))
log_msg(sprintf("Contains underscore: %d", has_underscore))
log_msg(sprintf("ALL CAPS format: %d", is_allcaps))

# Search each marker thoroughly
for (mk in all_markers) {
  exact <- mk %in% all_genes_raw
  ci <- sum(toupper(all_genes_raw) == toupper(mk)) > 0
  # Partial: gene name CONTAINS the marker
  partial <- grep(mk, all_genes_raw, ignore.case=TRUE, value=TRUE)
  partial_exact <- grep(paste0("^", mk, "$"), all_genes_raw, ignore.case=TRUE, value=TRUE)
  # Check if marker appears with a suffix
  with_suffix <- grep(paste0("^", mk, "[-\\.\\|]"), all_genes_raw, value=TRUE)

  marker_search_raw <- rbind(marker_search_raw, data.frame(
    marker = mk,
    exact = exact,
    case_insensitive = ci,
    partial_count = length(partial),
    partial_examples = paste(head(partial, 3), collapse="; "),
    with_suffix = paste(head(with_suffix, 3), collapse="; "),
    stringsAsFactors = FALSE
  ))
}

write.csv(marker_search_raw, "tables/single_cell_diagnosis/GSE138709_marker_search_in_raw_files.csv", row.names=FALSE)
log_msg(sprintf("Markers exact match: %d/%d", sum(marker_search_raw$exact), nrow(marker_search_raw)))
log_msg(sprintf("Markers case-insensitive match: %d/%d", sum(marker_search_raw$case_insensitive), nrow(marker_search_raw)))
log_msg(sprintf("Markers with any partial match: %d/%d", sum(marker_search_raw$partial_count > 0), nrow(marker_search_raw)))

# ============================================================
# Step 3: 检查 Seurat object features
# ============================================================
log_msg("\n=== Step 3: 检查 Seurat object features ===")

merged <- readRDS("data/processed/GSE138709/GSE138709_seurat_processed.rds")

ncells <- ncol(merged)
ngenes <- nrow(merged)
log_msg(sprintf("Seurat: %d cells x %d genes", ncells, ngenes))

obj_genes <- rownames(merged)
log_msg(sprintf("First 5 features: %s", paste(head(obj_genes, 5), collapse=", ")))
log_msg(sprintf("Last 5 features: %s", paste(tail(obj_genes, 5), collapse=", ")))

# Check overlap with raw file genes
common_genes <- intersect(obj_genes, all_genes_raw)
log_msg(sprintf("Genes in both Seurat and raw: %d / %d (%.1f%%)",
  length(common_genes), length(obj_genes), length(common_genes)/length(obj_genes)*100))

# Feature format analysis
seurat_has_ensg <- sum(grepl("^ENSG", obj_genes))
seurat_has_dot <- sum(grepl("\\.\\d+$", obj_genes))
seurat_has_dash <- sum(grepl("-", obj_genes))
seurat_allcaps <- sum(grepl("^[A-Z][A-Z0-9]+$", obj_genes))

log_msg(sprintf("Features: ENSG=%d, dot_version=%d, dash=%d, ALLCAPS=%d",
  seurat_has_ensg, seurat_has_dot, seurat_has_dash, seurat_allcaps))

# Show features whose names changed (raw → Seurat)
raw_set <- all_genes_raw
seurat_set <- obj_genes
only_raw <- setdiff(raw_set, seurat_set)
only_seurat <- setdiff(seurat_set, raw_set)
log_msg(sprintf("Genes only in raw: %d, only in Seurat: %d", length(only_raw), length(only_seurat)))
if (length(only_raw) > 0) log_msg(sprintf("  Raw-only examples: %s", paste(head(only_raw, 10), collapse=", ")))
if (length(only_seurat) > 0) log_msg(sprintf("  Seurat-only examples: %s", paste(head(only_seurat, 10), collapse=", ")))

# Feature example table
feature_examples <- data.frame(
  index = 1:min(200, ngenes),
  feature_name = obj_genes[1:min(200, ngenes)],
  stringsAsFactors = FALSE
)
feature_examples$has_dot <- grepl("\\.\\d+$", feature_examples$feature_name)
feature_examples$has_dash <- grepl("-", feature_examples$feature_name)
feature_examples$allcaps <- grepl("^[A-Z][A-Z0-9]+$", feature_examples$feature_name)
write.csv(feature_examples, "tables/single_cell_diagnosis/GSE138709_feature_name_examples.csv", row.names=FALSE)

# ============================================================
# Step 4: 在 Seurat object 中搜索 marker
# ============================================================
log_msg("\n=== Step 4: 在 Seurat object 中搜索 marker ===")

marker_search_seurat <- data.frame(stringsAsFactors=FALSE)

for (mk in all_markers) {
  exact <- mk %in% obj_genes
  ci <- sum(toupper(obj_genes) == toupper(mk)) > 0
  partial <- grep(mk, obj_genes, ignore.case=TRUE, value=TRUE)

  # Also check in cell barcodes (colnames)
  in_colnames <- mk %in% colnames(merged)
  in_meta <- mk %in% colnames(merged@meta.data)

  marker_search_seurat <- rbind(marker_search_seurat, data.frame(
    marker = mk,
    exact_in_features = exact,
    case_insensitive = ci,
    partial_count = length(partial),
    partial_examples = paste(head(partial, 3), collapse="; "),
    in_cell_barcodes = in_colnames,
    in_metadata = in_meta,
    stringsAsFactors = FALSE
  ))
}

write.csv(marker_search_seurat, "tables/single_cell_diagnosis/GSE138709_marker_search_in_Seurat_object.csv", row.names=FALSE)
log_msg(sprintf("Seurat exact match: %d/%d", sum(marker_search_seurat$exact_in_features), nrow(marker_search_seurat)))
log_msg(sprintf("Seurat partial match: %d/%d", sum(marker_search_seurat$partial_count > 0), nrow(marker_search_seurat)))

# ============================================================
# Step 5: 矩阵方向诊断
# ============================================================
log_msg("\n=== Step 5: 矩阵方向诊断 ===")

# Transpose test: what if we swap rows and columns?
transpose_hypothesis <- data.frame(
  test = c(
    "Cells count (expected ~300-10000 per sample)",
    "Genes count (expected ~15000-25000)",
    "Cells/genes ratio (should be << 1 for scRNA-seq)",
    "Markers in colnames (would indicate matrix is transposed)",
    "Markers in features (expected location)",
    "Colnames look like cell barcodes",
    "Rownames look like genes",
    "Orientation from raw file"
  ),
  result = c(
    as.character(ncells),
    as.character(ngenes),
    sprintf("%.3f", ncells/ngenes),
    as.character(sum(all_markers %in% colnames(merged))),
    as.character(sum(all_markers %in% rownames(merged))),
    paste(head(colnames(merged), 3), collapse="; "),
    paste(head(rownames(merged), 3), collapse="; "),
    raw_structure$orientation[1]
  ),
  interpretation = c(
    if (ncells > 50000) "WARNING: unusually many cells" else "OK",
    if (ngenes < 5000) "WARNING: too few genes — might be cells" else "OK",
    if (ncells/ngenes > 1) "WARNING: more cells than genes — possible transpose" else "OK",
    if (sum(all_markers %in% colnames(merged)) > 5) "WARNING: markers in column names!" else "OK",
    if (sum(all_markers %in% rownames(merged)) > 5) "OK: markers found" else "WARNING: no markers",
    "Cell barcodes: should be long strings",
    "Gene symbols: should be short names",
    raw_structure$orientation[1]
  ),
  stringsAsFactors = FALSE
)

write.csv(transpose_hypothesis, "tables/single_cell_diagnosis/GSE138709_orientation_diagnosis.csv", row.names=FALSE)
log_msg("Saved orientation diagnosis")

# ============================================================
# Step 6: 综合诊断结论
# ============================================================
log_msg("\n=== Step 6: 综合诊断 ===")

# Determine primary issue
nm_exact_raw <- sum(marker_search_raw$exact)
nm_partial_raw <- sum(marker_search_raw$partial_count > 0)
nm_exact_seurat <- sum(marker_search_seurat$exact_in_features)
nm_partial_seurat <- sum(marker_search_seurat$partial_count > 0)

has_markers_raw <- nm_exact_raw > 5 || nm_partial_raw > 10
has_markers_seurat <- nm_exact_seurat > 5 || nm_partial_seurat > 10
ratio_suspicious <- ncells/ngenes > 1
is_transposed <- ratio_suspicious && (sum(all_markers %in% colnames(merged)) > 5)

# Diagnosis
if (!has_markers_raw) {
  diag_issue <- "GENE_ID_MISMATCH"
  diag_detail <- "Markers not found in raw CSV files. Gene identifiers are not standard HGNC symbols."
  diag_solution <- "Need gene ID mapping (e.g., Gencode transcript ID → HGNC symbol via biomaRt)"
} else if (!has_markers_seurat && has_markers_raw) {
  diag_issue <- "SEURAT_IMPORT_CORRUPTION"
  diag_detail <- "Markers found in raw files but NOT in Seurat object. Gene names were modified during import."
  diag_solution <- "Re-check CreateSeuratObject logic; possible gene name mangling by Seurat v5 layer merging"
} else if (is_transposed) {
  diag_issue <- "MATRIX_TRANSPOSED"
  diag_detail <- "Matrix orientation appears reversed (cells as rows, genes as columns)"
  diag_solution <- "Transpose matrix before CreateSeuratObject"
} else {
  diag_issue <- "NEED_DEEPER_INVESTIGATION"
  diag_detail <- sprintf("Markers in raw: %d exact, %d partial. Markers in Seurat: %d exact, %d partial.", nm_exact_raw, nm_partial_raw, nm_exact_seurat, nm_partial_seurat)
  diag_solution <- "Review gene name format carefully"
}

cat(sprintf("\n--- DIAGNOSIS ---\n"))
cat(sprintf("Issue: %s\n", diag_issue))
cat(sprintf("Detail: %s\n", diag_detail))
cat(sprintf("Solution: %s\n", diag_solution))

conclusions <- data.frame(
  question = c(
    "1. 矩阵方向错误?",
    "2. Gene symbol 被破坏?",
    "3. 需要 Ensembl ID 映射?",
    "4. 应该重新读取 raw CSV?",
    "5. 应该继续使用 GSE138709?",
    "6. 需要切换 GSE125449/GSE151530?"
  ),
  answer = c(
    if (is_transposed) "YES — matrix is transposed" else "No — orientation appears correct",
    if (!has_markers_seurat && has_markers_raw) "YES — markers in raw but not in Seurat" else if (!has_markers_raw) "No — markers absent from raw data too, not a Seurat issue" else "No — markers preserved",
    if (!has_markers_raw && sum(marker_search_raw$partial_count > 0) < 5) "YES — raw gene IDs are non-standard, need mapping" else "No — gene symbols are usable",
    if (!has_markers_seurat && has_markers_raw) "YES — re-read with fixed import" else "No — import is fine",
    if (!has_markers_raw) "NO — gene IDs cannot be used without mapping; consider switching" else "YES — but fix the specific issue first",
    if (!has_markers_raw && !has_markers_seurat) "Recommended — GSE138709 uses non-standard IDs" else "Not yet — fix GSE138709 first"
  ),
  stringsAsFactors = FALSE
)

write.csv(conclusions, "tables/single_cell_diagnosis/GSE138709_gene_identifier_diagnosis_summary.csv", row.names=FALSE)

cat(sprintf("\n========================================\n"))
cat(sprintf("诊断完成\n"))
cat(sprintf("========================================\n"))
print(conclusions)
cat(sprintf("========================================\n"))
