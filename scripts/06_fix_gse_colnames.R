# Fix X-prefixed column names in GSE107943 CSVs
fix_csv <- function(path) {
  x <- read.csv(path, check.names = FALSE)
  cn <- colnames(x)
  for (i in seq_along(cn)) {
    if (cn[i] != "gene_symbol" && grepl("^X[0-9]", cn[i])) {
      cn[i] <- sub("^X", "", cn[i])
    }
  }
  colnames(x) <- cn
  write.csv(x, path, row.names = FALSE)
  cat(sprintf("Fixed: %s (%d cols)\n", path, length(cn)))
}
fix_csv("data/processed/GSE107943/GSE107943_log2CPM_gene_symbol.csv")
fix_csv("data/processed/GSE107943/GSE107943_counts_gene_symbol.csv")
cat("Done\n")
