
# Differential expression: DCM Cluster 2 vs Cluster 1
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(limma)
  library(here)
})

dir.create(here("results"), showWarnings = FALSE, recursive = TRUE)

# 1) Load expression matrix

expr <- fread(here("data", "MAGNET_GeneExpressionData_CPM_19112020.txt")) %>%
  as.data.frame()

colnames(expr) <- gsub('"', "", colnames(expr))
expr$EnsemblGeneID <- gsub('"', "", expr$EnsemblGeneID)

rownames(expr) <- expr$EnsemblGeneID
expr <- expr[, -1]

expr_mat <- as.matrix(expr)
mode(expr_mat) <- "numeric"

# 2) Load cluster assignments

clust <- fread(here("results", "dcm_cluster_assignments.csv")) %>%
  as.data.frame()

clust$cluster <- factor(clust$cluster, levels = c("C1", "C2"))

# Reorder samples to match expression columns
expr_mat <- expr_mat[, clust$sample_name]

# 3) limma design

design <- model.matrix(~ cluster, data = clust)
colnames(design) <- c("Intercept", "C2_vs_C1")

fit <- lmFit(expr_mat, design)
fit <- eBayes(fit)

# 4) Extract DE genes

de_genes <- topTable(
  fit,
  coef = "C2_vs_C1",
  number = Inf,
  sort.by = "P"
) %>%
  rownames_to_column("gene_symbol")

write.csv(
  de_genes,
  here("results", "DE_genes_C2_vs_C1.csv"),
  row.names = FALSE
)

message("Saved: results/DE_genes_C2_vs_C1.csv")