
# Integrate DE miRNAs with cluster-specific DE genes

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(multiMiR)
  library(here)
})

dir.create(here("results"), showWarnings = FALSE, recursive = TRUE)

# 1) Load DE table (DCM vs control) with gene_biotype + q.value

mirna_de <- fread(here("data", "Table2_Differential_Expression_Results_miRNAs.csv")) %>%
  as.data.frame()

colnames(mirna_de) <- tolower(colnames(mirna_de))

# Significant miRNAs (based on biotype + q.value)
sig_mirna <- mirna_de %>%
  filter(tolower(gene_biotype) == "mirna", q.value < 0.05) %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  distinct(gene_symbol, .keep_all = TRUE)

message("Significant miRNAs (q<0.05): ", nrow(sig_mirna))

# 2) Load cluster-based DE genes (C2 vs C1)

gene_de <- fread(here("results", "DE_genes_C2_vs_C1.csv")) %>%
  as.data.frame()

colnames(gene_de) <- tolower(colnames(gene_de))

# Significant cluster-DE genes (adj.p.val < 0.05)
sig_genes <- gene_de %>%
  filter(adj.p.val < 0.05) %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  distinct(gene_symbol, .keep_all = TRUE)

message("Significant cluster-DE genes (adj.p.val<0.05): ", nrow(sig_genes))

# 3) Prepare miRNA names for multiMiR

mirna_query_raw <- sig_mirna$gene_symbol %>%
  as.character() %>%
  trimws() %>%
  .[. != ""] %>%
  unique()

mirna_query <- ifelse(
  grepl("^MIRLET", mirna_query_raw),
  paste0("hsa-let-", tolower(sub("^MIRLET", "", mirna_query_raw))),
  paste0("hsa-miR-", tolower(sub("^MIR", "", mirna_query_raw)))
)

message("Example miRNA query names: ", paste(head(mirna_query, 5), collapse = ", "))

# 4) Query multiMiR (predicted targets)

mm <- get_multimir(
  mirna = mirna_query,
  table = "predicted",
  legacy.out = TRUE
)

mirna_targets <- mm$predicted

if (is.null(mirna_targets) || nrow(mirna_targets) == 0) {
  stop("multiMiR returned 0 predicted targets. Check miRNA naming / internet / package installation.")
}

message("multiMiR predicted targets returned: ", nrow(mirna_targets), " rows")

# 5) Identify the target gene symbol column and filter to cluster-DE genes

target_col_candidates <- c("target_symbol", "target", "target_gene", "targetgene", "target.gene")
target_col <- target_col_candidates[target_col_candidates %in% colnames(mirna_targets)][1]

if (is.na(target_col)) {
  stop("Could not find a target gene column in mirna_targets. Columns are: ",
       paste(colnames(mirna_targets), collapse = ", "))
}

mirna_targets_filt <- mirna_targets %>%
  filter(.data[[target_col]] %in% sig_genes$gene_symbol)

message("Edges after intersecting with cluster-DE genes: ", nrow(mirna_targets_filt))

# 6) Export network edge list for Cytoscape

mirna_col_candidates <- c("mirna", "mature_mirna_id", "mirna_name", "mature_mirna", "mirna_id")
mirna_col <- mirna_col_candidates[mirna_col_candidates %in% colnames(mirna_targets_filt)][1]

if (is.na(mirna_col)) {
  stop("Could not find a miRNA column in mirna_targets_filt. Columns are: ",
       paste(colnames(mirna_targets_filt), collapse = ", "))
}
edge_list <- mirna_targets_filt %>%
  transmute(
    mirna = .data[[mirna_col]],
    target_symbol = .data[[target_col]]
  ) %>%
  distinct()

write.csv(
  edge_list,
  here("results", "miRNA_gene_network_edges.csv"),
  row.names = FALSE
)

message("Exported: results/miRNA_gene_network_edges.csv")
message("Done.")