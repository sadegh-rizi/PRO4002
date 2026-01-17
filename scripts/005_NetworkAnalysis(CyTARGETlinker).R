# ==============================================================================
# Script: 005_miRNA_Cluster_Analysis.R
# Purpose: Check if General DCM miRNAs are specific to any Cluster
# ==============================================================================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("multiMiR")

# ==============================================================================
# Script: 005_miRNA_Target_Integration.R
# Purpose: Link 19 DCM miRNAs to Cluster Endotypes via Target Prediction
# ==============================================================================

source(here::here("scripts", "000_setup.R"))
library(multiMiR)
library(dplyr)

message("\n--- Starting Target Integration ---")

# ------------------------------------------------------------------------------
# 1. LOAD YOUR LISTS
# ------------------------------------------------------------------------------
message("(I) Loading Lists")

# A. Load the 19 Significant miRNAs (From your Table 2)
# We need the "gene_symbol" column (e.g., "MIR27B")

setwd(data_path)
mirna_file <- file.data(data_path, "Table2_Differential_Expression_Results_miRNAs.csv")
mirna_table <- read.csv("Table2_Differential_Expression_Results_miRNAs.csv")
sig_mirnas <- subset(mirna_table, comparison == "DCM_vs_Control" & q.value < 0.05)

# Clean names for Database (remove "MIR" prefix often helps, or keep standard)
# multiMiR expects "hsa-miR-27b". Your table has "MIR27B". We need to convert.
# Quick/Dirty conversion: "MIR27B" -> "hsa-miR-27b"
# (Note: This is an estimation. For publication, check exact ID mapping).
candidate_mirnas <- paste0("hsa-miR-", tolower(sub("MIR", "", sig_mirnas$gene_symbol)))
# Handle special cases (e.g., "MIR378A" -> "hsa-miR-378a")
candidate_mirnas <- gsub("mir([0-9]+)([a-z]*)", "mir-\\1\\2", candidate_mirnas)

message(paste("Scanning targets for", length(candidate_mirnas), "miRNAs..."))

# B. Load Cluster Marker Genes (The CSVs you generated in 002)
# We will check Cluster 1 (The Fibrotic/Mitochondrial one) first as proof of concept.
c1_genes <- read.csv(file.path(tables_path, "Cluster_1_Top_Genes.csv"))
c2_genes <- read.csv(file.path(tables_path, "Cluster_2_Top_Genes.csv"))
c3_genes <- read.csv(file.path(tables_path, "Cluster_3_Top_Genes.csv"))

view(c1_genes)

# Combine into one target list with cluster tags
all_targets <- rbind(
  c1_genes %>% mutate(Cluster = "1"),
  c2_genes %>% mutate(Cluster = "2"),
  c3_genes %>% mutate(Cluster = "3")
)


# ------------------------------------------------------------------------------
# 2. RUN PREDICTION (multiMiR)
# ------------------------------------------------------------------------------
message("(II) Querying Databases (TargetScan, miRTarBase, etc.)")
# This might take 1-2 minutes

results <- get_multimir(
  mirna = candidate_mirnas,
  target = all_targets$Symbol, # Check if your miRNAs hit these genes
  table = "validated",         # Use "validated" for strong evidence (experimental)
  summary = TRUE
)

# If 'validated' yields 0 results, switch table to "predicted" (TargetScan)
if(nrow(results@data) == 0) {
  message("No experimental hits found. Switching to Predicted targets (TargetScan)...")
  results <- get_multimir(
    mirna = candidate_mirnas,
    target = all_targets$Symbol,
    table = "predicted",
    summary = TRUE
  )
}

interactions <- results@data


#TRYING
# 1. Get Validated (Experimental Proof)
results_val <- get_multimir(
  mirna = candidate_mirnas,
  target = all_targets$Symbol,
  table = "validated",
  summary = TRUE
)

# 2. Get Predicted (TargetScan - Mathematical Guesses)
# We run this REGARDLESS of whether we found validated hits
results_pred <- get_multimir(
  mirna = candidate_mirnas,
  target = all_targets$Symbol,
  table = "predicted",
  summary = TRUE
)

# 3. Combine Them
# Use bind_rows or simple rbind if columns match. 
# multiMiR structures are complex, so we extract @data first.
interactions <- dplyr::bind_rows(results_val@data, results_pred@data)

message(paste("  Validated Hits:", nrow(results_val@data)))
message(paste("  Predicted Hits:", nrow(results_pred@data)))

# ------------------------------------------------------------------------------
# 3. FILTER & MERGE
# ------------------------------------------------------------------------------
message("(III) Processing Interactions")

# We only care about "Inverse Correlations" typically
# (Upregulated miRNA -> Downregulated Gene)
# Let's merge the Interaction data with your Cluster Data

# Clean miRNA names back to match (optional)
interactions$Symbol <- interactions$target_symbol

# Merge with your Cluster Gene info (to see logFC of the gene)
final_network <- merge(interactions, all_targets, by="Symbol")

# Select useful columns
final_network <- final_network %>%
  dplyr::select(miRNA = mature_mirna_id, 
                Target_Gene = Symbol, 
                Target_Cluster = Cluster, 
                Gene_logFC = logFC, 
                Database = database, 
                Evidence = experiment) %>%
  distinct()

# ------------------------------------------------------------------------------
# 4. EXPORT
# ------------------------------------------------------------------------------
outfile <- file.path(tables_path, "miRNA_Cluster_Network.csv")
write.csv(final_network, outfile)

message(paste("Found", nrow(final_network), "connections."))
message("Check 'miRNA_Cluster_Network.csv' to see which miRNAs are driving your Clusters.")

# ------------------------------------------------------------------------------
# 5. QUICK VIEW
# ------------------------------------------------------------------------------
# Show hits for Cluster 1 (Mitochondrial/Fibrosis)
print("--- TOP HITS FOR CLUSTER 1 ---")
c1_hits <- final_network %>% filter(Target_Cluster == "1")
print(head(c1_hits))

c2_hits <- final_network %>% filter(Target_Cluster == "2")
print(head(c2_hits))

c3_hits <- final_network %>% filter(Target_Cluster == "3")
print(head(c3_hits))

view(c1_hits)
view(c2_hits)
view(c3_hits)
