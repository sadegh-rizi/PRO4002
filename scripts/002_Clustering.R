 # ==============================================================================
# Script: 002_DCM_Consensus_Clustering.R
# Purpose: Identify DCM Endotypes clusters
# Inputs: 'geneExpressionData_DCM_meanFiltered.rds' (from QC script)
# ==============================================================================


analysis_packages <- c( 
  "tidyverse", "tidyr", "dplyr", "gridExtra", "pcaMethods",
  "data.table", "tableone", "kableExtra", "rmarkdown",
  "readr", "readxl", "gprofiler2", "knitr", "Hmisc", "table1",
  "gt", "gtsummary", "ggsci", "DESeq2", "edgeR"
)

bio_packages <- c("limma", "qvalue", "biomaRt", "Biocmanager")

for (pkg in analysis_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message(paste("Installing Analysis package:", pkg))
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in bio_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message(paste("Installing Bioconductor package:", pkg))
    BiocManager::install(pkg)
  }
  library(pkg, character.only = TRUE)
}

search()

#-----------------------------------------------------------------------------#
# STYLE SETUP
#-----------------------------------------------------------------------------#

# Prepare Color Values to Match Nature Color Scheme
npg_colors <- pal_npg("nrc", alpha = 1)(10)
npg_additional_colors <- colorRampPalette(npg_colors)(20)
continuous_npg_colors <- colorRampPalette(c(npg_colors[2], "white"))(100)

# GGplot Plotting Settings 
center_title <- theme(plot.title = element_text(hjust = 0.5, vjust = 1))

target_font <- "sans" 
my_style <- theme(
  text = element_text(family = target_font),        
  plot.title = element_text(hjust = 0.5, face="bold"),
  axis.title = element_text(face = "bold"),    
  legend.position = "right"                    
)


#-----------------------------------------------------------------------------#
# SETUP
#-----------------------------------------------------------------------------#

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) { 
  script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(script_path)
} else {
  warning("Not running in RStudio. Please set working directory manually for portability.") 
}

data_path <- file.path(dirname(getwd()), "data")
if(!dir.exists(data_path)) { dir.create(data_path) }
message(paste("Data Directory:", data_path))

plots_path <- file.path(dirname(getwd()), "plots")
if(!dir.exists(plots_path)) { dir.create(plots_path) }
message(paste("Plots Directory:", plots_path))

results_path <- file.path(dirname(getwd()), "results")
if(!dir.exists(results_path)) { dir.create(results_path) }
message(paste("Results Directory:", results_path))

tables_path <- file.path(dirname(getwd()), "tables")
if(!dir.exists(tables_path)) { dir.create(tables_path) }
message(paste("Tables Directory:", tables_path))

cache_path <- file.path(dirname(getwd()), "cache")
if(!dir.exists(cache_path)) { dir.create(cache_path) }
message(paste("Cache Directory:", cache_path))

# ------------------------------------------------------------------------------
# 1. CONFIRM LIBRARIES
# ------------------------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE)) { install.packages("BiocManager") }
BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)
BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(tidyverse)
library(ConsensusClusterPlus) # The engine for "Robust K-means"
library(pheatmap)             # For the final look
library(factoextra)           # For PCA visualization
library(edgeR)                # For Log2 normalization
library(readxl)               # To read metadata
library(ggpubr)               # For clinical stats


------------------------------------------------------------------------------
  # 2. LOAD DATA (Output from 001_QualityControl.R)
  # ------------------------------------------------------------------------------
# Load the filtered DCM expression data (Mean Threshold Filter applied)
dcm_cpm <- readRDS("geneExpressionData_DCM_meanFiltered.rds")

# Load Raw Metadata
sample_data_raw <- read_excel("MAGNET_SampleData_18112022_WithDescriptions.xlsx", sheet=1)

# Clean Metadata (Match 001_QualityControl logic)
dcm_meta <- sample_data_raw %>%
  filter(etiology == "DCM") %>%
  mutate(
    # FIX 1: Handle European decimals (comma) or symbols (%)
    LVEF_clean = gsub(",", ".", LVEF),        # Swap 30,5 to 30.5
    LVEF_clean = gsub("%", "", LVEF_clean),   # Remove % sign
    
    # FIX 2: Convert to numeric safely
    LVEF = suppressWarnings(as.numeric(LVEF_clean)),
    
    age = as.numeric(age),
    gender = as.factor(gender)
  ) %>%
  # Remove the temporary cleaning column
  dplyr::select(-LVEF_clean) %>% 
  data.frame()

# Align Rownames for safety
rownames(dcm_meta) <- dcm_meta$sample_name

# ------------------------------------------------------------------------------
# 3. ALIGNMENT CHECK (Safety Step)
# ------------------------------------------------------------------------------
# Ensure expression columns match metadata rows
common_ids <- intersect(colnames(dcm_cpm), rownames(dcm_meta))
dcm_cpm <- dcm_cpm[, common_ids]
dcm_meta <- dcm_meta[common_ids, ]

print(paste("Analysis ready for:", ncol(dcm_cpm), "DCM Patients"))



# ------------------------------------------------------------------------------
# 4. NORMALIZATION & FEATURE SELECTION
# ------------------------------------------------------------------------------

# A. Log2 Transformation
# Essential for Pearson Correlation and Heatmaps
# We use log2(CPM + 1) to handle zeros
dcm_log2 <- dcm_cpm #data is already transform, this is just a check 

# B. Variance Filtering (Feature Selection)
# We want the genes that *change* between patients.
# Calculate Variance for every gene
gene_vars <- apply(dcm_log2, 1, var)

# Select Top 2000 Genes (Consensus between 1000 and 5000)
top_n <- 2000
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:top_n]

dcm_variable <- dcm_log2[top_genes, ]

# C. Z-Score Scaling
# Centers genes so "High" is Red and "Low" is Blue.
# Transpose so Scale works on columns (genes), then transpose back.
dcm_scaled <- t(scale(t(dcm_variable)))

print(paste("Clustering Matrix Dimensions:", nrow(dcm_scaled), "Genes x", ncol(dcm_scaled), "Patients"))

# ------------------------------------------------------------------------------
# 5. CONSENSUS CLUSTERING 
# ------------------------------------------------------------------------------
# Method: K-Means
# Distance: Pearson Correlation (1 - r) #not possible with consensus, as Z is check we are trying with euclidean for now
# Reps: 1000 (Very robust)
# ------------------------------------------------------------------------------

# 2. Count is there is any NA or infinite value as a check
num_bad <- sum(is_bad_row)
print(paste("Found", num_bad, "genes that failed scaling (Zero Variance or NA)."))

if(num_bad > 0) {
  print("REMOVING these genes to prevent crashes...")
  
  # 3. Subset to keep only the "Good" rows
  dcm_scaled_clean <- dcm_scaled[!is_bad_row, ]
  
} else {
  print("Matrix is clean. No removal needed.")
  dcm_scaled_clean <- dcm_scaled
}

# 4. Final Safety Check (Dimensions)
print(paste("Final Matrix for Clustering:", nrow(dcm_scaled_clean), "Genes x", ncol(dcm_scaled_clean), "Patients"))
title_dir <- "Consensus_Plots"

results <- ConsensusClusterPlus(
  dcm_scaled_clean,
  maxK = 6,              # Check k=2 to k=6
  reps = 1000,           # Resample 1000 times
  pItem = 0.8,           # Use 80% of patients per run
  pFeature = 1,          # Use 100% of the top 2000 genes
  title = title_dir,
  clusterAlg = "hc",     # "km" = K-MEANS
  distance = "pearson",  # "pearson" = CORRELATION
  seed = 123,
  plot = "png"
)


# ------------------------------------------------------------------------------
# 6. EXTRACT CLUSTERS (Decision Time)
# ------------------------------------------------------------------------------

final_k <- 3
dcm_meta$Endotype <- as.factor(results[[final_k]][["consensusClass"]])

print("Final Cluster Sizes (k=3):")
print(table(dcm_meta$Endotype))

# ------------------------------------------------------------------------------
# 7. VALIDATION 1: STABILITY Y PCA PLOT
# ------------------------------------------------------------------------------
library(cluster)
library(factoextra)
# ------------------------------------------------------------------------------
# A. CLUSTER STABILITY (Silhouette Plot)
# ------------------------------------------------------------------------------
# Rationale: Measures how similar a patient is to its own cluster vs. others.
# Range: -1 (Wrong cluster) to +1 (Perfect fit).

print("--- Calculating Cluster Stability (Silhouette Width) ---")

# 1. Calculate Pearson Distance Matrix (Must match your clustering logic)
# We use (1 - correlation) as the distance
dist_matrix <- as.dist(1 - cor(dcm_scaled_clean, method = "pearson"))

# 2. Calculate Silhouette Scores
sil <- silhouette(as.numeric(dcm_meta$Endotype), dist_matrix)

# 3. Visualize
# This plot shows which patients are "core" members (high bars) 
# and which are ambiguous (short/negative bars).
sil_plot <- fviz_silhouette(sil, 
                            palette = "jco", 
                            ggtheme = theme_minimal(),
                            main = paste0("Cluster Stability (k=", final_k, "): Silhouette Plot"))

# 4. Save
ggsave("plots/Stability_Silhouette.pdf", sil_plot, width = 8, height = 6)

# 5. Print Summary Stats
avg_sil_width <- summary(sil)$avg.width
print(paste("Average Silhouette Width:", round(avg_sil_width, 3)))
print("Interpretation: > 0.25 is decent. > 0.5 is strong structure.")

# ------------------------------------------------------------------------------
# B. PCA PLOT (Visual Separation)
# ------------------------------------------------------------------------------


pca_plot <- fviz_cluster(list(data = t(dcm_scaled_clean), cluster = dcm_meta$Endotype),
                         geom = "point", ellipse.type = "convex", palette = "jco",
                         ggtheme = theme_minimal(), main = "PCA: DCM Endotypes (No Double-Log)")
ggsave("plots/PCA_Validation_Corrected.pdf", pca_plot, width = 6, height = 5)

# ------------------------------------------------------------------------------
# C. CLINICAL PHENOCOPY CHECK
# ------------------------------------------------------------------------------
p_lvef <- ggplot(dcm_meta, aes(x = Endotype, y = LVEF, fill = Endotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(method = "kruskal.test", label.y = 50) +
  labs(title = "Clinical Severity (LVEF)", y = "LVEF (%)") +
  theme_minimal() + scale_fill_jco()
ggsave("plots/Clinical_Validation_Corrected.pdf", p_lvef, width = 5, height = 4)

# ------------------------------------------------------------------------------
# D. HEATMAP
# ------------------------------------------------------------------------------
ann_df <- data.frame(Endotype = dcm_meta$Endotype)
rownames(ann_df) <- rownames(dcm_meta)
pheatmap(dcm_scaled_clean, show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = ann_df, cutree_cols = final_k,
         clustering_distance_cols = "correlation", 
         clustering_method = "ward.D2", scale = "none",
         main = "DCM Molecular Endotypes",
         filename = "plots/Final_Heatmap_Corrected.pdf", width = 8, height = 10)

# ------------------------------------------------------------------------------
# 8. SAVE DATA
# ------------------------------------------------------------------------------
saveRDS(dcm_meta, "dcm_meta_with_clusters.rds")
print("SUCCESS: Pipeline Corrected.")
dcm_cluster <- readRDS("dcm_meta_with_clusters.rds")

# ==============================================================================
#_Differential_Expression_K3.R
# Purpose: Identify Top Marker Genes for 3 Consensus Endotypes
# ==============================================================================


