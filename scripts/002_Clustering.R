###############################################################################
###############################################################################
# PRO4002 Research Project                                                    #
#																	                                            #
# File: 002_Clustering.R                                                      #
# Date: Jan 22, 2026											                                    #
# Author: Sadegh, Matas, Nur, Arlin & Oona                                    #  
###############################################################################
###############################################################################

#-----------------------------------------------------------------------------#
# SETUP
#-----------------------------------------------------------------------------#

source(here::here("scripts", "000_setup.R"))

message("\n--- Starting Clustering ---")

#-----------------------------------------------------------------------------#
# DATA IMPORT
#-----------------------------------------------------------------------------#

message("(I) Data Import ")

geneExpressionData.CPM.meanFiltered.DCM <- readRDS(file.path(cache_path, "geneExpressionData_DCM_meanFiltered.rds"))
sampleData.DCM <- readRDS(file.path(cache_path, "sampleData_DCM.rds"))
geneListInfo <- readRDS(file.path(cache_path, "geneListInfo.rds"))

#-----------------------------------------------------------------------------#
# PARAMETER SETUP
#-----------------------------------------------------------------------------#

seed <- 123
num_genes <- 1000    # Number of Genes to Use for Clustering
distance_measure <- "euclidean"
clustering_algorithm <- "km"

#-----------------------------------------------------------------------------#
# FEATURE SELECTION - Variance-Based Filtering (Highly Variable Genes)
#-----------------------------------------------------------------------------#

message("(II) Feature Selection ")


# Impute the one RIN value that is missing
sampleData.DCM$rin[is.na(sampleData.DCM$rin)] <- median(sampleData.DCM$rin, na.rm = TRUE)

# Remove dwares with weird BMIs
sampleData.DCM <- sampleData.DCM %>% filter(height>100)
rownames(sampleData.DCM) <- sampleData.DCM$sample_name

geneExpressionData.CPM.meanFiltered.DCM <- geneExpressionData.CPM.meanFiltered.DCM[,rownames(sampleData.DCM)]

geneData.variance <- apply(geneExpressionData.CPM.meanFiltered.DCM, 1, var)
geneData.variance.selected <- names(sort(geneData.variance, decreasing = TRUE))[1:num_genes]
geneExpressionData.variance.selected <- geneExpressionData.CPM.meanFiltered.DCM[geneData.variance.selected, ]
geneExpressionData.variance.selected.centered <- t(scale(t(geneExpressionData.variance.selected)))

#-----------------------------------------------------------------------------#
# CLUSTERING SETUP
#-----------------------------------------------------------------------------#

message("(III) Clustering ")

consensus_cluster_path <- file.path(plots_path, "ConsensusClustering")
if(!dir.exists(consensus_cluster_path)) dir.create(consensus_cluster_path, recursive = TRUE)

results <- ConsensusClusterPlus(
  geneExpressionData.variance.selected.centered,
  maxK = 6,              # Test k=2 to k=6 clusters
  reps = 1000,           # Resample 1000 times (High robustness)
  pItem = 0.8,           # Use 80% of samples per iteration
  pFeature = 1,          # Use 100% of the top genes
  title = consensus_cluster_path,      # Where to save the PDF plots
  clusterAlg = clustering_algorithm,     # Hierarchical clustering
  distance = distance_measure,
  seed = seed,
  plot = "png"
)

#-----------------------------------------------------------------------------#
# CLUSTERING ANALYSIS
#-----------------------------------------------------------------------------#

message("(III) Clustering - Analysis")

cluster_path <- file.path(plots_path, "ClusteringAnalysis")
if(!dir.exists(cluster_path)) dir.create(cluster_path, recursive = TRUE)

k_choice <- 3
clusterLabels <- results[[k_choice]][["consensusClass"]]
clusterLabels.formatted <- paste0("Cluster_", clusterLabels)
all(sampleData.DCM$sample_name == names(clusterLabels.formatted))
sampleData.DCM$subtype <- clusterLabels.formatted



table(sampleData.DCM$subtype)

geneExpressionData.variance.selected.centered <- as.matrix(geneExpressionData.variance.selected.centered)


sampleData.annotated <- sampleData.DCM %>%
  dplyr::select(subtype, gender, lvef, race) %>%
  mutate(lvef = lvef * 100) %>%
  as.data.frame()
rownames(sampleData.annotated) <- sampleData.DCM$sample_name

sampleData.annotated$subtype <- factor(sampleData.annotated$subtype)
sampleData.annotated$gender  <- factor(sampleData.annotated$gender)
sampleData.annotated$race  <- factor(sampleData.annotated$race)
sampleData.annotated$lvef    <- as.numeric(as.character(sampleData.annotated$lvef))



colors.annotated <- list(
  subtype = c(
    Cluster_1 = npg_colors[[1]],  
    Cluster_2 = npg_colors[[2]],  
    Cluster_3 = npg_colors[[3]]   
  ),
  gender = c(Male = "grey30", Female = "grey80"),
  race = c(Caucasian = "grey30", "African American" = "grey80"),
  lvef = colorRampPalette(c("white", npg_colors[[1]]))(100) # Continuous color scale for heart function
)

geneExpressionData.variance.selected.centered.ordered <- geneExpressionData.variance.selected.centered[,order(sampleData.annotated$subtype)]

cor_mat <- cor(
  geneExpressionData.variance.selected.centered.ordered,
  method = "pearson",
  use = "everything"
)
cor_mat_heat <- cov2cor(cor_mat)

pheatmap(
  cor_mat_heat,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  clustering_method = "complete",
  main = "Patient–Patient Gene Expression Correlation",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,      # Too many genes to show names
  show_colnames = FALSE,   
  annotation_col = sampleData.annotated,
  annotation_colors = colors.annotated,
  border_color = NA,
  filename = file.path(cluster_path, "subtypeHeatmap_DCM.jpg"),
)

pheatmap(
  geneExpressionData.variance.selected.centered.ordered,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c(npg_colors[[4]], "white", npg_colors[[1]]))(100),
  breaks = seq(-2, 2, length.out = 101), # Cap scale at z-score +/- 2
  show_rownames = TRUE,      # Too many genes to show names
  show_colnames = TRUE,      # Hiding sample names for cleaner look
  annotation_col = sampleData.annotated,
  annotation_colors = colors.annotated,
  na_col = "grey30",
  main = paste0("DCM Molecular Subtypes (Top ", num_genes, " Variable Genes)"),
  filename = file.path(cluster_path, "subtypeHeatmap_DCM.jpg"),
  width = 10, height = 8
)

#-----------------------------------------------------------------------------#
# CLUSTERING STABILITY 
#-----------------------------------------------------------------------------#

message("(IV) Clustering - Stability Analysis")

geneData.bootstrap <- clusterboot(
  t(geneExpressionData.variance.selected.centered), 
  B = 100,
  bootmethod = "boot", 
  clustermethod = kmeansCBI,
  k = k_choice,
  count = FALSE
) 

geneData.stabilityScores <- geneData.bootstrap$bootmean
message(paste("  Stability Scores:", round(geneData.stabilityScores, 3)))

jpeg(filename = file.path(cluster_path, "ClusterStability.jpg"), 
     width = 600, height = 500, quality = 90)
barplot(geneData.stabilityScores, 
        ylim = c(0, 1), 
        col = ifelse(geneData.stabilityScores > 0.6, "steelblue", "firebrick"),
        main = "Cluster Stability (Bootstrap)",
        ylab = "Jaccard Similarity Score",
        # This adds the "Cluster 1, Cluster 2..." labels automatically
        names.arg = paste("Cluster", 1:length(geneData.stabilityScores)))

abline(h = 0.6, lty = 2, col = "black")
dev.off()

print("SUCCESS: stability check complete.")

#-----------------------------------------------------------------------------#
# CLUSTERING VISUALIZATION
#-----------------------------------------------------------------------------#

message("(V) Clustering - Visualization")

# PCA PLOT 

geneExpressionDataPCAPlot <- fviz_cluster(list(data = t(geneExpressionData.variance.selected.centered), 
  cluster = sampleData.annotated$subtype), 
  geom = "point", ellipse.type = "convex", palette = npg_colors, 
  ggtheme = my_style, main = "PCA: DCM Endotypes (No Double-Log)") 
ggsave(file.path(cluster_path, "subtypePCAPlot.jpg"), geneExpressionDataPCAPlot, width = 6, height = 5)

# CLINICAL PHENOTYPE CHECK 

sampleDataLVEFPlot <- ggplot(sampleData.annotated, aes(x = subtype, y = lvef, fill = subtype)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + 
  geom_jitter(width = 0.2, alpha = 0.5) + 
  stat_compare_means(method = "kruskal.test", label.y = 50) + 
  labs(title = "Clinical Severity (LVEF)", y = "LVEF (%)") + 
  scale_color_npg() +
  scale_fill_npg() +
  center_title + my_style
ggsave(file.path(cluster_path, "subtypeLVEFPlot.jpg"), sampleDataLVEFPlot, width = 5, height = 4) 

#-----------------------------------------------------------------------------#
# DATA EXPORT
#-----------------------------------------------------------------------------#

message("(IV) Data Export")

saveRDS(sampleData.annotated, file.path(cache_path, "sampleData_DCM_subtypes.rds"))

#-----------------------------------------------------------------------------#
# CLUSTER DIFFERENTIAL GENE EXPRESSION ANALYSIS
#----------------------------------------------------------------------------#
message("(V) CLUSTER DGE Analysis")

message("  Cluster Counts:")
print(table(sampleData.annotated$subtype))

message(paste("  Analyzing", nrow(geneExpressionData.variance.selected.centered), "Genes across", ncol(geneExpressionData.variance.selected.centered), "Patients."))



print(paste("Patients in Expression Data:", ncol(geneExpressionData.variance.selected.centered)))
print(paste("Patients in Design Matrix:", nrow(design)))

# Check for NAs in your metadata (The likely killer)
print("--- Checking for NAs in Covariates ---")




subtypes <- factor(sampleData.DCM$subtype, levels = c("Cluster_1", "Cluster_2", "Cluster_3"))
design <- model.matrix(~ 0 + subtype+rin+age+gender+Library.Pool, data=sampleData.DCM)
colnames(design) <- gsub("subtypeCluster_", "C", colnames(design))# B. Fit Linear Model
fit <- lmFit(geneExpressionData.CPM.meanFiltered.DCM, design)

# C. Create Contrasts (One-vs-All for K=3)
# Logic: Compare Cluster X against the Average of the other TWO.
contrast_matrix <- makeContrasts(
  C1_Unique = C1 - (C2 + C3)/2,
  C2_Unique = C2 - (C1 + C3)/2,
  C3_Unique = C3 - (C1 + C2)/2,
  levels = design
)

# D. Apply Contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

export_top_genes <- function(contrast_name, cluster_id) {
  
  # Get Top 50 Genes
  top_table <- topTable(fit2, coef = contrast_name, number = 50, adjust.method = "fdr")
  
  # Annotate
  top_table$EnsemblID <- rownames(top_table)
  top_table$Symbol <- geneListInfo$hgnc_symbol[match(top_table$EnsemblID, geneListInfo$ensembl_gene_id)]
  top_table$Description <- geneListInfo$description[match(top_table$EnsemblID, geneListInfo$ensembl_gene_id)]
  
  # Clean & Sort
  final_table <- top_table %>%
    dplyr::select(Symbol, EnsemblID, logFC, adj.P.Val, Description) %>%
    arrange(adj.P.Val)
  
  # Print Preview
  print(paste("--- TOP MARKERS FOR CLUSTER", cluster_id, "---"))
  print(head(final_table, 5))
  
    # Save
  filename <- file.path(tables_path, paste0("Cluster_", cluster_id, "_Top_Genes.csv"))
  write.csv(final_table, filename, row.names = FALSE)
}

export_top_genes("C1_Unique", "1")
export_top_genes("C2_Unique", "2")
export_top_genes("C3_Unique", "3")

#-----------------------------------------------------------------------------#
# COMPLETE
#-----------------------------------------------------------------------------#

message("--- Finished Clustering ---")













