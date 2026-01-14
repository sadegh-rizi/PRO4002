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

geneExpressionData.CPM.meanFiltered.DCM <- readRDS(file.path(cache_path, "geneExpressionData_DCM_meanFiltered.RDS"))
sampleData.DCM <- readRDS(file.path(cache_path, "sampleData_DCM.RDS"))
geneListInfo <- readRDS(file.path(cache_path, "geneListInfo.RDS"))

#-----------------------------------------------------------------------------#
# PARAMETER SETUP
#-----------------------------------------------------------------------------#

seed <- 123
num_genes <- 1000    # Number of Genes to Use for Clustering
distance_measure <- "pearson"
clustering_algorithm <- "pam"

#-----------------------------------------------------------------------------#
# FEATURE SELECTION - Variance-Based Filtering (Highly Variable Genes)
#-----------------------------------------------------------------------------#

message("(II) Feature Selection ")

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

sampleData.annotated$subtype <- factor(sampleData.annotated$subtype)
sampleData.annotated$gender  <- factor(sampleData.annotated$gender)
sampleData.annotated$race  <- factor(sampleData.annotated$race)
sampleData.annotated$lvef    <- as.numeric(as.character(sampleData.annotated$lvef))

rownames(sampleData.annotated) <- sampleData.DCM$sample_name

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

pheatmap(
  geneExpressionData.variance.selected.centered,
  color = colorRampPalette(c(npg_colors[[4]], "white", npg_colors[[1]]))(100),
  breaks = seq(-2, 2, length.out = 101), # Cap scale at z-score +/- 2
  show_rownames = FALSE,      # Too many genes to show names
  show_colnames = FALSE,      # Hiding sample names for cleaner look
  annotation_col = sampleData.annotated,
  annotation_colors = colors.annotated,
  clustering_method = "ward.D2",
  cutree_cols = k_choice,
  na_col = "grey30",
  main = paste0("DCM Molecular Subtypes (Top ", num_genes, " Variable Genes)"),
  filename = file.path(cluster_path, "subtypeHeatmap_DCM.jpg"),
  width = 10, height = 8
)

#-----------------------------------------------------------------------------#
# CLUSTERING STABILITY 
#-----------------------------------------------------------------------------#

message("(IV) Clustering - Stability Analysis")

geneData.distance <- as.dist(1 - cor(geneExpressionData.variance.selected.centered, method = "pearson"))
silhouetteData <- silhouette(as.numeric(sampleData.annotated$subtype), geneData.distance) 
silhouetteDataPlot <- fviz_silhouette(silhouetteData,  
                            palette = npg_colors,  
                            ggtheme = my_style, 
                            main = paste0("Cluster Stability (k=", k_choice, "): Silhouette Plot")) 
ggsave(file.path(cluster_path, "silhouetteDataPlot.jpg"), silhouetteDataPlot, width = 8, height = 6) 

silhouetteData.width <- summary(silhouetteData)$avg.width 
message(paste("Average Silhouette Width:", round(silhouetteData.width, 3))) 

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
# COMPLETE
#-----------------------------------------------------------------------------#

message("--- Finished Clustering ---")
