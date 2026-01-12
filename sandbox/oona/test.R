###############################################################################
###############################################################################
# PRO4002 Research Project                                                    #
#																	                                       		  #
# Date: Jan 11, 2026											                                    #
# Author: Oona Kintscher                                                      #  
###############################################################################
###############################################################################

#=============================================================================#
#=============================================================================#
# START
#=============================================================================#
#=============================================================================#

#-----------------------------------------------------------------------------#
# Library Installation 
#-----------------------------------------------------------------------------#

if (!require("ggsci", quietly = TRUE)) { install.packages("ggsci") }
library(ggsci)

if (!require("rstudioapi", quietly = TRUE)) { install.packages("rstudioapi") }
library(rstudioapi)

if (!require("openxlsx", quietly = TRUE)) { install.packages("openxlsx") }
library(openxlsx)

if (!require("readxl", quietly = TRUE)) { install.packages("readxl") }
library(readxl)

if (!require("BiocManager", quietly = TRUE)) { install.packages("BiocManager") }
BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)

if (!require("biomaRt", quietly = TRUE)) { install.packages("biomaRt") }
library(biomaRt)

if (!require("pcaMethods", quietly = TRUE)) { install.packages("pcaMethods") }
library(pcaMethods)

if (!require("ggplot2", quietly = TRUE)) { install.packages("ggplot2") }
library(ggplot2)

if (!require("dplyr", quietly = TRUE)) { install.packages("dplyr") }
library(dplyr)

if (!require("sva", quietly = TRUE)) { install.packages("sva") }
library(sva)

#-----------------------------------------------------------------------------#
# SETUP
#-----------------------------------------------------------------------------#

# Set Working Directory to the Location of this Script
setwd("C:/Users/oonak/Desktop/PRO4002/") 

# Verify Working Directory
getwd()

# Setup Seed 
set.seed(1262118388)

#-----------------------------------------------------------------------------#
# DATA IMPORT
#-----------------------------------------------------------------------------#

# FILE DATA IMPORT 
# Import Various Data Files (txt & xlsx) - specify NA values string, sheet, and header attributes
exonLengthsData <- read.delim("data/MAGNET_exonLengths.txt", header=TRUE, row.names = 1, as.is = T)
geneExpressionData <- read.delim("data/MAGNET_GeneExpressionData_CPM_19112020.txt", header=TRUE, row.names = 1, as.is = T)
sampleDescriptionData <- read_excel("data/MAGNET_SampleData_18112022_WithDescriptions.xlsx", sheet=2, na="NA")
sampleDataRaw <- read_excel("data/MAGNET_SampleData_18112022_WithDescriptions.xlsx", sheet=1, na="NA")

# ENSEMBLE DATA IMPORT 
# Get Gene Information from ENSEMBL
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
searchAttributes(mart = ensembl)

# Extract Gene Information of Genes in the Expression Data
geneListInfo <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description", "chromosome_name", "gene_biotype", "name_1006", "definition_1006", 
                 "namespace_1003"),
  filters = "ensembl_gene_id",,
  values = rownames(geneExpressionData),
  mart = ensembl
)

#-----------------------------------------------------------------------------#
# DATA PREPROCESSING 
#-----------------------------------------------------------------------------#

# Remove Batch Effect 
# mod = model.matrix(~as.factor(etiology), data=sampleDataRaw)
# mod0 = model.matrix(~1,data=sampleDataRaw)
# n.sv = num.sv(geneExpressionData,mod,method="leek")
# svobj = sva(geneExpressionData,mod,mod0,n.sv=n.sv)

# Get DCM Patients 
DCM_idx <- sampleDataRaw$etiology == "DCM"
geneExpressionData.DCM <- geneExpressionData[,DCM_idx]
sampleDataRaw.DCM <- sampleDataRaw[DCM_idx,]

# Get Category-Specific Genes 
immune.go_id <- "GO:0002376"
geneListInfo.immuneIDs <- getBM(
  attributes = c("ensembl_gene_id"),
  filters = "go",,
  values = immune.go_id,
  mart = ensembl
)

# Function to Convert logCPM Values into FPMK Values
cpm2fpkm <- function(x, y) {
  .t <- 2^(x) * 1E3 / y[, 1]
}

# Function to Center Around the Median 
center_by_median <- function(x) {
  x <- as.matrix(x)
  sweep(x, 1, matrixStats::rowMedians(x, na.rm = TRUE))
}

#-----------------------------------------------------------------------------#
# CLUSTERING: COMPLETE
#-----------------------------------------------------------------------------#

# Convert Gene Expression Data from logCPM Values to FPMK Values
geneExpressionDataFPMK.DCM <- cpm2fpkm(geneExpressionData.DCM, exonLengthsData)

# Using only the top 5000 most variable genes
geneExpressionData.DCM.mads <- apply(geneExpressionDataFPMK.DCM,1,mad)
geneExpressionDataFPMK.DCM.selected <- as.matrix(geneExpressionDataFPMK.DCM[rev(order(geneExpressionData.DCM.mads))[1:5000],])

# Centering Gene Expression to Gene Median
geneExpressionDataFPMK.DCM.centered <- sweep(geneExpressionDataFPMK.DCM.selected, 1, matrixStats::rowMedians(geneExpressionDataFPMK.DCM.selected, na.rm = TRUE))

# Clustering
title = file.path(dirname(getActiveDocumentContext()$path),"output","ConsensusClusterPlus_Complete")
dir.create(title, recursive = TRUE, showWarnings = FALSE)
results = ConsensusClusterPlus(geneExpressionDataFPMK.DCM.centered, maxK=6, reps=10, pItem=0.8,
                               pFeature=1, title=title, clusterAlg="hc", distance="pearson",
                               seed=1262118388.71279, plot="png")
results$consensusClass

pcaGeneExpressionDataFPMK.DCM.centered <- pca(t(geneExpressionDataFPMK.DCM.centered), nPcs = 10)
summary(pcaGeneExpressionDataFPMK.DCM.centered)

pcaSampleInfoDF.DCM <- cbind(scores(pcaGeneExpressionDataFPMK.DCM.centered), sampleDataRaw.DCM)
pcaSampleInfoDF.DCM$cluster <- as.factor(results[[2]]$consensusClass)

# Plot PCA Results (PC1 vs PC2)
scatterVariablePlotPCA <- ggplot(pcaSampleInfoDF.DCM, aes(PC1, PC2, color=cluster)) +
  geom_point(alpha=0.8) +
  xlab(paste("PC1 (", pcaGeneExpressionData@R2[1] * 100, "% of the variance)")) +
  ylab(paste("PC2 (", pcaGeneExpressionData@R2[2] * 100, "% of the variance)")) +
  llabs(title="PC1 vs. PC2", color="Cluster Nr.") +
  theme_minimal()
scatterVariablePlotPCA

#-----------------------------------------------------------------------------#
# CLUSTERING: GO-ID-SPECIFIC GENES
#-----------------------------------------------------------------------------#

immune_genes <- geneListInfo.immuneIDs$ensembl_gene_id 
geneExpressionDataFPMK.DCM.immune <- geneExpressionDataFPMK.DCM[rownames(geneExpressionDataFPMK.DCM) %in% immune_genes,]

# Centering Gene Expression to Gene Median
rownames(geneExpressionDataFPMK.DCM.immune) <- geneExpressionDataFPMK.DCM.immune[,1]
geneExpressionDataFPMK.DCM.immune <- geneExpressionDataFPMK.DCM.immune[,-1]
geneExpressionDataFPMK.DCM.immune.centered <- sweep(geneExpressionDataFPMK.DCM.immune, 1, 
                                                    apply(geneExpressionDataFPMK.DCM.immune, 1, median, na.rm=T))
geneExpressionDataFPMK.DCM.immune.centered <- as.matrix(geneExpressionDataFPMK.DCM.immune.centered )

# Clustering by Category (Immune System)
title = file.path(dirname(getActiveDocumentContext()$path),"output","ConsensusClusterPlus_Immune")
dir.create(title, recursive = TRUE, showWarnings = FALSE)
results.immune = ConsensusClusterPlus(geneExpressionDataFPMK.DCM.immune.centered, maxK=6, reps=10, pItem=0.8,
                                      pFeature=1, title=title, clusterAlg="hc", distance="pearson",
                                      seed=1262118388.71279, plot="png")


geneExpressionData.immune <- geneExpressionData[rownames(geneExpressionData) %in% immune_genes]
pcaGeneExpressionData.immune <- pca(t(geneExpressionData.immune), nPcs = 10)
summary(pcaGeneExpressionData.immune)
pcaSampleInfoDF.immune <- cbind(scores(pcaGeneExpressionData.immune), sampleDataRaw)
pcaSampleInfoDF.immune.DCM <- subset(pcaSampleInfoDF, pcaSampleInfoDF$etiology == "DCM")
pcaSampleInfoDF.immune.DCM <- cbind(pcaSampleInfoDF.DCM, results.immune[[2]]$consensusClass)

# Plot PCA Results (PC1 vs PC2)
scatterVariablePlotPCA.immune <- ggplot(pcaSampleInfoDF.immune.DCM, aes(PC1, PC2, color=as.factor(`results[[2]]$consensusClass`))) +
  geom_point(alpha=0.8) +
  xlab(paste("PC1 (", pcaGeneExpressionData@R2[1] * 100, "% of the variance)")) +
  ylab(paste("PC2 (", pcaGeneExpressionData@R2[2] * 100, "% of the variance)")) +
  labs(title="PC1 vs. PC2", color="Cluster Nr.") +
  theme_minimal()
scatterVariablePlotPCA.immune
