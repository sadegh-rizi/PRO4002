###############################################################################
###############################################################################
# MSB1005_FinalAssignment_2025.R                                              #
#																	                                       		  #
# Date: Dec 12, 2025											                                    #
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

if (!require("openxlsx", quietly = TRUE)) { install.packages("openxlsx") }
library(openxlsx)

if (!require("rstudioapi", quietly = TRUE)) { install.packages("rstudioapi") }
library(rstudioapi)

if (!require("readxl", quietly = TRUE)) { install.packages("readxl") }
library(readxl)

if (!require("Hmisc", quietly = TRUE)) { install.packages("Hmisc") }
library(Hmisc)

if (!require("table1", quietly = TRUE)) { install.packages("table1") }
library(table1)

if (!require("tidyr", quietly = TRUE)) { install.packages("tidyr") }
library(tidyr)

if (!require("dplyr", quietly = TRUE)) { install.packages("dplyr") }
library(dplyr)

if (!require("ggplot2", quietly = TRUE)) { install.packages("ggplot2") }
library(ggplot2)

if (!require("pcaMethods", quietly = TRUE)) { install.packages("pcaMethods") }
library(pcaMethods)

if (!require("ggpubr", quietly = TRUE)) { install.packages("ggpubr") }
library(ggpubr)

if (!require("gridExtra", quietly = TRUE)) { install.packages("gridExtra") }
library(gridExtra)

if (!require("grid", quietly = TRUE)) { install.packages("grid") }
library(grid)

if (!require("pheatmap", quietly = TRUE)) { install.packages("pheatmap") }
library(pheatmap)

if (!require("tidyverse", quietly = TRUE)) { install.packages("tidyverse") }
library(tidyverse)

if (!require("limma", quietly = TRUE)) { install.packages("limma") }
library(limma)

if (!require("EnhancedVolcano", quietly = TRUE)) { install.packages("EnhancedVolcano") }
library(EnhancedVolcano)

if (!require("biomaRt", quietly = TRUE)) { install.packages("biomaRt") }
library(biomaRt)

if (!require("pheatmap", quietly = TRUE)) { install.packages("pheatmap") }
library(pheatmap)

if (!require("VennDiagram", quietly = TRUE)) { install.packages("VennDiagram") }
library(VennDiagram) 

if (!require("readr", quietly = TRUE)) { install.packages("readr") }
library(readr)

#-----------------------------------------------------------------------------#
# Design & Export Setup
#-----------------------------------------------------------------------------#

# Prepare Color Values to Match Nature Color Scheme
npgColors <- pal_npg("nrc", alpha = 1)(10)
npgAdditionalColors <- colorRampPalette(npgColors)(20)
continuousNPGColors <- colorRampPalette(c(npgColors[2], "white"))(100)

# Create a Excel Workbook to Add to Throughout The Script 
workbook <- createWorkbook()

# Create Folder to Store All Images & Graphs 
graphs_folder <- "graphs"
if(!dir.exists(graphs_folder)) {
  dir.create(graphs_folder)
}

#=============================================================================#
#=============================================================================#
# Part 1: DATA IMPORT                                                         #
#=============================================================================#
#=============================================================================#

# Set Working Directory to the Location of this Script
setwd(dirname(getActiveDocumentContext()$path)) 

# Verify Working Directory
getwd()

#-----------------------------------------------------------------------------#
# The Working Directory Should Have the Following Setup
#
# Working Directory (wd)
# ├─ MAGNET_exonLengths.txt
# ├─ MAGNET_GeneExpressionData_CPM_19112020.txt
# ├─ MAGNET_SampleData_18112022.csv
# ├─ MAGNET_SampleData_18112022_WithDescriptions.xlsx
# └─ MSB1005_FinalAssignment_2025.R 
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Part 1a: Import all the data files                                          
#-----------------------------------------------------------------------------#

# Import Various Data Files (txt & xlsx) - specify NA values string, sheet, and header attributes
exonLengthsData <- read.delim("MAGNET_exonLengths.txt", header=TRUE, row.names = 1, as.is = T)
geneExpressionData <- read.delim("MAGNET_GeneExpressionData_CPM_19112020.txt", header=TRUE, row.names = 1, as.is = T)
sampleDescriptionData <- read_excel("MAGNET_SampleData_18112022_WithDescriptions.xlsx", sheet=2, na="NA")
sampleDataRaw <- read_excel("MAGNET_SampleData_18112022_WithDescriptions.xlsx", sheet=1, na="NA")

# View Imported Data As Simple Quality Check 
View(exonLengthsData)
View(geneExpressionData)
View(sampleDescriptionData)
View(sampleDataRaw)

# Inspect Dimensions, Storage and NAs in Preparation of Preprocessing 
Hmisc::contents(exonLengthsData)
Hmisc::contents(geneExpressionData)
Hmisc::contents(sampleDataRaw)

# Check for Data Conformance - do Patient-IDs and Gene-IDs match across tables
all(rownames(geneExpressionData) == rownames(exonLengthsData))
all(colnames(geneExpressionData) == sampleDataRaw$sample_name)

# Data Preprocessing - factoring of specific columns for better understanding
sampleData <- sampleDataRaw
sampleData$race <- factor(sampleData$race, levels=c("AA", "Caucasian"), 
                          labels = c("African American", "Caucasian"))
sampleData$afib <- factor(sampleData$afib, levels = c("Yes", "No"), 
                          labels = c("Presence", "Absent"))
sampleData$VTVF <- factor(sampleData$VTVF, levels = c("Yes", "No"), 
                          labels = c("Presence", "Absent"))
sampleData$Diabetes <- factor(sampleData$Diabetes, levels = c("Yes", "No"), 
                              labels = c("Presence", "Absent"))
sampleData$Hypertension <- factor(sampleData$Hypertension, levels = c("Yes", "No"), 
                                  labels = c("Presence", "Absent"))
sampleData$etiology <- factor(sampleData$etiology, levels = c("NF", "DCM", "HCM", "PPCM"), 
                              labels = c("Control", "DCM", "HCM", "PPCM"))

#-----------------------------------------------------------------------------#
# Part 1b: Participant Characteristics                                        #
#-----------------------------------------------------------------------------#

# Setup Table
demoTableLabels <- list(variables = list(age="Age (years)",
                                         weight="Weight (kg)",
                                         height="Height (cm)",
                                         gender="Gender",
                                         race="Race",
                                         hw="Heart Weight (g)",
                                         lv_mass="LVᵃ Mass (g)",
                                         afib="Atrial Fibrillation",
                                         VTVF="VTVFᵇ",
                                         Diabetes="Diabetes",
                                         Hypertension="Hypertension",
                                         LVEF="LVᵃ Ejection Fraction" ),
                        groups=list("DCM", "HCM", "PPCM", "NF"))
demoTableFootnote  <- c("ᵃ Left Ventricle", "ᵇ Ventricular Tachycardia/Ventricular Fibrillation")
demoTableStrata <- c(list(Total=sampleData), split(sampleData, sampleData$etiology))
 
demographicTable <- table1(demoTableStrata, demoTableLabels, 
                           footnote=demoTableFootnote, topclass="Rtable1-zebra")
demographicTable

# Export Table
addWorksheet(workbook, "Participant Characteristics")
writeData(workbook, "Participant Characteristics", demographicTable)

#=============================================================================#
#=============================================================================#
# Part 2: DIAGNOSTIC PLOTS                                                    #
#=============================================================================#
#=============================================================================#

#-----------------------------------------------------------------------------#
# Part 2a: Data Distribution                                                  #
#-----------------------------------------------------------------------------#

# Extending Gene Expression Data Table - one row represents a unique combination of gene and patient
extendedGeneExpressionData <- geneExpressionData %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols=-gene, names_to="patient", values_to="gene_expression_level") %>%
  left_join(sampleData, by=c('patient'='sample_name'))

# Data Distribution (Box Plots)
boxPlotDataDistribution.etiology <- ggplot(data = extendedGeneExpressionData, 
                                           aes(x = patient, 
                                               y = gene_expression_level, 
                                               color = etiology)) +
  geom_boxplot() + 
  labs(title = "Patients vs. Gene Expression Level", x = "Patient ID",
       y = "Gene Expression Level (logCPM)", color = 'Etiology') +
  scale_color_npg() +
  scale_fill_npg() +
  facet_wrap(~etiology, scales="free_x") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust=1))
boxPlotDataDistribution.etiology

# Data Distribution (Density Plots per Etiology)
densityPlotDataDistribution.etiology <- ggplot(data = extendedGeneExpressionData, 
                                               aes(x = gene_expression_level, 
                                                   fill = etiology, color = etiology)) +
  geom_density(alpha=0.7) +
  labs(title="Gene Expression Level Density Plot", x="Gene Expression Level (logCPM)",
       y="Density", color="Etiology", fill="Etiology") +
  scale_color_npg() +
  scale_fill_npg() +
  facet_wrap(~etiology, scales="free_x")
densityPlotDataDistribution.etiology

# Data Distribution (Density Plots) - allows to investigate the overlap of different etiologies 
densityPlotDataDistribution.overlap <- ggplot(data = extendedGeneExpressionData, 
                                               aes(x = gene_expression_level, 
                                                   fill = etiology, color = etiology)) +
  geom_density(alpha=0.3) +
  labs(title  ="Gene Expression Level Density Plot", x="Gene Expression Level (logCPM)",
       y = "Density", color = "Etiology", fill = "Etiology") +
  scale_color_npg() +
  scale_fill_npg()
densityPlotDataDistribution.overlap

#-----------------------------------------------------------------------------#
# Part 2b: PCA                                                                #
#-----------------------------------------------------------------------------#
#                                                                             #
# PCA (Principle Component Analysis) is a technique that reduces the          #
# dimensionality of a given dataset. It can be used to check for population   #
# stratification by examining the PCA plot for clusters and other patterns.   #
#                                                                             #
#-----------------------------------------------------------------------------#

 

# Performing a PCA on the Gene Expression Data - specified to compute 10 principle components
pcaGeneExpressionData <- pca(t(geneExpressionData), nPcs = 10)
summary(pcaGeneExpressionData)

# Create Data Frame with Summary of PCA 
pcaSummaryDF <- data.frame(PC = paste0("PC",factor(1:pcaGeneExpressionData@nPcs)), 
                                     var = pcaGeneExpressionData@R2, var_cum=pcaGeneExpressionData@R2cum)
pcaScoresDF <- scores(pcaGeneExpressionData)
pcaLoadingsDF <- loadings(pcaGeneExpressionData)

# Add PCA Data Frame To Excel File 
addWorksheet(workbook, "PCA Summary")
writeData(workbook, "PCA Summary", pcaSummaryDF)
addWorksheet(workbook, "PCA Scores")
writeData(workbook, "PCA Scores", pcaScoresDF, rowNames = TRUE)
addWorksheet(workbook, "PCA Loadings")
writeData(workbook, "PCA Loadings", pcaLoadingsDF, rowNames = TRUE)

# Plot Principle Components vs Variance in Data Explained by Component & Sum of Variance
pcaVarianceExplainedPlot <- ggplot(pcaSummaryDF, aes(x = PC)) +
  geom_col(aes(y = var*100), fill = npgColors[2], alpha = 0.7, color = npgColors[2]) +   
  geom_line(aes(y = var_cum*100, group = 1)) + 
  geom_point(aes(y = var_cum*100)) +
  labs(title="Principle Component vs. Variance Explained",
       x = "Principal Component", y = "Percentage of Variance Explained (%)") +
  theme_minimal()
pcaVarianceExplainedPlot

# Combining Principle Components with the Corresponding Sample Data
pcaSampleInfoDF <- cbind(scores(pcaGeneExpressionData), sampleData)

# Plot PCA Results (PC1 vs PC2 & PC3 vs PC4)
scatterVariablePlotPCA.PC1PC2 <- ggplot(pcaSampleInfoDF, aes(PC1, PC2, color=etiology)) +
  geom_point(alpha=0.8) +
  xlab(paste("PC1 (", pcaGeneExpressionData@R2[1] * 100, "% of the variance)")) +
  ylab(paste("PC2 (", pcaGeneExpressionData@R2[2] * 100, "% of the variance)")) +
  labs(title="PC1 vs. PC2", color='Etiology') +
  scale_color_npg() +
  scale_fill_npg() +
  theme_minimal()
scatterVariablePlotPCA.PC1PC2

scatterVariablePlotPCA.PC3PC4 <- ggplot(pcaSampleInfoDF, aes(PC3, PC4, color = etiology)) +
  geom_point(alpha=0.8) +
  xlab(paste("PC3 (", pcaGeneExpressionData@R2[3] * 100, "% of the variance)")) +
  ylab(paste("PC4 (", pcaGeneExpressionData@R2[4] * 100, "% of the variance)")) +
  labs(title="PC3 vs. PC4", color='Etiology') +
  scale_color_npg() +
  scale_fill_npg() +
  theme_minimal()
scatterVariablePlotPCA.PC3PC4

# Combine Both PCA Plots into Single Image
sharedLegendPlotPCA <- get_legend(scatterVariablePlotPCA.PC1PC2)
scatterVariablePlotPCA.PC1PC2 <- scatterVariablePlotPCA.PC1PC2 + theme(legend.position = "none")
scatterVariablePlotPCA.PC3PC4 <- scatterVariablePlotPCA.PC3PC4 + theme(legend.position = "none")
scatterVariablePlotPCA.combined <- grid.arrange(arrangeGrob(scatterVariablePlotPCA.PC1PC2, 
                                                            scatterVariablePlotPCA.PC3PC4, 
                                                            nrow = 1), 
                                                sharedLegendPlotPCA,
                                                ncol = 2, widths = c(8,1),
                                                top = "Principle Component Analysis")

# Plot PCA Results (PC1 vs PC2) vs Gender, Age, Race & Library.Pool
scatterVariablePlotPCA.PC1PC2.gender <- ggplot(pcaSampleInfoDF, aes(PC1, PC2, color = gender)) +
  geom_point(alpha=0.8) +
  xlab(paste("PC1 (", pcaGeneExpressionData@R2[1] * 100, "% of the variance)")) +
  ylab(paste("PC2 (", pcaGeneExpressionData@R2[2] * 100, "% of the variance)")) +
  labs(color='Gender') +
  scale_color_manual(values = npgAdditionalColors) +
  theme_minimal()
scatterVariablePlotPCA.PC1PC2.gender

scatterVariablePlotPCA.PC1PC2.age <- ggplot(pcaSampleInfoDF, aes(PC1, PC2, color = age)) +
  geom_point(alpha=0.8) +
  xlab(paste("PC1 (", pcaGeneExpressionData@R2[1] * 100, "% of the variance)")) +
  ylab(paste("PC2 (", pcaGeneExpressionData@R2[2] * 100, "% of the variance)")) +
  labs(color='Age') +
  scale_color_gradientn(colors = continuousNPGColors) +
  theme_minimal()
scatterVariablePlotPCA.PC1PC2.age

scatterVariablePlotPCA.PC1PC2.race <- ggplot(pcaSampleInfoDF, aes(PC1, PC2, color=race)) +
  geom_point(alpha=0.8) +
  xlab(paste("PC1 (", pcaGeneExpressionData@R2[1] * 100, "% of the variance)")) +
  ylab(paste("PC2 (", pcaGeneExpressionData@R2[2] * 100, "% of the variance)")) +
  labs(color='Race') +
  scale_color_manual(values = npgAdditionalColors) +
  theme_minimal()
scatterVariablePlotPCA.PC1PC2.race

scatterVariablePlotPCA.PC1PC2.libraryPool <- ggplot(pcaSampleInfoDF, aes(PC1, PC2, color=Library.Pool)) +
  geom_point(alpha=0.8) +
  xlab(paste("PC1 (", pcaGeneExpressionData@R2[1] * 100, "% of the variance)")) +
  ylab(paste("PC2 (", pcaGeneExpressionData@R2[2] * 100, "% of the variance)")) +
  labs(color='Library.Pool') +
  scale_color_manual(values = npgAdditionalColors) +
  stat_ellipse(aes(color = Library.Pool)) +
  theme_minimal()
scatterVariablePlotPCA.PC1PC2.libraryPool

# Combine Four PC1vsPC2 Plots into Single Image
scatterVariablePlotPCA.variables <- grid.arrange(scatterVariablePlotPCA.PC1PC2.gender, 
                                                 scatterVariablePlotPCA.PC1PC2.age,
                                                 scatterVariablePlotPCA.PC1PC2.race, 
                                                 scatterVariablePlotPCA.PC1PC2.libraryPool,
                                                 top = "Principle Component Analysis (PC1vsPC2)")

# Plot Interpretation --------------------------------------------------------#
#
# scatterVariablePlotPCA.combined
# The PCA (PC1 vs. PC2) shows clustering based on the etiology. 
# In this case, we expect to see differences between these groups, 
# therefore this clustering does not indicate problematic population stratification.  
#
# scatterVariablePlotPCA.variables
# These plots were used to assess whether clustering occurse based on gender, age, 
# race and Library.Pool. The plots show no clustering with respect to gender, age 
# or race. However, some clustering can be observed based on Library.Pool (ellipses), 
# which should be corrected for in DGEA. 
#
# ----------------------------------------------------------------------------#

# Function to Plot PCA of a Specific Etiology Type 
pcaPlot <- function(type) { 
  pcaGeneExpressionData <- pca(t(geneExpressionData[,sampleData$etiology==type]))
  summary(pcaGeneExpressionData)
  pcaGeneExpressionDataDf <- cbind(scores(pcaGeneExpressionData), 
                                   sampleData[sampleData$etiology==type,])
  scatterVariablePlotPCA <- ggplot(pcaGeneExpressionDataDf, aes(PC1, PC2)) +
    geom_point(alpha=0.8, color = npgColors[2]) +
    xlab(paste("PC1 (", pcaGeneExpressionData@R2[1] * 100, "% of the variance)")) +
    ylab(paste("PC2 (", pcaGeneExpressionData@R2[2] * 100, "% of the variance)")) +
    labs(title=type) +
    scale_color_manual(values = npgAdditionalColors) +
    theme_minimal()
  return(scatterVariablePlotPCA)
}

pcaPlot.DCM <- pcaPlot("DCM")
pcaPlot.HCM <- pcaPlot("HCM") 
pcaPlot.PPCM <- pcaPlot("PPCM") 
pcaPlot.Control <- pcaPlot("Control") 

# Combine PCA Plots of each Etiology Together
scatterVariablePlotPCA.etiology <- grid.arrange(pcaPlot.Control, pcaPlot.DCM, 
                                                pcaPlot.HCM, pcaPlot.PPCM,
                                                nrow=2, ncol = 2,
                                                top = "Principle Component Analysis")

# Plot Interpretation -------------------------------------------------------#
# 
# scatterVariablePlotPCA.etiology
# A separate PCA was performed for each group to check for population 
# stratification within the groups. The plots show no clustering which tells 
# us there are no hidden subgroups or co-founding variables. 
#
# ----------------------------------------------------------------------------#

# ----------------------------------------------------------------------------#
# Part 2c: Plot Export
# ----------------------------------------------------------------------------#

ggsave(file.path(graphs_folder, "data_distribution_box_plot.jpg"), 
       plot = boxPlotDataDistribution.etiology, width = 8, height = 6, dpi = 300)
ggsave(file.path(graphs_folder, "data_distribution_density_plot_etiology.jpg"), 
       plot = densityPlotDataDistribution.etiology, width = 8, height = 6, dpi = 300)
ggsave(file.path(graphs_folder, "data_distribution_density_plot_overlap.jpg"), 
       plot = densityPlotDataDistribution.overlap, width = 8, height = 6, dpi = 300)

ggsave(file.path(graphs_folder, "pca_plot_components.jpg"), 
       plot = scatterVariablePlotPCA.combined, width = 12, height = 6, dpi = 300)
ggsave(file.path(graphs_folder, "pca_plot_variables.jpg"), 
       plot = scatterVariablePlotPCA.variables, width = 12, height = 12, dpi = 300)
ggsave(file.path(graphs_folder, "pca_plot_etiology.jpg"), 
       plot = scatterVariablePlotPCA.etiology, width = 12, height = 12, dpi = 300)

#=============================================================================#
#=============================================================================#
# Part 3: STATISTICAL ANALYSIS
#=============================================================================#
#=============================================================================#

# Setting P-Value and Log2 Fold-Change Cutoff
log2FC.cutoff <- 0.585 
pval.cutoff <- 0.05 

# Check for Covariates and Principle Componenet Correlation/Association
pcNames <- colnames(pcaSampleInfoDF[1:4])
variateNames <- colnames(sampleData[3:17])
pcaVariateDF <- matrix(NA, nrow = length(pcNames), ncol = length(variateNames),
                       dimnames = list(pcNames, variateNames))
pcaVariateDF <- as.data.frame(pcaVariateDF)

for (variate in variateNames) {
  for (pc in pcNames) {
    if (is.numeric(pcaSampleInfoDF[[variate]])){
      pcValues <- pcaSampleInfoDF[[pc]]
      variateValues <- pcaSampleInfoDF[[variate]]
      fit <- cor.test(pcValues, variateValues)
      p_value = fit$p.value
    } else {
      formula <- as.formula(paste(pc, "~", variate))
      fit <- aov(formula, data = pcaSampleInfoDF)
      p_value <- summary(fit)[[1]][["Pr(>F)"]][1]
    }
    pcaVariateDF[pc,variate] = p_value
  }
}

pcaVariateDF <- as.matrix(sapply(pcaVariateDF, as.numeric))

pcaVariateDFSigLevels <- ifelse(pcaVariateDF < 0.001, "***",
                                 ifelse(pcaVariateDF < 0.01, "**",
                                        ifelse(pcaVariateDF < 0.05, "*", "")))

heatmapPCVariateSignifcance <- pheatmap(t(pcaVariateDF),
                                        color = continuousNPGColors,
                                        main = "PC & Covariate Statistical Significance",
                                        display_numbers = t(pcaVariateDFSigLevels),
                                        labels_col = pcNames,
                                        number_color = "black",
                                        border_color = "black",
                                        cluster_rows = FALSE,
                                        cluster_cols = FALSE)

# Plot Interpretation --------------------------------------------------------#
#
# heatmapPCVariateSignifcance
# This heatmap allows us to investigate the correlation between each variable 
# and principle component. It shows that etiology is statistically significantly
# correlated with the first three PCs. PC4 shows only correlation with height 
# and Library.Pool, while all other correlation are statistically insignificant.
#
# ----------------------------------------------------------------------------#

# DGEA Formula Design --------------------------------------------------------#
#
# The PCA did not reveal significant population stratification. However, to correct 
# for some sources of unwanted variance in the dataset, certain co-variates were 
# added to the design formula. This includes standard demographic variables such 
# as gender, age, and race. The 4th principal component was also included to account 
# for other sources of variation. This ensures that the clustering observed in 
# Library.Pool is accounted for, as its effect is represented in the 4th component, 
# as shown in the heatmap above.
#
# ----------------------------------------------------------------------------#

# Perform DGEA Function
performDGEA <- function(gene_expression, sample_data, contrast,
                        design_formula = ~ 0 + etiology + gender + age + race + PC4) {
  design <- model.matrix(design_formula, data = sample_data)
  fit <- lmFit(gene_expression, design)
  cont.matrix <- makeContrasts(contrasts = contrast, levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  ebFit <- eBayes(fit2, trend = TRUE)
  resultDGEA <- topTable(ebFit, coef = 1, number = nrow(gene_expression))
  return(resultDGEA)
}

# Plot DGEA using Enhanced Volcano Plot Function 
plotEnhancedVolcano <- function(results_dgea, title, pval_cutoff = pval.cutoff, 
                                log2Fc_cutoff = log2FC.cutoff) {
  plot <- EnhancedVolcano(results_dgea,  title = title, labSize = 3, 
                          x = 'logFC', y = 'P.Value', 
                          lab = row.names(results_dgea),
                          pCutoff = pval_cutoff, 
                          FCcutoff = log2Fc_cutoff)
  return(plot)
}

# DCM vs Control
resultDGEA.DCM <- performDGEA(geneExpressionData, pcaSampleInfoDF,
                             contrast = "etiologyDCM - etiologyControl")

# HCM vs Control
resultDGEA.HCM <- performDGEA(geneExpressionData, pcaSampleInfoDF,
                             contrast = "etiologyHCM - etiologyControl")

# PPCM vs Control
resultDGEA.PPCM <- performDGEA(geneExpressionData, pcaSampleInfoDF,
                              contrast = "etiologyPPCM - etiologyControl")

# Create Enhanced Volcano Plots for each DGE 
enhancedVolcanoPlot.DCM <- plotEnhancedVolcano(resultDGEA.DCM, "DCM vs. Control")
sharedLegendEnhancedVolcanoPlot <- get_legend(enhancedVolcanoPlot.DCM)
enhancedVolcanoPlot.DCM <- enhancedVolcanoPlot.DCM + theme(legend.position = "none")
enhancedVolcanoPlot.HCM <- plotEnhancedVolcano(resultDGEA.HCM, "HCM vs. Control") + 
  theme(legend.position = "none")
enhancedVolcanoPlot.PPCM <- plotEnhancedVolcano(resultDGEA.PPCM, "PPCM vs. Control") + 
  theme(legend.position = "none")

# Combine the Three Enhanced Volcano Plots to Single Image
enhancedVolcanoPlot.combined <- grid.arrange(sharedLegendEnhancedVolcanoPlot,
                                             arrangeGrob(enhancedVolcanoPlot.DCM, 
                                                         enhancedVolcanoPlot.HCM, 
                                                         enhancedVolcanoPlot.PPCM,
                                                         nrow=1),
                                             nrow = 2, heights = c(1,8),
                                             top = "Enhanced Volcano Plots")
enhancedVolcanoPlot.combined

# Summary Statistics of Differential Expressed Genes for Each Group
deg.DCM <- resultDGEA.DCM[resultDGEA.DCM$P.Value < pval.cutoff & 
                            abs(resultDGEA.DCM$logFC) > log2FC.cutoff, 1:5]
deg.HCM <- resultDGEA.HCM[resultDGEA.HCM$P.Value < pval.cutoff & 
                            abs(resultDGEA.HCM$logFC) > log2FC.cutoff, 1:5]
deg.PPCM <- resultDGEA.PPCM[resultDGEA.PPCM$P.Value < pval.cutoff 
                            & abs(resultDGEA.PPCM$logFC) > log2FC.cutoff, 1:5]

cat("Number of Differentially Expressed Genes (DEGs):",
    "\n - DCM:", nrow(deg.DCM),
    "\n - HCM:", nrow(deg.HCM),
    "\n - PPCM:", nrow(deg.PPCM), "\n")

# Venn Diagram of DGEA Comparing Each Group
venn.diagram(
  x = list(
    row.names(deg.DCM),
    row.names(deg.HCM),
    row.names(deg.PPCM)
  ),
  category.names = c("DCM", "HCM", "PPCM"),
  main = 'Differential Gene Expression Analysis',
  filename = file.path(graphs_folder, "DGEA_venn_diagram.jpg"),
  output = FALSE,
  col = npgColors[2:4],
  fill = c(
    alpha(npgColors[2], 0.3), 
    alpha(npgColors[3], 0.3),
    alpha(npgColors[4], 0.3)
  ),
  cex = 1.5
)

# Save Graphs Created
ggsave(file.path(graphs_folder, "PC_variate_signifcance_heatmap.jpg"), 
       plot = heatmapPCVariateSignifcance, width = 8, height = 6, dpi = 300)
ggsave(file.path(graphs_folder, "DGEA_enhanced_volcano_plots.jpg"),
       plot = enhancedVolcanoPlot.combined, width = 16, height = 6, dpi = 300)


#=============================================================================#
#=============================================================================#
# Part 4: ADDITIONAL GENE ANNOTATION
#=============================================================================#
#=============================================================================#

#-----------------------------------------------------------------------------#
# Part 4a: Retrieve Gene Symbols and Gene Names                               #
#-----------------------------------------------------------------------------#

# Import Human Ensembl Dataset
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Extract Gene Information of Genes in the Expression Data
geneListInfo <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description", "chromosome_name"),
  filters = "ensembl_gene_id",,
  values = rownames(geneExpressionData),
  mart = ensembl
)

#-----------------------------------------------------------------------------#
# Part 4b: Merge Annotations with the Gene Expression Data                    #
#-----------------------------------------------------------------------------#

annotatedGeneExpressionData <- merge(geneExpressionData, geneListInfo, by.x="row.names", by="ensembl_gene_id")

#=============================================================================#
#=============================================================================#
# Part 5: RELATIVE GENE EXPRESSION
#=============================================================================#
#=============================================================================#

# Function to Convert logCPM Values into FPMK Values
cpm2fpkm <- function(x, y) {
  .t <- 2^(x) * 1E3 / y[, 1]
}

# Convert Gene Expression Data from logCPM Values to FPMK Values
geneExpressionDataFPMK <- cpm2fpkm(geneExpressionData, exonLengthsData)

# Y-Chromosome Gene Expression in Female Patients ----------------------------#
#
# - extract expression data of genes on the Y chromosome of female patients
# - genes of Y chromosome should not be expressed in a female patient
# - expression level of these genes in female patients can be used to threshold noise
#
#-----------------------------------------------------------------------------#
geneExpressionDataFPMK.femaleY <- geneExpressionDataFPMK[annotatedGeneExpressionData$chromosome_name == "Y",
                                                    sampleData$gender == "Female", drop=FALSE]

# Create Numerical List of Gene Expression Data of Female Patient Y-Chromosome Genes
backgroundExpressionData.femaleY <- as.numeric(as.matrix(geneExpressionDataFPMK.femaleY))

# 95% Quantile Threshold -----------------------------------------------------#
#                                                             
# - use 95% quantile gene expression data of female patient Y-chromosome genes as 
#   a threshold for background expression levels
# - quantiles are preferred over the mean and max because they provide a better 
#   measure of the upper range background noise 
# - quantiles are widely used in RNA-seq quality control 
#
#-----------------------------------------------------------------------------#
backgroundThreshold95PExpressionData <- quantile(backgroundExpressionData.femaleY, 
                                                 probs=0.95, na.rm=TRUE)

# Mean Threshold -------------------------------------------------------------#
#                                                             
# - use mean gene expression data of female patient Y-chromosome genes as 
#   a threshold for background expression levels
#
#-----------------------------------------------------------------------------#
backgroundThresholdMeanExpressionData <- mean(backgroundExpressionData.femaleY)

# Reshaping Data Frame To Have Gene ID, Patient and Expression Level in Each Row
extendedGeneExpressionDataFPMK.femaleY <- geneExpressionDataFPMK.femaleY %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  pivot_longer(cols=-ensembl_gene_id, names_to="patient", values_to="gene_expression_level") %>%
  left_join(sampleData, by=c('patient'='sample_name'))

extendedGeneExpressionDataFPMK <- geneExpressionDataFPMK %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  pivot_longer(cols=-ensembl_gene_id, names_to="patient", values_to="gene_expression_level") %>%
  left_join(sampleData, by=c('patient'='sample_name'))

# Plot Density Expression Level of Female Patient Y-Chromosome Genes with Threshold 
noiseExpressionDensityPlot <- ggplot(data=extendedGeneExpressionDataFPMK.femaleY,
                                     aes(x=log10(gene_expression_level))) +
  geom_density(fill = npgColors[2], color = npgColors[2], alpha = 0.7) +
  geom_vline(xintercept = log10(backgroundThreshold95PExpressionData), 
             color = npgColors[1], linetype = "dashed", linewidth=0.8) +
  annotate("text", x = log10(backgroundThreshold95PExpressionData), y = 0, 
           label = "95% Threshold", color = npgColors[1], angle =90, vjust = -0.5, hjust=-0.3, size =4) +
  geom_vline(xintercept = log10(backgroundThresholdMeanExpressionData), 
             color = npgColors[4], linetype = "dashed", linewidth=0.8) +
  annotate("text", x = log10(backgroundThresholdMeanExpressionData), y = 0, 
           label = "Mean Threshold", color = npgColors[4], angle =90, vjust = -0.5, hjust=-0.3, size =4) +
  labs(title="Female - Y Chromosome Gene Expression with Threshold", x="Gene Expression Level (FPKM)", y="Density") 
noiseExpressionDensityPlot

densityPlotDataDistribution.threshold <- ggplot(data = extendedGeneExpressionDataFPMK, 
                                              aes(x = log10(gene_expression_level))) +
  geom_density(fill = npgColors[2], color = npgColors[2], alpha = 0.7) +
  geom_vline(xintercept=log10(backgroundThreshold95PExpressionData), 
             color = npgColors[1], linetype = "dashed", linewidth=0.8) +
  annotate("text", x = log10(backgroundThreshold95PExpressionData), y = 0, 
           label = "95% Threshold", color = npgColors[1], angle =90, vjust = -0.5, hjust=-0.3, size =4) +
  geom_vline(xintercept = log10(backgroundThresholdMeanExpressionData), 
             color = npgColors[4], linetype = "dashed", linewidth=0.8) +
  annotate("text", x = log10(backgroundThresholdMeanExpressionData), y = 0, 
           label = "Mean Threshold", color = npgColors[4], angle =90, vjust = -0.5, hjust=-0.3, size =4) +
  labs(title  ="Gene Expression Level Density Plot with Threshold", x="Gene Expression Level (FPMK)",y = "Density")
densityPlotDataDistribution.threshold

# Compute Mean Expression Level of Each Gene Across All Samples
meanGeneExpressionData <- rowMeans(geneExpressionDataFPMK, na.rm=TRUE)

# Above/Below Threshold ------------------------------------------------------#
#
# - check if mean expression level of each gene is above background threshold 
# - true: gene is biologically expressed 
# - false: gene expression in part of background noise
#
#-----------------------------------------------------------------------------#
aboveBackground95PThresholdGenes <- meanGeneExpressionData > backgroundThreshold95PExpressionData
aboveBackgroundMeanThresholdGenes <- meanGeneExpressionData > backgroundThresholdMeanExpressionData

cat("95% Percentile Threshold: ", backgroundThreshold95PExpressionData, "FPKM",
    "\nNumber of Genes with Expression Level Above or Below the Threshold",
    "\n - above:", sum(aboveBackground95PThresholdGenes),
    "\n - below:", sum(!aboveBackground95PThresholdGenes))

cat("Mean Threshold: ", backgroundThresholdMeanExpressionData, "FPKM",
    "\nNumber of Genes with Expression Level Above or Below the Threshold",
    "\n - above:", sum(aboveBackgroundMeanThresholdGenes),
    "\n - below:", sum(!aboveBackgroundMeanThresholdGenes))

# Combine All Information Into Single Data Frame
geneExpressionDataInfo <- data.frame(
  ensembl_gene_id = rownames(geneExpressionDataFPMK),
  mean_expression_fpkm = meanGeneExpressionData,
  mean_dcm_expression_fpkm = rowMeans(geneExpressionDataFPMK[,sampleData$etiology == "DCM",], na.rm=TRUE),
  mean_control_expression_fpkm = rowMeans(geneExpressionDataFPMK[,sampleData$etiology == "Control"], na.rm=TRUE),
  expressed_above_background_95P_threshold = aboveBackground95PThresholdGenes,
  expressed_above_background_mean_threshold = aboveBackgroundMeanThresholdGenes
)

# Save Graphs Created
ggsave(file.path(graphs_folder, "background_data_distribution_density_plot_with_threshold.jpg"), plot = noiseExpressionDensityPlot, width = 8, height = 6, dpi = 300)
ggsave(file.path(graphs_folder, "data_distribution_density_plot_with_threshold.jpg"), plot = densityPlotDataDistribution.threshold, width = 8, height = 6, dpi = 300)

#=============================================================================#
#=============================================================================#
# Part 6: EXPORT RESULTS
#=============================================================================#
#=============================================================================#

# Prepare Individual Table for Merging 
tableResultDGEA.DCM <- resultDGEA.DCM %>%
  rename_with(~ paste0("DCMvsControl_", .)) %>%
  tibble::rownames_to_column("ensembl_gene_id") 

tableResultDGEA.HCM <- resultDGEA.HCM %>%
  rename_with(~ paste0("HCMvsControl_", .)) %>%
  tibble::rownames_to_column("ensembl_gene_id") 

tableResultDGEA.PPCM <- resultDGEA.PPCM %>%
  rename_with(~ paste0("PPCMvsControl_", .)) %>%
  tibble::rownames_to_column("ensembl_gene_id") 

# Merge All Results Into Single Data Frame 
dataCollection <- data.frame(ensembl_gene_id = row.names(geneExpressionData))
tableResultDGEA.list <- list(tableResultDGEA.DCM, tableResultDGEA.HCM, tableResultDGEA.PPCM)
dataCollection <- Reduce(function(x,y) left_join(x,y, by="ensembl_gene_id"), 
                         c(list(dataCollection), tableResultDGEA.list, list(geneListInfo), list(geneExpressionDataInfo)))

# Export Table as Tab-Separated Text File
write_tsv(dataCollection, "dge_results_collection.txt")

#-----------------------------------------------------------------------------#
# Exporting Excel File 
#-----------------------------------------------------------------------------#

# Add Result Data Frame to Workbook
addWorksheet(workbook, "DGEA Results")
writeData(workbook, "DGEA Results", dataCollection)

# Save Excel File 
saveWorkbook(workbook, "study_data_collection.xlsx", overwrite = TRUE)

#=============================================================================#
#=============================================================================#
# END
#=============================================================================#
#=============================================================================#




