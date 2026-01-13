###############################################################################
###############################################################################
# PRO4002 Research Project                                                    #
#																	                                            #
# Date: Jan 22, 2026											                                    #
# Author: Sadegh, Matas, Nur, Arlin & Oona                                    #  
###############################################################################
###############################################################################

#-----------------------------------------------------------------------------#
# LIBRARY/PACKAGE INSTALLATION
#-----------------------------------------------------------------------------#

analysis_packages <- c( 
  "tidyverse", "tidyr", "dplyr", "gridExtra", "pcaMethods",
  "data.table", "tableone", "kableExtra", "rmarkdown",
  "readr", "readxl", "gprofiler2", "knitr", "Hmisc", "table1",
  "gt", "gtsummary", "ggsci"
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
  text = element_text(family = target_font),        # Forces the font
  plot.title = element_text(hjust = 0.5, face="bold"), # Centers & bolds title
  axis.title = element_text(face = "bold"),         # Bolds axis titles
  legend.position = "right"                         # Default legend position
)


#-----------------------------------------------------------------------------#
# SETUP
#-----------------------------------------------------------------------------#

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) { 
  script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
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

#-----------------------------------------------------------------------------#
# DATA IMPORT
#-----------------------------------------------------------------------------#
setwd(data_path)

exonLengthsData <- read.delim("MAGNET_exonLengths.txt", header=TRUE, row.names = 1, as.is = T)
geneExpressionData <- read.delim("MAGNET_GeneExpressionData_CPM_19112020.txt", header=TRUE, row.names = 1, as.is = T)
sampleDescriptionData <- read_excel("MAGNET_SampleData_18112022_WithDescriptions.xlsx", sheet=2, na="NA")
sampleDataRaw <- read_excel("MAGNET_SampleData_18112022_WithDescriptions.xlsx", sheet=1, na="NA")

setwd(script_path)

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

#-----------------------------------------------------------------------------#
# PARTICIPANT INFO DATA - PREPROCESSING
#-----------------------------------------------------------------------------#

sampleData <- sampleDataRaw

sampleData.colnames <- colnames(sampleData)
sampleData.colnames <- gsub("LVEF", "lvef", sampleData.colnames)
sampleData.colnames <- gsub("RIN", "rin", sampleData.colnames)
sampleData.colnames <- gsub("TIN.median.", "tin_median", sampleData.colnames, fixed = TRUE) 
colnames(sampleData) <- sampleData.colnames

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

sampleData <- sampleData %>% mutate(bmi = weight/(height*0.01)^2)

#-----------------------------------------------------------------------------#
# PARTICIPANT INFO DATA - DEMOGRAPHIC TABLE
#-----------------------------------------------------------------------------#

demoTableLabels <- list(variables = list(age="Age (years)",
                                         weight="Weight (kg)",
                                         height="Height (cm)",
                                         bmi="BMI (kg/m²)",
                                         gender="Gender",
                                         race="Ethnicity",
                                         hw="Heart Weight (g)",
                                         lv_mass="LVᵃ Mass (g)",
                                         afib="Atrial Fibrillation",
                                         VTVF="VTVFᵇ",
                                         Diabetes="Diabetes Mellitus",
                                         Hypertension="Hypertension",
                                         lvef="LVᵃ Ejection Fraction" ),
                        groups=list("DCM", "HCM", "PPCM", "NF"))
demoTableFootnote  <- c("ᵃ Left Ventricle", "ᵇ Ventricular Tachycardia/Ventricular Fibrillation")
demoTableStrata <- c(list(Total=sampleData), split(sampleData, sampleData$etiology))
 
demographicTable <- table1(demoTableStrata, demoTableLabels, 
                           footnote=demoTableFootnote, topclass="Rtable1-zebra")

#-----------------------------------------------------------------------------#
# PARTICIPANT INFO DATA - CLINICAL CHARACTERISTICS BOXPLOT
#-----------------------------------------------------------------------------#

clinicalSampleData <- sampleData %>%
  dplyr::select(etiology, age, weight, height, lvef, rin) %>%
  tidyr::pivot_longer(
    cols = -etiology,         
    names_to = "characteristic",
    values_to = "value"        
  ) %>%
  dplyr::mutate(
    characteristic = str_to_title(characteristic),
    characteristic = str_replace_all(
      characteristic,
      c(
        "Age"    = "Age (years)",
        "Weight" = "Weight (kg)",
        "Height" = "Height (cm)",
        "Lvef"   = "LVEF (%)",
        "Rin"    = "RIN (mg/dL)"
      )
    )
  )

clinicalSampleDataPlot <- ggplot(clinicalSampleData, aes(x = etiology, y = value, fill = etiology)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white") +
  facet_wrap(~ characteristic, scales = "free_y", nrow = 2) + 
  scale_color_npg() +
  scale_fill_npg() +
  labs(
    title = "Patient Characteristics by Etiology",
    x = "Disease Group",
    y = "Value (scales vary)",
    fill = "Etiology"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none", # Legend is redundant since x-axis is labeled
    strip.background = element_rect(fill = "lightgrey"), # Gray headers for panels
    strip.text = element_text(face = "bold")
  ) +
  center_title + my_style

ggsave(file.path(plots_path, "clinicalSampleDataPlot.png"), clinicalSampleDataPlot, width = 10, height = 8)

# Violoing Plots -------------------------------------------------------------#

# Age Distribution (Violin Plot per Etiology) 
violinPlotAgeDistribution <- ggplot(sampleData, aes(x = etiology, y = age, fill = etiology)) +
  geom_violin() +
  geom_jitter(width=0.2, size=1, alpha=0.5)+ 
  labs(title = "Distribution of Sample Ages per Etiology", x = "Etiology",
       y = "Age (years)") +
  scale_color_npg() +
  scale_fill_npg() + 
  center_title + my_style +
  theme(legend.position="none")  

ggsave(file.path(plots_path, "violinPlotAgeDistribution.jpg"), 
       plot = violinPlotAgeDistribution, 
       width = 6, height = 4)

# RIN Distribution (Violin Plot per Etiology) 
violinPlotRINDistribution <- ggplot(sampleData, aes(x = etiology, y = rin, fill = etiology)) +
  geom_violin() +
  geom_jitter(width=0.2, size=1, alpha=0.5)+ 
  labs(title = "Distribution of Sample RIN (RNA Integrity) per Etiology", x = "Etiology",
       y = "RIN (mg/dL)") +
  scale_color_npg() +
  scale_fill_npg() + 
  center_title + my_style +
  theme(legend.position="none")  

ggsave(file.path(plots_path, "violinPlotRINDistribution.jpg"), 
       plot = violinPlotRINDistribution, 
       width = 6, height = 4)

#-----------------------------------------------------------------------------#
# GENE EXPRESSION DATA - DATA DISTRIBUTION 
#-----------------------------------------------------------------------------#

extendedGeneExpressionData <- geneExpressionData %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols=-gene, names_to="patient", values_to="gene_expression_level") %>%
  left_join(sampleData, by=c('patient'='sample_name'))

# Box Plots ------------------------------------------------------------------#

# Data Distribution (Box Plots per Etiology - Grid)
boxPlotDataDistribution.etiology.grid <- ggplot(data = extendedGeneExpressionData, 
                                           aes(x = patient, 
                                               y = gene_expression_level, 
                                               color = etiology)) +
  geom_boxplot() + 
  labs(title = "Patients vs. Gene Expression Level", x = "Patient ID",
       y = "Gene Expression Level (logCPM)", color = 'Etiology') +
  scale_color_npg() +
  scale_fill_npg() +
  facet_wrap(~etiology, scales="free_x") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(angle = 90, hjust = 0.5)
    ) +
  center_title + my_style

ggsave(file.path(plots_path, "boxPlotDataDistribution_etiology_grid.jpg"), 
       plot = boxPlotDataDistribution.etiology.grid, 
       width = 8, height = 6, dpi = 300)

# Data Distribution (Box Plots per Etiology - Inline)
boxPlotDataDistribution.etiology.inline <- ggplot(data = extendedGeneExpressionData, 
                                           aes(x = patient, 
                                               y = gene_expression_level, 
                                               color = etiology)) +
  geom_boxplot() + 
  labs(title = "Patients vs. Gene Expression Level", x = "Patient ID",
       y = "Gene Expression Level (logCPM)", color = 'Etiology') +
  scale_color_npg() +
  scale_fill_npg() +
  facet_grid(~etiology, scales = "free_x", space="free_x") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(angle = 90, hjust = 0.5)
    ) +
  center_title + my_style

ggsave(file.path(plots_path, "boxPlotDataDistribution_etiology_inline.jpg"), 
       plot = boxPlotDataDistribution.etiology.inline,
       width = 20, height = 8 )

# Density Plots -------------------------------------------------------------#

# Data Distribution (Density Plots per Etiology - Grid)
densityPlotDataDistribution.etiology.grid <- ggplot(data = extendedGeneExpressionData, 
                                               aes(x = gene_expression_level, 
                                                   fill = etiology, color = etiology)) +
  geom_density(alpha=0.7) +
  labs(title="Gene Expression Level Density Plot", x="Gene Expression Level (logCPM)",
       y="Density", color="Etiology", fill="Etiology") +
  scale_color_npg() +
  scale_fill_npg() +
  facet_wrap(~etiology, scales="free_x") +
  center_title + my_style

ggsave(file.path(plots_path, "densityPlotDataDistribution_etiology_grid.jpg"), 
       plot = densityPlotDataDistribution.etiology.grid, 
       width = 8, height = 6, dpi = 300 )

# Data Distribution (Density Plots per Etiology - Overlap) 
densityPlotDataDistribution.etiology.overlap <- ggplot(data = extendedGeneExpressionData, 
                                               aes(x = gene_expression_level, 
                                                   fill = etiology, color = etiology)) +
  geom_density(alpha=0.3) +
  labs(title  ="Gene Expression Level Density Plot", x="Gene Expression Level (logCPM)",
       y = "Density", color = "Etiology", fill = "Etiology") +
  scale_color_npg() +
  scale_fill_npg() +
  center_title + my_style

ggsave(file.path(plots_path, "densityPlotDataDistribution_etiology_overlap.jpg"), 
       plot = densityPlotDataDistribution.etiology.overlap, 
       width = 8, height = 6, dpi = 300 )

#-----------------------------------------------------------------------------#
# GENE EXPRESSION DATA - PRINCIPLE COMPONENT ANALYSIS 
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
