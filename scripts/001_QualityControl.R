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

#-----------------------------------------------------------------------------#
# DATA IMPORT
#-----------------------------------------------------------------------------#
setwd(data_path)

exonLengthsData <- read.delim("MAGNET_exonLengths.txt", header=TRUE, row.names = 1, as.is = T)
geneExpressionData.CPM <- read.delim("MAGNET_GeneExpressionData_CPM_19112020.txt", header=TRUE, row.names = 1, as.is = T)
sampleDescriptionData <- read_excel("MAGNET_SampleData_18112022_WithDescriptions.xlsx", sheet=2, na="NA")
sampleData.RAW <- read_excel("MAGNET_SampleData_18112022_WithDescriptions.xlsx", sheet=1, na="NA")

setwd(script_path)

# View Imported Data As Simple Quality Check 
View(exonLengthsData)
View(geneExpressionData.CPM)
View(sampleDescriptionData)
View(sampleData.RAW)

# Inspect Dimensions, Storage and NAs in Preparation of Preprocessing 
Hmisc::contents(exonLengthsData)
Hmisc::contents(geneExpressionData.CPM)
Hmisc::contents(sampleData.RAW)

# Check for Data Conformance - do Patient-IDs and Gene-IDs match across tables
all(rownames(geneExpressionData.CPM) == rownames(exonLengthsData))
all(colnames(geneExpressionData.CPM) == sampleData.RAW$sample_name)

#-----------------------------------------------------------------------------#
# IMPORT ENSEMBLE DATA
#-----------------------------------------------------------------------------#

# Import Human Ensembl Dataset
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Extract Gene Information of Genes in the Expression Data
geneListInfo <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description", "chromosome_name"),
  filters = "ensembl_gene_id",,
  values = rownames(geneExpressionData.CPM),
  mart = ensembl
)

#-----------------------------------------------------------------------------#
# PARTICIPANT INFO DATA - PREPROCESSING
#-----------------------------------------------------------------------------#

sampleData <- sampleData.RAW

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

# Violine Plots -------------------------------------------------------------#

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
# CPM TO FPKM CONVERSION 
#-----------------------------------------------------------------------------#

cpm2fpkm <- function(x, y) {
  .t <- 2^(x) * 1E3 / y[, 1]
}

# Convert Gene Expression Data from logCPM Values to FPKM Values
geneExpressionData.FPKM <- cpm2fpkm(geneExpressionData.CPM, exonLengthsData)

#-----------------------------------------------------------------------------#
# GENE EXPRESSION DATA - BACKGROUND NOISE REMOVAL
#-----------------------------------------------------------------------------#

annotatedGeneExpressionData.FPKM <- merge(geneExpressionData.FPKM, geneListInfo, by.x="row.names", by="ensembl_gene_id")

geneExpressionData.FPKM.femaleY <- geneExpressionData.FPKM[annotatedGeneExpressionData.FPKM$chromosome_name == "Y",
                                                    sampleData$gender == "Female", drop=FALSE]

# Create Numerical List of Gene Expression Data of Female Patient Y-Chromosome Genes
backgroundExpressionData.femaleY <- as.numeric(as.matrix(geneExpressionData.FPKM.femaleY))

backgroundThreshold95PExpressionData <- quantile(backgroundExpressionData.femaleY, 
                                                 probs=0.95, na.rm=TRUE)
message(paste("Background Noise 95% Quantile Threshold (FPKM):", round(backgroundThreshold95PExpressionData, 4)))

backgroundThresholdMeanExpressionData <- mean(backgroundExpressionData.femaleY)
message(paste("Background Noise Mean Threshold (FPKM):", round(backgroundThresholdMeanExpressionData, 4)))

extendedGeneExpressionData.FPKM.femaleY <- geneExpressionData.FPKM.femaleY %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  pivot_longer(cols=-ensembl_gene_id, names_to="patient", values_to="gene_expression_level") %>%
  left_join(sampleData, by=c('patient'='sample_name'))

extendedGeneExpressionData.FPKM <- geneExpressionData.FPKM %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  pivot_longer(cols=-ensembl_gene_id, names_to="patient", values_to="gene_expression_level") %>%
  left_join(sampleData, by=c('patient'='sample_name'))

# Plot Density Expression Level of Female Patient Y-Chromosome Genes with Threshold 
noiseExpressionDensityPlot <- ggplot(data=extendedGeneExpressionData.FPKM.femaleY,
                                     aes(x=log10(gene_expression_level))) +
  geom_density(fill = npg_colors[2], color = npg_colors[2], alpha = 0.7) +
  geom_vline(xintercept = log10(backgroundThreshold95PExpressionData), 
             color = npg_colors[1], linetype = "dashed", linewidth=0.8) +
  annotate("text", x = log10(backgroundThreshold95PExpressionData), y = 0, 
           label = "95% Threshold", color = npg_colors[1], angle =90, vjust = -0.5, hjust=-0.3, size =4) +
  geom_vline(xintercept = log10(backgroundThresholdMeanExpressionData), 
             color = npg_colors[4], linetype = "dashed", linewidth=0.8) +
  annotate("text", x = log10(backgroundThresholdMeanExpressionData), y = 0, 
           label = "Mean Threshold", color = npg_colors[4], angle =90, vjust = -0.5, hjust=-0.3, size =4) +
  labs(title="Female - Y Chromosome Gene Expression with Threshold", x="Gene Expression Level (FPKM)", y="Density") 
noiseExpressionDensityPlot

densityPlotDataDistribution.threshold <- ggplot(data = extendedGeneExpressionData.FPKM, 
                                              aes(x = log10(gene_expression_level))) +
  geom_density(fill = npg_colors[2], color = npg_colors[2], alpha = 0.7) +
  geom_vline(xintercept=log10(backgroundThreshold95PExpressionData), 
             color = npg_colors[1], linetype = "dashed", linewidth=0.8) +
  annotate("text", x = log10(backgroundThreshold95PExpressionData), y = 0, 
           label = "95% Threshold", color = npg_colors[1], angle =90, vjust = -0.5, hjust=-0.3, size =4) +
  geom_vline(xintercept = log10(backgroundThresholdMeanExpressionData), 
             color = npg_colors[4], linetype = "dashed", linewidth=0.8) +
  annotate("text", x = log10(backgroundThresholdMeanExpressionData), y = 0, 
           label = "Mean Threshold", color = npg_colors[4], angle =90, vjust = -0.5, hjust=-0.3, size =4) +
  labs(title  ="Gene Expression Level Density Plot with Threshold", x="Gene Expression Level (FPKM)",y = "Density")
densityPlotDataDistribution.threshold

meanGeneExpressionData <- rowMeans(geneExpressionData.FPKM, na.rm=TRUE)

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

geneExpressionData.CPM.95PFiltered <- geneExpressionData.CPM[aboveBackground95PThresholdGenes,]
geneExpressionData.CPM.meanFiltered <- geneExpressionData.CPM[aboveBackgroundMeanThresholdGenes,]

#-----------------------------------------------------------------------------#
# DCM PATIENT GENE SET
#-----------------------------------------------------------------------------#

sampleData.DCM <- subset(sampleData, sampleData$etiology == "DCM")
geneExpressionData.CPM.DCM <- geneExpressionData.CPM[, colnames(geneExpressionData.CPM) %in% sampleData.DCM$sample_name]
geneExpressionData.CPM.95PFiltered.DCM <- geneExpressionData.CPM.95PFiltered[, 
  colnames(geneExpressionData.CPM.95PFiltered) %in% sampleData.DCM$sample_name] 
geneExpressionData.CPM.meanFiltered.DCM <- geneExpressionData.CPM.meanFiltered[, 
  colnames(geneExpressionData.CPM.meanFiltered) %in% sampleData.DCM$sample_name] 

#-----------------------------------------------------------------------------#
# EXPORT DCM PATIENT GENE SETS
#-----------------------------------------------------------------------------#

saveRDS(geneExpressionData.CPM.DCM, file.path(cache_path, "geneExpressionData_DCM.rds"))
saveRDS(geneExpressionData.CPM.95PFiltered.DCM, file.path(cache_path, "geneExpressionData_DCM_95PFiltered.rds"))
saveRDS(geneExpressionData.CPM.meanFiltered.DCM, file.path(cache_path, "geneExpressionData_DCM_meanFiltered.rds"))

#-----------------------------------------------------------------------------#
# DCM PATIENT GENE SET - DENSITY PLOTs
#-----------------------------------------------------------------------------#

extendedGeneExpressionData.CPM.meanFiltered.DCM <- geneExpressionData.CPM.meanFiltered.DCM %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols=-gene, names_to="patient", values_to="gene_expression_level") %>%
  left_join(sampleData, by=c('patient'='sample_name'))

# Data Distribution 
densityPlotDataDistribution.meanFiltered.DCM <- ggplot(data = extendedGeneExpressionData.CPM.meanFiltered.DCM, 
  aes(x = gene_expression_level)) +
  geom_density(alpha=0.7, fill=npg_colors[[1]]) +
  labs(title="Gene Expression Level Density Plot", x="Gene Expression Level (logCPM)", y="Density") +
  scale_color_npg() +
  scale_fill_npg() +
  center_title + my_style

ggsave(file.path(plots_path, "densityPlotDataDistribution_meanFiltered_DCM.jpg"), 
       plot = densityPlotDataDistribution.meanFiltered.DCM, 
       width = 8, height = 6, dpi = 300 )

extendedGeneExpressionData.CPM.95PFiltered.DCM <- geneExpressionData.CPM.95PFiltered.DCM %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols=-gene, names_to="patient", values_to="gene_expression_level") %>%
  left_join(sampleData, by=c('patient'='sample_name'))

# Data Distribution 
densityPlotDataDistribution.95PFiltered.DCM <- ggplot(data = extendedGeneExpressionData.CPM.95PFiltered.DCM, 
  aes(x = gene_expression_level)) +
  geom_density(alpha=0.7, fill=npg_colors[[1]]) +
  labs(title="Gene Expression Level Density Plot", x="Gene Expression Level (logCPM)", y="Density") +
  scale_color_npg() +
  scale_fill_npg() +
  center_title + my_style

ggsave(file.path(plots_path, "densityPlotDataDistribution_95PFiltered_DCM.jpg"), 
       plot = densityPlotDataDistribution.95PFiltered.DCM, 
       width = 8, height = 6, dpi = 300 )

#-----------------------------------------------------------------------------#
# SYSTEM-SPECIFIC GENE SET - IMMUNE
#-----------------------------------------------------------------------------#

innateDB.RAW <- read.csv(file.path(data_path, "InnateDB_genes.csv"), header=TRUE, na="NA")
innateDB <- subset(innateDB.RAW, innateDB.RAW$species == "Homo sapiens")

geneExpressionData.CPM.DCM.INNATE <- geneExpressionData.CPM.DCM[rownames(geneExpressionData.CPM.DCM) %in% innateDB$ensembl, ]
geneExpressionData.CPM.meanFiltered.DCM.INNATE <- geneExpressionData.CPM.meanFiltered[
  rownames(geneExpressionData.CPM.meanFiltered) %in% innateDB$ensembl, ]

saveRDS(geneExpressionData.CPM.DCM.INNATE, file.path(cache_path, "geneExpressionData_DCM_INNATE.rds"))
saveRDS(geneExpressionData.CPM.meanFiltered.DCM.INNATE, file.path(cache_path, "geneExpressionData_DCM_meanFiltered_INNATE.rds"))
