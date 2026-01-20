###############################################################################
###############################################################################
# PRO4002 Research Project                                                    #
#																	          #
# File: 003_PhenotypeValidation.R                                             #
# Date: Jan 22, 2026											              #
# Author: Sadegh, Matas, Nur, Arlin & Oona                                    #  
###############################################################################
###############################################################################

#-----------------------------------------------------------------------------#
# SETUP
#-----------------------------------------------------------------------------#

source(here::here("scripts", "000_setup.R"))

phenotype_validation_path <- file.path(plots_path, "PhenotypeValidation")
if(!dir.exists(phenotype_validation_path)) dir.create(phenotype_validation_path, recursive = TRUE)

message("\n--- Starting Phenotype Validation ---")

#-----------------------------------------------------------------------------#
# DATA IMPORT
#-----------------------------------------------------------------------------#

message("(I) Data Import ")

# Load Expression Data (Same as 002)
geneExpressionData <- readRDS(file.path(cache_path, "geneExpressionData_DCM_meanFiltered.rds"))
# Load Full Metadata (For RIN, Pool, etc.)
sampleData <- readRDS(file.path(cache_path, "sampleData_DCM.rds"))

# Load Your New Clusters (From 002)
sampleData.subtypes <- readRDS(file.path(cache_path, "sampleData_DCM_subtypes.rds"))

# Merge: Add 'subtype' to the full metadata
# We match by row names (sample names)
sampleData$subtype <- sampleData.subtypes$subtype[match(sampleData$sample_name, rownames(sampleData.subtypes))]
sampleData$rin[is.na(sampleData$rin)] <- median(sampleData$rin, na.rm = TRUE)

# Remove dwares with weird BMIs
sampleData <- sampleData %>% filter(height>100)
rownames(sampleData) <- sampleData$sample_name
# Filter: Ensure we only keep samples that have a subtype
sampleData <- sampleData %>% filter(!is.na(subtype))

#-----------------------------------------------------------------------------#
# DATA VISUALIZATION
#-----------------------------------------------------------------------------#

message("(II) Data Visualization ")

#-----------------------------------------------------------------------------#
# DATA VISUALIZATION - Disease vs Subtype 
#-----------------------------------------------------------------------------#

table_data <- table(sampleData$subtype, sampleData$Diabetes)
chi_test <- chisq.test(table_data)
chi_test$p.value
p_label <- ifelse(chi_test$p.value < 0.001, "p < 0.001", paste0("p = ", signif(chi_test$p.value, 3)))

sampleDataDiabetesSubtypeBarPlot <- ggplot(sampleData, aes(x = subtype, fill = Diabetes, color = Diabetes)) +
  geom_bar(position = "dodge", alpha=0.7) + 
  labs(title = "Diabetes", x ="", y = "Count", fill = "Subtype", color = "Subtype") +
  scale_fill_manual(values = c(standard_color, accent_color)) +
  scale_color_manual(values = c(standard_color, accent_color)) +
  annotate("text", x = 3, y = 70, label = p_label, size = 8) +
  scale_y_continuous(limits = c(0, 75)) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataDiabetesSubtypeBarPlot.png"), plot = sampleDataDiabetesSubtypeBarPlot,
  width = 10, height = 8)

table_data <- table(sampleData$subtype, sampleData$afib)
chi_test <- chisq.test(table_data)
chi_test$p.value
p_label <- ifelse(chi_test$p.value < 0.001, "p < 0.001", paste0("p = ", signif(chi_test$p.value, 3)))

sampleDataAFibSubtypeBarPlot <- ggplot(sampleData, aes(x = subtype, fill = afib, color = afib)) +
  geom_bar(position = "dodge", alpha=0.7) + 
  labs(title = "Atrial Fibrillation", x ="", y = "Count", fill = "Subtype", color = "Subtype") +
  scale_fill_manual(values = c(standard_color, accent_color)) +
  scale_color_manual(values = c(standard_color, accent_color)) +
  annotate("text", x = 3, y = 70, label = p_label, size = 8) +
  scale_y_continuous(limits = c(0, 75)) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataAFibSubtypeBarPlot.png"), plot = sampleDataAFibSubtypeBarPlot,
  width = 10, height = 8)

table_data <- table(sampleData$subtype, sampleData$VTVF)
chi_test <- chisq.test(table_data)
chi_test$p.value
p_label <- ifelse(chi_test$p.value < 0.001, "p < 0.001", paste0("p = ", signif(chi_test$p.value, 3)))

sampleDataVTVFSubtypeBarPlot <- ggplot(sampleData, aes(x = subtype, fill = VTVF, color = VTVF)) +
  geom_bar(position = "dodge", alpha=0.7) + 
  labs(title = "VTVF", x ="", y = "Count", fill = "Subtype", color = "Subtype") +
  scale_fill_manual(values = c(standard_color, accent_color)) +
  scale_color_manual(values = c(standard_color, accent_color)) +
  annotate("text", x = 3, y = 70, label = p_label, size = 8) +
  scale_y_continuous(limits = c(0, 75)) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataVTVFSubtypeBarPlot.png"), plot = sampleDataVTVFSubtypeBarPlot,
  width = 10, height = 8)

table_data <- table(sampleData$subtype, sampleData$Hypertension)
chi_test <- chisq.test(table_data)
chi_test$p.value
p_label <- ifelse(chi_test$p.value < 0.001, "p < 0.001", paste0("p = ", signif(chi_test$p.value, 3)))

sampleDataHypertensionSubtypeBarPlot <- ggplot(sampleData, aes(x = subtype, fill = Hypertension, color = Hypertension)) +
  geom_bar(position = "dodge", alpha=0.7) + 
  labs(title = "Hypertension", x ="", y = "Count", fill = "Subtype", color = "Subtype") +
  scale_fill_manual(values = c(standard_color, accent_color)) +
  scale_color_manual(values = c(standard_color, accent_color)) +
  annotate("text", x = 3, y = 70, label = p_label, size = 8) +
  scale_y_continuous(limits = c(0, 75)) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataHypertensionSubtypeBarPlot.png"), plot = sampleDataHypertensionSubtypeBarPlot,
  width = 10, height = 8)

# sampleDataDiabetesSubtypeBarPlot2 <- ggplot(sampleData, aes(x = subtype, fill = subtype, color=subtype)) +
#  geom_bar(alpha=0.7, show.legend = TRUE) +
#  labs(title = "Diabetes", x ="", y = "Count", fill = "Subtype", color = "Subtype") +
#  scale_fill_manual(values = npg_colors[1:3]) +
#  scale_color_manual(values = npg_colors[1:3]) +
#  facet_wrap( ~ Diabetes)+
#  center_title + my_style
# ggsave(filename = file.path(phenotype_validation_path, "sampleDataDiabetesSubtypeBarPlot2.png"), plot = sampleDataDiabetesSubtypeBarPlot2)

#-----------------------------------------------------------------------------#
# DATA VISUALIZATION - Clinical Data
#-----------------------------------------------------------------------------#

sampleDataLVEFSubtypeViolinPlot <- ggplot(sampleData, aes(x = subtype, y = (lvef*100), fill = subtype, color = subtype)) +
  geom_violin(alpha= 0.7) +
  geom_jitter(width = 0.2, color="black", size = 1) + 
  stat_compare_means(method = "kruskal.test", label.y = 95, size = 8) +
  labs(title = "LVEF", y = "LVEF (%)", fill = "Subtype", color = "Subtype", x = "") +
  scale_fill_manual(values = npg_colors[1:3]) +
  scale_color_manual(values = npg_colors[1:3]) +
  scale_y_continuous(limits = c(0, 100)) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataLVEFSubtypeViolinPlot.png"), plot = sampleDataLVEFSubtypeViolinPlot,
  width = 10, height = 8)

sampleDataAgeSubtypeViolinPlot <- ggplot(sampleData, aes(x = subtype, y = age, fill = subtype, color = subtype)) +
  geom_violin(alpha= 0.7) +
  geom_jitter(width = 0.2, color="black", size = 1) + 
  stat_compare_means(method = "kruskal.test", label.y = 75, size = 8) +
  labs(title = "Age", y = "Age (Years)", fill = "Subtype", color = "Subtype", x = "") +
  scale_fill_manual(values = npg_colors[1:3]) +
  scale_color_manual(values = npg_colors[1:3]) +
  scale_y_continuous(limits = c(15, 80)) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataAgeSubtypeViolinPlot.png"), plot = sampleDataAgeSubtypeViolinPlot,
  width = 10, height = 8)

sampleDataBMISubtypeViolinPlot <- ggplot(sampleData, aes(x = subtype, y = bmi, fill = subtype, color = subtype)) +
  geom_violin(alpha= 0.7) +
  geom_jitter(width = 0.2, color="black", size = 1) + 
   stat_compare_means(method = "kruskal.test", label.y = 75, size = 8) +
  labs(title = "Body Mass Index", y = "BMI (kg/m²)", fill = "Subtype", color = "Subtype", x = "") +
  scale_fill_manual(values = npg_colors[1:3]) +
  scale_color_manual(values = npg_colors[1:3]) +
  scale_y_continuous(limits = c(15, 80)) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataBMISubtypeViolinPlot.png"), plot = sampleDataBMISubtypeViolinPlot,
  width = 10, height = 8)

sampleDataTINSubtypeViolinPlot <- ggplot(sampleData, aes(x = subtype, y = tin_median, fill = subtype, color = subtype)) +
  geom_violin(alpha= 0.7) +
  geom_jitter(width = 0.2, color="black", size = 1) + 
  stat_compare_means(method = "kruskal.test", label.y = 85, size = 8) +
  labs(title = "Transcript Integrity Number ", y = "TIN Median", fill = "Subtype", color = "Subtype", x = "") +
  scale_fill_manual(values = npg_colors[1:3]) +
  scale_color_manual(values = npg_colors[1:3]) +
  scale_y_continuous(limits = c(20, 90)) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataTINSubtypeViolinPlot.png"), plot = sampleDataTINSubtypeViolinPlot,
  width = 10, height = 8)

sampleDataRINSubtypeViolinPlot <- ggplot(sampleData, aes(x = subtype, y = rin, fill = subtype, color = subtype)) +
  geom_violin(alpha= 0.7) +
  geom_jitter(width = 0.2, color="black", size = 1) + 
  stat_compare_means(method = "kruskal.test", label.y = 10.5, size = 8) +
  labs(title = "RNA Integrity Number", y = "RIN", fill = "Subtype", color = "Subtype", x = "") +
  scale_fill_manual(values = npg_colors[1:3]) +
  scale_color_manual(values = npg_colors[1:3]) + 
  scale_y_continuous(limits = c(5, 11)) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataRINSubtypeViolinPlot.png"), plot = sampleDataRINSubtypeViolinPlot,
  width = 10, height = 8)

sampleDataHWSubtypeViolinPlot <- ggplot(sampleData, aes(x = subtype, y = hw, fill = subtype, color = subtype)) +
  geom_violin(alpha= 0.7) +
  geom_jitter(width = 0.2, color="black", size = 1) + 
  labs(title = "Heart Weight", y = "Heart Weight (g)", fill = "Subtype", color = "Subtype", x = "") +
  stat_compare_means(method = "kruskal.test", label.y = 90, size = 8) +
  scale_fill_manual(values = npg_colors[1:3]) +
  scale_color_manual(values = npg_colors[1:3]) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataHWSubtypeViolinPlot.png"), plot = sampleDataHWSubtypeViolinPlot,
  width = 10, height = 8)

sampleDataLVSubtypeViolinPlot <- ggplot(sampleData, aes(x = subtype, y = lv_mass, fill = subtype, color = subtype)) +
  geom_violin(alpha= 0.7) +
  geom_jitter(width = 0.2, color="black", size = 1) + 
  labs(title = "LV Mass", y = "LV Mass (g)", fill = "Subtype", color = "Subtype", x = "") +
  stat_compare_means(method = "kruskal.test", label.y = 150, size = 8) +
  scale_fill_manual(values = npg_colors[1:3]) +
  scale_color_manual(values = npg_colors[1:3]) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataLVSubtypeViolinPlot.png"), plot = sampleDataLVSubtypeViolinPlot,
  width = 10, height = 8)

sampleDataLPSubtypeBarPlot <- ggplot(sampleData, aes(x = subtype, fill = Library.Pool, color = Library.Pool)) +
  geom_bar(position = "dodge", alpha = 0.7) + 
  labs(title = "Library Pools", x = "", y = "Count", fill = "Library.Pool") +
  scale_fill_manual(values = tail(npg_additional_colors,-3)) +
  scale_color_manual(values = tail(npg_additional_colors,-3)) +
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataLPSubtypeBarPlot.png"), plot = sampleDataLPSubtypeBarPlot,
  width = 10, height = 8)

#-----------------------------------------------------------------------------#
# COMPLETE
#-----------------------------------------------------------------------------#

message("--- Finished Phenotype Validation  ---")