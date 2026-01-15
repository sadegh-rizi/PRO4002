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

sampleData.DCM.annotated <- readRDS(file.path(cache_path, "sampleData_DCM_subtypes.rds"))
sampleData.DCM <- readRDS(file.path(cache_path, "sampleData_DCM.rds"))
sampleData.DCM.complete <- sampleData.DCM
sampleData.DCM.complete$subtype <- sampleData.DCM.annotated$subtype

#-----------------------------------------------------------------------------#
# DATA VISUALIZATION
#-----------------------------------------------------------------------------#

message("(II) Data Visualization ")


sampleDataLVEFSubtypeBarPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, y = lvef, fill = subtype)) +
  geom_bar(stat = "identity") +
  stat_compare_means(method = "kruskal.test", label.y = 50) + 
  labs(title = "Is Cluster 4 the 'Severe' Group?", y = "LVEF (%)") +
  scale_color_npg() +
  scale_fill_npg() +
  scale_y_continuous(limits = c(0, 20)) +
  center_title + my_style

sampleDataLVEFSubtypeViolinPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, y = lvef, fill = subtype)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.5) + # Add dots to see the 9 patients
  stat_compare_means(method = "kruskal.test", label.y = 50) + 
  labs(title = "Is Cluster 4 the 'Severe' Group?", y = "LVEF (%)") +
  scale_color_npg() +
  scale_fill_npg() +
  scale_y_continuous(limits = c(0, 1)) + 
  center_title + my_style

sampleDataAgeSubtypeViolinPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, y = age, fill = subtype)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(method = "kruskal.test", label.y = 90) +
  labs(title = "Violin plot of Age Distribution", y = "Age (Years)") +
  scale_color_npg() +
  scale_fill_npg() +
  center_title + my_style

sampleDataBMISubtypeViolinPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, y = bmi, fill = subtype)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(method = "kruskal.test", label.y = 90) +
  labs(title = "Violin plot of BMI Distribution", y = "BMI - Weight(kg) / Height(m)^2") +
  scale_color_npg() +
  scale_fill_npg() +
  center_title + my_style
  # scale_y_continuous(limits = c(0, 100)) 
  # with scale_y_continuous it excludes the outliers of BMI, which greatly increases the p value 
  # from kruskal test. This means the outliers are impactful data, enough to alter the p value

sampleDataBMISubtypeBoxPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, y = bmi, fill = subtype)) +
  geom_boxplot(outlier.shape = 2) + # "2" is used as the shape - triangle to visualise the outliers
  geom_jitter(width = 0.2, alpha = 0.5) + # Add dots to see the 9 patients
  stat_compare_means(method = "kruskal.test", label.y = 50) + 
  labs(title = "Box plot of BMI Distribution", y = "BMI - Weight(kg) / Height(m)^2") +
  scale_color_npg() +
  scale_fill_npg() +
  center_title + my_style

sampleDataDiabetesSubtypeBarPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, fill = Diabetes)) +
  geom_bar(position = "dodge") + 
  labs(title = "Bar plot of Diabetes for each Endotype", y = "Diabetes", fill = "Diabetes") +
  scale_color_npg() +
  scale_fill_npg() +
  center_title + my_style

sampleDataHWSubtypeViolinPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, y = hw, fill = subtype)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(method = "kruskal.test", label.y = 90) +
  labs(title = "Violin plot of Heart Weight Distribution", y = "Heart Weight(g)") +
  scale_color_npg() +
  scale_fill_npg() +
  center_title + my_style
sampleDataHWSubtypeViolinPlot

sampleDataTINSubtypeViolinPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, y = tin_median, fill = subtype)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(method = "kruskal.test", label.y = 90) +
  labs(title = "Violin plot of Transcript Integrity Number (TIN) Distribution", y = "TIN(median)") +
  scale_color_npg() +
  scale_fill_npg() +
  center_title + my_style


sampleDataRINSubtypeViolinPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, y = RIN, fill = subtype)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(method = "kruskal.test", label.y = 90) +
  labs(title = "Violin plot of RNA Integrity Number (RIN) Distribution", y = "RIN") +
  scale_color_npg() +
  scale_fill_npg() +
  center_title + my_style
  scale_y_continuous(limits = c(0, 20))

sampleDataLPSubtypeBarPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, fill = Library.Pool)) +
  geom_bar(position = "dodge") + 
  labs(title = "Bar plot of Library Pools for each Endotype", y = "Samples", fill = "Library.Pool") +
  scale_color_npg() +
  scale_fill_manual(values = npg_additional_colors) +
  center_title + my_style

ggsave(filename = file.path(phenotype_validation_path, "LVEF_Bar_Plot.png"), plot = sampleDataLVEFSubtypeBarPlot)
ggsave(filename = file.path(phenotype_validation_path, "LVEF_Violin_Plot.png"), plot = sampleDataLVEFSubtypeViolinPlot)
ggsave(filename = file.path(phenotype_validation_path, "Age_Violin_Plot.png"), plot = sampleDataAgeSubtypeViolinPlot)
ggsave(filename = file.path(phenotype_validation_path, "BMI_Violin_Plot.png"), plot = sampleDataBMISubtypeViolinPlot)
ggsave(filename = file.path(phenotype_validation_path, "BMI_Box_Plot.png"), plot = sampleDataBMISubtypeBoxPlot)
ggsave(filename = file.path(phenotype_validation_path, "Diabetes_Bar_Plot.png"), plot = sampleDataDiabetesSubtypeBarPlot)
ggsave(filename = file.path(phenotype_validation_path, "HW_Violin_Plot.png"), plot = sampleDataHWSubtypeViolinPlot)
ggsave(filename = file.path(phenotype_validation_path, "TIN_Violin_Plot.png"), plot = sampleDataTINSubtypeViolinPlot)
ggsave(filename = file.path(phenotype_validation_path, "RIN_Violin_plot.png"), plot = sampleDataRINSubtypeViolinPlot)
ggsave(filename = file.path(phenotype_validation_path, "Library_Pool_Bar_Plot.png"), plot = sampleDataLPSubtypeBarPlot)

#-----------------------------------------------------------------------------#
# COMPLETE
#-----------------------------------------------------------------------------#

message("--- Finished Phenotype Validation  ---")