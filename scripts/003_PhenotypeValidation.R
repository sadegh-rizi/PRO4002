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

sampleData.DCM.complete <- sampleData.DCM %>%
  filter(height > 100)

rownames(sampleData.DCM.complete) <- sampleData.DCM.complete$sample_name
sampleData.DCM.complete$subtype <- sampleData.DCM.annotated$subtype[match(sampleData.DCM.complete$sample_name, 
  rownames(sampleData.DCM.annotated))]
sampleData.DCM.complete <- sampleData.DCM.complete %>% 
  filter(!is.na(subtype))

sampleData.DCM.complete$rin[is.na(sampleData.DCM.complete$rin)] <- median(sampleData.DCM.complete$rin, na.rm = TRUE)

#-----------------------------------------------------------------------------#
# DATA VISUALIZATION
#-----------------------------------------------------------------------------#

message("(II) Data Visualization ")

#-----------------------------------------------------------------------------#
# DATA VISUALIZATION - Disease vs Subtype 
#-----------------------------------------------------------------------------#

table_data <- table(sampleData.DCM.complete$subtype, sampleData.DCM.complete$Diabetes)
chi_test <- chisq.test(table_data)
chi_test$p.value
p_label <- ifelse(chi_test$p.value < 0.001, "p < 0.001", paste0("p = ", signif(chi_test$p.value, 3)))

sampleDataDiabetesSubtypeBarPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, fill = Diabetes, color = Diabetes)) +
  geom_bar(position = "dodge", alpha=0.7) + 
  labs(title = "Diabetes", x ="", y = "Count", fill = "Subtype", color = "Subtype") +
  scale_fill_manual(values = c(standard_color, accent_color)) +
  scale_color_manual(values = c(standard_color, accent_color)) +
  annotate("text", x = 2, y = 70, label = p_label, size = 5) +
  scale_y_continuous(limits = c(0, 75)) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataDiabetesSubtypeBarPlot.png"), plot = sampleDataDiabetesSubtypeBarPlot)

table_data <- table(sampleData.DCM.complete$subtype, sampleData.DCM.complete$afib)
chi_test <- chisq.test(table_data)
chi_test$p.value
p_label <- ifelse(chi_test$p.value < 0.001, "p < 0.001", paste0("p = ", signif(chi_test$p.value, 3)))

sampleDataAFibSubtypeBarPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, fill = afib, color = afib)) +
  geom_bar(position = "dodge", alpha=0.7) + 
  labs(title = "Atrial Fibrillation", x ="", y = "Count", fill = "Subtype", color = "Subtype") +
  scale_fill_manual(values = c(standard_color, accent_color)) +
  scale_color_manual(values = c(standard_color, accent_color)) +
  annotate("text", x = 2, y = 70, label = p_label, size = 5) +
  scale_y_continuous(limits = c(0, 75)) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataAFibSubtypeBarPlot.png"), plot = sampleDataAFibSubtypeBarPlot)

table_data <- table(sampleData.DCM.complete$subtype, sampleData.DCM.complete$VTVF)
chi_test <- chisq.test(table_data)
chi_test$p.value
p_label <- ifelse(chi_test$p.value < 0.001, "p < 0.001", paste0("p = ", signif(chi_test$p.value, 3)))

sampleDataVTVFSubtypeBarPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, fill = VTVF, color = VTVF)) +
  geom_bar(position = "dodge", alpha=0.7) + 
  labs(title = "VTVF", x ="", y = "Count", fill = "Subtype", color = "Subtype") +
  scale_fill_manual(values = c(standard_color, accent_color)) +
  scale_color_manual(values = c(standard_color, accent_color)) +
  annotate("text", x = 2, y = 70, label = p_label, size = 5) +
  scale_y_continuous(limits = c(0, 75)) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataVTVFSubtypeBarPlot.png"), plot = sampleDataVTVFSubtypeBarPlot)

table_data <- table(sampleData.DCM.complete$subtype, sampleData.DCM.complete$Hypertension)
chi_test <- chisq.test(table_data)
chi_test$p.value
p_label <- ifelse(chi_test$p.value < 0.001, "p < 0.001", paste0("p = ", signif(chi_test$p.value, 3)))

sampleDataHypertensionSubtypeBarPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, fill = Hypertension, color = Hypertension)) +
  geom_bar(position = "dodge", alpha=0.7) + 
  labs(title = "Hypertension", x ="", y = "Count", fill = "Subtype", color = "Subtype") +
  scale_fill_manual(values = c(standard_color, accent_color)) +
  scale_color_manual(values = c(standard_color, accent_color)) +
  annotate("text", x = 2, y = 70, label = p_label, size = 5) +
  scale_y_continuous(limits = c(0, 75)) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataHypertensionSubtypeBarPlot.png"), plot = sampleDataHypertensionSubtypeBarPlot)

# sampleDataDiabetesSubtypeBarPlot2 <- ggplot(sampleData.DCM.complete, aes(x = subtype, fill = subtype, color=subtype)) +
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

sampleDataLVEFSubtypeViolinPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, y = (lvef*100), fill = subtype, color = subtype)) +
  geom_violin(alpha= 0.7) +
  geom_jitter(width = 0.2, color="black", size = 1) + 
   stat_compare_means(method = "kruskal.test", label.y = 95) +
  labs(title = "", y = "LVEF (%)", fill = "Subtype", color = "Subtype", x = "") +
  scale_fill_manual(values = npg_colors[1:3]) +
  scale_color_manual(values = npg_colors[1:3]) +
  scale_y_continuous(limits = c(0, 100)) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataLVEFSubtypeViolinPlot.png"), plot = sampleDataLVEFSubtypeViolinPlot)

sampleDataAgeSubtypeViolinPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, y = age, fill = subtype, color = subtype)) +
  geom_violin(alpha= 0.7) +
  geom_jitter(width = 0.2, color="black", size = 1) + 
   stat_compare_means(method = "kruskal.test", label.y = 75) +
  labs(title = "", y = "Age (Years)", fill = "Subtype", color = "Subtype", x = "") +
  scale_fill_manual(values = npg_colors[1:3]) +
  scale_color_manual(values = npg_colors[1:3]) +
  scale_y_continuous(limits = c(15, 80)) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataAgeSubtypeViolinPlot.png"), plot = sampleDataAgeSubtypeViolinPlot)

sampleDataBMISubtypeViolinPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, y = bmi, fill = subtype, color = subtype)) +
  geom_violin(alpha= 0.7) +
  geom_jitter(width = 0.2, color="black", size = 1) + 
   stat_compare_means(method = "kruskal.test", label.y = 75) +
  labs(title = "", y = "BMI (kg/m²)", fill = "Subtype", color = "Subtype", x = "") +
  scale_fill_manual(values = npg_colors[1:3]) +
  scale_color_manual(values = npg_colors[1:3]) +
  scale_y_continuous(limits = c(15, 80)) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataBMISubtypeViolinPlot.png"), plot = sampleDataBMISubtypeViolinPlot)

sampleDataTINSubtypeViolinPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, y = tin_median, fill = subtype, color = subtype)) +
  geom_violin(alpha= 0.7) +
  geom_jitter(width = 0.2, color="black", size = 1) + 
  stat_compare_means(method = "kruskal.test", label.y = 85) +
  labs(title = "", y = "TIN Median", fill = "Subtype", color = "Subtype", x = "") +
  scale_fill_manual(values = npg_colors[1:3]) +
  scale_color_manual(values = npg_colors[1:3]) +
  scale_y_continuous(limits = c(20, 90)) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataTINSubtypeViolinPlot.png"), plot = sampleDataTINSubtypeViolinPlot)

sampleDataRINSubtypeViolinPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, y = rin, fill = subtype, color = subtype)) +
  geom_violin(alpha= 0.7) +
  geom_jitter(width = 0.2, color="black", size = 1) + 
  stat_compare_means(method = "kruskal.test", label.y = 10.5) +
  labs(title = "", y = "RIN", fill = "Subtype", color = "Subtype", x = "") +
  scale_fill_manual(values = npg_colors[1:3]) +
  scale_color_manual(values = npg_colors[1:3]) + 
  scale_y_continuous(limits = c(5, 11)) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataRINSubtypeViolinPlot.png"), plot = sampleDataRINSubtypeViolinPlot)

sampleDataHWSubtypeViolinPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, y = hw, fill = subtype, color = subtype)) +
  geom_violin(alpha= 0.7) +
  geom_jitter(width = 0.2, color="black", size = 1) + 
  labs(title = "", y = "Heart Weight (g)", fill = "Subtype", color = "Subtype", x = "") +
  stat_compare_means(method = "kruskal.test", label.y = 90) +
  scale_fill_manual(values = npg_colors[1:3]) +
  scale_color_manual(values = npg_colors[1:3]) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataHWSubtypeViolinPlot.png"), plot = sampleDataHWSubtypeViolinPlot)

sampleDataLVSubtypeViolinPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, y = lv_mass, fill = subtype, color = subtype)) +
  geom_violin(alpha= 0.7) +
  geom_jitter(width = 0.2, color="black", size = 1) + 
  labs(title = "", y = "LV Mass (g)", fill = "Subtype", color = "Subtype", x = "") +
  stat_compare_means(method = "kruskal.test", label.y = 150) +
  scale_fill_manual(values = npg_colors[1:3]) +
  scale_color_manual(values = npg_colors[1:3]) + 
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataLVSubtypeViolinPlot.png"), plot = sampleDataLVSubtypeViolinPlot)

sampleDataLPSubtypeBarPlot <- ggplot(sampleData.DCM.complete, aes(x = subtype, fill = Library.Pool, color = Library.Pool)) +
  geom_bar(position = "dodge", alpha = 0.7) + 
  labs(title = "Bar plot of Library Pools for each Endotype", x = "", y = "Samples", fill = "Library.Pool") +
  scale_fill_manual(values = tail(npg_additional_colors,-3)) +
  scale_color_manual(values = tail(npg_additional_colors,-3)) +
  center_title + my_style
ggsave(filename = file.path(phenotype_validation_path, "sampleDataLPSubtypeBarPlot.png"), plot = sampleDataLPSubtypeBarPlot)

#-----------------------------------------------------------------------------#
# COMPLETE
#-----------------------------------------------------------------------------#

message("--- Finished Phenotype Validation  ---")