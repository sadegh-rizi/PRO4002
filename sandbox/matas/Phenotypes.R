lvef_vs_endo_bar <- ggplot(dcm_meta, aes(x = Endotype, y = LVEF, fill = Endotype)) +
  geom_bar(stat = "identity") +
  stat_compare_means(method = "kruskal.test", label.y = 50) + 
  labs(title = "Is Cluster 4 the 'Severe' Group?", y = "LVEF (%)") +
  theme_minimal() + 
  scale_fill_jco() + 
  scale_y_continuous(limits = c(0, 20))
lvef_vs_endo_bar

lvef_vs_endo_violin <- ggplot(dcm_meta, aes(x = Endotype, y = LVEF, fill = Endotype)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.5) + # Add dots to see the 9 patients
  stat_compare_means(method = "kruskal.test", label.y = 50) + 
  labs(title = "Is Cluster 4 the 'Severe' Group?", y = "LVEF (%)") +
  theme_minimal() +
  scale_fill_jco() + 
  scale_y_continuous(limits = c(0, 1))
lvef_vs_endo_violin

age_vs_endo_violin <- ggplot(dcm_meta, aes(x = Endotype, y = age, fill = Endotype)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(method = "kruskal.test", label.y = 90) +
  labs(title = "Violin plot of Age Distribution", y = "Age (Years)") +
  theme_minimal() +
  scale_fill_jco()
age_vs_endo_violin

dcm_meta$BMI <- round(dcm_meta$weight / ((dcm_meta$height) / 100) ^ 2, 1)

bmi_vs_endo_violin <- ggplot(dcm_meta, aes(x = Endotype, y = BMI, fill = Endotype)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(method = "kruskal.test", label.y = 90) +
  labs(title = "Violin plot of BMI Distribution", y = "BMI - Weight(kg) / Height(m)^2") +
  theme_minimal() +
  scale_fill_jco() + 
  # scale_y_continuous(limits = c(0, 100)) 
  # with scale_y_continuous it excludes the outliers of BMI, which greatly increases the p value 
  # from kruskal test. This means the outliers are impactful data, enough to alter the p value
  bmi_vs_endo_violin

bmi_vs_endo_box <- ggplot(dcm_meta, aes(x = Endotype, y = BMI, fill = Endotype)) +
  geom_boxplot(outlier.shape = 2) + # "2" is used as the shape - triangle to visualise the outliers
  geom_jitter(width = 0.2, alpha = 0.5) + # Add dots to see the 9 patients
  stat_compare_means(method = "kruskal.test", label.y = 50) + 
  labs(title = "Box plot of BMI Distribution", y = "BMI - Weight(kg) / Height(m)^2") +
  theme_minimal() +
  scale_fill_jco()
bmi_vs_endo_box # not such a good visual

diabetes_vs_endo_bar <- ggplot(dcm_meta, aes(x = Endotype, fill = Diabetes)) +
  geom_bar(position = "dodge") + 
  labs(title = "Bar plot of Diabetes for each Endotype", y = "Diabetes", fill = "Diabetes") +
  theme_minimal()
diabetes_vs_endo_bar

hw_vs_endo_violin <- ggplot(dcm_meta, aes(x = Endotype, y = hw, fill = Endotype)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(method = "kruskal.test", label.y = 90) +
  labs(title = "Violin plot of Heart Weight Distribution", y = "Heart Weight(g)") +
  theme_minimal() +
  scale_fill_jco()
hw_vs_endo_violin

tin_vs_endo_violin <- ggplot(dcm_meta, aes(x = Endotype, y = TIN.median., fill = Endotype)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(method = "kruskal.test", label.y = 90) +
  labs(title = "Violin plot of Transcript Integrity Number (TIN) Distribution", y = "TIN(median)") +
  theme_minimal() +
  scale_fill_jco()
tin_vs_endo_violin

rin_vs_endo_violin <- ggplot(dcm_meta, aes(x = Endotype, y = RIN, fill = Endotype)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(method = "kruskal.test", label.y = 90) +
  labs(title = "Violin plot of RNA Integrity Number (RIN) Distribution", y = "RIN") +
  theme_minimal() +
  scale_fill_jco() + 
  scale_y_continuous(limits = c(0, 20))
rin_vs_endo_violin

library_pool_vs_endo_bar <- ggplot(dcm_meta, aes(x = Endotype, fill = Library.Pool)) +
  geom_bar(position = "dodge") + 
  labs(title = "Bar plot of Library Pools for each Endotype", y = "Samples", fill = "Library.Pool") +
  theme_minimal()
library_pool_vs_endo_bar # heatmap will be better, but for validation

ggsave(filename = "LVEF_Bar_plot.pdf", 
       plot = lvef_vs_endo_bar,
       path = "/Users/matasm/Downloads/Ms System Bio/Maastricht Uni/Y1/P3/Research Project 1/Phenotype Results")
ggsave(filename = "LVEF_Violin_plot.pdf", 
       plot = lvef_vs_endo_violin,
       path = "/Users/matasm/Downloads/Ms System Bio/Maastricht Uni/Y1/P3/Research Project 1/Phenotype Results")
ggsave(filename = "Age_Violin_plot.pdf", 
       plot = age_vs_endo_violin,
       path = "/Users/matasm/Downloads/Ms System Bio/Maastricht Uni/Y1/P3/Research Project 1/Phenotype Results")
ggsave(filename = "BMI_Violin_plot.pdf", 
       plot = bmi_vs_endo_violin,
       path = "/Users/matasm/Downloads/Ms System Bio/Maastricht Uni/Y1/P3/Research Project 1/Phenotype Results")
ggsave(filename = "Diabetes_Bar_plot.pdf", 
       plot = diabetes_vs_endo_bar,
       path = "/Users/matasm/Downloads/Ms System Bio/Maastricht Uni/Y1/P3/Research Project 1/Phenotype Results")
ggsave(filename = "HW_Violin_plot.pdf", 
       plot = hw_vs_endo_violin,
       path = "/Users/matasm/Downloads/Ms System Bio/Maastricht Uni/Y1/P3/Research Project 1/Phenotype Results")
ggsave(filename = "TIN_Violin_plot.pdf", 
       plot = tin_vs_endo_violin,
       path = "/Users/matasm/Downloads/Ms System Bio/Maastricht Uni/Y1/P3/Research Project 1/Phenotype Results")
ggsave(filename = "RIN_Violin_plot.pdf", 
       plot = rin_vs_endo_violin,
       path = "/Users/matasm/Downloads/Ms System Bio/Maastricht Uni/Y1/P3/Research Project 1/Phenotype Results")
ggsave(filename = "Library_Pool_Bar_plot.pdf", 
       plot = library_pool_vs_endo_bar,
       path = "/Users/matasm/Downloads/Ms System Bio/Maastricht Uni/Y1/P3/Research Project 1/Phenotype Results")