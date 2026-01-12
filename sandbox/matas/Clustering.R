library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(gtsummary)
library(pcaMethods)
library(gridExtra)
library(biomaRt)
library(tibble)
library(janitor)
library(limma)
library(BiocManager)
library(conflicted)

setwd("/Users/matasm/Downloads/Ms System Bio/Maastricht Uni/Y1/P2/Experimental Design & Data Management/Assignment/MAGNET_GX_Assignment_MSB1005")
sample_data <- read.csv("MAGNET_SampleData_18112022.csv", row.names = 1)
gene_exp_data <- read.delim("MAGNET_GeneExpressionData_CPM_19112020.txt", row.names = 1)
exon_length <- read.delim("MAGNET_exonLengths.txt", row.names = 1)

DCM_samples_exp <- sample_data %>%
  dplyr::select(etiology, gender, race, age, weight, height, hw, LVEF, TIN.median.) %>%
  mutate(race = case_match(race, "AA" ~ "African American", .default = race)) %>%
  filter(etiology == "DCM") %>%
  group_by(gender) %>%
  tbl_summary(by = gender, label = list(etiology = "Etiology", 
                                        gender = "Gender", 
                                        race = "Race",
                                        age = "Age",
                                     weight = "Weight (kg)", 
                                     height = "Height (cm)", 
                                     hw = "Heart Weight (g)")) %>%
  modify_header(label = "**DCM Characteristics**") %>%
  bold_labels()
DCM_samples_exp

DCM_samples <- sample_data %>%
  filter(etiology == "DCM")
DCM_samples

samples_numeric <- sample_data[sapply(sample_data, is.numeric)]
samples_pca <- pca((samples_numeric), nPcs = 10)
samples_plot_pca <- cbind(data.frame(samples_pca@scores), sample_data)
ggplot(samples_plot_pca, aes(x = PC1, y = PC2)) + geom_point(aes(size = sample_data$age, col = sample_data$gender, alpha = 0.4)) + 
  facet_wrap(~ sample_data$gender, scales = "fixed")

DCM_numeric <- DCM_samples[sapply(DCM_samples, is.numeric)]
DCM_pca <- pca((DCM_numeric), nPcs = 10)
DCM_plot_pca <- cbind(data.frame(DCM_pca@scores), DCM_samples)

ggplot(DCM_plot_pca, aes(x = PC1, y = PC2)) + geom_point(aes(size = DCM_samples$age, col = DCM_samples$gender, alpha = 0.4)) + 
  facet_wrap(~ DCM_samples$gender, scales = "fixed")

ggplot(DCM_plot_pca, aes(x = PC1, y = PC2)) + geom_violin(aes(col = gender)) + 
  facet_wrap(~ gender)

# filter samples and create data frames
dcm_samples <- rownames(sample_data[sample_data$etiology == "DCM", ])
gene_exp_dcm <- gene_exp_data[, dcm_samples]
data_dcm <- sample_data[dcm_samples, ]

# normalise and filter the data, log2 transform data to compress it, filters out highly expressed genes
# +1 removes log0
gene_exp_log <- log2(gene_exp_dcm + 1)
# gene variance - reduce noise, keeps genes that only vary and not flat values, 0.5 is imperical
gene_exp_var <- gene_exp_log[apply(gene_exp_log, 1, var) > 0.5, ]
# scale the genes so there is no gene dominating euclidean distance, scale only works on columns
gene_exp_scaled <- t(scale(t(gene_exp_var)))

#explore the data now with pca
gene_pca <- prcomp(t(gene_exp_scaled), center = FALSE, scale. = FALSE)
gene_pca_data_frame <- data.frame(PC1 = gene_pca$x[, 1], PC2 = gene_pca$x[, 2], data_dcm)
ggplot(gene_pca_data_frame, aes(x = PC1, y = PC2)) + geom_point(aes(color = gender), size = 2) + theme_classic()

# it seems the samples separate by gender

# umap can be used to reduce the data dimension so it is better to cluster lower volume data
library(umap)
set.seed(500)
umap_dim <- umap(t(gene_exp_scaled))
umap_data_frame <- data.frame(UMAP1 = umap_dim$layout[, 1], UMAP2 = umap_dim$layout[, 2])
plot(umap_data_frame, pch = 16)

library(ConsensusClusterPlus)
# repeated cluster to see if there are any patterns emerging, if the same samples cluster together each time = structure
cluster_res_k_6 <- ConsensusClusterPlus(as.matrix(gene_exp_scaled), maxK = 6, reps = 1000, pItem = 0.8, pFeature = 1, clusterAlg = "pam", distance = "euclidean", seed = 1)
# pam - partitioning around mediods, stronger than k means, can also used k means, pFeature - uses all genes, pItem - each iteration use 80% of samples for robustness,
# if cluster does not work when 20% are removed then something is wrong

# in case of batches appearing, before umap
# library(sva)
# gene_exp_bc <- ComBat(dat = gene_exp_scaled, batch = data_dcm$gender, mod = NULL)
# umap_res <- umap(t(gene_exp_bc))
# umap_df <- data.frame(UMAP1 = umap_res$layout[, 1], UMAP2 = umap_res$layout[, 2])
# plot(umap_df, pch = 16)

cluster_res_k_3 <- ConsensusClusterPlus(as.matrix(gene_exp_scaled), maxK = 3, reps = 1000, pItem = 0.8, pFeature = 1, clusterAlg = "pam", distance = "euclidean", seed = 1)
# I believe not enough info provided for k = 3

# create the subtypes found and assign them to variable subtype, then add to main data frames
# choosing k = 4, seems the best since 5, and 6 start to create overlapping clusters of subtypes
subtype <- cluster_res_k_6[[4]]$consensusClass
data_dcm$Subtype <- factor(subtype)
umap_data_frame$Subtype <- data_dcm$Subtype
ggplot(umap_data_frame, aes(x = UMAP1, y = UMAP2, color = Subtype)) + geom_point(size = 2) + theme_classic()

ggplot(data_dcm, aes(Subtype, LVEF, color = Subtype)) + geom_boxplot() + geom_jitter(width = 0.3, alpha = 0.3) + theme_classic()

# since data is not necessarily producing a specific distribution we can use kruskal.test or wilcox.test
#   kruskal is better for more than 2 independent groups to test, while wilcoxon better for 2 groups (paired)
# otherwise if data follow a certain distribution say normal, then we can use ANOVA test

kruskal.test(LVEF ~ Subtype, data = data_dcm) # check formula = response ~ group/explanatory
# p-val = 0.9416, would claim the null hypothesis, it is not less than 0.05

# extract categorical data and make a table, chisq.test specifically test categorical variables, to see if there is any significance between them
table(data_dcm$Subtype, data_dcm$gender)
chisq.test(data_dcm$Subtype, data_dcm$gender)
# p-val is extremely small, reject null hypothesis

# DEA, see the contrasts between the subtypes
design1 <- model.matrix(~ 0 + Subtype, data_dcm)
colnames(design1)[1:4] <- c("S1", "S2", "S3", "S4")
design1
fit1 <- lmFit(gene_exp_scaled, design1)
con_mat <- makeContrasts(Subtype_1_vs_Subtype_2 = S1 - S2,
                         Subtype_1_vs_Subtype_3 = S1 - S3,
                         Subtype_1_vs_Subtype_4 = S1 - S4,
                         Subtype_2_vs_Subtype_3 = S2 - S3,
                         Subtype_2_vs_Subtype_4 = S2 - S4,
                         Subtype_3_vs_Subtype_4 = S3 - S4,
                         levels = design1)
fit2 <- contrasts.fit(fit1, con_mat)
fit_Bayes <- eBayes(fit2, trend = TRUE)
# extract genes
# BH by default
gene_res_S1_S2 <- topTable(fit_Bayes, coef = 'Subtype_1_vs_Subtype_2', number = Inf, adjust = "BH")
gene_res_S1_S3 <- topTable(fit_Bayes, coef = 'Subtype_1_vs_Subtype_3', number = Inf, adjust = "BH")
gene_res_S1_S4 <- topTable(fit_Bayes, coef = 'Subtype_1_vs_Subtype_4', number = Inf, adjust = "BH")
gene_res_S2_S3 <- topTable(fit_Bayes, coef = 'Subtype_2_vs_Subtype_3', number = Inf, adjust = "BH")
gene_res_S2_S4 <- topTable(fit_Bayes, coef = 'Subtype_2_vs_Subtype_4', number = Inf, adjust = "BH")
gene_res_S3_S4 <- topTable(fit_Bayes, coef = 'Subtype_3_vs_Subtype_4', number = Inf, adjust = "BH")

ggplot(gene_res_S1_S2, aes(x = logFC, y = -log10(adj.P.Val), color = adj.P.Val < 0.05)) + geom_point() + geom_smooth()
ggplot(gene_res_S1_S3, aes(x = logFC, y = -log10(adj.P.Val), color = adj.P.Val < 0.05)) + geom_point() + geom_smooth()

library(pheatmap)
pheatmap(gene_exp_scaled, annotation_col = data_dcm["Subtype"])



