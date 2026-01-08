#=============================================================================#
# RNAseq.R                                          #
#																	                                       		  #
# Version: 1.0   															                                #
# Date: Dec 12, 2025											                                    #
# Author: Sadegh Rizi             #
# History:																	                                  #
#  1.0: Creation
#  1.1: Update the publication ready table with gt package
#  1.2: Fix bugs
#  1.3: Create a function for PCA plots to remove redundancy
#  1.4: Reorder steps so that filtering is before statistical analysis
# Description: Differential Expression Analysis of MAGNET dataset

#=============================================================================#


# This script serves as a scaffold for RNA-seq anaylsis of MAGNET dataset. 

# MAGNET dataset overview: 
# WORKFLOW 

# step0: library loading
# step1: Data import
# step2: Diagnostic plots
# step3: Statistical Analysis
# step4: Gene Annotation
# step5: Relative expression level
# step6: Export the results


### Step 0: Importing Libraries and setting working directories
# Before that we check for package installations

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bio_pkgs <- c("limma", "edgeR", "qvalue", "biomaRt", "pcaMethods","EnhancedVolcano")

for (pkg in bio_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    message(paste("Installing Bioconductor package:", pkg))
    BiocManager::install(pkg, ask = FALSE)
    library(pkg, character.only = TRUE)
  }
}

cran_pkgs <- c("tidyverse", "gt", "gtsummary", "here")

for (pkg in cran_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    message(paste("Installing CRAN package:", pkg))
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

if (!require("here", quietly = TRUE)) install.packages("here")
library(here)

#Load necessary packages
library(tidyverse) # Covers ggplot2, dplyr, readr, tidyr

# For differential analysis
library(limma)
library(edgeR)
library(qvalue)

# For gene annotation
library(biomaRt)


# Libraries for Exporting participants table 
#install.packages("gt")
#install.packages("webshot2")

# For exporting publication-ready tables
library(gt)
library(gtsummary)


# PCA
library(pcaMethods)


# Set working directory to the location of this script
setwd(here()) 

# Verify working directory
print(paste("Working directory set to:", getwd()))



if (!dir.exists("plots")) {
  dir.create("plots")
}

if (!dir.exists("tables")) {
  dir.create("tables")
}





# GGplot Plotting Settings 
center_title <- theme(plot.title = element_text(hjust = 0.5, vjust = 1))

# To harmonize font and style
center_title <- theme(plot.title = element_text(hjust = 0.5, vjust = 1))

# To harmonize font and style
target_font <- "sans" 
my_style <- theme(
  text = element_text(family = target_font),        # Forces the font
  plot.title = element_text(hjust = 0.5, face="bold"), # Centers & bolds title
  axis.title = element_text(face = "bold"),         # Bolds axis titles
  legend.position = "right"                         # Default legend position
)


### Step 1: Data import 
## 1.a : Importing the data
expression_input_file <- "data/MAGNET_GeneExpressionData_CPM_19112020.txt"
sample_input_file <- "data/MAGNET_SampleData_18112022.csv"
gene_length_input <- "data/MAGNET_exonLengths.txt"


# Control check to see if files exist in the data directory
#otherwise stop
for (file_path in c(expression_input_file,sample_input_file,gene_length_input)){
  if (!file.exists(file_path)) stop
}

expression <- as.matrix(read.delim(expression_input_file,
                                   row.names=1))
sampleInfo <- read.csv(sample_input_file, row.names=1)
gene_lengths <- read.delim(gene_length_input, row.names=1)

# Sanity check to see if rows of sample metadata table match the columns of gene
# expression table
(all(colnames(expression) == rownames(sampleInfo))) #TRUE, just to make sure


## 1.b : Exporting a publication-ready table of patient characteristics

sampleInfo <-sampleInfo %>% mutate(bmi= weight/(height*0.01)^2)


# For race change AA to African American to make it more readable
sampleInfo<- sampleInfo %>% mutate(race=case_match(race, "AA" ~ "African American",.default=race))
sampleInfo

# We use the gtsummary package, that automatically reads which variables are continuous and which ones are categorical and even performs statistical tests on them
#Here non-parametric tests have been chosen. (Fisher's exact test for numerical and Kruskal wallis for categorical)
# This could also be done by dplyr select, mutate and summarize functions

participant_table <- sampleInfo %>% 
  select(everything()) %>%
  mutate(etiology = factor(etiology, 
                           levels = c("NF", "PPCM", "HCM", "DCM"),
                           labels = c("Non-Failing (NF)", 
                                      "Peripartum Cardiomyopathy (PPCM)", 
                                      "Hypertrophic CM (HCM)", 
                                      "Dilated CM (DCM)"))) %>%
  tbl_summary(
    by= etiology,
    
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),
    
    label= list(
      age ~ "Age (years)",
      gender ~ "Sex",
      race ~ "Ethnicity",
      weight ~ "Weight (kg)",
      height ~ "Height (Cm)",
      bmi ~ "BMI (kg/m2)",
      Diabetes ~ "Diabetes Mellitus"
    ),
    missing = "ifany" 
  )  %>% 
  add_p(
    pvalue_fun = ~style_pvalue(.x,digits=3) 
    
  ) %>%
  add_overall() %>%
  modify_caption("**Table 1. Participant characteristics by etiology**") %>% 
  bold_labels()

participant_table %>%
  as_gt() %>%             # Convert gtsummary object to gt object
  gtsave(filename = "tables/Table1-Patientcharacteristics.pdf")  # Save to your working directory



### Step 2: Diagnostic plots

## a. At least one data distribution figure (e.g. boxplots, density plots) that enables 
# comparing samples.
#here we want to see if the disribution of gene expression for each sample is 
#comparable to the other samples

# Make the expression matrix into a long format, so it can be merged with 
# sample metadata and ease of plotting
expression_long <- expression %>% as.data.frame() %>% pivot_longer(cols=-1,names_to = "sample_name",
                                                                   values_to="expression"
)

sampleInfo$sample_name <- rownames(sampleInfo)
#merging sample metadata with gene expression matrix (long-format)
merged_df <- expression_long %>% left_join(sampleInfo,by="sample_name")

# ggplot for ploting the boxplots of each sample, grouped by etiology
sample_dist_boxplot <- ggplot(merged_df,aes(x=sample_name,y=expression,fill=etiology))+
  geom_boxplot()+
  facet_grid(
    ~ etiology,   # Creates columns based on etiology
    scales = "free_x",
    space="free_x"# Each panel only shows its *own* samples
    # Panel width is proportional to sample count
  ) +
  
  labs(
    title = "Gene Expression Distribution by Sample",
    x = "Sample",
    y = "Gene Expression Level (logCPM)"
  ) + 
  theme_minimal()+
  theme(
    axis.text.x = element_blank(), # Still hide sample names
    axis.ticks.x = element_blank(),
    legend.position = "none" ,
    strip.text.x = element_text(angle = 90, hjust = 0.5) ) +
  center_title + my_style

sample_dist_boxplot

ggsave("plots/expression_boxplots_across_samples.pdf", plot = sample_dist_boxplot, 
       width = 20, height = 8)
ggsave("plots/expression_boxplots_across_samples.png", plot = sample_dist_boxplot,
       width = 20, height = 8)





dist_age <- ggplot(sampleInfo, aes(x = etiology, y = age, fill = etiology)) +
  geom_violin() +
  geom_jitter(width=0.2, size=1, alpha=0.5)+ 
  labs(title = "Distribution of sample ages") + 
  center_title + my_style
dist_age
ggsave("plots/violin_age_by_etiology.pdf", plot = dist_age, width = 6, height = 4)
ggsave("plots/violin_age_by_etiology.png", plot = dist_age, width = 6, height = 4)


# Check if RIN (Quality) differs
dist_rin <- ggplot(sampleInfo, aes(x = etiology, y = RIN, fill = etiology)) +
  geom_violin() +
  geom_jitter(width=0.2, size=1, alpha=0.5)+ 
  
  labs(title = "Distribution of RIN (RNA integrity)")+
  center_title + my_style
dist_rin
ggsave("plots/violin_RIN_by_etiology.pdf", plot = dist_rin, width = 6, height = 4)
ggsave("plots/violin_RIN_by_etiology.png", plot = dist_rin, width = 6, height = 4)


## b. At least one PCA figure showing the sample clustering colored by relevant co
#variates.  



# Perform pca using pcaMethods package, and using svd method(singlular-value decomposition )
pca_results <- pca(t(expression),method="svd")
# save the scores
pca_scores <- as.data.frame(scores(pca_results))
var_expl <- round(pca_results@R2 * 100, 1) # Extract R2 slot for variance

colnames(pca_scores) <- c("PC1","PC2")
# merge scores with sample metadata
plotting_data <- cbind(pca_scores,sampleInfo)


# define a theme for PCA plots to remove redundancy
pca_theme <-   theme_bw()+
  theme(panel.grid = element_blank())+ # to remove the rectangles inside
  center_title # center the title

plot_pca <- function(data, color_col,v1,v2) {
    ggplot(data, aes_string(x = "PC1", y = "PC2", color = color_col))+
    geom_point(size = 3, alpha = 0.8) +
    scale_color_brewer(palette = "Set1") + # Nice distinct colors
    labs(
      title = paste("PCA by", str_to_title(color_col)),
      x = paste0("PC1 (", v1, "% Var)"),
      y = paste0("PC2 (", v2, "% Var)")
    ) + my_style
}


pca_gender <- plot_pca(data=plotting_data,color_col="gender",var_expl[1], var_expl[2])

pca_gender 
ggsave("plots/PCA_gender.pdf", plot = pca_gender, width = 10, height = 6)
ggsave("plots/PCA_gender.png", plot = pca_gender, width = 10, height = 6)



pca_race <- plot_pca(data=plotting_data,color_col="race",var_expl[1], var_expl[2])
pca_race
ggsave("plots/PCA_race.pdf", plot = pca_race, width = 10, height = 6)
ggsave("plots/PCA_race.png", plot = pca_race, width = 10, height = 6)



pca_etiology <- plot_pca(data=plotting_data,color_col="etiology",var_expl[1], var_expl[2])
pca_etiology
ggsave("plots/PCA_etiology.pdf", plot = pca_etiology, width = 6, height = 4)
ggsave("plots/PCA_etiology.png", plot = pca_etiology, width = 6, height = 4)


## Step 3: Gene annotation 
# This step is necessary before filtering the gene expressions
# Because we need to know which genes are on Y chromosome


# use homo sapiens gene ensembl dataset
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
ensembl_ids <- rownames(expression)


# mapping: INPUT: ensembl_gene_id
# OUTPUT: gene symbol, location of gene, and description
gene_map <- getBM(attributes = c("ensembl_gene_id","external_gene_name","description","chromosome_name", 
                                 "start_position", 
                                 "end_position"),
                  filters = "ensembl_gene_id",
                  values=ensembl_ids,
                  mart=ensembl)
head(gene_map)

# change the column  name of to gene_symbol for easier readability
gene_map<-gene_map %>% rename("gene_symbol"="external_gene_name")








## Step 4: relative expression levels

# This is performed before Step 3 (which is the statistical analysis)! 
# Before doing the differential expression we need to filter the genes
# that we assume are not expressed.
# We use genes expressed in chromsome Y of female patients as the noise or 
# background value


all(rownames(gene_lengths) == rownames(expression)) # TRUE (just a check)


# Function to convert logCPM to fpkm value, from skill sessions
cpm2fpkm <- function(x) {
  .t <- 2^(x) * 1E3 / gene_lengths[, 1]
}


expression_fpkm <- cpm2fpkm(expression)

# To define baseline value, we take Y chromosome gene expressions in female samples

female_samples <- sampleInfo %>% filter(gender=="Female") %>% pull(sample_name)
ychromosome_genes<- gene_map %>% filter(chromosome_name=="Y") %>% pull(ensembl_gene_id)
expression_noise<- expression_fpkm[ychromosome_genes,female_samples]
background_threshold <- mean(expression_noise,na.rm=TRUE)

print(paste("Background Noise Threshold (FPKM):", round(background_threshold, 4)))

# Now if the average expression of a gene is lower than this baseline we consider
# it noise. 



gene_means <- rowMeans(expression_fpkm)
is_expressed <- gene_means>background_threshold
expressed_filtered <- expression_fpkm[is_expressed,]


mean_expression_dataframe <- data.frame(mean_expression=gene_means)
avg_exp_plot<- ggplot(mean_expression_dataframe,
                      aes(x=mean_expression))+
  geom_histogram(aes(fill=mean_expression>background_threshold),
                 bins=50,color="white",alpha=0.7)+
  geom_vline(xintercept = background_threshold,color="red",size=1,linetype="dashed")+
  #use log scale because FPKM is not log transformed
  scale_x_log10()+
  scale_fill_manual(values = c("gray", "cornflowerblue"), 
                    labels = c("Noise", "Expressed"), 
                    name = "Status") +
  labs(title = "Distribution of Average expressions vs Background Threshold",
       subtitle = paste("Threshold:", round(background_threshold, 4), "FPKM"),
       x = "Mean Expression (FPKM)",
       y = "Count of Genes") +
  theme_minimal() + my_style

ggsave("plots/expression_filtering.pdf", plot = avg_exp_plot, width = 6, height = 4)
ggsave("plots/expression_filtering.png", plot = avg_exp_plot, width = 6, height = 4)



# How many genes are selected? 15,585/20,781 (around 75%)
table(is_expressed)

background_df <- data.frame(
  EnsemblGeneID = names(is_expressed), 
  Above_Background_Threshold = is_expressed
)

expression <- expression[is_expressed,]

## Step 5: Statistical analysis


# The covariates to be adjusted 
# gender: to remove sex-specific differences
# age: to remove age-specific differences
# RIN: to remove effects caused by RNA integrity
# Library.pool: to remove technical batch effects 
# Since the sample size is large enough we automatically correct for these 
# covariates without prior statistical testing like ANOVA 
# These co-variates are known to influence RNA-seq analysis, so it is a good
# practice to correct for them regardless of the significance of ANOVA or kruskal
# wallis test. The result of those tests are available in Patient characteristic 
# table.




print(all(rownames(sampleInfo) == colnames(expression)))

# model.matrxi removes rows with missing values
# for covariates so we need to account for that 
cat("Number of missing RIN values:", sum(is.na(sampleInfo$RIN)), "\n")

complete_samples <- sampleInfo %>% 
  filter(!is.na(RIN) & !is.na(etiology) & !is.na(age) & !is.na(gender) & !is.na(Library.Pool)) %>%
  rownames()


cat("Number of Remaining Samples after filtering for missing values:"
    , length(complete_samples),"/",nrow(sampleInfo), "\n")


sampleInfo <- sampleInfo[complete_samples, ]
expression <- expression[, complete_samples]

design <- model.matrix(~0+etiology+RIN+age+gender+Library.Pool, data = sampleInfo)
colnames(design) <- gsub("etiology", "", colnames(design))


# These are the groups that are compared based on the description
# Alternatively we can also compare different disease cases with each other
# Or we can also create interaction models to investigate the sex-specific 
# differences.
cm <- makeContrasts(
  DCM_vs_Control = DCM - NF,
  HCM_vs_Control = HCM - NF,
  PPCM_vs_Control = PPCM - NF,
  levels = design
)

# Fit this generalized linear model from limma package
fit <- lmFit(expression, design)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrasts = cm)

# Calculate the t-statistics for the contrasts
# Trend is sett to true to check for mean-variance trend
fit2 <- eBayes(fit2,trend=TRUE)

# Summarize results
results <- decideTests(fit2,p.value=0.05,adjust.method="fdr")
summary(results)


# Define a function to get the q-values and statistics
get_results <- function(fit_obj, coef_name) {
  
  # n=Inf ensures we get all genes 
  res <- topTable(fit_obj, coef = coef_name, number = Inf, sort.by = "none") 
  
  # Calculate q-values using the qvalue package
  qobj <- qvalue(p = res$P.Value)
  
  # Add q-values
  res %>%
    rownames_to_column(var = "EnsemblGeneID") %>%
    mutate(
      q.value = qobj$qvalues,    # Storey's q-value
      comparison = coef_name     # Label the comparison
    ) %>%
    arrange(P.Value)             # Sort by significance
}

# Run for all comparisons and combine into one table
all_results <- bind_rows(
  get_results(fit2, "DCM_vs_Control"),
  get_results(fit2, "HCM_vs_Control"),
  get_results(fit2, "PPCM_vs_Control")
)

head(all_results)

# Count significant genes (e.g., q-value < 0.05)
table(all_results$comparison, all_results$q.value < 0.05)




## Step 6: Export
# Here we merge 3 datasets
# 1. average expression of genes for each etiology
# 2. Gene annotation from part 5 
# 3. Results (fold-changes, pvalues and qvalues from part 4)


# 1. average expression of genes for each etiology
avg_expression <- expression %>%
  as.data.frame() %>%
  rownames_to_column(var = "EnsemblGeneID") %>%
  # longer pivot for joining
  pivot_longer(-EnsemblGeneID, names_to = "sample_name", values_to = "logCPM") %>% 
  # joining with sample metadata
  left_join(sampleInfo %>% select(sample_name, etiology), by = "sample_name") %>%
  # grouping by genes and etiology 
  group_by(EnsemblGeneID, etiology) %>%
  # for each gene,etiology pair take the average of logCPM value (or FPKM alternatively)
  summarise(Mean_logCPM = mean(logCPM), .groups = "drop") %>%
  # Convert back to wide format 
  pivot_wider(
    names_from = etiology, 
    values_from = Mean_logCPM, 
    names_prefix = "Avg_Expr_"
  )




# Merging the three tables 
final_table <- all_results %>% #DGE analysis results
  # Add Annotation
  left_join(gene_map, by = c("EnsemblGeneID"="ensembl_gene_id")) %>%
  # Add Group Averages
  left_join(avg_expression, by = "EnsemblGeneID") %>%
  # Reorder columns
  select(
    EnsemblGeneID,gene_symbol, chromosome_name,start_position, end_position,comparison,
    logFC, P.Value, q.value,
    Avg_Expr_NF, Avg_Expr_DCM, Avg_Expr_HCM, Avg_Expr_PPCM
  ) %>%
  # Sort by most significant q-value
  arrange(q.value)



final_table <- final_table %>%
 
  
  #  Pivot the statistics to Wide format
  # I pivot logFC, P.value, and q.value so columns become "logFC_DCM_vs_Control", etc.
  # similar to the input file for skills session week6
  pivot_wider(
    id_cols = c(EnsemblGeneID, gene_symbol, starts_with("Avg_Expr")), # Keep these fixed
    names_from = comparison, 
    values_from = c(logFC, P.Value, q.value),
    names_glue = "{.value}_{comparison}" # Custom column naming pattern
  ) %>%
  
  # Rename columns to be more readable
  rename_with(~ str_replace_all(., "Avg_Expr_", "Mean Expression (logCPM): "), starts_with("Avg_Expr")) %>%
  rename_with(~ str_replace_all(., "logFC_", "Log2FC: "), starts_with("logFC")) %>%
  rename_with(~ str_replace_all(., "q.value_", "FDR: "), starts_with("q.value")) %>%
  
  # Add if it is expressed above background level or not
  left_join(background_df, by = "EnsemblGeneID")   %>% 


  # Final ordering (gene id and name first and then the statistics)
  select(
    EnsemblGeneID, gene_symbol,  Above_Background_Threshold,
    starts_with("Mean Expression"), 
    everything()
  )

write_tsv(final_table, file = "tables/Table2_Differential_Expression_Results.txt")



# Additional Step: Volcano Plot



dcm_results <- all_results %>%
  filter(comparison == "DCM_vs_Control") %>%
  left_join(gene_map, by = c("EnsemblGeneID" = "ensembl_gene_id"))


# This is based on Skills session week 6
log2FC.cutoff <- 0.585  # ~1.5x fold change
pval.cutoff <- 0.05     # Significance threshold (can also use q.value if preferred)

volcano_plot <- EnhancedVolcano(
  dcm_results,
  lab = dcm_results$gene_symbol,
  x = 'logFC',
  y = 'P.Value',
  
  # Title and subtitles
  title = 'DCM vs Control',
  subtitle = 'Differential expression adjusted for Age, Sex, RIN, and Batch',
  
  # Thresholds
  pCutoff = pval.cutoff,
  FCcutoff = log2FC.cutoff,
  
  # Visual Customization
  labSize = 3,

) + center_title

volcano_plot
#  Save the plot
ggsave("plots/Volcano_DCM_vs_Control.pdf", plot = volcano_plot, width = 10, height = 10)
ggsave("plots/Volcano_DCM_vs_Control.png", plot = volcano_plot, width = 10, height = 8)
