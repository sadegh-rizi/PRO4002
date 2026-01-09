# ==============================================================================
# PROJECT: GENE EXPRESSION ANALYSIS (DCM, HCM, PPCM vs NF)
# STYLE: Tidyverse (http://style.tidyverse.org/)
# ==============================================================================

# ------------------------------------------------------------------------------
# 01. SETUP & PACKAGE INSTALLATION
# ------------------------------------------------------------------------------


# 01 Installing packages --------------------------------------------------------
analysis_packages <- c(
  "tidyverse", "tidyr", "dplyr", "gridExtra", "pcaMethods", 
  "data.table", "tableone", "kableExtra", "rmarkdown", 
  "readr", "readxl", "gprofiler2", "knitr"
)

bio_packages <- c("limma", "qvalue", "biomaRt", "Biocmanager")




# Loading packages, and checking warnings --------------------------------------

#This section will confirm to load all necessary packages
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



#Confirm all packages are installed

search()

# # ------------------------------------------------------------------------------
# 02. DATA IMPORT & CLEANING
# --------------------------------------------------------------------------------
#Confirm necessary packages for this section
library(readr)
require(readr)
library(readxl)
require(readxl)


#The next code allows you to check your directory file, c=make sure to have the 
#data in the same file, as a result you have set the script location

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(current_dir)
  message(paste("Working directory set to script location:", getwd()))
} else {
  # Fallback: If not in RStudio, the user must ensure they are running the script
  # from the correct location, or you can add a generic error message/stop.
  warning("Not running in RStudio. Please set working directory manually for portability.")
}

#Importing data to r 

#clear variable names for the main files where sample_data have the phenotype information
#related to the patients , gene expression containsthe actual sequencing results.
# Unit: CPM (Counts Per Million). Finally, exon lenght contains the Total Exon Length for every gene


sample_data <- read.csv("MAGNET_SampleData_18112022.csv", as.is = T, 
                               row.names = 1)

gene_expression <- read.delim("MAGNET_GeneExpressionData_CPM_19112020.txt", as.is = T,
                              row.names = 1)


exons_lenght <-  read.delim("MAGNET_exonLengths.txt", as.is = T, 
                           row.names = 1) 


data_description <- read_excel("MAGNET_SampleData_18112022_WithDescriptions.xlsx")

View(data_description)



# Section 2 ---------------------------------------------------------------
#Associate with the data structure
str(sample_data)
head(sample_data)
tail(sample_data)





# ==============================================================================
# 3. PARTICIPANT CHARACTERISTICS TABLE (Assignment 1b)
# ==============================================================================

# Confirm the pakages that you need

library(data.table)
library(tableone) #help in the exporting, and format of tables

# ==============================================================================
# 3. PARTICIPANT CHARACTERISTICS TABLE 
# ==============================================================================

# 3a. Define the variables for the table, cleaning conflictive variable names

# Continuous variables (will be summarized as Mean(SD) or Median[IQR])
current_names <- colnames(sample_data)


# Step 2: Create a mapping to rename the colunms

new_names <- gsub("LVEF", "lvef", current_names)
new_names <- gsub("RIN", "rin", new_names)
new_names <- gsub("TIN.median.", "tin_median", new_names, fixed = TRUE) 


# Step 3: Apply the new names
colnames(sample_data) <- new_names

# Continuous variables:
cont_vars <- c("age", "weight", "height", "lvef", "rin", "tin_median") 

# Categorical variables
# NOTE: 'etiology' is the grouping variable, not a variable to summarize!
cat_vars <- c("gender", "race", "afib", "VTVF", "Diabetes", "Hypertension")

  
# Defining the first table about table of participants

table_one <- tableone::CreateTableOne(
  vars = c(cont_vars, cat_vars), 
  strata = "etiology", 
  data = sample_data,
  factorVars = cat_vars, # Ensure categorical variables are treated as factors
  test = TRUE # This computes the p-values for comparison across strata
)

print(table_one)


# 3c. Format the table for publication


# a. Convert the table object to a data frame 
table_one_formatted <- print(
  table_one, 
  printToggle = FALSE,         
  noSpaces = TRUE,              
  # (uses Kruskal-Wallis instead of ANOVA for p-value)
  contDigits = 1,               # Number of decimal places for continuous variables
  pDigits = 3                   # Number of decimal places for p-values
)

# Clean up column names in the exported data frame
table_one_formatted <- table_one_formatted %>%
  as.data.frame() %>%
  rownames_to_column(var = "Variable")


# Display the first few rows to show the user the output structure
print(head(table_one_formatted))
message("Table 1 calculated. (Export will happen in exporting section )")



# ==============================================================================
# TASK 2. DIAGNOSTIC PLOTS
# ==============================================================================

#Quality control of rows and colunms in the data sets: 
#phenotype patient data set
pheno_qc <-  sample_data %>%
  tibble::rownames_to_column(var = "sample_id") %>%
  dplyr::select(sample_id, etiology)
    
expr_long <- gene_expression %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Ensembl_ID") %>%
  tidyr::pivot_longer(
    cols = -Ensembl_ID, 
    names_to = "sample_id", 
    values_to = "expression_value"
  )
qc_data <- expr_long %>%
  dplyr::left_join(pheno_qc, by = "sample_id") %>%
  # FIX: Use tidyr::drop_na explicitly to avoid the error
  tidyr::drop_na(expression_value, etiology)
boxplot_faceted <- ggplot(qc_data, aes(x = sample_id, y = expression_value, fill = etiology)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) +
  scale_fill_brewer(palette = "Set1") +
  # Split the plot into panels based on 'etiology'
  facet_wrap(~ etiology, scales = "free_x", ncol = 2) + 
  labs(
    title = "Gene Expression Distribution by Etiology",
    subtitle = "Separate boxplots for Healthy (NF), DCM, HCM, and PPCM",
    x = "Sample ID",
    y = expression("log"[2] * "(CPM) Expression Value"),
    fill = "Etiology"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    # Rotate x-axis text 90 degrees so sample IDs are readable
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    plot.title = element_text(face = "bold"),
    legend.position = "none" # 
  )

print(boxplot_faceted)
message("Diagnostic plot created. (Export will happen in Section 10)")


# Analyzing individually the density plots
etiology1 <- sample_data %>%
  dplyr::select(etiology) %>% 
  table()

etiology2 <- sample_data %>%
  count(etiology)

etiology3 <- sample_data %>%
  group_by(etiology) %>% 
  tally()

avg_etiology_by_age <- sample_data %>%
  group_by(etiology) %>%
  summarise(
    n = n(),
    age = mean(age, na.rm = TRUE)
  )
avg_etiology_by_age

names_et <- etiology2$etiology

# Get the names of the samples per etiology
for (i in names_et){
  # Define the name of the dataframes with the list of samples per etiology
  var_name <- paste(i,"columns",sep = "_")
  # Get the dataframes with the list of sample names per etiology
  df <- rownames_to_column(sample_data, var = "sample") %>% 
    dplyr::filter(etiology == i) %>% 
    dplyr::select(sample) %>% 
    pull(sample)
  # Assign the names to the dataframes
  assign(var_name, df)
}

NF_data <- gene_expression %>%
  dplyr::select(NF_columns)

NF_data_group <- gather(NF_data) # melt also works, but it is deprecated

ggplot(NF_data_group, aes(x = key, y = value)) +
  geom_boxplot()

# ==============================================================================
#  Clinical Characteristics Boxplots 
# ==============================================================================

# 1. Select the variables you want to plot
# We select 'etiology' (for grouping) and the numeric variables of interest.
clinical_data_long <- sample_data %>%
  dplyr::select(etiology, age, weight, height, lvef, rin) %>%
  # 2. Convert to long format: piles all values into one column, with a label column
  tidyr::pivot_longer(
    cols = -etiology,          # Keep etiology as the ID
    names_to = "characteristic", # Column for variable names (age, weight...)
    values_to = "value"        # Column for the actual numbers
  )

# 3. Create the Faceted Plot
clinical_plot <- ggplot(clinical_data_long, aes(x = etiology, y = value, fill = etiology)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white") + # Nice looking outliers
  scale_fill_brewer(palette = "Set1") +
  # 4. facet_wrap creates a separate panel for each 'characteristic'
  facet_wrap(~ characteristic, scales = "free_y", nrow = 2) + 
  labs(
    title = "Patient Characteristics by Etiology",
    x = "Disease Group",
    y = "Value (scales vary)",
    fill = "Etiology"
  ) +
  theme_bw(base_size = 12) + # A clean black-and-white theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none", # Legend is redundant since x-axis is labeled
    strip.background = element_rect(fill = "lightgrey"), # Gray headers for panels
    strip.text = element_text(face = "bold")
  )

# Export
ggsave("clinical_characteristics_panel.png", clinical_plot, width = 10, height = 8)

message("Clinical characteristics plot exported.")

print(clinical_plot)



# PCA ANALYSIS ------------------------------------------------------------
library(ggplot2)

t_gxData <- t(gene_expression)

# Calculate pca

pca_hf <- pca(t_gxData, method = "svd")

summary(pca_hf)

plot(pca_hf)

df_hf <- merge(scores(pca_hf), sample_data, by = 0)

plotPcs(pca_hf)

var_pc1 <- round(pca_hf@R2[1] * 100, 1)
var_pc2 <- round(pca_hf@R2[2] * 100, 1)

pca_plot_final <- ggplot(df_hf, aes(x = PC1, y = PC2, color = etiology)) + 
  stat_ellipse(aes(fill = etiology), 
               geom = "polygon", alpha = 0.1, color = NA, show.legend = FALSE) +
  
  # Shape is defined ONLY here, for the points
  geom_point(aes(shape = gender), size = 3, alpha = 0.8) +
  
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  
  xlab(paste0("PC1 (", var_pc1, "% Variance)")) +
  ylab(paste0("PC2 (", var_pc2, "% Variance)")) +
  
  labs(title = "PCA: Gene Expression by Etiology",
       color = "Etiology",
       shape = "Gender") +
  
  theme(aspect.ratio = 1)
# 3. Save the Plot as an Image
ggsave("PCA_Analysis_Plot.png", pca_plot_final, width = 8, height = 8, dpi = 300)




# 3. Statistical analysis -------------------------------------------------

library(limma)

# ==============================================================================
# TASK 3a: SIMPLE ANALYSIS (NO COVARIATES)
# ==============================================================================

#Expressing new variables to avoid overwriting

#We use new names so we don't overwrite data
raw_counts_stats <- gene_expression
raw_meta_stats   <- sample_data

# 2. Fix Orientation (Genes MUST be Rows for Limma)
if (nrow(raw_counts_stats) < ncol(raw_counts_stats)) {
  expr_final <- t(raw_counts_stats)
} else {
  expr_final <- raw_counts_stats
}

# 3. Match the Samples (Intersection)
common_ids <- intersect(colnames(expr_final), rownames(raw_meta_stats))

if (length(common_ids) == 0) {
  stop("CRITICAL ERROR: No matching sample IDs found. Check your names!")
}

# 4. Create the FINAL Clean Datasets for Statistics
# We use distinct names (_final) so you don't confuse them with your plotting data
expr_final <- expr_final[, common_ids]
meta_final <- raw_meta_stats[common_ids, ]

# 1. Matrix design 
# We only use 'etiology'. We ignore age, gender, and RIN for now.
design_simple <- model.matrix(~ 0 + etiology, data = meta_final)


# Keeping column names clean:
colnames(design_simple) <- gsub("etiology", "", colnames(design_simple))


if (ncol(expr_final) != nrow(design_simple)) {
  message("MISMATCH DETECTED! Checking if transposition fixes it...")
  
  # Check if transposing makes the numbers match
  if (nrow(expr_final) == nrow(design_simple)) {
    message("Fix found: Transposing expression data...")
    expr_final <- t(expr_final)
  } else {
    # If numbers still don't match, we strictly subset to common samples again
    common_ids <- intersect(colnames(expr_final), rownames(design_simple))
    if (length(common_ids) == 0) {
      common_ids <- intersect(rownames(expr_final), rownames(design_simple))
      if (length(common_ids) > 0) {
        expr_final <- t(expr_final) # Flip it because IDs were on rows
      }
    }
    
    if (length(common_ids) == 0) stop("CRITICAL ERROR: Sample IDs strictly do not match.")
    
    # Apply the strict subset to both
    expr_final <- expr_final[, common_ids]
    design_simple <- design[common_ids, ]
    message(paste("Force-aligned to", length(common_ids), "samples."))
  }
} else {
  message("Dimensions look correct. Proceeding...")
}

print(paste("Final Check - Expr Cols:", ncol(expr_final), " Design Rows:", nrow(design)))




# 2. FIT THE MODEL
# We use the same aligned data 'expr_final'
fit_simple <- lmFit(expr_final, design_simple)

# 3. DEFINE CONTRASTS
contrast_matrix_simple <- makeContrasts(
  DCM_vs_NF  = DCM - NF,
  HCM_vs_NF  = HCM - NF,
  PPCM_vs_NF = PPCM - NF,
  levels = design_simple
)

# 4. COMPUTE STATISTICS
fit2_simple <- contrasts.fit(fit_simple, contrast_matrix_simple)
fit2_simple <- eBayes(fit2_simple)

# 5. EXTRACT RESULTS (Uncorrected)
res_DCM_simple  <- topTable(fit2_simple, coef = "DCM_vs_NF", number = Inf, adjust.method = "BH")
res_HCM_simple  <- topTable(fit2_simple, coef = "HCM_vs_NF", number = Inf, adjust.method = "BH")
res_PPCM_simple <- topTable(fit2_simple, coef = "PPCM_vs_NF", number = Inf, adjust.method = "BH")
library(ggplot2)
library(gridExtra)

# 1. Volcano fucntion

create_volcano <- function(df, title_text) {
  
  # Create classification column
  df$diffexpressed <- "NO"
  df$diffexpressed[df$logFC > 1 & df$adj.P.Val < 0.05] <- "UP"
  df$diffexpressed[df$logFC < -1 & df$adj.P.Val < 0.05] <- "DOWN"
  
  # Plot
  ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = diffexpressed)) +
    geom_point(alpha = 0.6, size = 1.5) +
    theme_bw() +
    scale_color_manual(values = c("blue", "grey", "red")) +
    
    # Threshold lines
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    
    # Labels (Updated for 3a)
    labs(title = title_text,
         subtitle = "Statistical analysis
         (Task 3a)", 
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-Value") +
    
    theme(legend.position = "none", 
          plot.title = element_text(face = "bold", size = 11))
}

# 2. Generate the 3 Uncorrected Plots
# We use the '_simple' result objects created in the previous step
p_dcm_simple  <- create_volcano(res_DCM_simple, "DCM vs NF")
p_hcm_simple  <- create_volcano(res_HCM_simple, "HCM vs NF")
p_ppcm_simple <- create_volcano(res_PPCM_simple, "PPCM vs NF")

# 3. Arrange them side-by-side
grid.arrange(p_dcm_simple, p_hcm_simple, p_ppcm_simple, ncol = 3)




# 5. Verify data for our covariates 

print(paste("Matching samples:", length(common_ids)))
print(paste("Expr Dimensions (Rows x Cols):", nrow(expr_final), "x", ncol(expr_final)))
print(paste("Meta Dimensions (Rows x Cols):", nrow(meta_final), "x", ncol(meta_final)))

# --- STEP 2: BUILD DESIGN MATRIX ---
# We build the design first so we have a target to match against.
design <- model.matrix(~ 0 + etiology + gender + age + rin, data = meta_final)
colnames(design) <- gsub("etiology", "", colnames(design))
colnames(design) <- gsub("gender", "", colnames(design))

if (ncol(expr_final) != nrow(design)) {
  message("MISMATCH DETECTED! Checking if transposition fixes it...")
  
  # Check if transposing makes the numbers match
  if (nrow(expr_final) == nrow(design)) {
    message("Fix found: Transposing expression data...")
    expr_final <- t(expr_final)
  } else {
    # If numbers still don't match, we strictly subset to common samples again
    common_ids <- intersect(colnames(expr_final), rownames(design))
    if (length(common_ids) == 0) {
      common_ids <- intersect(rownames(expr_final), rownames(design))
      if (length(common_ids) > 0) {
        expr_final <- t(expr_final) # Flip it because IDs were on rows
      }
    }
    
    if (length(common_ids) == 0) stop("CRITICAL ERROR: Sample IDs strictly do not match.")
    
    # Apply the strict subset to both
    expr_final <- expr_final[, common_ids]
    design <- design[common_ids, ]
    message(paste("Force-aligned to", length(common_ids), "samples."))
  }
} else {
  message("Dimensions look correct. Proceeding...")
}

print(paste("Final Check - Expr Cols:", ncol(expr_final), " Design Rows:", nrow(design)))

fit <- lmFit(expr_final, design)

contrast_matrix <- makeContrasts(
  DCM_vs_NF  = DCM - NF,
  HCM_vs_NF  = HCM - NF,
  PPCM_vs_NF = PPCM - NF,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Check results
print("Success! Top 5 DCM Genes:")
print(head(topTable(fit2, coef = "DCM_vs_NF")))

# ==============================================================================
# 4. EXTRACT RESULTS (DCM, HCM, PPCM vs NF)
# ==============================================================================

# 1. Extract DCM vs NF
# number = Inf ensures we get ALL genes, not just the top 10
res_DCM <- topTable(fit2, coef = "DCM_vs_NF", number = Inf, adjust.method = "BH")

# 2. Extract HCM vs NF
res_HCM <- topTable(fit2, coef = "HCM_vs_NF", number = Inf, adjust.method = "BH")

# 3. Extract PPCM vs NF
res_PPCM <- topTable(fit2, coef = "PPCM_vs_NF", number = Inf, adjust.method = "BH")

# --- VIEW THE TOP GENES ---
message("--- Top 5 Genes: DCM vs NF (Corrected) ---")
print(head(res_DCM))

message("--- Top 5 Genes: HCM vs NF (Corrected) ---")
print(head(res_HCM))

message("--- Top 5 Genes: PPCM vs NF (Corrected) ---")
print(head(res_PPCM))

# Save to CSV 
write.csv(res_DCM, "Results_DCM_vs_NF_Corrected.csv")
write.csv(res_HCM, "Results_HCM_vs_NF_Corrected.csv")
write.csv(res_PPCM, "Results_PPCM_vs_NF_Corrected.csv")

# ==============================================================================
# 5. VISUALIZATION (VOLCANO PLOTS)
# ==============================================================================
library(ggplot2)
library(gridExtra)

# Helper function to make the plots consistent
create_volcano <- function(df, title_text) {
  
  # Create a column to color-code Significant Genes
  # Criteria: Log2 Fold Change > 1 (or < -1) AND Adj.P.Val < 0.05
  df$diffexpressed <- "NO"
  df$diffexpressed[df$logFC > 1 & df$adj.P.Val < 0.05] <- "UP"
  df$diffexpressed[df$logFC < -1 & df$adj.P.Val < 0.05] <- "DOWN"
  
  # Plot
  p <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = diffexpressed)) +
    geom_point(alpha = 0.6, size = 1.5) +
    theme_bw() +
    scale_color_manual(values = c("blue", "grey", "red")) +
    
    # Add threshold lines
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    
    # Labels
    labs(title = title_text,
         subtitle = "Corrected for Age, Gender, RIN",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-Value") +
    
    theme(legend.position = "none", 
          plot.title = element_text(face = "bold", size = 11))
  
  return(p)
}

# Generate the 3 plots
p1 <- create_volcano(res_DCM, "DCM vs NF")
p2 <- create_volcano(res_HCM, "HCM vs NF")
p3 <- create_volcano(res_PPCM, "PPCM vs NF")

print(p1)
print(p2)
print(p3)

# --- COMPARISON: 3a vs 3b ---
# This is great for your report comments! 
# Let's count how many significant genes we found in each method.

# Helper function to count significant genes (Adj.P < 0.05)
count_sig <- function(df) { sum(df$adj.P.Val < 0.05) }

print("--- Comparison of Significant Genes (Adj.P < 0.05) ---")
print(paste("DCM vs NF (Simple 3a):", count_sig(res_DCM_simple)))
print(paste("DCM vs NF (Corrected 3b):", count_sig(res_DCM))) # Uses your result from previous step

print(paste("PPCM vs NF (Simple 3a):", count_sig(res_PPCM_simple)))
print(paste("PPCM vs NF (Corrected 3b):", count_sig(res_PPCM)))

library(knitr) # For nice table formatting (optional)

# 1. Define a helper function to count significant genes
# (Threshold: Adjusted P-value < 0.05)
get_sig_count <- function(df) {
  return(sum(df$adj.P.Val < 0.05, na.rm = TRUE))
}

# 2. Create the Data Frame
# We assume you have the result objects from previous steps:
# 3a objects: res_DCM_simple, res_HCM_simple, res_PPCM_simple
# 3b objects: res_DCM, res_HCM, res_PPCM

comparison_table <- data.frame(
  Contrast = c("DCM vs NF", "HCM vs NF", "PPCM vs NF"),
  
  # Column 1: Uncorrected Counts (Task 3a)
  Sig_Genes_Uncorrected = c(
    get_sig_count(res_DCM_simple),
    get_sig_count(res_HCM_simple),
    get_sig_count(res_PPCM_simple)
  ),
  
  # Column 2: Corrected Counts (Task 3b)
  Sig_Genes_Corrected = c(
    get_sig_count(res_DCM),
    get_sig_count(res_HCM),
    get_sig_count(res_PPCM)
  )
)

# 3. Add a "Difference" Column
# Negative number = Correction removed noise (likely false positives)
# Positive number = Correction found NEW genes (increased power)
comparison_table$Difference <- comparison_table$Sig_Genes_Corrected - comparison_table$Sig_Genes_Uncorrected

# 4. View the Table
print(comparison_table)




# 4. gene annotation ------------------------------------------------------

library(biomaRt)

# 1. Connect to the Ensembl database
# We use the 'hsapiens_gene_ensembl' dataset (Human genes)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 2. Get the list of Ensembl IDs from your data
# We use rownames(expr_final) because it contains ALL genes in your analysis.
my_ensembl_ids <- rownames(expr_final)

# 3. Retrieve the attributes
# - ensembl_gene_id: The key to match our data
# - external_gene_name: The symbol (e.g., TP53)
# - description: The full name (e.g., Tumor Protein P53)
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description"),
  filters = "ensembl_gene_id",
  values = my_ensembl_ids,
  mart = ensembl
)

# 4. Clean up the descriptions
# (Ensembl descriptions often have extra text like "[Source:HGNC...]" that we can remove)
gene_info$description <- gsub(" \\[Source:.*\\]", "", gene_info$description)

# Check what we got
head(gene_info)

# Function to annotate a result table
annotate_results <- function(res_table, annotation_data) {
  
  # 1. Convert result to data frame and move Ensembl IDs to a proper column
  res_df <- as.data.frame(res_table)
  res_df$ensembl_gene_id <- rownames(res_df)
  
  # 2. Merge with annotation data
  # all.x = TRUE ensures we keep all rows even if no name was found
  merged_df <- merge(res_df, annotation_data, by = "ensembl_gene_id", all.x = TRUE)
  
  # 3. Sort by significance (Adjusted P-value) again, as merge scrambles order
  merged_df <- merged_df[order(merged_df$adj.P.Val), ]
  
  # 4. Reorder columns so Gene Name is at the front (easier to read)
  # Moves 'external_gene_name' to the 2nd column
  cols <- c("ensembl_gene_id", "external_gene_name", "description", 
            setdiff(names(merged_df), c("ensembl_gene_id", "external_gene_name", "description")))
  merged_df <- merged_df[, cols]
  
  return(merged_df)
}

# Apply this to your tables from Task 3b (Corrected)
res_DCM_annotated  <- annotate_results(res_DCM, gene_info)
res_HCM_annotated  <- annotate_results(res_HCM, gene_info)
res_PPCM_annotated <- annotate_results(res_PPCM, gene_info)

# View the final human-readable table
print("--- Top Annotated Genes (DCM) ---")
head(res_DCM_annotated[, 1:5]) # Printing just the first 5 columns to fit screen

# Check column names of your annotated result
colnames(res_DCM_annotated)


library(kableExtra)
library(dplyr)

# 1. Prepare the Data Frame
top_genes_table <- res_DCM_annotated %>%
  head(50) %>%
  dplyr::select(external_gene_name, description, logFC, adj.P.Val, ensembl_gene_id) %>%
  # RENAME COLUMNS
  dplyr::rename(
    "Gene Symbol" = external_gene_name,
    "Description" = description,
    "Log2 Fold Change" = logFC,
    "Adj. P-Value" = adj.P.Val,
    "Ensembl ID" = ensembl_gene_id
  ) %>%
  # --- FIX: MUTATE BEFORE KBL ---
  # Format P-values to scientific notation (e.g., 1.2e-10) here
  mutate(`Adj. P-Value` = formatC(`Adj. P-Value`, format = "e", digits = 2)) %>%
  mutate(`Log2 Fold Change` = round(`Log2 Fold Change`, 2))

# 2. Create and Style the HTML Table
html_table <- top_genes_table %>%
  kbl(caption = "<b>Table 1: Top 50 Differentially Expressed Genes (DCM vs NF)</b>") %>%
  
  # Add professional styling
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "responsive"), 
    full_width = TRUE, 
    position = "center",
    font_size = 12
  ) %>%
  
  # Add scroll box
  scroll_box(width = "100%", height = "500px") %>%
  
  # Color code the Fold Change column (Column 3)
  # Note: We must check the original numeric values, but since we rounded them, 
  # it's safer to use the column index 3 directly.
  column_spec(3, color = ifelse(top_genes_table$`Log2 Fold Change` > 0, "green", "red"), bold = TRUE)

# 3. Save to file
save_kable(html_table, file = "Top_DCM_Genes_Annotated.html")

message("Success! Check your folder for 'Top_DCM_Genes_Annotated.html'")


# ==============================================================================
# MERGE ANNOTATION WITH FULL EXPRESSION DATA
# ==============================================================================

# 1. Convert your final expression matrix to a data frame
expr_df <- as.data.frame(expr_final)

# 2. Add the Ensembl ID as a column (currently it's just row names)
expr_df$ensembl_gene_id <- rownames(expr_df)

# 3. Merge with the 'gene_info' we got from biomaRt
# 'all.x = TRUE' keeps all genes, even if they don't have a symbol
full_expression_annotated <- merge(gene_info, expr_df, by = "ensembl_gene_id", all.y = TRUE)

# 4. View the result
# You will see: Ensembl ID | Gene Symbol | Description | Sample_101 | Sample_102 ...
head(full_expression_annotated[, 1:6]) 

# 5. (Optional) Export this master table
write.csv(full_expression_annotated, "Full_Expression_Data_Annotated.csv", row.names = FALSE)


# FPKM --------------------------------------------------------------------

common_genes <- intersect(rownames(expr_final), rownames(exons_lenght))
if (length(common_genes) == 0) {
  stop("Error: No matching gene IDs between expression data and exon lengths!")
}

expr_for_fpkm <- expr_final[common_genes, ]
lengths_clean <- exons_lenght[common_genes, , drop = FALSE]
if (!all(rownames(expr_for_fpkm) == rownames(lengths_clean))) stop("Alignment failed!")

cpm2fpkm <- function(x) {
  t <- 2^(x) * 1E3 / lengths_clean[, 1]
  return(t)
}
# 4. Apply the function
gxData_fpkm <- cpm2fpkm(expr_for_fpkm)

# 5. Check the result
message("Transformation successful!")
print("FPKM Data Dimensions:")
print(dim(gxData_fpkm))
print("Top 5 rows/cols of FPKM data:")
print(gxData_fpkm[1:5, 1:5])


# 5B; Gene ontology -------------------------------------------------------


install.packages("gprofiler2")
library(gprofiler2)

# 1. Filter for significant genes only (Adj P < 0.05)
sig_genes_DCM <- subset(res_DCM_annotated, adj.P.Val < 0.05)

# 2. Separate them into UP and DOWN regulated
# This helps us see if "Inflammation" is being turned ON (Up) or OFF (Down)
up_genes   <- subset(sig_genes_DCM, logFC > 0)$external_gene_name
down_genes <- subset(sig_genes_DCM, logFC < 0)$external_gene_name

# 3. Combine them into a named list for the tool
query_list <- list(
  "DCM_Up" = up_genes,
  "DCM_Down" = down_genes
)


# Check how many we have
print(paste("UP genes:", length(up_genes)))
print(paste("DOWN genes:", length(down_genes)))

# Run Enrichment Analysis
gost_res <- gost(query = query_list, 
                 organism = "hsapiens", # Human
                 sources = c("GO:BP", "KEGG"), # GO Biological Process + KEGG Pathways
                 user_threshold = 0.05)

# Check if we found anything
print(head(gost_res$result))

# Interactive Plot (Great for exploring)
p <- gostplot(gost_res, interactive = TRUE)


# Static Plot (For saving to your report)
p_static <- gostplot(gost_res, interactive = FALSE)
publish_gostplot(p_static, highlight_terms = head(gost_res$result$term_id, 10))

top_pathways <- gost_res$result %>%
  dplyr::filter(source %in% c("KEGG", "REAC")) %>%
  dplyr::select(term_name, source, p_value, term_size, intersection_size) %>%
  dplyr::arrange(p_value) %>%
  head(10)

print(top_pathways)



# EXPORTING EVERYTHING ----------------------------------------------------


# ==============================================================================
# 3e. EXPORT AS HTML TABLE (Publication Format)
# ==============================================================================

table_one_clean <- table_one_formatted %>%
  dplyr::select(1:6) 

# 4. Define the Description Text
desc_table1 <- paste(
  "<div style='text-align: left; margin-bottom: 10px;'>",
  "<strong>Table 1: Participant Characteristics by Etiology.</strong><br>",
  "This table summarizes the demographic and clinical variables for Healthy Donors (NF) ",
  "and Heart Failure patients (DCM, HCM, PPCM).<br>",
  "<em>Note: P-values calculate comparison across all strata.</em>",
  "</div>"
)

# 5. Create the Styled HTML Table
html_table1 <- table_one_clean %>%
  kbl(caption = desc_table1, align = "c") %>%
  
  # Professional Styling
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "responsive"),
    full_width = FALSE, 
    position = "center",
    font_size = 14
  ) %>%
  
  # Add Group Header (1 empty + 4 Groups + 1 Stat = 6 Columns)
  add_header_above(c(" " = 1, "Patient Groups" = 4, "Statistics" = 1)) %>%
  
  # Styling specific columns
  column_spec(1, bold = TRUE, width = "5cm") %>%
  row_spec(0, bold = TRUE, color = "white", background = "#2c3e50")

# 6. Save Files
save_kable(html_table1, file = "Participant_Characteristics_Table.html")
write_delim(table_one_formatted, "Participant_Characteristics_Table.txt", delim = "\t") # Saves full version to TXT

message("Table 1 Exported Successfully.")

# Print the current working directory to the console
print(paste("The HTML file will be saved in:", getwd()))

# ==============================================================================
# 10. TASK 6: FINAL EXPORT (ALL FILES)
# ==============================================================================

message("Generating Final Export Files...")

# ... [Keep Part A (Table 1) and Part B (Master Results) as they were] ...

# ------------------------------------------------------------------------------
# C. DIAGNOSTIC PLOTS (PNG & HTML)
# ------------------------------------------------------------------------------

# 1. Save the Plot as an Image
ggsave("QC_Boxplots.png", boxplot_faceted, width = 10, height = 8, dpi = 300)

# 2. Define the Description Text
desc_plots <- paste(
  "<div style='font-family: Arial; padding: 20px;'>",
  "<h2>Diagnostic Plots: Gene Expression Distribution</h2>",
  "<p>The following boxplots display the distribution of normalized gene expression values ",
  "(log2 CPM) for every sample, grouped by etiology.</p>",
  "<ul>",
  "<li><strong>Purpose:</strong> To visually inspect data quality and normalization consistency.</li>",
  "<li><strong>Interpretation:</strong> Medians (center lines) should be roughly aligned across samples. ",
  "The 'notches' indicate the confidence interval around the median.</li>",
  "</ul>",
  "<br>",
  # Embed the image directly into the HTML
  "<img src='QC_Boxplots.png' alt='QC Boxplots' style='width: 100%; max-width: 1000px; border: 1px solid #ddd;'/>",
  "</div>"
)

# 3. Save as HTML File
# We use writeLines because we are writing raw HTML/Image tags, not a table
fileConn <- file("Diagnostic_Plots.html")
writeLines(desc_plots, fileConn)
close(fileConn)

message("Diagnostic Plots Exported: 'QC_Boxplots.png' and 'Diagnostic_Plots.html'")

# D. CLINICAL DIAGNOSTIC PLOTS (PNG & HTML)
# ------------------------------------------------------------------------------

# 1. Save the Clinical Plot as an Image
ggsave("QC_Clinical_Data.png", clinical_plot, width = 10, height = 8, dpi = 300)

# 2. Define the HTML Content for Clinical Data
desc_clinical <- paste(
  "<html>",
  "<head>",
  "<style>",
  "  body { font-family: Arial, sans-serif; padding: 40px; background-color: #f9f9f9; }",
  "  .container { background-color: white; padding: 30px; border-radius: 8px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }",
  "  h2 { color: #2c3e50; border-bottom: 2px solid #eee; padding-bottom: 10px; }",
  "  p { color: #555; line-height: 1.6; }",
  "  ul { color: #555; line-height: 1.6; }",
  "  img { width: 100%; max-width: 1000px; border: 1px solid #ddd; display: block; margin: 20px auto; }",
  "</style>",
  "</head>",
  "<body>",
  
  "<div class='container'>",
  "  <h2>Clinical Characteristics by Etiology</h2>",
  "  <p>The following faceted boxplots visualize the distribution of key clinical variables across the different disease groups (DCM, HCM, PPCM, and NF).</p>",
  "  <ul>",
  "    <li><strong>Variables Plotted:</strong> Age, Weight, Height, LVEF (Left Ventricular Ejection Fraction), and RIN (RNA Integrity Number).</li>",
  "    <li><strong>Interpretation:</strong> These plots help identify potential confounding factors (e.g., if one group is significantly older) and verify disease phenotypes (e.g., lower LVEF is expected in Heart Failure groups).</li>",
  "  </ul>",
  "  <br>",
  "  <img src='QC_Clinical_Data.png' alt='Clinical Characteristics Boxplots'/>",
  "</div>",
  
  "</body>",
  "</html>"
)

# 3. Save as a new HTML File
fileConn <- file("Clinical_Diagnostic_Plots.html")
writeLines(desc_clinical, fileConn)
close(fileConn)


message("Clinical Export Complete: 'QC_Clinical_Data.png' and 'Clinical_Diagnostic_Plots.html'")


# 4. Define the HTML Content with Interpretation
desc_pca <- paste(
  "<html>",
  "<head>",
  "<style>",
  "  body { font-family: Arial, sans-serif; padding: 40px; background-color: #f9f9f9; }",
  "  .container { background-color: white; padding: 30px; border-radius: 8px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }",
  "  h2 { color: #2c3e50; border-bottom: 2px solid #eee; padding-bottom: 10px; }",
  "  h3 { color: #e67e22; margin-top: 20px; }",
  "  p { color: #555; line-height: 1.6; }",
  "  li { margin-bottom: 10px; color: #555; }",
  "  img { width: 100%; max-width: 800px; border: 1px solid #ddd; display: block; margin: 20px auto; }",
  "</style>",
  "</head>",
  "<body>",
  
  "<div class='container'>",
  "  <h2>Principal Component Analysis (PCA)</h2>",
  "  <p>This plot reduces the dimensionality of the gene expression dataset to visualize sample similarity. It uses the Singular Value Decomposition (SVD) method.</p>",
  
  "  <img src='PCA_Analysis_Plot.png' alt='PCA Plot'/>",
  
  "  <h3>Interpretation of Results:</h3>",
  "  <ul>",
  "    <li><strong>Healthy vs. Disease Separation (PC2):</strong> The plot reveals a distinct separation along the second principal component (PC2). The <strong>Non-Failing (NF)</strong> samples (teal) cluster separately from the disease groups, indicating a clear transcriptional difference between healthy and failing hearts.</li>",
  "    <li><strong>Disease Group Overlap:</strong> The disease etiologies (DCM, HCM, PPCM) show significant overlap with each other. This suggests that despite different clinical causes, the end-stage heart failure phenotype shares a common underlying gene expression profile.</li>",
  "    <li><strong>Gender Distribution:</strong> The shape of the points indicates gender. PPCM samples (pink) are exclusively female (circles), which is consistent with the disease pathology (Peripartum Cardiomyopathy).</li>",
  "  </ul>",
  "</div>",
  
  "</body>",
  "</html>"
)

# 5. Save as HTML File
fileConn <- file("PCA_Analysis_Report.html")
writeLines(desc_pca, fileConn)
close(fileConn)

message("PCA Export Complete: 'PCA_Analysis_Plot.png' and 'PCA_Analysis_Report.html'")

# ==============================================================================
# 4. EXTRACT RESULTS (DCM, HCM, PPCM vs NF)
# ==============================================================================

# 1. Extract DCM vs NF
# number = Inf ensures we get ALL genes, not just the top 10
res_DCM <- topTable(fit2, coef = "DCM_vs_NF", number = Inf, adjust.method = "BH")

# 2. Extract HCM vs NF
res_HCM <- topTable(fit2, coef = "HCM_vs_NF", number = Inf, adjust.method = "BH")

# 3. Extract PPCM vs NF
res_PPCM <- topTable(fit2, coef = "PPCM_vs_NF", number = Inf, adjust.method = "BH")

# --- VIEW THE TOP GENES ---
message("--- Top 5 Genes: DCM vs NF (Corrected) ---")
print(head(res_DCM))

message("--- Top 5 Genes: HCM vs NF (Corrected) ---")
print(head(res_HCM))

message("--- Top 5 Genes: PPCM vs NF (Corrected) ---")
print(head(res_PPCM))

# Optional: Save to CSV if you need to share them
write.csv(res_DCM, "Results_DCM_vs_NF_Corrected.csv")
write.csv(res_HCM, "Results_HCM_vs_NF_Corrected.csv")
write.csv(res_PPCM, "Results_PPCM_vs_NF_Corrected.csv")


# ------------------------------------------------------------------------------
# F. VOLCANO PLOTS & COMPARISON TABLE REPORT (PNG & HTML)
# ------------------------------------------------------------------------------

# 1. Save the Corrected Volcano Plots as Images
ggsave("Volcano_DCM_Corrected.png", p1, width = 6, height = 6, dpi = 300)
ggsave("Volcano_HCM_Corrected.png", p2, width = 6, height = 6, dpi = 300)
ggsave("Volcano_PPCM_Corrected.png", p3, width = 6, height = 6, dpi = 300)

# 2. Convert Comparison Table to HTML
# We use kable to make a simple HTML string of the table to embed below
table_html <- comparison_table %>%
  kbl(align = "c", caption = "Impact of Covariate Correction on Gene Discovery") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE) %>%
  row_spec(0, bold = TRUE, color = "white", background = "#2c3e50") %>%
  column_spec(4, bold = TRUE, color = ifelse(comparison_table$Difference > 0, "green", "red"))

# 3. Define the Full HTML Content
desc_volcano <- paste(
  "<html>",
  "<head>",
  "<style>",
  "  body { font-family: Arial, sans-serif; padding: 40px; background-color: #f9f9f9; }",
  "  .container { background-color: white; padding: 30px; border-radius: 8px; box-shadow: 0 0 10px rgba(0,0,0,0.1); margin-bottom: 30px; }",
  "  h2 { color: #2c3e50; border-bottom: 2px solid #eee; padding-bottom: 10px; }",
  "  h3 { color: #e67e22; margin-top: 20px; }",
  "  p, ul { color: #555; line-height: 1.6; }",
  "  img { width: 100%; max-width: 600px; border: 1px solid #ddd; display: block; margin: 20px auto; }",
  "  table { margin: 20px auto; border-collapse: collapse; }", # Centers the table
  "</style>",
  "</head>",
  "<body>",
  
  "<h1>Differential Expression Results</h1>",
  "<p>The following report compares the Uncorrected Analysis (Task 3a) with the Corrected Analysis (Task 3b), followed by the visual results (Volcano Plots).</p>",
  
  # --- SECTION 1: COMPARISON TABLE ---
  "<div class='container'>",
  "  <h2>1. Impact of Covariate Correction (3a vs 3b)</h2>",
  "  <p>This table shows how the number of significant genes (Adj. P < 0.05) changed after correcting for Age, Gender, and RIN.</p>",
  "  <ul>",
  "    <li><strong>Negative Difference (<span style='color:red'>Red</span>):</strong> Indicates that the correction removed false positives (noise) driven by confounders like Gender.</li>",
  "    <li><strong>Positive Difference (<span style='color:green'>Green</span>):</strong> Indicates that the correction increased statistical power, revealing new genes that were previously hidden by noise.</li>",
  "  </ul>",
  "  <br>",
  # Embed the Kable HTML string directly here
  as.character(table_html),
  "</div>",
  
  # --- SECTION 2: VOLCANO PLOTS ---
  "<div class='container'>",
  "  <h2>2. Corrected Volcano Plots (Task 3b)</h2>",
  "  <p>Visualization of the final, corrected results. <strong>Red</strong> dots are Upregulated, <strong>Blue</strong> dots are Downregulated.</p>",
  
  "  <h3>DCM vs Healthy (NF)</h3>",
  "  <img src='Volcano_DCM_Corrected.png' alt='Volcano DCM'/>",
  
  "  <h3>HCM vs Healthy (NF)</h3>",
  "  <img src='Volcano_HCM_Corrected.png' alt='Volcano HCM'/>",
  
  "  <h3>PPCM vs Healthy (NF)</h3>",
  "  <p><em>Note: PPCM results are now corrected for the female-only bias.</em></p>",
  "  <img src='Volcano_PPCM_Corrected.png' alt='Volcano PPCM'/>",
  "</div>",
  
  "</body>",
  "</html>"
)

# 4. Save as HTML File
fileConn <- file("Volcano_Plots_Report.html")
writeLines(desc_volcano, fileConn)
close(fileConn)

message("Volcano Export Complete: 'Volcano_Plots_Report.html' (Includes Comparison Table)")

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# G. ENRICHMENT ANALYSIS & MANHATTAN PLOT (PNG & HTML)
# ------------------------------------------------------------------------------

# 1. Save the Static Manhattan Plot as an Image


png("Gost_Manhattan_Plot.png", width = 1000, height = 600)
publish_gostplot(p_static, highlight_terms = head(gost_res$result$term_id, 10))
dev.off() # Close the device to save the file

# 2. Convert Top Pathways Table to HTML
pathway_html <- top_pathways %>%
  kbl(caption = "Top 10 Enriched Pathways (KEGG & Reactome)", align = "c") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE) %>%
  row_spec(0, bold = TRUE, color = "white", background = "#e67e22") %>% # Orange header matching the plot
  column_spec(3, bold = TRUE) # Bold the P-value column

# 3. Define the Full HTML Content
desc_enrichment <- paste(
  "<html>",
  "<head>",
  "<style>",
  "  body { font-family: Arial, sans-serif; padding: 40px; background-color: #f9f9f9; }",
  "  .container { background-color: white; padding: 30px; border-radius: 8px; box-shadow: 0 0 10px rgba(0,0,0,0.1); margin-bottom: 30px; }",
  "  h2 { color: #2c3e50; border-bottom: 2px solid #eee; padding-bottom: 10px; }",
  "  h3 { color: #d35400; margin-top: 20px; }",
  "  p, ul { color: #555; line-height: 1.6; }",
  "  img { width: 100%; max-width: 900px; border: 1px solid #ddd; display: block; margin: 20px auto; }",
  "  table { margin: 20px auto; }",
  "</style>",
  "</head>",
  "<body>",
  
  "<h1>Functional Enrichment Analysis (GO & Pathways)</h1>",
  "<p>This report details the biological functions associated with the differentially expressed genes in DCM. The analysis was performed using g:Profiler.</p>",
  
  # --- SECTION 1: MANHATTAN PLOT ---
  "<div class='container'>",
  "  <h2>1. Manhattan Plot of Functional Terms</h2>",
  "  <p>The plot below visualizes significant biological terms found in the gene list.</p>",
  "  <ul>",
  "    <li><strong>X-Axis (Data Sources):</strong> The terms are grouped by database. <span style='color:orange'><strong>GO:BP</strong></span> (Gene Ontology Biological Process) shows broad functions. <span style='color:deeppink'><strong>KEGG/REAC</strong></span> (Right side) shows specific molecular pathways.</li>",
  "    <li><strong>Y-Axis (Significance):</strong> The height represents the significance (-log10 P-value). Higher dots are more significant.</li>",
  "    <li><strong>Split Results:</strong> The plot is split into <strong>DCM_Up</strong> (Top panel) and <strong>DCM_Down</strong> (Bottom panel) to distinguish processes being activated vs. suppressed.</li>",
  "  </ul>",
  "  <img src='Gost_Manhattan_Plot.png' alt='Manhattan Plot'/>",
  "</div>",
  
  # --- SECTION 2: PATHWAY TABLE ---
  "<div class='container'>",
  "  <h2>2. Top Enriched Pathways</h2>",
  "  <p>The table below lists the top 10 most significant signaling pathways (KEGG & Reactome) identified in the analysis.</p>",
  "  <br>",
  # Embed the Pathway Table HTML string
  as.character(pathway_html),
  "</div>",
  
  "</body>",
  "</html>"
)

# 4. Save as HTML File
fileConn <- file("Enrichment_Analysis_Report.html")
writeLines(desc_enrichment, fileConn)
close(fileConn)

message("Enrichment Export Complete: 'Enrichment_Analysis_Report.html' and 'Gost_Manhattan_Plot.png'")

 
# ------------------------------------------------------------------------------
# H. MERGE EVERYTHING: CREATE MASTER REPORT (HTML & TXT)
# ------------------------------------------------------------------------------

message("Creating Merged Master Reports...")

# --- 1. CREATE MASTER HTML REPORT ---

# We concatenate all the HTML strings we created in previous steps.
# We add a Table of Contents (TOC) at the top for navigation.
master_html <- paste(
  "<html><head>",
  "<style>",
  "  body { font-family: Arial, sans-serif; padding: 40px; background-color: #f4f4f9; color: #333; }",
  "  .section { background: white; padding: 40px; margin-bottom: 40px; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }",
  "  h1 { color: #2c3e50; text-align: center; border-bottom: 3px solid #3498db; padding-bottom: 20px; }",
  "  h2 { color: #2c3e50; margin-top: 0; }",
  "  img { max-width: 90%; border: 1px solid #ddd; display: block; margin: 20px auto; }",
  "  .toc { background: #e8f6f3; padding: 20px; border-radius: 8px; margin-bottom: 40px; border: 1px solid #a2d9ce; }",
  "  .toc a { text-decoration: none; color: #16a085; font-weight: bold; font-size: 18px; display: block; margin: 10px 0; }",
  "  .toc a:hover { text-decoration: underline; }",
  "</style>",
  "</head><body>",
  
  "<h1>Full Project Report: Gene Expression Analysis</h1>",
  
  # Table of Contents
  "<div class='toc'>",
  "  <h2>Table of Contents</h2>",
  "  <a href='#sec1'>1. Participant Characteristics</a>",
  "  <a href='#sec2'>2. Diagnostic Plots (Gene Expression)</a>",
  "  <a href='#sec3'>3. Diagnostic Plots (Clinical Data)</a>",
  "  <a href='#sec4'>4. PCA Analysis</a>",
  "  <a href='#sec5'>5. Differential Expression (Volcano Plots)</a>",
  "  <a href='#sec6'>6. Functional Enrichment</a>",
  "  <a href='#sec7'>7. Top 50 Genes List</a>",
  "</div>",
  
  # Section 1: Table 1 (We use the kable output 'html_table1')
  "<div id='sec1' class='section'>", desc_table1, as.character(html_table1), "</div>",
  
  # Section 2: QC Plots
  "<div id='sec2' class='section'>", desc_plots, "</div>",
  
  # Section 3: Clinical Plots
  "<div id='sec3' class='section'>", desc_clinical, "</div>",
  
  # Section 4: PCA
  "<div id='sec4' class='section'>", desc_pca, "</div>",
  
  # Section 5: Volcano Plots (Includes Comparison Table)
  "<div id='sec5' class='section'>", desc_volcano, "</div>",
  
  # Section 6: Enrichment
  "<div id='sec6' class='section'>", desc_enrichment, "</div>",
  
  # Section 7: Top Genes (We use the kable output 'html_genes')
  "<div id='sec7' class='section'>", 
  "<h2>7. Top 50 Differentially Expressed Genes (DCM vs NF)</h2>",
  "This table lists the most significant genes identified in the corrected analysis.",
  as.character(pathway_html), 
  "</div>",
  
  "</body></html>"
)

# Save Master HTML
fileConn <- file("Project_Full_Report.html")
writeLines(master_html, fileConn)
close(fileConn)


# --- 2. CREATE MASTER TXT REPORT ---

# Merged
clean_text <- function(html_string) {
  # Remove <style> blocks content first
  txt <- gsub("<style>.*?</style>", "", html_string)
  # Remove tags but keep content
  txt <- gsub("<[^>]+>", "", txt)
  # Fix spacing
  txt <- gsub("\\s+", " ", txt)
  return(txt)
}

# Define the Text Content
# We manually reconstruct the layout: TITLE -> DESCRIPTION -> DATA TABLE
master_txt <- paste0(
  "==============================================================================\n",
  "                     FULL PROJECT REPORT: GENE EXPRESSION ANALYSIS            \n",
  "==============================================================================\n\n",
  
  "SECTION 1: PARTICIPANT CHARACTERISTICS\n",
  "------------------------------------------------------------------------------\n",
  clean_text(desc_table1), "\n\n",
  # Add the actual data table using knitr::kable for simple text formatting
  paste(knitr::kable(table_one_formatted, format="simple"), collapse="\n"), "\n\n",
  
  "SECTION 2: DIAGNOSTIC PLOTS (GENE EXPRESSION)\n",
  "------------------------------------------------------------------------------\n",
  clean_text(desc_plots), "\n",
  "(See QC_Boxplots.png)\n\n",
  
  "SECTION 3: DIAGNOSTIC PLOTS (CLINICAL DATA)\n",
  "------------------------------------------------------------------------------\n",
  clean_text(desc_clinical), "\n",
  "(See QC_Clinical_Data.png)\n\n",
  
  "SECTION 4: PCA ANALYSIS\n",
  "------------------------------------------------------------------------------\n",
  clean_text(desc_pca), "\n",
  "(See PCA_Analysis_Plot.png)\n\n",
  
  "SECTION 5: DIFFERENTIAL EXPRESSION (VOLCANO PLOTS)\n",
  "------------------------------------------------------------------------------\n",
  clean_text(desc_volcano), "\n",
  "(See Volcano_*.png images)\n\n",
  "Comparison of Correction Methods (3a vs 3b):\n",
  paste(knitr::kable(comparison_table, format="simple"), collapse="\n"), "\n\n",
  
  "SECTION 6: FUNCTIONAL ENRICHMENT\n",
  "------------------------------------------------------------------------------\n",
  clean_text(desc_enrichment), "\n",
  "(See Gost_Manhattan_Plot.png)\n\n",
  "Top 10 Enriched Pathways:\n",
  paste(knitr::kable(top_pathways, format="simple"), collapse="\n"), "\n\n",
  
  "SECTION 7: TOP 50 GENES (DCM vs NF)\n",
  "------------------------------------------------------------------------------\n",
  paste(knitr::kable(head(top_pathways, 50), format="simple"), collapse="\n"), "\n"
)

# Save Master TXT
fileConn <- file("Project_Full_Report.txt")
writeLines(master_txt, fileConn)
close(fileConn)

message("MERGE COMPLETE.")
message(" - Project_Full_Report.html")
message(" - Project_Full_Report.txt (Plain Text Summary)")

