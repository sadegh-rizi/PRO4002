###############################################################################
###############################################################################
# PRO4002 Research Project                                                    #
#                                                                             #
# File: 003_PathwayEnrichment.R                                               #
# Date: Jan 22, 2026                                                          #
# Author: Sadegh, Matas, Nur, Arlin & Oona                                    #  
###############################################################################
###############################################################################

#-----------------------------------------------------------------------------#
# SETUP
#-----------------------------------------------------------------------------#

source(here::here("scripts", "000_setup.R"))

enrichment_plots_path <- file.path(plots_path, "PathwayEnrichment")
if(!dir.exists(enrichment_plots_path)) dir.create(enrichment_plots_path, recursive = TRUE)

library(tidyverse)
library(limma)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

message("\n--- Starting Pathway Enrichment ---")

#-----------------------------------------------------------------------------#
# 1. DATA IMPORT & MERGE
#-----------------------------------------------------------------------------#
message("(I) Data Import")

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



geneExpressionData <- geneExpressionData[, sampleData$sample_name]

message(paste("  Analyzing", ncol(geneExpressionData), "samples with defined subtypes."))

#-----------------------------------------------------------------------------#
# 2. LIMMA ANALYSIS (Get Gene Lists)
#-----------------------------------------------------------------------------#
message("(II) Differential Expression for Pathways")

# A. Design Matrix (Including Covariates)
# Ensure factors are set correctly
group <- factor(sampleData$subtype, levels = c("Cluster_1", "Cluster_2", "Cluster_3"))
gender <- factor(sampleData$gender)
pool <- factor(sampleData$Library.Pool)
age <- as.numeric(sampleData$age)
rin <- as.numeric(sampleData$rin)

# Check for NAs in covariates
if(any(is.na(rin))) {
  message("  Warning: NAs found in RIN. Imputing with median.")
  rin[is.na(rin)] <- median(rin, na.rm=TRUE)
}

design <- model.matrix(~ 0 + group + age + gender + rin + pool)
colnames(design) <- gsub("group", "", colnames(design))

# B. Fit Model
fit <- lmFit(geneExpressionData, design)

# C. Define Contrasts (Matching 002 Logic)
# We look for the "Unique Identity" of each cluster
cm <- makeContrasts(
  C1_Identity = Cluster_1 - (Cluster_2 + Cluster_3)/2,
  C2_Identity = Cluster_2 - (Cluster_1 + Cluster_3)/2,
  C3_Identity = Cluster_3 - (Cluster_1 + Cluster_2)/2,
  levels = design
)

fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2, trend = TRUE)

#-----------------------------------------------------------------------------#
# 3. EXTRACT UP/DOWN GENES
#-----------------------------------------------------------------------------#
message("(III) Extracting Marker Genes")

gene_lists <- list()
comparisons <- colnames(cm)

# Helper function to map Ensembl to Entrez (Better for GO)
get_entrez <- function(ensembl_ids) {
  bitr(ensembl_ids, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
}

all_genes_ensembl <- rownames(geneExpressionData)

# Convert Universe to Entrez (Just like we did for the gene lists)
all_genes_entrez <- bitr(all_genes_ensembl, 
                         fromType = "ENSEMBL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)$ENTREZID

message(paste("  Defined Universe with", length(all_genes_entrez), "genes."))

for(comp in comparisons) {
  
  # Get all results
  res <- topTable(fit2, coef = comp, number = Inf) %>% 
    rownames_to_column("EnsemblID")
  
  # Filter UP (Higher in this cluster)
  up_ids <- res %>% filter(adj.P.Val < 0.05 & logFC > 0.5) %>% pull(EnsemblID)
  
  # Filter DOWN (Lower in this cluster)
  down_ids <- res %>% filter(adj.P.Val < 0.05 & logFC < -0.5) %>% pull(EnsemblID)
  
  # Save to list (Using Entrez IDs for better GO matching)
  # If you prefer Ensembl, remove the get_entrez() wrapper
  if(length(up_ids) > 0) gene_lists[[paste0(comp, "__UP")]] <- get_entrez(up_ids)
  if(length(down_ids) > 0) gene_lists[[paste0(comp, "__DOWN")]] <- get_entrez(down_ids)
  
  message(paste0("  ", comp, ": ", length(up_ids), " UP, ", length(down_ids), " DOWN genes."))
}

#-----------------------------------------------------------------------------#
# 4. RUN GO ENRICHMENT (compareCluster)
#-----------------------------------------------------------------------------#
message("(IV) Running GO Enrichment")


ck <- compareCluster(
  geneCluster = gene_lists, 
  universe=all_genes_entrez,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",    # Changed to ENTREZID
  ont = "BP",              # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE          # Converts Entrez back to Gene Symbols in output
)

# Save the raw table
write.csv(ck@compareClusterResult, file.path(tables_path, "GO_Enrichment_Results.csv"))













#-----------------------------------------------------------------------------#
# 4b. SIMPLIFY & CATEGORIZE (The "Facet" Strategy)
#-----------------------------------------------------------------------------#
message("(IV-b) Categorizing Pathways for Faceted Plot")

# 1. Reduce Redundancy (Crucial!)
# This merges terms like "cardiac muscle contraction" and "muscle contraction"
# so you don't have 50 dots for the same thing.
ck_simplified <- simplify(ck, cutoff = 0.6, by = "p.adjust", select_fun = min)

# 2. Extract Data for Manual Categorization
df_plot <- ck_simplified@compareClusterResult %>% 
  separate(Cluster, into = c("Contrast", "Direction"), sep = "__") %>%
  mutate(
    Direction = factor(Direction, levels = c("DOWN", "UP")),
    logP = -log10(p.adjust)
  )

# 3. Clean Names (Same as before)
df_plot$Contrast <- gsub("_Identity", "", df_plot$Contrast)
df_plot$Contrast <- gsub("C1", "Cluster 1", df_plot$Contrast)
df_plot$Contrast <- gsub("C2", "Cluster 2", df_plot$Contrast)
df_plot$Contrast <- gsub("C3", "Cluster 3", df_plot$Contrast)

# 4. ASSIGN CATEGORIES (The Magic Step)
# We use keywords to force terms into the groups seen in your image.
# You can add/change keywords here based on what you find!

df_plot <- df_plot %>% mutate(Category = case_when(
  # Group 1: Mitochondria & Energy (Common in DCM)
  str_detect(Description, "(?i)mitochon|respiration|ATP|oxidative|electron|energy|metabolic|fatty") ~ "Mitochondria & Metabolism",
  
  # Group 2: Heart Structure & Contraction
  str_detect(Description, "(?i)muscle|cardiac|contract|sarcomere|myofibril|heart|actin|z-disc") ~ "Heart Function & Structure",
  
  # Group 3: Extracellular Matrix (Fibrosis)
  str_detect(Description, "(?i)matrix|collagen|adhesion|junction|fibro") ~ "Extracellular Matrix & Fibrosis",
  
  # Group 4: Immune Response
  str_detect(Description, "(?i)immune|viral|defense|b cell|t cell|leukocyte|cytokine|inflam") ~ "Immune Response",
  
  # Group 5: Translation/Protein (High turnover)
  str_detect(Description, "(?i)translat|ribosom|protein target|folding") ~ "Protein Synthesis & Processing",
  
  # Catch-all for everything else
  TRUE ~ "Signaling & Regulation"
))

# 5. FILTER (Optional but Recommended)
# Keep only the top 5 most significant terms per Category per Cluster
# This prevents the "Signaling" group from being huge.
df_filtered <- df_plot %>% 
  group_by(Contrast, Category) %>%
  slice_min(p.adjust, n = 5) %>%
  ungroup()

#-----------------------------------------------------------------------------#
# 5. VISUALIZATION (Faceted Plot)
#-----------------------------------------------------------------------------#
message("(V) Generating Faceted Dotplot")

p_faceted <- ggplot(df_filtered, aes(x = Contrast, y = Description)) +
  
  # The Dots (Dodged so Blue/Red don't overlap)
  geom_point(aes(size = logP, color = Direction), 
  ) +
  
  # THE FACETING (This creates the blocks on the right side)
  # scales="free" ensures each block only shows its own genes
  # space="free" allows big blocks to be taller than small blocks
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  
  # Colors
  scale_color_manual(values = c("DOWN" = "#377EB8", "UP" = "#E41A1C")) +
  
  # Formatting
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title = element_blank(),
    panel.grid.minor = element_blank(),
    
    # Beautify the Facet Labels (The grey boxes on the right)
    strip.text.y = element_text(angle = 0, face = "bold", size = 10),
    strip.background = element_rect(fill = "grey90", color = NA)
  ) +
  labs(
    title = "Functional Landscape of DCM Subtypes",
    size = "-log10(p-adj)", 
    color = "Regulation"
  )

# Save
ggsave(file.path(enrichment_plots_path, "GO_Enrichment_Faceted.pdf"), p_faceted, width = 14, height = 12)