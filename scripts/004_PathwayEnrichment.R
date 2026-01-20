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
library(rstudioapi)        # For setting working directory
library(readxl)            # For reading Excel files
library(EnhancedVolcano)   # For creating volcano plots
library(VennDiagram)       # For creating Venn diagrams
library(rWikiPathways)     # For accessing WikiPathways
library(dplyr)             # For data manipulation
library(RCy3)              # For Cytoscape integration
library(ragg)
message("\n--- Starting Pathway Enrichment ---")

#-----------------------------------------------------------------------------#
# 1. DATA IMPORT & MERGE
#-----------------------------------------------------------------------------#
message("(I) Data Import")

# Load Expression Data (Same as 002)
geneExpressionData <- readRDS(file.path(cache_path, "geneExpressionData_DCM_meanFiltered.rds"))
# Load Full Metadata (For RIN, Pool, etc.)
sampleData <- readRDS(file.path(cache_path, "sampleData_DCM.rds"))
geneListInfo <- readRDS(file.path(cache_path, "geneListInfo.rds"))

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
tin_median <- as.numeric()
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
ggsave(file.path(enrichment_plots_path, "GO_Enrichment_Faceted_presentation.png"), p_faceted,dpi=1200, width = 10, height = 8)


ggsave(file.path(enrichment_plots_path, "GO_Enrichment_Faceted.pdf"), p_faceted, width = 14, height = 12)


message("(III) Preparing Matrix for ComplexHeatmap")




if(!exists("ck")) ck <- readRDS(file.path(cache_path, "ck_result.rds")) # Or re-run Step 4 of 003
if(!exists("geneExpressionData")) geneExpressionData <- readRDS(file.path(cache_path, "geneExpressionData_DCM_meanFiltered.rds"))
if(!exists("sampleData")) sampleData <- readRDS(file.path(cache_path, "sampleData_DCM.rds"))

#-----------------------------------------------------------------------------#
# 2. EXTRACT GENES FROM PATHWAY CATEGORIES
#-----------------------------------------------------------------------------#
message("(I) Parsing Genes from Pathway Terms")

# We use the same dataframe 'df_filtered' you created for the Dotplot
# If it's not in memory, we recreate the categorization logic briefly here:
df_raw <- ck@compareClusterResult %>% 
  mutate(geneID = str_split(geneID, "/")) # Split "GENEA/GENEB" into lists

# Re-apply your categorization logic (Must match your Dotplot code!)
df_annotated <- df_raw %>% mutate(Category = case_when(
  str_detect(Description, "(?i)mitochon|respiration|ATP|oxidative|electron|energy|metabolic|fatty") ~ "Mitochondria & Metabolism",
  str_detect(Description, "(?i)muscle|cardiac|contract|sarcomere|myofibril|heart|actin|z-disc") ~ "Heart Function & Structure",
  str_detect(Description, "(?i)matrix|collagen|adhesion|junction|fibro") ~ "Extracellular Matrix & Fibrosis",
  str_detect(Description, "(?i)immune|viral|defense|b cell|t cell|leukocyte|cytokine|inflam") ~ "Immune Response",
  str_detect(Description, "(?i)translat|ribosom|protein target|folding") ~ "Protein Synthesis & Processing",
  TRUE ~ "Signaling & Regulation"
))

# Filter to keep only the TOP terms you visualized (e.g. Top 3 per category)
# This prevents the heatmap from having 5000 genes
top_terms <- df_annotated %>%
  group_by(Category) %>%
  slice_min(p.adjust, n = 3) %>% # Keep genes from top 3 pathways per group
  pull(Description)

df_subset <- df_annotated %>% filter(Description %in% top_terms)

#-----------------------------------------------------------------------------#
# 3. BUILD THE GENE LIST
#-----------------------------------------------------------------------------#
message("(II) Building Unique Gene List per Category")

# We need a clean table of: GeneSymbol | Category
gene_category_map <- df_subset %>%
  dplyr::select(Category, geneID) %>%
  unnest(geneID) %>%            # Expand the list so each gene is a row
  distinct(geneID, .keep_all = TRUE) # Remove duplicates (if gene is in 2 pathways, keep 1)

# Get the list of genes
target_genes <- gene_category_map$geneID

# Check if these symbols exist in your expression matrix
# (Your matrix has EnsemblIDs? Or Symbols? This assumes Rownames = Symbols)
# IF your matrix has EnsemblIDs, we need to map names first (See Note below)
valid_genes <- intersect(target_genes, rownames(geneExpressionData))

message(paste("Found", length(valid_genes), "unique genes across", length(unique(gene_category_map$Category)), "categories."))
# Map Symbols -> Ensembl (Only if needed)
symbol_to_ensembl <- geneListInfo$ensembl_gene_id[match(target_genes, geneListInfo$hgnc_symbol)]
valid_genes <- symbol_to_ensembl[!is.na(symbol_to_ensembl)]
heatmap_matrix <- geneExpressionData[valid_genes, ]
# (Then rename rows back to Symbols for the plot)
rownames(heatmap_matrix) <- geneListInfo$hgnc_symbol[match(rownames(heatmap_matrix), geneListInfo$ensembl_gene_id)]





# 1. Subset & Scale Expression Data
heatmap_matrix <- geneExpressionData[valid_genes, ]
heatmap_scaled <- t(scale(t(heatmap_matrix)))

# 2. Prepare Match Vector for Rows (Category)
# We need a vector that matches the order of rows in the matrix exactly
row_categories <- gene_category_map$Category[match(rownames(heatmap_scaled), gene_category_map$geneID)]

# 3. Define Colors (ComplexHeatmap format)
# We use the 'circlize' package for the gradient
library(ComplexHeatmap)
library(circlize)

# Heatmap Body Colors (Blue -> White -> Red)
col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

# Annotation Colors
# Note: Ensure npg_colors are defined (from your setup script)
category_colors <- c(
  "Mitochondria & Metabolism" = "#E64B35",
  "Heart Function & Structure" = "#4DBBD5",
  "Extracellular Matrix & Fibrosis" = "#00A087",
  "Immune Response" = "#3C5488",
  "Protein Synthesis & Processing" = "#F39B7F",
  "Signaling & Regulation" = "#8491B4"
)

subtype_colors <- c(
  "Cluster_1" = npg_colors[[1]], 
  "Cluster_2" = npg_colors[[2]], 
  "Cluster_3" = npg_colors[[3]]
)
#-----------------------------------------------------------------------------#
# 4. PREPARE HEATMAP DATA (ROBUST MAPPING)
#-----------------------------------------------------------------------------#
message("(III) Preparing Matrix with Correct ID Mapping")

# A. Create a Master Table: Symbol | Ensembl | Category
# We start with the symbols from your enrichment result
master_map <- gene_category_map %>%
  dplyr::rename(Symbol = geneID) %>%
  # Add Ensembl IDs from your geneListInfo
  left_join(geneListInfo, by = c("Symbol" = "hgnc_symbol")) %>%
  # Filter: Keep only genes that actually exist in your Expression Data
  filter(ensembl_gene_id %in% rownames(geneExpressionData)) %>%
  # Deduplicate: If a gene maps to multiple categories, keep the first one (or most relevant)
  distinct(ensembl_gene_id, .keep_all = TRUE)

message(paste("Identified", nrow(master_map), "genes matching between Pathway terms and Expression Data."))

# B. Subset Expression Matrix using Ensembl IDs (The stable ID)
heatmap_matrix <- geneExpressionData[master_map$ensembl_gene_id, ]

# C. Set Rownames to Symbols (for the plot)
# We use make.unique to handle duplicates like "GeneA.1" without breaking the map
plot_symbols <- make.unique(master_map$Symbol)
rownames(heatmap_matrix) <- plot_symbols

# D. Create the Category Vector (The Fix)
# Now we rely on the master_map order, which matches the matrix exactly.
row_categories <- master_map$Category

# E. Scale the Data
heatmap_scaled <- t(scale(t(heatmap_matrix)))

#-----------------------------------------------------------------------------#
# 5. GENERATE COMPLEX HEATMAP
#-----------------------------------------------------------------------------#
message("(IV) Generating ComplexHeatmap")

library(ComplexHeatmap)
library(circlize)

# Colors
col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

category_colors <- c(
  "Mitochondria & Metabolism" = "#E64B35",
  "Heart Function & Structure" = "#4DBBD5",
  "Extracellular Matrix & Fibrosis" = "#00A087",
  "Immune Response" = "#3C5488",
  "Protein Synthesis & Processing" = "#F39B7F",
  "Signaling & Regulation" = "#8491B4"
)

subtype_colors <- c(
  "Cluster_1" = npg_colors[[1]], 
  "Cluster_2" = npg_colors[[2]], 
  "Cluster_3" = npg_colors[[3]]
)

# Annotations
ha_top <- HeatmapAnnotation(
  Subtype = sampleData$subtype,
  col = list(Subtype = subtype_colors),
  simple_anno_size = unit(0.5, "cm"),
  show_annotation_name = TRUE
)

ha_left <- rowAnnotation(
  Category = row_categories,
  col = list(Category = category_colors),
  show_legend = FALSE, 
  width = unit(5, "mm")
)

# Draw Plot
pdf(file.path(enrichment_plots_path, "Pathway_Genes_Heatmap_Complex.pdf"), width = 12, height = 14)

ht <- Heatmap(
  heatmap_scaled,
  name = "Z-Score",
  col = col_fun,
  
  # Attach Annotations
  top_annotation = ha_top,
  left_annotation = ha_left,
  
  # SPLIT: This groups the heatmap blocks by Category
  row_split = row_categories,
  column_split = sampleData$subtype,
  
  # Settings
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  
  # Fix Text Cutting
  row_title_rot = 0,             # Horizontal Titles
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_names_gp = gpar(fontsize = 8),
  
  # Gaps
  row_gap = unit(2, "mm"),
  column_gap = unit(1, "mm")
)

draw(ht, 
     merge_legend = TRUE, 
     column_title = "Gene Expression Landscape of Enriched Pathways", 
     column_title_gp = gpar(fontsize = 16, fontface = "bold")
)
dev.off()

message("Saved: Pathway_Genes_Heatmap_Complex.pdf")

message("Saving PNG version...")

# Define the filename

png_filename <- file.path(enrichment_plots_path, "Pathway_Genes_Heatmap_Complex.png")

# Open PNG device
# width/height are in inches (same as PDF)
# res=300 is standard for publication quality (prevents blurry text)

agg_png(png_filename, 
        width = 8, height = 7, units = "in", 
        res = 1200, 
        scaling = 1) # scaling=1 ensures text size stays consistent
# Draw the exact same heatmap object 'ht'
draw(ht, 
     merge_legend = TRUE, 
     column_title = "Gene Expression Landscape of Enriched Pathways", 
     column_title_gp = gpar(fontsize = 16, fontface = "bold")
)

# Close the device
dev.off()

message("Saved: Pathway_Genes_Heatmap_Complex.png")




png_filename <- file.path(enrichment_plots_path, "Pathway_Genes_Heatmap_Complex_Presentation.png")

agg_png(png_filename, 
        width = 8, height = 7, units = "in", 
        res = 1200, 
        scaling = 1) # scaling=1 ensures text size stays consistent
# Draw the exact same heatmap object 'ht'
draw(ht, 
     merge_legend = TRUE, 
     column_title = "Gene Expression Landscape of Enriched Pathways", 
     column_title_gp = gpar(fontsize = 16, fontface = "bold")
)

# Close the device
dev.off()

message("Saved: Pathway_Genes_Heatmap_Complex_Presentation.png")









#-----------------------------------------------------------------------------#
# Upset plots
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
# UPSET PLOTS: VISUALIZE GENE OVERLAPS
#-----------------------------------------------------------------------------#
message("(III-b) Generating UpSet Plots")

library(UpSetR)

# 1. PREPARE DATA
# We need to separate the 'gene_lists' into UP and DOWN sets
# and rename them for the plot (e.g., "C1_Identity__UP" -> "Cluster 1")

# Extract UP sets
up_list <- gene_lists[grep("__UP", names(gene_lists))]
names(up_list) <- gsub("_Identity__UP", "", names(up_list))
names(up_list) <- gsub("C1", "Cluster 1", names(up_list))
names(up_list) <- gsub("C2", "Cluster 2", names(up_list))
names(up_list) <- gsub("C3", "Cluster 3", names(up_list))

# Extract DOWN sets
down_list <- gene_lists[grep("__DOWN", names(gene_lists))]
names(down_list) <- gsub("_Identity__DOWN", "", names(down_list))
names(down_list) <- gsub("C1", "Cluster 1", names(down_list))
names(down_list) <- gsub("C2", "Cluster 2", names(down_list))
names(down_list) <- gsub("C3", "Cluster 3", names(down_list))

#-----------------------------------------------------------------------------#
# 2. PLOT UPREGULATED INTERSECTIONS
#-----------------------------------------------------------------------------#

pdf(file.path(plots_path, "UpSet_Upregulated_Genes.pdf"), width = 8, height = 6, onefile=FALSE)
upset(fromList(up_list), 
      nsets = 3, 
      order.by = "freq", 
      empty.intersections = "on",
      mainbar.y.label = "Number of Upregulated Genes",
      sets.x.label = "Total Genes per Cluster",
      text.scale = c(1.5, 1.2, 1.2, 1, 1.5, 1.3), 
      main.bar.color = "#E41A1C", 
      sets.bar.color = "#E41A1C",
      queries = list(
        list(query = intersects, params = list("Cluster 1"), color = npg_colors[[1]], active = T),
        list(query = intersects, params = list("Cluster 2"), color = npg_colors[[2]], active = T),
        list(query = intersects, params = list("Cluster 3"), color = npg_colors[[3]], active = T)
      )
)
dev.off() # Close PDF

# --- Save PNG ---
png(file.path(plots_path, "UpSet_Upregulated_Genes.png"), width = 6, height = 6, units = "in", res = 600)
upset(fromList(up_list), 
      nsets = 3, 
      order.by = "freq", 
      empty.intersections = "on",
      mainbar.y.label = "Number of Upregulated Genes",
      sets.x.label = "Total Genes per Cluster",
      text.scale = c(1.5, 1.2, 1.2, 1, 1.5, 1.3), 
      main.bar.color = "#E41A1C", 
      sets.bar.color = "#E41A1C",
      queries = list(
        list(query = intersects, params = list("Cluster 1"), color = npg_colors[[1]], active = T),
        list(query = intersects, params = list("Cluster 2"), color = npg_colors[[2]], active = T),
        list(query = intersects, params = list("Cluster 3"), color = npg_colors[[3]], active = T)
      )
)
dev.off() # Close PNG

# -----------------------------------------------------------------------------
# 3. PLOT DOWNREGULATED INTERSECTIONS
# -----------------------------------------------------------------------------

# --- Save PDF ---
pdf(file.path(plots_path, "UpSet_Downregulated_Genes.pdf"), width = 8, height = 6, onefile=FALSE)
upset(fromList(down_list), 
      nsets = 3, 
      order.by = "freq",
      empty.intersections = "on",
      mainbar.y.label = "Number of Downregulated Genes",
      sets.x.label = "Total Genes per Cluster",
      text.scale = c(1.5, 1.2, 1.2, 1, 1.5, 1.3),
      main.bar.color = "#377EB8", 
      sets.bar.color = "#377EB8",
      queries = list(
        list(query = intersects, params = list("Cluster 1"), color = npg_colors[[1]], active = T),
        list(query = intersects, params = list("Cluster 2"), color = npg_colors[[2]], active = T),
        list(query = intersects, params = list("Cluster 3"), color = npg_colors[[3]], active = T)
      )
)
dev.off() # Close PDF

# --- Save PNG ---
png(file.path(plots_path, "UpSet_Downregulated_Genes.png"), width = 6, height = 6, units = "in", res = 600)
upset(fromList(down_list), 
      nsets = 3, 
      order.by = "freq",
      empty.intersections = "on",
      mainbar.y.label = "Number of Downregulated Genes",
      sets.x.label = "Total Genes per Cluster",
      text.scale = c(1.5, 1.2, 1.2, 1, 1.5, 1.3),
      main.bar.color = "#377EB8", 
      sets.bar.color = "#377EB8",
      queries = list(
        list(query = intersects, params = list("Cluster 1"), color = npg_colors[[1]], active = T),
        list(query = intersects, params = list("Cluster 2"), color = npg_colors[[2]], active = T),
        list(query = intersects, params = list("Cluster 3"), color = npg_colors[[3]], active = T)
      )
)
dev.off() # Close PNG

message("UpSet plots saved to plots/UpSet_*.pdf and *.png")



#

#-----------------------------------------------------------------------------#
# Cytoscape mapping
#-----------------------------------------------------------------------------#

# Test connection to Cytoscape
RCy3::cytoscapePing()

# Install required Cytoscape apps if not already installed
if (!grepl("status: Installed", RCy3::getAppStatus("wikipathways"))) {
  RCy3::installApp("wikipathways")
}

if (!grepl("status: Installed", RCy3::getAppStatus("enhancedGraphics"))) {
  RCy3::installApp("enhancedGraphics")
}

cat("Cytoscape is ready!\n")


message("(I) Preparing Differential Expression Data")

# We need a single table with: EnsemblID | LogFC_C1 | LogFC_C2 | LogFC_C3

# Extract logFC for each contrast
# Note: Ensure 'fit2' exists. If not, re-run the limma step from 003.
c1_data <- topTable(fit2, coef = "C1_Identity", number = Inf, sort.by = "none") %>% dplyr::select(logFC) %>% dplyr::rename(C1_logFC = logFC)
c2_data <- topTable(fit2, coef = "C2_Identity", number = Inf, sort.by = "none") %>% dplyr::select(logFC) %>% dplyr::rename(C2_logFC = logFC)
c3_data <- topTable(fit2, coef = "C3_Identity", number = Inf, sort.by = "none") %>% dplyr::select(logFC) %>% dplyr::rename(C3_logFC = logFC)

# Combine columns
cytoscape_data <- bind_cols(c1_data, c2_data, c3_data)
cytoscape_data$Ensembl_ID <- rownames(c1_data)

# CRITICAL: Clean Ensembl IDs (Remove version suffix .1, .2)
# Cytoscape usually expects pure IDs (ENSG000001...)
cytoscape_data$Ensembl_Clean <- gsub("\\..*", "", cytoscape_data$Ensembl_ID)

# View check
print(head(cytoscape_data))

#-----------------------------------------------------------------------------#
# II. IMPORT PATHWAY (Oxidative Phosphorylation)
#-----------------------------------------------------------------------------#
message("(II) Importing Pathway")

# We use WP111 (Electron Transport Chain) as discussed.
# You can swap this for "WP383" (Striated Muscle Contraction) if you prefer.
wp_id <- "WP111" 

commandsRun(paste0('wikipathways import-as-pathway id=', wp_id)) 

# Toggle graphics to make it look sharp
toggleGraphicsDetails()

#-----------------------------------------------------------------------------#
# III. LOAD DATA INTO CYTOSCAPE
#-----------------------------------------------------------------------------#
message("(III) Loading Data onto Nodes")

# We map our "Ensembl_Clean" column to the "Ensembl" column in the Pathway
loadTableData(
  cytoscape_data, 
  data.key.column = "Ensembl_Clean", # Our column
  table.key.column = "Ensembl"       # Cytoscape's column (Standard for WP)
)

#-----------------------------------------------------------------------------#
# IV. APPLY HEATMAP STYLE
#-----------------------------------------------------------------------------#
message("(IV) Applying Heatmap Style")

# 1. Create a new style based on the WikiPathways default
style_name <- "DCM_Subtype_Style"
copyVisualStyle("WikiPathways", style_name)
setVisualStyle(style_name)

# 2. Define the Columns to Visualize
# These must match the column names in 'cytoscape_data'
columns_to_map <- c("C1_logFC", "C2_logFC", "C3_logFC")

# 3. Create the Heatmap (Node Custom Graphics)
# Colors: Blue (-1) -> White (0) -> Red (1)
setNodeCustomHeatMapChart(
  columns = columns_to_map,
  slot = 2,                     # Slot 2 puts it right on the node body
  style.name = style_name,
  colors = c("#E41A1C", "#FFFFFF","#377EB8" ,"#E4E3E3"), # Matches your previous plots
  range = c(-1.5, 1.5)          # Cap the color scale at +/- 1.5 logFC
)

#-----------------------------------------------------------------------------#
# V. EXPORT
#-----------------------------------------------------------------------------#
message("(V) Exporting Image")

# Fit content to screen
fitContent()

# Save
output_filename <- file.path(enrichment_plots_path, paste0("Cytoscape_", wp_id, "_Subtypes.png"))
exportImage(output_filename, type = "PNG", zoom = 300)

message(paste("Saved visualization to:", output_filename))



