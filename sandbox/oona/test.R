###############################################################################
###############################################################################
# PRO4002 Research Project                                                    #
#																	                                       		  #
# Date: Jan 11, 2026											                                    #
# Author: Oona Kintscher                                                      #  
###############################################################################
###############################################################################

#=============================================================================#
#=============================================================================#
# START
#=============================================================================#
#=============================================================================#

#-----------------------------------------------------------------------------#
# Library Installation 
#-----------------------------------------------------------------------------#

if (!require("ggsci", quietly = TRUE)) { install.packages("ggsci") }
library(ggsci)

if (!require("rstudioapi", quietly = TRUE)) { install.packages("rstudioapi") }
library(rstudioapi)

if (!require("openxlsx", quietly = TRUE)) { install.packages("openxlsx") }
library(openxlsx)

if (!require("readxl", quietly = TRUE)) { install.packages("readxl") }
library(readxl)

if (!require("BiocManager", quietly = TRUE)) { install.packages("BiocManager") }
BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)
BiocManager::install("clusterProfiler")
library(clusterProfiler)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)  
BiocManager::install("RCy3")  
library(RCy3)  
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
BiocManager::install("rWikiPathways")
library(rWikiPathways)

if (!require("biomaRt", quietly = TRUE)) { install.packages("biomaRt") }
library(biomaRt)

if (!require("pcaMethods", quietly = TRUE)) { install.packages("pcaMethods") }
library(pcaMethods)

if (!require("ggplot2", quietly = TRUE)) { install.packages("ggplot2") }
library(ggplot2)

if (!require("dplyr", quietly = TRUE)) { install.packages("dplyr") }
library(dplyr)

if (!require("sva", quietly = TRUE)) { install.packages("sva") }
library(sva)

if (!require("mclust", quietly = TRUE)) { install.packages("mclust") }
library(mclust)

if (!require("limma", quietly = TRUE)) { install.packages("limma") }
library(limma)

if (!require("plyr", quietly = TRUE)) { install.packages("plyr") }
library(plyr)

if (!require("rWikiPathways", quietly = TRUE)) { install.packages("rWikiPathways") }
library(rWikiPathways)    

if (!require("RCy3", quietly = TRUE)) { install.packages("RCy3") }
library(RCy3)             

if (!require("enrichplot", quietly = TRUE)) { install.packages("enrichplot") }
library(enrichplot)        # For visualization of enrichment result

if (!require("EnhancedVolcano", quietly = TRUE)) { install.packages("EnhancedVolcano") }
library(EnhancedVolcano)

#-----------------------------------------------------------------------------#
# SETUP
#-----------------------------------------------------------------------------#

# Set Working Directory to the Location of this Script
setwd("C:/Users/oonak/Desktop/PRO4002/") 

# Verify Working Directory
getwd()

# Setup Folder System 
cache_path <- file.path(dirname(getActiveDocumentContext()$path), "cache")
dir.create(cache_path, recursive = TRUE, showWarnings = FALSE)
output_path = file.path(dirname(getActiveDocumentContext()$path),"output")
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# Setup Seed 
set.seed(1262118388)

# P Value Cutoff
pval.cutoff <- 0.05

# Log2 Fold-Change Cutoff
log2FC.cutoff <- 0.585 

#-----------------------------------------------------------------------------#
# DATA IMPORT
#-----------------------------------------------------------------------------#

# FILE DATA IMPORT 
# Import Various Data Files (txt & xlsx) - specify NA values string, sheet, and header attributes
exonLengthsData <- read.delim("data/MAGNET_exonLengths.txt", header=TRUE, row.names = 1, as.is = T)
geneExpressionData <- read.delim("data/MAGNET_GeneExpressionData_CPM_19112020.txt", header=TRUE, row.names = 1, as.is = T)
sampleDescriptionData <- read_excel("data/MAGNET_SampleData_18112022_WithDescriptions.xlsx", sheet=2, na="NA")
sampleDataRaw <- read_excel("data/MAGNET_SampleData_18112022_WithDescriptions.xlsx", sheet=1, na="NA")

# ENSEMBLE DATA IMPORT 
# Get Gene Information from ENSEMBL
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
searchAttributes(mart = ensembl)

# Extract Gene Information of Genes in the Expression Data
geneListInfo <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description", "chromosome_name", "gene_biotype", "name_1006", "definition_1006", 
                 "namespace_1003"),
  filters = "ensembl_gene_id",
  values = rownames(geneExpressionData),
  mart = ensembl
)

#-----------------------------------------------------------------------------#
# DATA PREPROCESSING 
#-----------------------------------------------------------------------------#

# Get DCM Patients 
DCM_idx <- sampleDataRaw$etiology == "DCM"
geneExpressionData.DCM <- geneExpressionData[,DCM_idx]
sampleDataRaw.DCM <- sampleDataRaw[DCM_idx,]

summary(as.vector(geneExpressionData.DCM))
geneExpressionData.DCM.CPM_filtered_idx <- rowSums(geneExpressionData.DCM > 1) >= (0.2 * ncol(geneExpressionData.DCM))
geneExpressionData.DCM.filtered <- geneExpressionData.DCM[geneExpressionData.DCM.CPM_filtered_idx, ]

cat("Genes before:", nrow(geneExpressionData.DCM), "\n")
cat("Genes after:", nrow(geneExpressionData.DCM.filtered), "\n")

table(sampleDataRaw.DCM$Library.Pool)

mod <- model.matrix(~ gender + age, data = sampleDataRaw.DCM)

geneExpressionData.DCM.combat <- ComBat(
  dat = as.matrix(geneExpressionData.DCM.filtered),
  batch = sampleDataRaw.DCM$Library.Pool,
  mod = mod,
  par.prior = TRUE
)

geneExpressionData.DCM.corrected <- limma::removeBatchEffect(
  geneExpressionData.DCM.combat,
  covariates = sampleDataRaw.DCM$age,
  design = model.matrix(~ gender, sampleDataRaw.DCM)
)

# Get Category-Specific Genes 
immune.go_id <- "GO:0002376"
geneListInfo.immuneIDs <- getBM(
  attributes = c("ensembl_gene_id"),
  filters = "go",
  values = immune.go_id,
  mart = ensembl
)

# Function to Convert logCPM Values into FPMK Values
cpm2fpkm <- function(x, y) {
  .t <- 2^(x) * 1E3 / y[, 1]
}

# Function to Center Around the Median 
center_by_median <- function(x) {
  x <- as.matrix(x)
  sweep(x, 1, matrixStats::rowMedians(x, na.rm = TRUE))
}

# Plot Enhanced Volcano Plot
plotEnhancedVolcano <- function(results_dgea, title, pval_cutoff = pval.cutoff, 
                                log2Fc_cutoff = log2FC.cutoff) {
  plot <- EnhancedVolcano(results_dgea,  title = title, labSize = 3, 
                          x = 'logFC', y = 'P.Value', 
                          lab = row.names(results_dgea),
                          pCutoff = pval_cutoff, 
                          FCcutoff = log2Fc_cutoff)
  return(plot)
}

#-----------------------------------------------------------------------------#
# CLUSTERING: COMPLETE
#-----------------------------------------------------------------------------#

# Using only the top 5000 most variable genes
geneExpressionData.DCM.mads <- apply(geneExpressionData.DCM.corrected,1,mad)
geneExpressionData.DCM.selected <- as.matrix(geneExpressionData.DCM[rev(order(geneExpressionData.DCM.mads))[1:5000],])

# Centering Gene Expression to Gene Median
geneExpressionData.DCM.selected.centered <- center_by_median(geneExpressionData.DCM.selected)

# Clustering
cluster_output_path = file.path(output_path, "ConsensusClusterPlus-Complete")
dir.create(cluster_output_path, recursive = TRUE, showWarnings = FALSE)
results = ConsensusClusterPlus(geneExpressionData.DCM.selected.centered, maxK=6, reps=100, pItem=0.8,
                               pFeature=1, title=cluster_output_path, clusterAlg="hc", distance="pearson", plot="png")


adjustedRandIndex(
  results[[2]]$consensusClass,
  results[[3]]$consensusClass
)


# PCA
pcaGeneExpressionData.DCM.selected.centered <- pca(t(geneExpressionData.DCM.selected.centered), nPcs = 10)
summary(pcaGeneExpressionData.DCM.selected.centered)

pcaSampleInfoDF.DCM <- cbind(scores(pcaGeneExpressionData.DCM.selected.centered), sampleDataRaw.DCM)
pcaSampleInfoDF.DCM$cluster <- as.factor(results[[2]]$consensusClass)

# Plot PCA Results (PC1 vs PC2)
scatterVariablePlotPCA <- ggplot(pcaSampleInfoDF.DCM, aes(PC1, PC2, color=cluster)) +
  geom_point(alpha=0.8) +
  xlab(paste("PC1 (", pcaGeneExpressionData.DCM.selected.centered@R2[1] * 100, "% of the variance)")) +
  ylab(paste("PC2 (", pcaGeneExpressionData.DCM.selected.centered@R2[2] * 100, "% of the variance)")) +
  labs(title="PC1 vs. PC2", color="Cluster Nr.") +
  theme_minimal()
scatterVariablePlotPCA

#-----------------------------------------------------------------------------#
# DEG: COMPLETE CLUSTERS 
#-----------------------------------------------------------------------------#

cluster <- revalue(as.factor(results[[2]]$consensusClass), c("1" = "c1","2" = "c2"))
design <- model.matrix(~ 0 + cluster )
colnames(design) <- levels(cluster)

fit <- lmFit(geneExpressionData.DCM.selected.centered, design)
contrast <- makeContrasts(
  C1vsC2 = c1 - c2,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

DEG.results <- topTable(fit2, adjust.method = "BH", number=Inf)

head(DEG.results)
hist(DEG.results$P.Value)

plotEnhancedVolcano(DEG.results, title="Cluster 1 vs 2")

DEG <- DEG.results[
  DEG$adj.P.Val < pval.cutoff & 
  abs(DEG$logFC) > log2FC.cutoff, 
]

DEG.up <- DEG.results[
  DEG$adj.P.Val < pval.cutoff & 
  DEG$logFC > log2FC.cutoff, 
]

DEG.down <- DEG.results[
  DEG$adj.P.Val < pval.cutoff & 
  DEG$logFC < log2FC.cutoff, 
]

cat("DEG total:", nrow(DEG), "\n")
cat(" - up regulated:", nrow(DEG.up), "\n")
cat(" - down regulated:", nrow(DEG.down), "\n")

#-----------------------------------------------------------------------------#
# PATHWAY ANALYSIS: COMPLETE CLUSTERS
#-----------------------------------------------------------------------------#


DEG.ranks <- DEG.results$logFC
names(DEG.ranks) <- rownames(DEG.results)
DEG.ranks <- sort(DEG.ranks, decreasing = TRUE)


gmt <- rWikiPathways::downloadPathwayArchive(
  organism = "Homo sapiens", 
  date = "20251110", 
  format = "gmt"
)

wp2gene <- readPathwayGMT(gmt)

wpid2gene <- wp2gene %>% dplyr::select(wpid, gene)  
wpid2name <- unique(wp2gene %>% dplyr::select(wpid, name)) 

bkgd.genes <- unique(rownames(geneExpressionData.DCM))

num.pathways <- length(unique(wpid2name$wpid))
num.genes <- length(unique(wpid2gene$gene))

bkgd.entrez <- bitr(
  bkgd.genes,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

DEG.entrez <- bitr(
  rownames(DEG),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

ewp <- clusterProfiler::enricher(
  gene = DEG.entrez$ENTREZID,
  universe = bkgd.entrez$ENTREZID,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.02,
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name
)

ewp.results <- as.data.frame(ewp)

num.pathways <- min(20, nrow(ewp.results))

if (num.pathways > 0) {
  print(
    ggplot(ewp[1:num.pathways], aes(x = reorder(Description, -pvalue), y = Count)) +
      geom_bar(stat = "identity", fill = "#BA8CD7") +
      coord_flip() +
      labs(x = "", y = "Gene count", title = "C1vsC2 DCM - Enriched Pathways") +
      theme_minimal()
  )
  
  ewp.sim <- enrichplot::pairwise_termsim(ewp)
  print(treeplot(ewp.sim, label_format = 40))
}

gene.map <- bitr(
  names(DEG.ranks),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

DEG.ranks.entrez <- DEG.ranks[gene.map$ENSEMBL]
names(DEG.ranks.entrez) <- gene.map$ENTREZID

rank.df <- data.frame(
  ENTREZID = names(DEG.ranks.entrez),
  logFC = as.numeric(DEG.ranks.entrez)
)

rank.df.collapsed <- rank.df %>%
  group_by(ENTREZID) %>%
  slice_max(order_by = abs(logFC), n = 1) %>%
  ungroup()

DEG.ranks.entrez.final <- rank.df.collapsed$logFC
names(DEG.ranks.entrez.final) <- rank.df.collapsed$ENTREZID

DEG.ranks.entrez.final <- sort(DEG.ranks.entrez.final, decreasing = TRUE)

length(DEG.ranks.entrez.final)
length(unique(names(DEG.ranks.entrez.final)))

DEG.ranks.entrez.final <- sort(DEG.ranks.entrez.final, decreasing = TRUE)

GSEA.WP <- GSEA(
  geneList = DEG.ranks.entrez.final,
  TERM2GENE = wpid2gene,
  pvalueCutoff = pval.cutoff
)

GSEA.WP@result$Description <- sapply(GSEA.WP@result$ID, function(x) {
  wp2gene %>% filter(wpid == x) %>% pull(name) %>% unique()
})

dotplot(GSEA.WP, showCategory = 20)





#-----------------------------------------------------------------------------#
# FUNCTIONAL ENRICHMENT: COMPLETE CLUSTERS 
#-----------------------------------------------------------------------------#

sig_genes <- rownames(DEG)[DEG$adj.P.Val < p_value.cutoff]

enrichment_go <- enrichGO(
  gene = sig_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH"
)

dotplot(enrichment_go)

library(glmnet)

train.idx <- createDataPartition(cluster, p = 0.7, list = FALSE)

x <- t(expr.resid[top.genes, train.idx])
y <- cluster[train.idx]

lasso <- cv.glmnet(x, y, family = "binomial", alpha = 1)

#-----------------------------------------------------------------------------#
# CLUSTERING: GO-ID-SPECIFIC GENES
#-----------------------------------------------------------------------------#

immune_genes <- geneListInfo.immuneIDs$ensembl_gene_id 
geneExpressionData.DCM.immune <- geneExpressionData.DCM[rownames(geneExpressionData.DCM) %in% immune_genes,]

# Centering Gene Expression to Gene Median
geneExpressionData.DCM.immune.centered <- center_by_median(geneExpressionData.DCM.immune)

# Clustering by Category (Immune System)
cluster_output_path = file.path(output_path,"ConsensusClusterPlus-Immune")
dir.create(cluster_output_path, recursive = TRUE, showWarnings = FALSE)
results.immune = ConsensusClusterPlus(geneExpressionData.DCM.immune.centered, maxK=6, reps=10, pItem=0.8,
                                      pFeature=1, title=cluster_output_path, clusterAlg="hc", distance="pearson",, plot="png")

#PCA
pcaGeneExpressionData.DCM.immune <- pca(t(geneExpressionData.DCM.immune), nPcs = 10)
summary(pcaGeneExpressionData.DCM.immune)
pcaSampleInfoDF.DCM.immune <- cbind(scores(pcaGeneExpressionData.DCM.immune), sampleDataRaw.DCM)
pcaSampleInfoDF.DCM.immune$cluster <- as.factor(results.immune[[2]]$consensusClass)

# Plot PCA Results (PC1 vs PC2)
scatterVariablePlotPCA.immune <- ggplot(pcaSampleInfoDF.DCM.immune, aes(PC1, PC2, color=cluster)) +
  geom_point(alpha=0.8) +
  xlab(paste("PC1 (", pcaGeneExpressionData.DCM.immune@R2[1] * 100, "% of the variance)")) +
  ylab(paste("PC2 (", pcaGeneExpressionData.DCM.immune@R2[2] * 100, "% of the variance)")) +
  labs(title="PC1 vs. PC2", color="Cluster Nr.") +
  theme_minimal()
scatterVariablePlotPCA.immune

#-----------------------------------------------------------------------------#
# CACHING
#-----------------------------------------------------------------------------#

saveRDS(geneExpressionData.DCM, file.path(cache_path, "geneExpressionData_DCM_processed.rds"))
saveRDS(geneExpressionData.DCM.immune, file.path(cache_path, "geneExpressionData_DCM_immune_processed.rds"))
