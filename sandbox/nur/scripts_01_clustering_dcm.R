
# MAGNet: DCM clustering starter script
# - Reads metadata + logCPM expression matrix
# - Subsets to DCM samples
# - Filters genes by variance
# - PCA -> clustering (k-means on PCs)
# - K using silhouette + gap statistic

#Setup
install.packages(c(
  "data.table",
  "tidyverse",
  "cluster",
  "factoextra",
  "here"
))
suppressPackageStartupMessages({
  library(data.table)   
  library(tidyverse)    
  library(matrixStats)  
  library(cluster)      
  library(factoextra)   
  library(here)         
})

suppressPackageStartupMessages({
  library(here)
})

data_dir_candidates <- c(
  here::here("data"),
  file.path(path.expand("~"),
            "Library", "Mobile Documents", "com~apple~CloudDocs",
            "Desktop", basename(getwd()), "data")
)

data_dir <- NA_character_
for (d in data_dir_candidates) {
  if (dir.exists(d) && length(list.files(d)) > 0) {
    data_dir <- d
    break
  }
}

if (is.na(data_dir)) {
  stop(
    paste0(
      "Could not find MAGNet files.\n\n",
      "Tried these data folders:\n- ",
      paste(data_dir_candidates, collapse = "\n- "),
      "\n\nFix:\n",
      "1) Make sure your files are in <project>/data/\n",
      "2) Or move the whole project out of iCloud Desktop (recommended).\n",
      "3) Then reopen the .Rproj and re-run.\n"
    )
  )
}

message("Using data directory: ", data_dir)
message("Files found: ", paste(list.files(data_dir), collapse = ", "))


meta_path <- file.path(data_dir, "MAGNET_SampleData_18112022.csv")
expr_path <- file.path(data_dir, "MAGNET_GeneExpressionData_CPM_19112020.txt")
exon_path <- file.path(data_dir, "MAGNET_exonLengths.txt")  

if (!file.exists(meta_path)) stop("Missing metadata file: ", meta_path)
if (!file.exists(expr_path)) stop("Missing expression file: ", expr_path)


dir.create(here("results"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("figures"), showWarnings = FALSE, recursive = TRUE)

set.seed(1005)  # reproducibility for k-means / gap statistic

# 1) Load data

meta_path <- here("data", "MAGNET_SampleData_18112022.csv")
expr_path <- here("data", "MAGNET_GeneExpressionData_CPM_19112020.txt")

metadata <- fread(meta_path) %>% as.data.frame()

# Expression: genes x samples; first column is EnsemblGeneID
expr_dt <- fread(expr_path)
expr <- expr_dt %>% as.data.frame()

# Ensure first column name is consistent
stopifnot(colnames(expr)[1] %in% c("EnsemblGeneID", '"EnsemblGeneID"'))

expr$EnsemblGeneID <- gsub('"', "", expr$EnsemblGeneID)
expr$EnsemblGeneID <- sub("\\.\\d+$", "", expr$EnsemblGeneID)

# 2) Align samples between metadata and expression

stopifnot("sample_name" %in% colnames(metadata))

expr_sample_ids <- colnames(expr)[-1]
expr_sample_ids <- gsub('"', "", expr_sample_ids)  
colnames(expr)[-1] <- expr_sample_ids

common_samples <- intersect(metadata$sample_name, expr_sample_ids)

if (length(common_samples) < 10) {
  stop("Too few overlapping samples between metadata and expression. Check sample IDs.")
}

metadata2 <- metadata %>%
  filter(sample_name %in% common_samples) %>%
  arrange(match(sample_name, common_samples))

expr2 <- expr %>%
  select(EnsemblGeneID, all_of(metadata2$sample_name))

# Convert expression to numeric matrix (genes as rows, samples as columns)
expr_mat <- as.matrix(expr2[, -1])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- expr2$EnsemblGeneID

stopifnot(ncol(expr_mat) == nrow(metadata2))

# 3) Subset to DCM

stopifnot("etiology" %in% colnames(metadata2))

meta_dcm <- metadata2 %>% filter(etiology == "DCM")
dcm_samples <- meta_dcm$sample_name

expr_dcm <- expr_mat[, dcm_samples, drop = FALSE]

message("DCM samples: ", ncol(expr_dcm))
message("Genes: ", nrow(expr_dcm))

# 4) Feature selection for clustering

n_var_genes <- 3000

gene_vars <- rowVars(expr_dcm)
names(gene_vars) <- rownames(expr_dcm)

# Remove genes with NA variance 
gene_vars <- gene_vars[!is.na(gene_vars)]

top_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(min(n_var_genes, length(gene_vars)))]
expr_top <- expr_dcm[top_genes, , drop = FALSE]

# Transpose for PCA/clustering: samples as rows, features as columns
X <- t(expr_top)

# Center/scale genes so each gene contributes comparably
X_scaled <- scale(X)

# 5) PCA
pca <- prcomp(X_scaled, center = TRUE, scale. = FALSE)

pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("sample_name") %>%
  left_join(meta_dcm, by = "sample_name")

# Percent variance explained
pve <- (pca$sdev^2) / sum(pca$sdev^2)

# PCA plot (PC1 vs PC2)
p_pca <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 2, alpha = 0.9) +
  labs(
    title = "DCM PCA (top variable genes, scaled)",
    subtitle = sprintf("PC1: %.1f%%, PC2: %.1f%% variance", 100*pve[1], 100*pve[2])
  ) +
  theme_minimal()

ggsave(here("figures", "dcm_pca_pc1_pc2.png"), p_pca, width = 7, height = 5, dpi = 300)

# Optional: check confounding (maybe?)
# ggplot(pca_df, aes(PC1, PC2, color = Library.Pool)) + geom_point(size=2) + theme_minimal()
# ggplot(pca_df, aes(PC1, PC2, color = tissue_source)) + geom_point(size=2) + theme_minimal()
# ggplot(pca_df, aes(PC1, PC2, color = RIN)) + geom_point(size=2) + theme_minimal()

# 6) Clustering on PCA space

n_pcs <- 20
pc_mat <- pca$x[, seq_len(min(n_pcs, ncol(pca$x))), drop = FALSE]

# 7) Choose number of clusters K (2..6)
k_grid <- 2:6

# 7a) Silhouette for k-means
sil_scores <- map_dbl(k_grid, function(k) {
  km <- kmeans(pc_mat, centers = k, nstart = 50)
  ss <- silhouette(km$cluster, dist(pc_mat))
  mean(ss[, 3])
})

sil_df <- tibble(k = k_grid, silhouette = sil_scores)

p_sil <- ggplot(sil_df, aes(k, silhouette)) +
  geom_line() +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "K selection (k-means on PCs)", subtitle = "Average silhouette width")

ggsave(here("figures", "dcm_k_selection_silhouette.png"), p_sil, width = 7, height = 5, dpi = 300)

# 7b) Gap statistic 

gap <- clusGap(pc_mat, FUN = kmeans, K.max = max(k_grid), B = 100)

# Extract gap values for k=1..K.max 
gap_df <- as.data.frame(gap$Tab) %>%
  rownames_to_column("k") %>%
  mutate(k = as.integer(k)) %>%
  filter(k %in% k_grid) %>%
  rename(gap = gap)

p_gap <- ggplot(gap_df, aes(k, gap)) +
  geom_line() +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "K selection (Gap statistic)", subtitle = "Higher is better (rule-of-thumb)")

ggsave(here("figures", "dcm_k_selection_gap.png"), p_gap, width = 7, height = 5, dpi = 300)

# Decide K:
best_k_sil <- sil_df$k[which.max(sil_df$silhouette)]
message("Best K by silhouette: ", best_k_sil)

# 8) Final k-means clustering using chosen K

K <- best_k_sil 
km_final <- kmeans(pc_mat, centers = K, nstart = 100)

cluster_df <- tibble(
  sample_name = rownames(pc_mat),
  cluster = paste0("C", km_final$cluster)
) %>%
  left_join(meta_dcm, by = "sample_name")

write.csv(cluster_df, here("results", "dcm_cluster_assignments.csv"), row.names = FALSE)

# Visualize clusters on PCA plot
p_pca_cluster <- ggplot(pca_df %>% left_join(cluster_df %>% select(sample_name, cluster), by="sample_name"),
                        aes(PC1, PC2, color = cluster)) +
  geom_point(size = 2, alpha = 0.9) +
  theme_minimal() +
  labs(
    title = "DCM clusters projected onto PCA",
    subtitle = paste0("k-means on first ", n_pcs, " PCs; K = ", K)
  )

ggsave(here("figures", "dcm_pca_clusters.png"), p_pca_cluster, width = 7, height = 5, dpi = 300)

# 9) Quick clinical association table 

summary_table <- cluster_df %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    age_mean = mean(age, na.rm = TRUE),
    LVEF_mean = mean(LVEF, na.rm = TRUE),
    RIN_mean = mean(RIN, na.rm = TRUE),
    pct_female = mean(gender == "F", na.rm = TRUE) * 100
  ) %>%
  arrange(cluster)

write.csv(summary_table, here("results", "dcm_cluster_clinical_summary.csv"), row.names = FALSE)

message("Done. Outputs written to /figures and /results.")
