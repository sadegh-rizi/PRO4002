###############################################################################
###############################################################################
# PRO4002 Research Project                                                    #
#																	                                            #
# File: 000_Setup.R                                                           #
# Date: Jan 22, 2026											                                    #
# Author: Sadegh, Matas, Nur, Arlin & Oona                                    #  
###############################################################################
###############################################################################

message("\n--- Starting Setup ---")

#-----------------------------------------------------------------------------#
# LIBRARY/PACKAGE INSTALLATION
#-----------------------------------------------------------------------------#

packages <- c( 
  "tidyverse", "tidyr", "dplyr", "gridExtra", "pcaMethods",
  "data.table", "tableone", "kableExtra", "rmarkdown",
  "readr", "readxl", "gprofiler2", "knitr", "Hmisc", "table1",
  "gt", "gtsummary", "ggsci", "edgeR", "htmltools",
  "ConsensusClusterPlus", "pheatmap", "limma", "qvalue", 
  "biomaRt", "DESeq2", "dplyr", "factoextra", "cluster", "ggpubr",
  "fpc"
)

message("(I) Loading Packages")

for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        warning(paste("Package", pkg, "is not installed."))
    } else {
        library(pkg, character.only = TRUE)
    }
}

full_str <- paste(packages, collapse = ", ")
split_lines <- strwrap(full_str, width = 80)
for (line in split_lines) {
    message("  ", line)
}

#-----------------------------------------------------------------------------#
# PROJECT PATH SETUP 
#-----------------------------------------------------------------------------#

message("(II) Setting Up Project Directory")

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) { 
  script_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(script_path)
} else {
  warning("Not running in RStudio. Please set working directory manually for portability.") 
}

base_path <- dirname(getwd())

data_path    <- file.path(base_path, "data")
plots_path   <- file.path(base_path, "plots")
results_path <- file.path(base_path, "results")
tables_path  <- file.path(base_path, "tables")
cache_path   <- file.path(base_path, "cache")

paths <- list(data_path, plots_path, results_path, tables_path, cache_path)
sapply(paths, function(p) if(!dir.exists(p)) dir.create(p, recursive = TRUE))

message(paste("  Data Directory:   ", data_path))
message(paste("  Plots Directory:  ", plots_path))
message(paste("  Results Directory:", results_path))
message(paste("  Tables Directory: ", tables_path))
message(paste("  Cache Directory:  ", cache_path))

#-----------------------------------------------------------------------------#
# STYLE SETUP
#-----------------------------------------------------------------------------#

message("(III) Setting Up Style")

# Prepare Color Values to Match Nature Color Scheme
npg_colors <- pal_npg("nrc", alpha = 1)(10)
npg_additional_colors <- colorRampPalette(npg_colors)(20)
continuous_npg_colors <- colorRampPalette(c(npg_colors[2], "white"))(100)

# GGplot Plotting Settings 
center_title <- theme(plot.title = element_text(hjust = 0.5, vjust = 1))

target_font <- "sans" 
my_style <- theme(
  text = element_text(family = target_font),        
  plot.title = element_text(hjust = 0.5, face="bold"),
  axis.title = element_text(face = "bold"),    
  legend.position = "right"                    
)

#-----------------------------------------------------------------------------#
# COMPLETE
#-----------------------------------------------------------------------------#

message("--- Finished Setup ---")
