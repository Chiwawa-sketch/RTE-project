# scripts/03_rte_definition.R

library(Seurat)
library(dplyr)
library(ggplot2)

# Load helper function
source("R/functions.R")

# 1. Load T Cell Data
# ==============================================================================
if (!file.exists("results/tcell_subset.rds")) {
  stop("Input file 'results/tcell_subset.rds' not found. Please run script 02 first.")
}
tcell <- readRDS("results/tcell_subset.rds")

# 2. Define RTE Markers
# ==============================================================================
# The 10 candidate genes from the methodology
target_genes <- c("SOX4", "IKZF2", "TOX", "TOX2", "NREP", 
                  "AUTS2", "PECAM1", "CD38", "ZMIZ1", "PTK7")

message("Calculating thresholds for RTE markers in Naive T cells...")

# Create directory for density plots
if(!dir.exists("results/density_plots")) dir.create("results/density_plots")

# 3. Loop and Process
# ==============================================================================
for (g in target_genes) {
  # This function saves the density plot to 'results/density_plots' 
  # and adds metadata to the object
  tcell <- define_threshold_by_peak(tcell, gene_symbol = g, output_dir = "results/density_plots")
}

# 4. Save Updated Object
# ==============================================================================
saveRDS(tcell, "results/tcell_with_rte.rds")
message("Step 03 Complete! Density plots saved to results/density_plots/")