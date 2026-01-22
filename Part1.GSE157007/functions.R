# R/functions.R

library(ggplot2)
library(dplyr)
library(Seurat)

#' Plot Density and Define Threshold by Peak
#'
#' This function visualizes the expression density of a specific gene in 
#' Naive T cells, identifies the peak density, and sets a threshold 
#' to classify cells as Positive or Negative.
#'
#' @param seurat_obj A Seurat object containing the subsetted data.
#' @param gene_symbol Character string. The name of the gene to analyze.
#' @param output_dir Character string. Directory to save the plot.
#' 
#' @return The updated Seurat object with new metadata columns for expression status.
define_threshold_by_peak <- function(seurat_obj, gene_symbol, output_dir = "results") {
  
  # 1. Check if gene exists in the dataset
  if (!gene_symbol %in% rownames(seurat_obj)) {
    warning(paste("Gene", gene_symbol, "not found in dataset. Skipping."))
    return(seurat_obj)
  }
  
  # 2. Extract expression data for Naive T cells
  # We use the 'data' slot (normalized counts)
  expr_vec <- as.numeric(GetAssayData(seurat_obj, assay = "RNA", layer = "data")[gene_symbol, ])
  
  meta_df <- seurat_obj@meta.data %>%
    mutate(expr = expr_vec) %>%
    filter(celltype == "Naive T", expr > 0)
  
  if (nrow(meta_df) == 0) {
    warning(paste("No non-zero expression for", gene_symbol, "in Naive T"))
    return(seurat_obj)
  }
  
  # 3. Calculate density and find the peak
  dens <- density(meta_df$expr)
  peak_x <- dens$x[which.max(dens$y)]
  
  # 4. Update Metadata (Classify cells)
  col_status <- gene_symbol # Column name will be the gene symbol (e.g., "SOX4")
  
  # Initialize with NA
  seurat_obj@meta.data[[col_status]] <- NA
  
  # Identify Naive T cells
  is_naive <- seurat_obj@meta.data$celltype == "Naive T"
  
  # Logic:
  # - Expression > Peak -> Positive
  # - Expression <= Peak -> Negative
  # - Zero expression -> Negative (implicitly included in <= Peak usually, but handled carefully here)
  
  # Set Positive
  seurat_obj@meta.data[[col_status]][is_naive & expr_vec > peak_x] <- paste0(gene_symbol, "+ Naive T")
  
  # Set Negative (includes low expression and zero expression within Naive T)
  seurat_obj@meta.data[[col_status]][is_naive & expr_vec <= peak_x] <- paste0(gene_symbol, "- Naive T")
  
  # 5. Generate and Save Plot
  annot_y <- max(dens$y) * 0.9
  p <- ggplot(meta_df, aes(x = expr)) +
    geom_density(fill = "steelblue", alpha = 0.4) +
    geom_vline(xintercept = peak_x, color = "red", linetype = "dashed") +
    annotate("text", x = peak_x, y = annot_y, 
             label = sprintf("Peak: %.2f", peak_x), color = "red", hjust = -0.1) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_minimal() +
    labs(
      title = paste0(gene_symbol, " Expression in Naive T (Threshold Definition)"),
      x = "Expression Level (LogNormalized)",
      y = "Density"
    )
  
  # Ensure output directory exists
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  ggsave(filename = file.path(output_dir, paste0("Density_", gene_symbol, ".pdf")), 
         plot = p, width = 6, height = 4)
  message(paste("Processed:", gene_symbol, "| Threshold set at:", round(peak_x, 2)))
  return(seurat_obj)
}