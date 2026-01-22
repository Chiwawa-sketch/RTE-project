# scripts/Part2_CrossDisease/05_rte_analysis.R
library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)
library(patchwork)
# ==============================================================================
# 1. Configuration
# ==============================================================================
BASE_DIR <- "/path/to/ADT4DISEASES" # UPDATE THIS
OUTPUT_DIR <- file.path(BASE_DIR, "results/part2_cross_disease")
PLOT_DIR <- file.path(OUTPUT_DIR, "NaiveTPlots")
if(!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)
# Load Annotated Data from Step 04
merged_all <- readRDS(file.path(OUTPUT_DIR, "merged_all_annotated.rds"))
# ==============================================================================
# 2. Define Thresholding Function (Valley Detection)
# ==============================================================================
get_valley_threshold <- function(vec) {
  dens <- density(vec)
  dy <- diff(dens$y)
  sign_change <- diff(sign(dy))
  # Find peaks (change from + to -) and valleys (change from - to +)
  # Here we look for local minima (valleys)
  valley_idx <- which(sign_change == 2) + 1
  
  if (length(valley_idx) >= 1) {
    # If valley exists, pick the one closest to the center of the distribution
    center_idx <- valley_idx[which.min(abs(valley_idx - length(dens$x) / 2))]
    return(list(threshold = dens$x[center_idx], dens = dens))
  } else {
    # Fallback: if no clear valley, might pick peak or specific quantile
    # (Here returning peak as fallback based on provided logic)
    return(list(threshold = dens$x[which.max(dens$y)], dens = dens))
  }
}

# ==============================================================================
# 3. Calculate Thresholds & Define +/- Populations
# ==============================================================================
message("Calculating thresholds for RTE markers...")

target_genes <- c("SOX4", "IKZF2", "TOX", "TOX2", "NREP", 
                  "AUTS2", "PECAM1", "CD38", "ZMIZ1", "PTK7")

meta <- merged_all@meta.data
rna_data <- GetAssayData(merged_all, assay = "RNA", slot = "data")

for (gene in target_genes) {
  # Check if gene exists
  if (!gene %in% rownames(rna_data)) next
  expr_vec <- as.numeric(rna_data[gene, ])
  
  for (ctype in c("CD4 Naive", "CD8 Naive")) {
    idx_all <- which(meta$cell_type == ctype)
    # Only use non-zero expression for density calculation
    idx_plot <- which(meta$cell_type == ctype & expr_vec > 0)
    expr_vals_plot <- expr_vec[idx_plot]
    
    if (length(expr_vals_plot) == 0) {
      warning(paste(gene, "-", ctype, ": No expression found. Marking all as negative."))
      threshold <- 0
      dens_y <- 1
    } else {
      res <- get_valley_threshold(expr_vals_plot)
      threshold <- res$threshold
      dens_y <- max(res$dens$y)
    }
    
    # Define Column Names
    expr_col  <- paste0("expr_", gene, "_", gsub(" ", "_", ctype))
    label_col <- paste0(gene, "_", gsub(" ", "_", ctype))
    
    # Store Expression
    meta[[expr_col]] <- NA
    meta[[expr_col]][idx_all] <- expr_vec[idx_all]
    
    # Store Label (+/-)
    meta[[label_col]] <- NA
    # Logic: > threshold is Positive
    meta[[label_col]][idx_all[expr_vec[idx_all] > threshold]]  <- paste0(gene, "+ ", ctype)
    # Logic: <= threshold is Negative
    meta[[label_col]][idx_all[expr_vec[idx_all] <= threshold]] <- paste0(gene, "- ", ctype)
    
    # Save Density Plot
    if (length(expr_vals_plot) > 0) {
      p <- ggplot(data.frame(expr = expr_vals_plot), aes(x = expr)) +
        geom_density(fill = "steelblue", alpha = 0.4) +
        geom_vline(xintercept = threshold, color = "red", linetype = "dashed") +
        annotate("text", x = threshold, y = dens_y * 0.9, 
                 label = sprintf("Thresh: %.2f", threshold), color = "red") +
        theme_minimal() +
        labs(title = paste0(gene, " in ", ctype), x = "Expression", y = "Density")
      
      ggsave(file.path(PLOT_DIR, paste0("Density_", gene, "_", gsub(" ", "_", ctype), ".pdf")), 
             plot = p, width = 6, height = 4)
    }
  }
}

# Update metadata in Seurat object
merged_all@meta.data <- meta

# ==============================================================================
# 4. Export Statistics (Proportions)
# ==============================================================================
message("Exporting proportion CSVs...")

# Filter for Naive T cells only
meta_naive <- meta %>% filter(cell_type %in% c("CD4 Naive", "CD8 Naive"))

for (ctype in c("CD4 Naive", "CD8 Naive")) {
  meta_sub <- meta_naive %>% filter(cell_type == ctype)
  short_type <- gsub(" ", "_", ctype)
  
  # A. Proportion by Status (Group Level)
  status_list <- lapply(target_genes, function(gene) {
    label_col <- paste0(gene, "_", short_type)
    if (!label_col %in% colnames(meta_sub)) return(NULL)
    
    meta_sub %>%
      group_by(Status, GeneStatus = .data[[label_col]]) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(Status) %>%
      mutate(prop = round(100 * n / sum(n), 2), Gene = gene)
  })
  write_csv(bind_rows(status_list), file.path(PLOT_DIR, paste0("Prop_Status_", short_type, ".csv")))
  
  # B. Proportion by Subject (Individual Level)
  subject_list <- lapply(target_genes, function(gene) {
    label_col <- paste0(gene, "_", short_type)
    if (!label_col %in% colnames(meta_sub)) return(NULL)
    
    meta_sub %>%
      group_by(Subject, GeneStatus = .data[[label_col]]) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(Subject) %>%
      mutate(prop = round(100 * n / sum(n), 2), Gene = gene)
  })
  write_csv(bind_rows(subject_list), file.path(PLOT_DIR, paste0("Prop_Subject_", short_type, ".csv")))
}

# ==============================================================================
# 5. Define RTE (SOX4-based) & Update Cell Types
# ==============================================================================
message("Applying RTE Definition (SOX4+)...")

# IMPORTANT: We explicitly rename cells based on SOX4 status
# CD4 Naive + SOX4 High -> CD4 RTE
# CD8 Naive + SOX4 High -> CD8 RTE

meta$cell_type[meta$cell_type == "CD4 Naive" & meta$SOX4_CD4_Naive == "SOX4+ CD4 Naive"] <- "CD4 RTE"
meta$cell_type[meta$cell_type == "CD8 Naive" & meta$SOX4_CD8_Naive == "SOX4+ CD8 Naive"] <- "CD8 RTE"

merged_all@meta.data <- meta
Idents(merged_all) <- "cell_type"

# ==============================================================================
# 6. Visualization 1: UMAP with RTE Overlay
# ==============================================================================
message("Generating RTE Overlay UMAP...")

umap_df <- Embeddings(merged_all, "umap") %>% as.data.frame()
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$cell <- rownames(umap_df)
umap_df$cell_type <- merged_all$cell_type
umap_df$Status <- merged_all$Status

# Define specific colors
rte_colors <- c(
  "CD4 Naive" = "#1f77b4",
  "CD8 Naive" = "#ff7f0e",
  "CD4 RTE"   = "#e41a1c",
  "CD8 RTE"   = "#33a02c"
)

# Highlight only these 4 types
target_types <- names(rte_colors)

# Plot loop
plot_list <- lapply(c("Healthy", "RA", "PSA", "SPA"), function(st) {
  df_sub <- umap_df %>% filter(Status == st)
  
  ggplot() +
    # Background (All other cells in grey)
    geom_point(data = df_sub %>% filter(!cell_type %in% target_types),
               aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 0.3, alpha = 0.2) +
    # Layer 1: Naive (Background context)
    geom_point(data = df_sub %>% filter(cell_type %in% c("CD4 Naive", "CD8 Naive")),
               aes(x = UMAP_1, y = UMAP_2, color = cell_type), size = 0.4, alpha = 0.3) +
    # Layer 2: RTE (Foreground highlight)
    geom_point(data = df_sub %>% filter(cell_type %in% c("CD4 RTE", "CD8 RTE")),
               aes(x = UMAP_1, y = UMAP_2, color = cell_type), size = 0.7, alpha = 0.9) +
    scale_color_manual(values = rte_colors) +
    ggtitle(st) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face="bold")) +
    NoLegend()
})
pdf(file.path(PLOT_DIR, "Fig3C_UMAP_RTE_Overlay.pdf"), width = 16, height = 4)
wrap_plots(plot_list, ncol = 4)
dev.off()
# ==============================================================================
# 7. Visualization 2: Naive/RTE Composition Barplot (100% Stacked)
# ==============================================================================
message("Generating 100% Stacked Proportion Plot per Individual...")
# We calculate the proportion of the 4 subsets WITHIN the Naive/RTE compartment
# Total = CD4 Naive + CD4 RTE + CD8 Naive + CD8 RTE
# This ensures they sum to 100% per individual
prop_df <- meta %>%
  filter(cell_type %in% target_types) %>% # Only keep the 4 types
  group_by(Subject, Status, cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Subject) %>%
  mutate(fraction = count / sum(count)) # Normalize to 100%

# Order Subjects by Status for clean plotting
subject_order <- prop_df %>% 
  select(Subject, Status) %>% 
  distinct() %>% 
  arrange(factor(Status, levels = c("Healthy", "RA", "PSA", "SPA"))) %>% 
  pull(Subject)

prop_df$Subject <- factor(prop_df$Subject, levels = subject_order)

p_stack <- ggplot(prop_df, aes(x = Subject, y = fraction, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill", width = 0.9) +
  scale_fill_manual(values = rte_colors) +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(~Status, scales = "free_x", space = "free") +
  theme_minimal(base_size = 14) +
  labs(title = "Naive T Cell Compartment Composition", x = "Subject", y = "Proportion") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.title = element_blank()
  )
pdf(file.path(PLOT_DIR, "Fig3D_RTE_Proportions_Stacked.pdf"), width = 12, height = 6)
print(p_stack)
dev.off()
# ==============================================================================
# 8. Save Final Object
# ==============================================================================
saveRDS(merged_all, file.path(OUTPUT_DIR, "merged_all_final_with_RTE.rds"))
message("Step 05 Complete! RTE analysis finished.")