# scripts/Part2_CrossDisease/04_annotation_plots.R
library(Seurat)
library(ggplot2)
library(dplyr)
library(plyr)
library(patchwork)
library(RColorBrewer)
BASE_DIR <- "/path/to/ADT4DISEASES" # UPDATE THIS
OUTPUT_DIR <- file.path(BASE_DIR, "results/part2_cross_disease")
PLOT_DIR <- file.path(OUTPUT_DIR, "MarkerPlots1011")
if(!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)
# Load Integrated Data
merged_all <- readRDS(file.path(OUTPUT_DIR, "merged_all_harmony.rds"))
# ==============================================================================
# 1. Marker Visualization
# ==============================================================================
marker_genes <- list(
  "Naive B" = c("MS4A1", "IGHD", "CD19"),
  "Memory B" = c("MS4A1", "CD27"),
  "Plasma B" = c("MZB1", "XBP1", "PRDM1"),
  "CD4 T" = c("CD3E", "CD4", "IL7R"),
  "CD8 T" = c("CD8A", "CD8B", "CD3E"),
  "T naive" = c("CCR7", "SELL", "LEF1", "TCF7"),
  "T central memory" = c("CCR7", "IL7R", "CD27", "CD28"),
  "T effector memory" = c("GZMK", "GZMB", "GNLY", "NKG7"),
  "T exhausted" = c("PDCD1", "LAG3", "HAVCR2"),
  "Regulatory T" = c("FOXP3", "IL2RA", "CTLA4"),
  "Classical Mono" = c("CD14", "LYZ", "S100A8", "S100A9"),
  "Non-classical Mono" = c("FCGR3A", "MS4A7", "CX3CR1"),
  "Intermediate Mono" = c("CD14", "FCGR3A", "HLA-DRA", "CD74"),
  "pDC" = c("IRF8", "CLEC4C", "TCF4", "LILRA4"),
  "NK" = c("NKG7", "GNLY", "GZMB", "KLRD1", "KLRF1"),
  "Platelets" = c("PPBP", "PF4", "ITGA2B", "GP1BA"),
  "Intermediate_Unknown" = c("TYROBP", "CD7", "IL7R", "FCGR3A", "TRDC", "ISG15", "MX1", "PTPRC")
)
# FeaturePlot
message("Generating FeaturePlots...")
for (group in names(marker_genes)) {
  safe_group_name <- gsub("[^A-Za-z0-9]", "_", group)
  pdf(file = file.path(PLOT_DIR, paste0("Feature_", safe_group_name, ".pdf")), width = 8, height = 6)
  print(
    FeaturePlot(merged_all, features = marker_genes[[group]], 
                cols = c("grey90", "red"), max.cutoff = 3, ncol = 2) + NoLegend()
  )
  dev.off()
}
# VlnPlot
message("Generating VlnPlots...")
# Note: Using default colors here as cell_type isn't assigned yet in the logic flow
for (group in names(marker_genes)) {
  safe_group_name <- gsub("[^A-Za-z0-9]", "_", group)
  pdf(file = file.path(PLOT_DIR, paste0("Vln_", safe_group_name, ".pdf")), width = 8, height = 6)
  print(
    VlnPlot(merged_all, features = marker_genes[[group]], pt.size = 0, ncol = 2) + NoLegend()
  )
  dev.off()
}

# UMAP with Clusters
pdf(file.path(PLOT_DIR, "UMAP_cluster_labeled.pdf"), width = 7, height = 5)
print(DimPlot(merged_all, reduction = "umap", label = TRUE, label.size = 5) + NoLegend())
dev.off()

# DotPlot All Markers
all_markers <- unique(unlist(marker_genes))
pdf(file.path(PLOT_DIR, "DotPlot_AllMarkers.pdf"), width = 12, height = 8)
print(
  DotPlot(object = merged_all, features = all_markers, cols = c("white", "red"), dot.scale = 6) +
    RotatedAxis() +
    scale_colour_gradientn(colours = c("white", "red"), limits = c(0, 2), na.value = "white") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(face = "bold"))
)
dev.off()
# ==============================================================================
# 2. Annotation
# ==============================================================================
message("Applying Cell Type Annotations...")
cluster_to_celltype <- c(
  "0" = "CD4 Naive", "11" = "CD4 Naive",
  "9" = "CD8 Naive",
  "12" = "CD8 Tcm", "29" = "CD8 Tcm",
  "3" = "CD4 Tcm", "15" = "CD4 Tcm",
  "6" = "CD4 Tem",
  "16" = "CD8 Tem", "4" = "CD8 Tem",
  "24" = "NKT Cells",
  "17" = "NK Cells (CD56high)",
  "2" = "NK Cells (CD56low)", "27" = "NK Cells (CD56low)",
  "8" = "B Naive", "18" = "B Naive",
  "10" = "B Memory",
  "20" = "pDC",
  "22" = "Plasma B",
  "19" = "Platelets",
  "23" = "Hybrid",
  "14" = "cDC",
  "1" = "Classical Mono", "5" = "Classical Mono", "13" = "Classical Mono",
  "21" = "Classical Mono", "28" = "Classical Mono", "25" = "Classical Mono", "26" = "Classical Mono",
  "7" = "Non-classical Mono"
)
merged_all$cell_type <- plyr::mapvalues(
  as.character(merged_all$seurat_clusters),
  from = names(cluster_to_celltype),
  to = unname(cluster_to_celltype)
)
# Define Colors
celltype_colors <- c(
  "CD4 Naive" = "#1f77b4", "CD8 Naive" = "#ff7f0e", "CD8 Tcm" = "#2ca02c",
  "CD4 Tcm" = "#d62728", "CD4 Tem" = "#9467bd", "CD8 Tem" = "#8c564b",
  "NKT Cells" = "#e377c2", "NK Cells (CD56high)" = "#7f7f7f", "NK Cells (CD56low)" = "#bcbd22",
  "B Naive" = "#17becf", "B Memory" = "#aec7e8", "pDC" = "#ff9896",
  "Plasma B" = "#98df8a", "Platelets" = "#c5b0d5", "Hybrid" = "#ffbb78",
  "cDC" = "#f7b6d2", "Classical Mono" = "#c49c94", "Non-classical Mono" = "#dbdb8d"
)

# ==============================================================================
# 3. Final UMAPs
# ==============================================================================
message("Generating Annotated UMAPs...")
# Overall UMAP
pdf(file.path(PLOT_DIR, "UMAP_CellType_All.pdf"), width = 10, height = 7)
print(
  DimPlot(merged_all, reduction = "umap", group.by = "cell_type", label = TRUE, cols = celltype_colors) +
    ggtitle("Cell Type Distribution") + theme(plot.title = element_text(hjust = 0.5))
)
dev.off()
# Split UMAP by Status
Idents(merged_all) <- "cell_type"
ordered_status <- c("Healthy", "RA", "PSA", "SPA")
existing_status <- ordered_status[ordered_status %in% unique(merged_all$Status)]

umap_list <- lapply(existing_status, function(status_name) {
  obj_subset <- subset(merged_all, subset = Status == status_name)
  DimPlot(obj_subset, reduction = "umap", group.by = "cell_type", cols = celltype_colors, 
          label = TRUE, label.size = 3) +
    ggtitle(status_name) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 14), plot.margin = margin(5, 5, 5, 5)) +
    NoLegend()
})
pdf(file.path(PLOT_DIR, "UMAP_CellType_By_Status.pdf"), width = 14, height = 4)
print(wrap_plots(umap_list, ncol = length(umap_list)))
dev.off()
# Save Annotated Object
saveRDS(merged_all, file.path(OUTPUT_DIR, "merged_all_annotated.rds"))
message("Step 04 Complete! All plots saved to: MarkerPlots1011/")