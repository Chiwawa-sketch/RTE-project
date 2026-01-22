# scripts/01_preprocessing.R
# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
# ==============================================================================
# 1. Configuration
# ==============================================================================
# UPDATE THIS PATH: Point this to the folder containing your GSM folders
DATA_DIR <- "/mnt/data/home/tyu-a1ouoisq/GSE157007" 
OUTPUT_DIR <- "results"
# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
# ==============================================================================
# 2. Sample Metadata Definition (Complete List)
# ==============================================================================
# Explicitly listing all 17 samples to ensure reproducibility
sample_info <- list(
  "GSM4750306" = "Cord_blood",
  "GSM4750307" = "Cord_blood",
  "GSM4750308" = "Cord_blood",
  "GSM4750303" = "Young",
  "GSM4750304" = "Young",
  "GSM4750305" = "Young",
  "GSM4750309" = "Old",
  "GSM4750310" = "Old",
  "GSM4750311" = "Old",
  "GSM5684306" = "Old",
  "GSM5684307" = "Old",
  "GSM5684308" = "Old",
  "GSM4750298" = "Frail",
  "GSM4750299" = "Frail",
  "GSM4750300" = "Frail",
  "GSM4750301" = "Frail",
  "GSM4750302" = "Frail"
)
# ==============================================================================
# 3. Data Loading & Object Creation
# ==============================================================================
obj_list <- list()
message("Starting data loading...")
for (sample_id in names(sample_info)) {
  # Construct full path
  sample_path <- file.path(DATA_DIR, sample_id)
  
  # Check if path exists
  if (!dir.exists(sample_path)) {
    warning(paste("Directory not found:", sample_path, "- Skipping."))
    next
  }
  # Read 10X data
  data <- Read10X(data.dir = sample_path)
  # Handle multi-assay objects (extract 'Gene Expression' if present)
  if (is.list(data)) {
    if ("Gene Expression" %in% names(data)) {
      data <- data[["Gene Expression"]]
    } else {
      # Fallback: take the first list element
      data <- data[[1]]
    }
  }
  # Create Seurat Object
  seu <- CreateSeuratObject(
    counts = data, 
    project = sample_id, 
    min.cells = 3, 
    min.features = 200
  )
  # Add Metadata: Sample ID and Group
  seu$sample_id <- sample_id
  seu$Group <- sample_info[[sample_id]]
  # Standard Pre-processing (Normalization & Feature Selection)
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj_list[[sample_id]] <- seu
  message(paste("Loaded:", sample_id, "| Group:", sample_info[[sample_id]]))
}
# Stop if no data was loaded
if (length(obj_list) == 0) {
  stop("No data loaded! Please check your DATA_DIR path.")
}
# ==============================================================================
# 4. Integration (CCA)
# ==============================================================================
message("Starting data integration (this may take a while)...")
# Find Integration Anchors
anchors <- FindIntegrationAnchors(object.list = obj_list, dims = 1:30, verbose = TRUE)
# Integrate Data
integrated <- IntegrateData(anchorset = anchors, dims = 1:30, verbose = TRUE)
# ==============================================================================
# 5. Preprocessing & Global Clustering
# ==============================================================================
message("Running dimension reduction and clustering...")
# Switch to integrated assay
DefaultAssay(integrated) <- "integrated"
# Calculate Mitochondrial percentage for regression
integrated[["percent.mt"]] <- PercentageFeatureSet(integrated, pattern = "^MT-")
# Scale, PCA, UMAP, and Clustering
integrated <- ScaleData(integrated, vars.to.regress = c("percent.mt", "nCount_RNA", "nFeature_RNA"), verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
integrated <- FindClusters(integrated, resolution = 0.4, verbose = FALSE)
# ==============================================================================
# 6. Visualization & Saving
# ==============================================================================
message("Saving results...")
# Fig 1a: Global UMAP
p1 <- DimPlot(integrated, group.by = "seurat_clusters", label = TRUE, label.size = 5) + 
  ggtitle("Global Clustering") + NoLegend()
p2 <- DimPlot(integrated, group.by = "Group", shuffle = TRUE) + 
  ggtitle("Group Distribution")
pdf(file.path(OUTPUT_DIR, "supFig1a_UMAP_Global.pdf"), width = 12, height = 6)
print(p1 + p2)
dev.off()
# Fig 1b: Canonical Marker DotPlot (Lineage Identification)
markers <- c("CD3E", "CD4", "CD8A",   # T cells
             "MS4A1", "CD79A",        # B cells
             "CD14", "LYZ", "FCGR3A", # Monocytes
             "GNLY", "NKG7",          # NK cells
             "PPBP")                  # Platelets

pdf(file.path(OUTPUT_DIR, "supFig1b_DotPlot_Global.pdf"), width = 10, height = 5)
print(
  DotPlot(integrated, features = markers, group.by = "seurat_clusters") + 
    RotatedAxis() + 
    scale_color_gradientn(colours = c("white", "red"))
)
dev.off()
# Save the full integrated object for the next steps
saveRDS(integrated, file = file.path(OUTPUT_DIR, "integrated_obj.rds"))
message("Step 01 Complete! Integrated object saved to results/integrated_obj.rds")