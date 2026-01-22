library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(Matrix)
library(plyr)
library(cowplot)
library(RColorBrewer)
library(pheatmap)
library(ggsignif)
library(ggpubr)
library(NMF)
library(ggalluvial)
library(CellChat)
library(qs)
library(future)
library(irGSEA)
library(fmsb)
library(sceasy)
library(reticulate)

# ==============================================================================
# 1. Configuration & Paths
# ==============================================================================
# Server Path
BASE_DIR <- "/mnt/data/home/tyu-5shmhzlu/Figure3T2DGSE255566"
setwd(BASE_DIR)

# Output directory
PLOT_DIR <- file.path(BASE_DIR, "T2D_Plots")
if(!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, showWarnings = FALSE)

# Custom Colors
cell_type_cols <- c(
  brewer.pal(9, "Set1"),
  "#FF34B3", "#BC8F8F", "#20B2AA", "#00F5FF", "#FFA500", "#ADFF2F", "#FF6A6A", "#7FFFD4",
  "#AB82FF", "#90EE90", "#00CD00", "#008B8B", "#6495ED", "#FFC1C1", "#CD5C5C", "#8B008B", "#FF3030",
  "#7CFC00", "#000000", "#708090", "#FF7F50", "#6A5ACD", "#3CB371", "#B22222", "#FF1493"
)

# ==============================================================================
# 2. Data Loading & Integration
# ==============================================================================
# Con = HC, Mod = T2D
sample_names <- c(
  "GSM8075321Con", "GSM8075322Con", "GSM8075323Con",
  "GSM8075324Mod", "GSM8075325Mod", "GSM8075326Mod"
)
sample_paths <- setNames(file.path(BASE_DIR, sample_names), sample_names)

seurat_list <- list()
message("Loading T2D/HC samples...")
for (sample in sample_names) {
  if(dir.exists(sample_paths[[sample]])) {
    data <- Read10X(data.dir = sample_paths[[sample]])
    seu <- CreateSeuratObject(counts = data, project = sample, min.cells = 3, min.features = 200)
    seu$sample <- sample
    # Define Group: Con -> HC, Mod -> T2D
    seu$group <- ifelse(grepl("Con", sample, ignore.case = TRUE), "HC", "T2D")
    seu$disease <- seu$group # Alias for downstream compatibility
    seu <- NormalizeData(seu, verbose = FALSE)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    seurat_list[[sample]] <- seu
  } else {
    warning(paste("Path not found:", sample_paths[[sample]]))
  }
}

message("Integrating data (CCA)...")
T2D_anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30, anchor.features = 3500)
T2D <- IntegrateData(anchorset = T2D_anchors, dims = 1:30)
DefaultAssay(T2D) <- "integrated"

# ==============================================================================
# 3. Dimensionality Reduction & Clustering
# ==============================================================================
message("Running PCA, UMAP and Clustering...")
T2D <- ScaleData(T2D, vars.to.regress = c("percent.mt", "nCount_RNA", "nFeature_RNA"), verbose = FALSE)
T2D <- RunPCA(T2D, npcs = 50, verbose = FALSE)
T2D <- FindNeighbors(T2D, reduction = "pca", dims = 1:50)
T2D <- FindClusters(T2D, resolution = 0.8) 
T2D <- RunUMAP(T2D, reduction = "pca", dims = 1:50)

# ==============================================================================
# 4. Marker Definitions & Automated Annotation
# ==============================================================================
message("Running Automated Cluster Annotation...")
DefaultAssay(T2D) <- "RNA"
T2D <- NormalizeData(T2D, verbose = FALSE)

# Detailed Marker List (Same as COPD pipeline)
marker_genes <- list(
  "Naive B" = c("MS4A1", "IGHD", "CD19"),
  "Memory B" = c("MS4A1", "CD27"),
  "Plasma B" = c("MZB1", "XBP1", "PRDM1"),
  "CD4 Naive" = c("CD4", "CCR7", "SELL", "LEF1", "TCF7"), # Adjusted key to match exact output needed
  "CD4 Tcm" = c("CD4", "IL7R", "CD28"),
  "CD8 Naive" = c("CD8A", "CD8B", "CCR7", "SELL"), # Adjusted key
  "CD8 Tem" = c("CD8A", "GZMK", "GZMB", "NKG7"),
  "DNT Cells" = c("CD3E", "TRAC", "TRBC1"),
  "Regulatory T" = c("FOXP3", "IL2RA", "CTLA4"),
  "MAIT" = c("SLC4A10", "KLRB1", "RORC"),
  "Classical Monocyte" = c("CD14", "LYZ", "S100A8", "S100A9"),
  "Non-classical Monocyte" = c("FCGR3A", "MS4A7", "CX3CR1"),
  "pDC" = c("IRF8", "CLEC4C", "LILRA4"),
  "NK Cells" = c("NKG7", "GNLY", "KLRD1"),
  "Platelets" = c("PPBP", "PF4")
)

# --- Automated Annotation Logic ---
# 1. Calculate Module Score for each cell type signature
T2D <- AddModuleScore(T2D, features = marker_genes, name = "Sig_")

# 2. Rename score metadata to match cell types
sig_cols <- paste0("Sig_", 1:length(marker_genes))
names(sig_cols) <- names(marker_genes)

# 3. Assign each cluster to the signature with the highest average score
cluster_annotations <- character()
cluster_ids <- levels(Idents(T2D))

for (id in cluster_ids) {
  cells_in_cluster <- WhichCells(T2D, idents = id)
  # Calculate mean score for each signature in this cluster
  mean_scores <- colMeans(T2D@meta.data[cells_in_cluster, sig_cols])
  # Pick max
  best_match <- names(sig_cols)[which.max(mean_scores)]
  cluster_annotations[id] <- best_match
}

message("Automated Cluster Assignments:")
print(cluster_annotations)

# 4. Apply Annotation
T2D <- RenameIdents(T2D, cluster_annotations)
T2D$celltype <- Idents(T2D)

# 5. Save Annotated UMAP
pdf(file.path(PLOT_DIR, "UMAP_Annotated_Auto.pdf"), width = 10, height = 7)
DimPlot(T2D, group.by = "celltype", label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
dev.off()

# Save Annotated Object
saveRDS(T2D, "T2D_Annotated.rds")

# ==============================================================================
# 5. RTE Definition (T2D Specific Thresholds)
# ==============================================================================
message("Defining RTE based on SOX4 expression...")
# Extract SOX4 data
expr_gene <- as.numeric(LayerData(T2D[["RNA"]], layer = "data")["SOX4", ])
T2D$expr_SOX4 <- expr_gene

# Apply User-Defined Thresholds for T2D
# CD4 > 1.09
# CD8 > 0.94
T2D$SOX4_status <- NA

# Logic: Find cells annotated as Naive and check SOX4 levels
# Note: Auto-annotation ensures we have "CD4 Naive" and "CD8 Naive" labels
T2D$SOX4_status[T2D$celltype == "CD4 Naive" & expr_gene > 1.09] <- "SOX4+ CD4 Naive"
T2D$SOX4_status[T2D$celltype == "CD4 Naive" & expr_gene <= 1.09] <- "SOX4- CD4 Naive"

T2D$SOX4_status[T2D$celltype == "CD8 Naive" & expr_gene > 0.94] <- "SOX4+ CD8 Naive"
T2D$SOX4_status[T2D$celltype == "CD8 Naive" & expr_gene <= 0.94] <- "SOX4- CD8 Naive"

# Update Cell Types with RTE
T2D$celltype <- as.character(T2D$celltype)
T2D$celltype[T2D$SOX4_status == "SOX4+ CD4 Naive"] <- "CD4 RTE"
T2D$celltype[T2D$SOX4_status == "SOX4- CD4 Naive"] <- "CD4 T Naive"
T2D$celltype[T2D$SOX4_status == "SOX4+ CD8 Naive"] <- "CD8 RTE"
T2D$celltype[T2D$SOX4_status == "SOX4- CD8 Naive"] <- "CD8 T Naive"

# Export Cell Counts
write.csv(table(T2D$sample, T2D$celltype), file.path(BASE_DIR, "Cell_Counts_Sample_RTE_T2D.csv"))

# ==============================================================================
# 6. CellChat Analysis
# ==============================================================================
message("Running CellChat Analysis...")

run_cellchat_subset <- function(seurat_obj, group_name) {
  data_input <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  meta <- data.frame(celltype = seurat_obj$celltype, row.names = colnames(seurat_obj))
  cc <- createCellChat(object = data_input, meta = meta, group.by = "celltype")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, raw.use = FALSE)
  cc <- filterCommunication(cc, min.cells = 5)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  return(cc)
}

# Run for HC and T2D
cc_T2D <- run_cellchat_subset(subset(T2D, group == "T2D"), "T2D")
cc_HC <- run_cellchat_subset(subset(T2D, group == "HC"), "HC")

# Save CellChat Objects
saveRDS(cc_T2D, file.path(BASE_DIR, "cellchat_T2D.rds"))
saveRDS(cc_HC, file.path(BASE_DIR, "cellchat_HC_T2Dcontrol.rds"))

# Comparison Plot
pdf(file.path(PLOT_DIR, "CellChat_Circle_Comparison_T2D.pdf"), width = 10, height = 5)
par(mfrow=c(1,2))
netVisual_circle(cc_T2D@net$weight, title.name = "T2D Strength")
netVisual_circle(cc_HC@net$weight, title.name = "HC Strength")
dev.off()

# ==============================================================================
# 7. irGSEA Analysis (Radar Plots)
# ==============================================================================
message("Running irGSEA Scoring...")
markers_func <- list(
  Quiescence = c("BTG2", "KLF2", "FOXP1"),
  Regulating = c("TNFRSF9", "TNFRSF18", "ENTPD1", "FOXP3"),
  Proliferation = c("MKI67", "TYMS"),
  Helper = c("TNF", "IL2", "IFNG"),
  Cytotoxicity = c("GZMB", "PRF1", "GNLY"),
  Exhaustion = c("PDCD1", "HAVCR2", "TOX"),
  Senescence = c("KLRG1", "CDKN2A")
)

T2D <- irGSEA.score(object = T2D, assay = "RNA", slot = "data", geneset = markers_func, method = "AUCell", kcdf=FALSE)

# Extract scores
aucell_data <- as.data.frame(t(GetAssayData(T2D, assay = "AUCell", slot = "data")))
T2D <- AddMetaData(T2D, aucell_data)
state_features <- names(markers_func)

plot_radar_compare <- function(cell_type) {
  if(!cell_type %in% T2D$celltype) return(NULL)
  sub_obj <- subset(T2D, celltype == cell_type)
  df_list <- split(sub_obj@meta.data[, state_features], sub_obj$group)
  
  if(length(df_list) < 2) return(NULL)
  
  means <- do.call(rbind, lapply(df_list, colMeans, na.rm=TRUE))
  data_radar <- rbind(max = rep(max(means)*1.2, length(state_features)), min = rep(0, length(state_features)), means)
  radarchart(as.data.frame(data_radar), pcol = c("red", "blue"), plwd = 2, title = paste(cell_type, "State"))
  legend("topright", legend = rownames(means), col = c("red", "blue"), lty = 1)
}

pdf(file.path(PLOT_DIR, "RTE_Radar_Plots_T2D.pdf"), width = 10, height = 5)
par(mfrow=c(1,2))
plot_radar_compare("CD4 RTE")
plot_radar_compare("CD8 RTE")
dev.off()

# ==============================================================================
# 8. Save and Convert to h5ad
# ==============================================================================
message("Saving final object and converting to h5ad...")
saveRDS(T2D, "T2D_final.rds")

tryCatch({
  sceasy::convertFormat(T2D, from="seurat", to="anndata", outFile="T2D_scored.h5ad")
  message("Successfully created T2D_scored.h5ad")
}, error = function(e) {
  message("Error converting to h5ad. Please check environment.")
})
