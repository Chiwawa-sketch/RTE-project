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
library(sceasy) # For h5ad conversion
library(reticulate)

# ==============================================================================
# 1. Configuration & Paths
# ==============================================================================
# User specified path
BASE_DIR <- "/mnt/data/home/tyu-5shmhzlu/Figure2COPDGSE249584"
setwd(BASE_DIR)

# Output directory for plots
PLOT_DIR <- file.path(BASE_DIR, "COPD_Plots")
if(!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, showWarnings = FALSE)

# Custom Colors
cell_type_cols <- c(
  brewer.pal(9, "Set1"),
  "#FF34B3", "#BC8F8F", "#20B2AA", "#00F5FF", "#FFA500", "#ADFF2F", "#FF6A6A", "#7FFFD4",
  "#AB82FF", "#90EE90", "#00CD00", "#008B8B", "#6495ED", "#FFC1C1", "#CD5C5C", "#8B008B", "#FF3030",
  "#7CFC00", "#000000", "#708090", "#FF7F50", "#6A5ACD", "#3CB371", "#B22222", "#FF1493",
  "#20B2AA", "#BDB76B", "#9932CC", "#FFD700", "#00FA9A", "#8B0000", "#FF00FF", "#40E0D0", "#1E90FF",
  "#A0522D", "#D2691E", "#C71585", "#5F9EA0", "#FF4500", "#DA70D6", "#9ACD32", "#8FBC8F"
)

# ==============================================================================
# 2. Data Loading & Integration
# ==============================================================================
sample_names <- c(
  "CONTROLGSM7950680", "CONTROLGSM7950681", "CONTROLGSM7950682", "CONTROLGSM7950683",
  "CONTROLGSM7950684", "CONTROLGSM7950687", "CONTROLGSM7950691",
  "COPDGSM7950677", "COPDGSM7950678", "COPDGSM7950679",
  "COPDGSM7950685", "COPDGSM7950686", "COPDGSM7950688",
  "COPDGSM7950689", "COPDGSM7950690"
)
sample_paths <- setNames(file.path(BASE_DIR, sample_names), sample_names)

seurat_list <- list()
message("Loading samples...")
for (sample in sample_names) {
  if(dir.exists(sample_paths[[sample]])) {
    data <- Read10X(data.dir = sample_paths[[sample]])
    seu <- CreateSeuratObject(counts = data, project = sample, min.cells = 3, min.features = 200)
    seu$sample <- sample
    seu$group <- ifelse(grepl("CONTROL", sample, ignore.case = TRUE), "Control", "COPD")
    seu <- NormalizeData(seu, verbose = FALSE)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    seurat_list[[sample]] <- seu
  } else {
    warning(paste("Path not found:", sample_paths[[sample]]))
  }
}

message("Integrating data (CCA)...")
COPD_anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30, anchor.features = 3500)
COPD <- IntegrateData(anchorset = COPD_anchors, dims = 1:30)
DefaultAssay(COPD) <- "integrated"

# ==============================================================================
# 3. Dimensionality Reduction & Clustering
# ==============================================================================
message("Running PCA, UMAP and Clustering...")
COPD <- ScaleData(COPD, vars.to.regress = c("percent.mt", "nCount_RNA", "nFeature_RNA"), verbose = FALSE)
COPD <- RunPCA(COPD, npcs = 50, verbose = FALSE)
COPD <- FindNeighbors(COPD, reduction = "pca", dims = 1:50)
COPD <- FindClusters(COPD, resolution = 0.9)
COPD <- RunUMAP(COPD, reduction = "pca", dims = 1:50)
COPD <- RunTSNE(COPD, reduction = "pca", dims = 1:50)

# Add Comparison Metadata
current.cluster.ids <- sample_names
new.cluster.ids <- ifelse(grepl("CONTROL", sample_names, ignore.case = TRUE), "Control", "COPD")
COPD$Comparison <- plyr::mapvalues(COPD$sample, from = current.cluster.ids, to = new.cluster.ids)

# Plot UMAP by Cluster (Unknown Clusters)
pdf(file.path(PLOT_DIR, "UMAP_Clusters_Initial.pdf"), width = 8, height = 6)
DimPlot(COPD, group.by = "seurat_clusters", label = TRUE, label.size = 5) + NoLegend()
dev.off()

# ==============================================================================
# 4. Marker Visualization (BEFORE Annotation)
# ==============================================================================
message("Generating Marker Plots for Manual Annotation...")
DefaultAssay(COPD) <- "RNA"
COPD <- NormalizeData(COPD, verbose = FALSE)
COPD <- ScaleData(COPD, verbose = FALSE)

marker_genes <- list(
  "Naive B" = c("MS4A1", "IGHD", "CD19"),
  "Memory B" = c("MS4A1", "CD27"),
  "Plasma B" = c("MZB1", "XBP1", "PRDM1"),
  "CD4 T" = c("CD3E", "CD4", "IL7R", "CCR7", "TCF7"),
  "CD8 T" = c("CD8A", "CD8B", "CD3E", "GZMK", "NKG7"),
  "DNT Cells" = c("CD3E", "CD3D", "CD3G", "TRAC", "TRBC1", "TRBC2"),
  "T naive" = c("CCR7", "SELL", "LEF1", "TCF7", "IL7R"),
  "T central memory" = c("CCR7", "IL7R", "CD27", "CD28", "LEF1", "TCF7"),
  "T effector memory" = c("GZMK", "GZMB", "GNLY", "NKG7", "PRF1"),
  "T exhausted" = c("PDCD1", "LAG3", "HAVCR2", "TOX", "TIGIT"),
  "Regulatory T" = c("FOXP3", "IL2RA", "CTLA4"),  # <--- Added here
  "GammaDelta T" = c("TRDC", "TRGC1", "TRGC2", "KLRB1", "S1PR1"),
  "MAIT" = c("SLC4A10", "TRAV1-2", "KLRB1", "RORC", "IL18RAP"),
  "Classical Mono" = c("CD14", "LYZ", "S100A8", "S100A9"),
  "Non-classical Mono" = c("FCGR3A", "MS4A7", "CX3CR1"),
  "Intermediate Mono" = c("CD14", "FCGR3A", "HLA-DRA", "CD74"),
  "pDC" = c("IRF8", "CLEC4C", "TCF4", "LILRA4"),
  "NK" = c("NKG7", "GNLY", "GZMB", "KLRD1", "KLRF1", "NCAM1", "PRF1"),
  "NKT" = c("ZBTB16", "TRAV10", "TRBV25-1", "CD3E", "KLRK1", "TYROBP"),
  "Platelets" = c("PPBP", "PF4", "ITGA2B", "GP1BA"),
  "Intermediate/Unknown" = c("TYROBP", "CD7", "IL7R", "FCGR3A", "TRDC", "ISG15", "MX1", "PTPRC")
)

# FeaturePlots
dir.create(file.path(PLOT_DIR, "MarkerPlots"), showWarnings = FALSE)
for (group in names(marker_genes)) {
  safe_name <- gsub("[^A-Za-z0-9]", "_", group)
  pdf(file.path(PLOT_DIR, "MarkerPlots", paste0("Feature_", safe_name, ".pdf")), width = 8, height = 6)
  print(FeaturePlot(COPD, features = marker_genes[[group]], cols = c("grey90", "red"), ncol = 2) + NoLegend())
  dev.off()
}

# ViolinPlots
for (group in names(marker_genes)) {
  safe_name <- gsub("[^A-Za-z0-9]", "_", group)
  pdf(file.path(PLOT_DIR, "MarkerPlots", paste0("Vln_", safe_name, ".pdf")), width = 10, height = 16)
  print(VlnPlot(COPD, features = marker_genes[[group]], pt.size = 0, ncol = 1) + NoLegend())
  dev.off()
}

# ==============================================================================
# 5. Annotation (Apply Defined Clusters)
# ==============================================================================
message("Applying Cluster Annotations...")
cluster2celltype <- c(
  "0" = "Classical Monocyte", "3" = "Classical Monocyte", "19" = "Classical Monocyte",
  "26" = "Classical Monocyte", "23" = "Classical Monocyte", "18" = "Classical Monocyte", "41" = "Classical Monocyte",
  "11" = "Non-classical Monocyte",
  "15" = "NK Cells (CD56 low)", "12" = "NK Cells (CD56 low)", "2" = "NK Cells (CD56 low)",
  "37" = "NK Cells (CD56 low)", "5" = "NK Cells (CD56 low)", "8" = "NK Cells (CD56 low)",
  "31" = "NK Cells (CD56 high)",
  "33" = "pDC", "28" = "MAIT", "24" = "Platelet",
  "1" = "DNT Cells", "7" = "DNT Cells",
  "10" = "CD4 Naive", "20" = "CD4 Naive", "25" = "CD4 Tcm",
  "4" = "CD8 Naive", "21" = "CD8 Tem", "6" = "CD8 Tem", "13" = "CD8 Tem",
  "30" = "CD8 Tem", "36" = "CD8 Tem", "39" = "CD8 Tem", "29" = "CD8 Tem",
  "27" = "CD8 Tem", "32" = "CD8 Tem", "22" = "CD8 Tcm",
  "16" = "NKT Cells", "35" = "NKT Cells",
  "14" = "Naive B Cells", "9" = "Naive B Cells", "17" = "Memory B Cells", "34" = "Plasma B Cells",
  "38" = "Hybrid", "40" = "Hybrid"
)

Idents(COPD) <- "seurat_clusters"
COPD$celltype <- plyr::mapvalues(as.character(Idents(COPD)), from = names(cluster2celltype), to = as.character(cluster2celltype))

pdf(file.path(PLOT_DIR, "UMAP_Annotated.pdf"), width = 8, height = 6)
DimPlot(COPD, group.by = "celltype", label = TRUE, label.size = 4) + NoLegend()
dev.off()

# Save Annotated Object
saveRDS(COPD, "COPD_Annotated.rds")

# ==============================================================================
# 6. RTE Definition (SOX4 Peak Detection)
# ==============================================================================
message("Defining RTE based on SOX4 expression...")
expr_gene <- as.numeric(LayerData(COPD[["RNA"]], layer = "data")["SOX4", ])
COPD$expr_SOX4 <- expr_gene

# Visualize Density
plot_data <- COPD@meta.data %>% filter(celltype %in% c("CD8 Naive", "CD4 Naive"), expr_SOX4 > 0)
pdf(file.path(PLOT_DIR, "SOX4_Density_Peaks.pdf"), width = 7, height = 5)
print(ggplot(plot_data, aes(x = expr_SOX4)) + geom_density(fill="steelblue", alpha=0.4) + facet_wrap(~celltype) + theme_minimal())
dev.off()

# Apply Thresholds (Hardcoded from your analysis)
COPD$SOX4_status <- NA
COPD$SOX4_status[COPD$celltype == "CD4 Naive" & expr_gene > 1.25] <- "SOX4+ CD4 Naive"
COPD$SOX4_status[COPD$celltype == "CD4 Naive" & expr_gene <= 1.25] <- "SOX4- CD4 Naive"
COPD$SOX4_status[COPD$celltype == "CD8 Naive" & expr_gene > 1.18] <- "SOX4+ CD8 Naive"
COPD$SOX4_status[COPD$celltype == "CD8 Naive" & expr_gene <= 1.18] <- "SOX4- CD8 Naive"

# Update Cell Types
COPD$celltype <- as.character(COPD$celltype)
COPD$celltype[COPD$SOX4_status == "SOX4+ CD4 Naive"] <- "CD4 RTE"
COPD$celltype[COPD$SOX4_status == "SOX4- CD4 Naive"] <- "CD4 T Naive"
COPD$celltype[COPD$SOX4_status == "SOX4+ CD8 Naive"] <- "CD8 RTE"
COPD$celltype[COPD$SOX4_status == "SOX4- CD8 Naive"] <- "CD8 T Naive"

# Add disease grouping metadata (HC vs COPD)
COPD$disease <- ifelse(grepl("CONTROL", COPD$sample, ignore.case = TRUE), "HC", "COPD")

# ==============================================================================
# 7. CellChat Analysis
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

# Run for both groups
cc_COPD <- run_cellchat_subset(subset(COPD, disease == "COPD"), "COPD")
cc_HC <- run_cellchat_subset(subset(COPD, disease == "HC"), "HC")

# Save CellChat Objects
saveRDS(cc_COPD, file.path(BASE_DIR, "cellchat_COPD.rds"))
saveRDS(cc_HC, file.path(BASE_DIR, "cellchat_HC.rds"))

# Visualization
pdf(file.path(PLOT_DIR, "CellChat_Circle_Comparison.pdf"), width = 10, height = 5)
par(mfrow=c(1,2))
netVisual_circle(cc_COPD@net$weight, title.name = "COPD Strength")
netVisual_circle(cc_HC@net$weight, title.name = "HC Strength")
dev.off()

# ==============================================================================
# 8. irGSEA Analysis (Radar Plots)
# ==============================================================================
message("Running irGSEA Scoring...")
markers <- list(
  Quiescence = c("BTG2", "KLF2", "FOXP1", "RUNX1", "MYC"),
  Regulating = c("TNFRSF9", "TNFRSF18", "ENTPD1", "FOXP3"),
  Proliferation = c("MKI67", "TYMS", "CDK1"),
  Helper = c("LTA", "TNF", "IL2", "IL18", "IFNG"),
  Cytotoxicity = c("GZMB", "GZMH", "GNLY", "PRF1"),
  Progenitor_exhaustion = c("PDCD1", "TCF7", "SLAMF6"),
  Terminal_exhaustion = c("HAVCR2", "TIGIT", "LAG3", "ENTPD1"),
  Senescence = c("KLRG1", "B3GAT1", "CDKN2A")
)

COPD <- irGSEA.score(object = COPD, assay = "RNA", slot = "data", geneset = markers, method = "AUCell", kcdf=FALSE)

# Extract scores and plot Radar
aucell_data <- as.data.frame(t(GetAssayData(COPD, assay = "AUCell", slot = "data")))
COPD <- AddMetaData(COPD, aucell_data)
state_features <- names(markers)

plot_radar_compare <- function(cell_type) {
  sub_obj <- subset(COPD, celltype == cell_type)
  df_list <- split(sub_obj@meta.data[, state_features], sub_obj$disease)
  means <- do.call(rbind, lapply(df_list, colMeans, na.rm=TRUE))
  
  data_radar <- rbind(max = rep(max(means)*1.2, length(state_features)), min = rep(0, length(state_features)), means)
  radarchart(as.data.frame(data_radar), pcol = c("red", "blue"), plwd = 2, title = paste(cell_type, "State"))
  legend("topright", legend = rownames(means), col = c("red", "blue"), lty = 1)
}

pdf(file.path(PLOT_DIR, "RTE_Radar_Plots.pdf"), width = 10, height = 5)
par(mfrow=c(1,2))
if("CD4 RTE" %in% COPD$celltype) plot_radar_compare("CD4 RTE")
if("CD8 RTE" %in% COPD$celltype) plot_radar_compare("CD8 RTE")
dev.off()

# ==============================================================================
# 9. Save and Convert to h5ad (For MEBOCOST in Script 02)
# ==============================================================================
message("Saving final object and converting to h5ad...")
saveRDS(COPD, "COPD_final.rds")

# Convert to h5ad for Python Script 02
tryCatch({
  sceasy::convertFormat(COPD, from="seurat", to="anndata", outFile="COPD_scored.h5ad")
  message("Successfully created COPD_scored.h5ad")
}, error = function(e) {
  message("Error converting to h5ad. Please check environment.")
})
