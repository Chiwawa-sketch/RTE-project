library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(Matrix)
library(cowplot)
library(RColorBrewer)
library(data.table)
library(future)
library(future.apply)
library(irGSEA)
library(fmsb)
library(CellChat)
library(sceasy)
library(reticulate)

# ==============================================================================
# 1. Configuration
# ==============================================================================
BASE_DIR <- "/mnt/data/home/tyu-a1ouoisq/GSE227835_RAW"
setwd(BASE_DIR)

PLOT_DIR <- file.path(BASE_DIR, "MG_Plots")
if(!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, showWarnings = FALSE)

# Parallel Setup
plan("multisession", workers = 10) 
options(future.globals.maxSize = 200 * 1024^3)

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
# 2. Data Loading
# ==============================================================================
sample_files <- list.files(BASE_DIR, pattern = "\\.txt\\.gz$", full.names = TRUE)
sample_names <- tools::file_path_sans_ext(basename(sample_files))
sample_names <- gsub("\\.txt$", "", sample_names)

message("Loading samples in parallel...")
seurat_list <- future_lapply(seq_along(sample_files), function(i) {
  sample_file <- sample_files[i]
  sample_name <- sample_names[i]
  
  dt <- fread(sample_file, header = TRUE, sep = "\t", data.table = FALSE, check.names = FALSE, showProgress = FALSE, nThread = 2)
  rownames(dt) <- dt[[1]]
  dt <- dt[,-1]
  data <- as.matrix(dt)
  if (ncol(data) > nrow(data)) data <- t(data)
  
  seu <- CreateSeuratObject(counts = data, project = sample_name, min.cells = 3, min.features = 200)
  seu$sample <- sample_name
  seu$Group <- ifelse(grepl("H", sample_name, ignore.case = TRUE), "HC",
                      ifelse(grepl("A", sample_name, ignore.case = TRUE), "MG", "Unknown"))
  
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & 
                nCount_RNA > 500 & nCount_RNA < 30000 & percent.mt < 15)
  
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  return(seu)
})
names(seurat_list) <- sample_names

# ==============================================================================
# 3. Integration (RPCA) & Clustering
# ==============================================================================
message("Running Integration (RPCA)...")
seurat_list <- lapply(seurat_list, function(x) {
  x <- ScaleData(x, verbose = FALSE)
  x <- RunPCA(x, npcs = 30, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30, anchor.features = 2000, reduction = "rpca")
GSE227835 <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(GSE227835) <- "integrated"

GSE227835 <- ScaleData(GSE227835, vars.to.regress = c("percent.mt", "nCount_RNA", "nFeature_RNA"), verbose = FALSE)
GSE227835 <- RunPCA(GSE227835, npcs = 50, verbose = FALSE)
GSE227835 <- FindNeighbors(GSE227835, reduction = "pca", dims = 1:50)
GSE227835 <- FindClusters(GSE227835, resolution = 0.9)
GSE227835 <- RunUMAP(GSE227835, reduction = "pca", dims = 1:50)

# Add disease metadata
GSE227835$disease <- GSE227835$Group

pdf(file.path(PLOT_DIR, "UMAP_Clusters_Initial.pdf"), width = 10, height = 8)
DimPlot(GSE227835, group.by = "seurat_clusters", label = TRUE, label.size = 5) + NoLegend()
dev.off()

# ==============================================================================
# 4. Marker Visualization (Critical Step Before Annotation)
# ==============================================================================
message("Generating Marker Plots...")
DefaultAssay(GSE227835) <- "RNA"
GSE227835 <- NormalizeData(GSE227835, verbose = FALSE)
GSE227835 <- ScaleData(GSE227835, verbose = FALSE)

marker_genes <- list(
  "Naive B" = c("MS4A1", "IGHD", "CD19"),
  "Memory B" = c("MS4A1", "CD27"),
  "Plasma B" = c("MZB1", "XBP1", "PRDM1"),
  "Pathogenic B" = c("TNFRSF13B", "TCL1A", "FCRL5", "CD72", "CD180"),
  "CD4 T" = c("CD3E", "CD4", "IL7R", "CCR7", "TCF7"),
  "CD8 T" = c("CD8A", "CD8B", "CD3E", "GZMK", "NKG7"),
  "T naive" = c("CCR7", "SELL", "LEF1", "TCF7", "IL7R"),
  "T central memory" = c("CCR7", "IL7R", "CD27", "CD28", "LEF1", "TCF7"),
  "T effector memory" = c("GZMK", "GZMB", "GNLY", "NKG7", "PRF1"),
  "T exhausted" = c("PDCD1", "LAG3", "HAVCR2", "TOX", "TIGIT"),
  "Regulatory T" = c("FOXP3", "IL2RA", "CTLA4", "IKZF2"),
  "Autoimmune Treg-like" = c("FOXP3", "CTLA4", "IKZF2", "IL2RA", "TOX", "PDCD1", "LAG3"),
  "Thymic-like T" = c("AIRE", "FOXN1", "KRT5", "KRT14"),
  "Effector T" = c("IFNG", "TNF", "GZMB", "PRF1", "CD69"),
  "Proliferating T" = c("MKI67", "TOP2A", "TYMS", "STMN1", "PCNA", "UBE2C"),
  "ISG-high" = c("ISG15", "MX1", "IFIT1", "IFI6", "IFIT3"),
  "Th1" = c("TBX21", "IFNG", "STAT1"),
  "Th17" = c("RORC", "IL17A", "CCR6"),
  "Tfh" = c("BCL6", "CXCR5", "PDCD1", "MAF", "ICOS"),
  "Trm" = c("ITGAE", "CD69", "ZNF683"),
  "GammaDelta T" = c("TRDC", "TRGC1", "TRGC2", "KLRB1", "S1PR1"),
  "MAIT" = c("SLC4A10", "TRAV1-2", "KLRB1", "RORC", "IL18RAP"),
  "NKT" = c("ZBTB16", "TRAV10", "CD3E", "KLRK1", "TYROBP"),
  "Classical Mono" = c("CD14", "LYZ", "S100A8", "S100A9"),
  "Non-classical Mono" = c("FCGR3A", "MS4A7", "CX3CR1"),
  "Intermediate Mono" = c("CD14", "FCGR3A", "HLA-DRA", "CD74"),
  "Myeloid Progenitor-like" = c("CD34", "CSF1R", "MPO", "ITGAM"),
  "pDC" = c("IRF8", "CLEC4C", "TCF4", "LILRA4"),
  "cDC1" = c("CLEC9A", "BATF3", "XCR1"),
  "cDC2" = c("CD1C", "FCER1A", "IRF4"),
  "NK" = c("NKG7", "GNLY", "GZMB", "KLRD1", "KLRF1", "NCAM1", "PRF1"),
  "Platelets" = c("PPBP", "PF4", "ITGA2B", "GP1BA"),
  "DNT Cells" = c("CD3E", "CD3D", "CD3G", "TRAC", "TRBC1", "TRBC2"),
  "Intermediate/Unknown" = c("TYROBP", "CD7", "IL7R", "FCGR3A", "TRDC", "ISG15", "MX1", "PTPRC")
)

marker_out_dir <- file.path(PLOT_DIR, "MarkerPlots")
if(!dir.exists(marker_out_dir)) dir.create(marker_out_dir, showWarnings = FALSE)

# Generate Plots Loop
for (group in names(marker_genes)) {
  genes <- marker_genes[[group]]
  genes <- genes[genes %in% rownames(GSE227835)]
  if (length(genes) == 0) next
  
  safe_name <- gsub("[^A-Za-z0-9]", "_", group)
  
  # FeaturePlot
  pdf(file.path(marker_out_dir, paste0("Feature_", safe_name, ".pdf")), width = 12, height = 16)
  print(FeaturePlot(GSE227835, features = genes, cols = c("grey90", "red"), ncol = 2) + NoLegend())
  dev.off()
  
  # VlnPlot
  pdf(file.path(marker_out_dir, paste0("Vln_", safe_name, ".pdf")), width = 12, height = 13)
  print(VlnPlot(GSE227835, features = genes, pt.size = 0, ncol = 2) + NoLegend())
  dev.off()
}

# DotPlot
all_markers <- unique(unlist(marker_genes))
all_markers <- all_markers[all_markers %in% rownames(GSE227835)]
if(length(all_markers) > 0) {
  pdf(file.path(marker_out_dir, "DotPlot_All_Markers.pdf"), width = 16, height = 20)
  print(DotPlot(GSE227835, features = all_markers, group.by = "seurat_clusters") + RotatedAxis())
  dev.off()
}

# ==============================================================================
# 5. Annotation (Manual Mapping based on Visualization)
# ==============================================================================
message("Applying Cluster Annotations...")
cluster2celltype <- c(
  "31" = "Naive B Cells", "34" = "Naive B Cells", "2" = "Naive B Cells",
  "11" = "Memory B Cells", "21" = "Memory B Cells",
  "28" = "Basophils", "5" = "Classical Monocyte", "8" = "Classical Monocyte", "26" = "Classical Monocyte",
  "17" = "Non-classical Monocyte", "24" = "cDC2", "29" = "pDC",
  "13" = "Platelet", "15" = "Platelet", "23" = "Neutrophils", "20" = "Erythroid",
  "35" = "Mast Cells", "33" = "Plasma B", "30" = "Proliferation T",
  "0" = "CD4 Tcm", "25" = "CD4 Tcm", "19" = "CD4 Tcm", "36" = "CD4 Tcm", "1" = "CD4 Tcm",
  "4" = "CD4 Naive", "14" = "Autoimmunne-related Treg",
  "6" = "CD8 Naive", "7" = "CD8 Tcm", "32" = "CD8 Tcm",
  "9" = "CD8 Tem", "12" = "CD8 Tem", "16" = "CD8 Tem",
  "18" = "NKT",
  "3" = "NK Cells (CD56 low)", "10" = "NK Cells (CD56 low)", "27" = "NK Cells (CD56 low)",
  "22" = "NK Cells (CD56 high)"
)

Idents(GSE227835) <- "seurat_clusters"
GSE227835$celltype <- plyr::mapvalues(as.character(Idents(GSE227835)), from = names(cluster2celltype), to = as.character(cluster2celltype))

pdf(file.path(PLOT_DIR, "UMAP_Annotated.pdf"), width = 10, height = 8)
DimPlot(GSE227835, group.by = "celltype", label = TRUE, label.size = 4) + NoLegend()
dev.off()

# ==============================================================================
# 6. RTE Definition (MG Thresholds)
# ==============================================================================
message("Defining RTE based on SOX4 expression...")
DefaultAssay(GSE227835) <- "RNA"
GSE227835 <- JoinLayers(GSE227835)

expr_gene <- FetchData(GSE227835, vars = "SOX4")[, 1]
GSE227835$expr_SOX4 <- expr_gene

# Visual check of density peaks
plot_data <- GSE227835@meta.data %>% filter(celltype %in% c("CD8 Naive", "CD4 Naive"), expr_SOX4 > 0)
pdf(file.path(PLOT_DIR, "SOX4_Density_Peaks.pdf"), width = 7, height = 5)
print(ggplot(plot_data, aes(x = expr_SOX4)) + geom_density(fill="steelblue", alpha=0.4) + facet_wrap(~celltype) + theme_minimal())
dev.off()

# Apply Thresholds
GSE227835$SOX4_status <- NA
GSE227835$SOX4_status[GSE227835$celltype == "CD4 Naive" & expr_gene > 0.98] <- "SOX4+ CD4 Naive"
GSE227835$SOX4_status[GSE227835$celltype == "CD4 Naive" & expr_gene <= 0.98] <- "SOX4- CD4 Naive"
GSE227835$SOX4_status[GSE227835$celltype == "CD8 Naive" & expr_gene > 0.93] <- "SOX4+ CD8 Naive"
GSE227835$SOX4_status[GSE227835$celltype == "CD8 Naive" & expr_gene <= 0.93] <- "SOX4- CD8 Naive"

GSE227835$celltype <- as.character(GSE227835$celltype)
GSE227835$celltype[GSE227835$SOX4_status == "SOX4+ CD4 Naive"] <- "CD4 RTE"
GSE227835$celltype[GSE227835$SOX4_status == "SOX4- CD4 Naive"] <- "CD4 T Naive"
GSE227835$celltype[GSE227835$SOX4_status == "SOX4+ CD8 Naive"] <- "CD8 RTE"
GSE227835$celltype[GSE227835$SOX4_status == "SOX4- CD8 Naive"] <- "CD8 T Naive"

write.csv(table(GSE227835$sample, GSE227835$celltype), file.path(BASE_DIR, "Cell_Counts_Sample_RTE_MG.csv"))

# ==============================================================================
# 7. CellChat Analysis
# ==============================================================================
message("Running CellChat Analysis...")

run_cellchat_subset <- function(seurat_obj, group_name) {
  data_input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
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

cc_MG <- run_cellchat_subset(subset(GSE227835, disease == "MG"), "MG")
cc_HC <- run_cellchat_subset(subset(GSE227835, disease == "HC"), "HC")

saveRDS(cc_MG, file.path(BASE_DIR, "cellchat_MG.rds"))
saveRDS(cc_HC, file.path(BASE_DIR, "cellchat_HC_MGcontrol.rds"))

pdf(file.path(PLOT_DIR, "CellChat_Circle_Comparison_MG.pdf"), width = 10, height = 5)
par(mfrow=c(1,2))
netVisual_circle(cc_MG@net$weight, title.name = "MG Strength")
netVisual_circle(cc_HC@net$weight, title.name = "HC Strength")
dev.off()

# ==============================================================================
# 8. irGSEA Analysis (Radar Plots)
# ==============================================================================
message("Running irGSEA Scoring...")
markers_func <- list(
  Quiescence = c("BTG2", "KLF2", "FOXP1", "RUNX1", "MYC"),
  Regulating = c("TNFRSF9", "TNFRSF18", "ENTPD1", "IKZF4", "FOXP3"),
  Proliferation = c("MKI67", "TYMS", "TTK", "CKS1B"),
  Helper = c("TNF", "IL2", "IFNG", "TBX21"),
  Cytotoxicity = c("GZMB", "PRF1", "GNLY", "GZMK"),
  Progenitor_exhaustion = c("PDCD1", "TCF7", "SLAMF6"),
  Terminal_exhaustion = c("HAVCR2", "TIGIT", "LAG3"),
  Senescence = c("KLRG1", "B3GAT1", "CDKN2A")
)

GSE227835 <- irGSEA.score(object = GSE227835, assay = "RNA", slot = "data", geneset = markers_func, method = "AUCell", kcdf=FALSE)

aucell_data <- as.data.frame(t(GetAssayData(GSE227835, assay = "AUCell", slot = "data")))
GSE227835 <- AddMetaData(GSE227835, aucell_data)
state_features <- names(markers_func)

plot_radar_compare <- function(cell_type) {
  if(!cell_type %in% GSE227835$celltype) return(NULL)
  sub_obj <- subset(GSE227835, celltype == cell_type)
  df_list <- split(sub_obj@meta.data[, state_features], sub_obj$disease)
  
  if(length(df_list) < 2) return(NULL)
  
  means <- do.call(rbind, lapply(df_list, colMeans, na.rm=TRUE))
  data_radar <- rbind(max = rep(max(means)*1.2, length(state_features)), min = rep(0, length(state_features)), means)
  radarchart(as.data.frame(data_radar), pcol = c("red", "blue"), plwd = 2, title = paste(cell_type, "State"))
  legend("topright", legend = rownames(means), col = c("red", "blue"), lty = 1)
}

pdf(file.path(PLOT_DIR, "RTE_Radar_Plots_MG.pdf"), width = 10, height = 5)
par(mfrow=c(1,2))
plot_radar_compare("CD4 RTE")
plot_radar_compare("CD8 RTE")
dev.off()

# ==============================================================================
# 9. Save and Convert to h5ad
# ==============================================================================
message("Saving final object and converting to h5ad...")
saveRDS(GSE227835, "GSE227835_final.rds")

tryCatch({
  sceasy::convertFormat(GSE227835, from="seurat", to="anndata", outFile="GSE227835_scored.h5ad")
  message("Successfully created GSE227835_scored.h5ad")
}, error = function(e) {
  message("Error converting to h5ad. Please check environment.")
})
