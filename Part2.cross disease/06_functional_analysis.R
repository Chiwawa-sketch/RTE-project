# scripts/Part2_CrossDisease/06_functional_analysis.R

# ==============================================================================
# 1. Setup and Configuration
# ==============================================================================
library(Seurat)
library(dplyr)
library(Matrix)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

BASE_DIR <- "/path/to/ADT4DISEASES" # UPDATE THIS
OUTPUT_DIR <- file.path(BASE_DIR, "results/part2_cross_disease")
PLOT_DIR <- file.path(OUTPUT_DIR, "Functional_Analysis")

if(!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)

# Load Final Object from Step 05
merged_all <- readRDS(file.path(OUTPUT_DIR, "merged_all_final_with_RTE.rds"))

# ==============================================================================
# 2. GO & KEGG Enrichment Analysis (Per Methodology)
# ==============================================================================
message("Running Differential Expression and Enrichment Analysis...")

# Define comparisons: RA vs Healthy, etc.
# Focus on CD4 RTE and CD8 RTE
target_cells <- c("CD4 RTE", "CD8 RTE")

for (cell in target_cells) {
  message(paste("Analyzing enrichment for:", cell))
  
  # Subset object
  sub_obj <- subset(merged_all, subset = cell_type == cell)
  Idents(sub_obj) <- "Status"
  
  # 1. Identify DEGs (RA vs Healthy as primary example from text)
  # You can expand this logic for other comparisons
  degs <- FindMarkers(sub_obj, ident.1 = "RA", ident.2 = "Healthy", 
                      logfc.threshold = 0.25, min.pct = 0.1)
  
  # Filter significant genes
  sig_genes <- rownames(degs[degs$p_val_adj < 0.05, ])
  
  if (length(sig_genes) > 0) {
    # Convert Symbol to Entrez ID
    gene_ids <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    
    # 2. GO Enrichment (BP)
    ego <- enrichGO(gene = gene_ids$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE)
    
    if (!is.null(ego)) {
      pdf(file.path(PLOT_DIR, paste0("GO_Enrichment_", gsub(" ", "_", cell), "_RAvsHealthy.pdf")), width = 10, height = 8)
      print(dotplot(ego, showCategory = 20) + ggtitle(paste("GO:BP -", cell, "(RA vs Healthy)")))
      dev.off()
      write.csv(as.data.frame(ego), file.path(PLOT_DIR, paste0("GO_Table_", gsub(" ", "_", cell), ".csv")))
    }
    
    # 3. KEGG Enrichment
    ekegg <- enrichKEGG(gene = gene_ids$ENTREZID,
                        organism = 'hsa',
                        pvalueCutoff = 0.05)
    
    if (!is.null(ekegg)) {
      pdf(file.path(PLOT_DIR, paste0("KEGG_Enrichment_", gsub(" ", "_", cell), "_RAvsHealthy.pdf")), width = 10, height = 8)
      print(dotplot(ekegg, showCategory = 20) + ggtitle(paste("KEGG -", cell, "(RA vs Healthy)")))
      dev.off()
    }
  }
}

# ==============================================================================
# 3. Metabolic Gene Set Definition
# ==============================================================================
metabolic_genes <- list(
  "Fatty_Acid_Metabolism" = c("ACAA1", "ACAA2", "ACACA", "ACACB", "ACADM", "ACADS", "ACADSB", "ACADVL", "ACAT1", "ACAT2",
                              "ACLY", "ACSL3", "ACSL4", "ACSL5", "ACSL6", "ACSM3", "ACSS1", "ACSS2", "ACSS3", "CD36", "CIDEC",
                              "CPT1A", "CPT1B", "CPT2", "CRCT2", "CROT", "DGAT1", "DGAT2", "ECHS1", "ELOVL1", "ELOVL3",
                              "ELOVL4", "ELOVL5", "ELOVL6", "ELOVL7", "FABP4", "FABP5", "FASN", "FFAR2", "FFAR3", "FOXO1",
                              "GPR109A", "GPR84", "HADH", "HADHA", "HADHB", "HNF4A", "ME1", "ME2", "ME3", "OXSM", "PC", "PCK1",
                              "PCK2", "PRKAA1", "PUN1", "SCD", "SCD5", "SIRT1", "SLC25A20", "SLC27A2", "SREBF1", "SREBF2"),
  "Glycolysis" = c("AKT1", "AKT2", "ALDOA", "ALDOB", "ALDOC", "BPGM", "EGFR", "GAPDH", "GPI", "HIF1A", "HK1", "HK2", "IGF1R",
                   "INSR", "LDHA", "LDHB", "MTOR", "MYC", "MYCN", "PFKFB3", "PFKL", "PFKM", "PFKP", "PGK1", "PGK2", "PIK3CD",
                   "PIK3IP1", "PIK3R1", "PIK3R2", "PIK3R3", "PKM1", "PKM2", "PRKAA1", "PRKAB1", "PRKAG1", "PTEN", "RPTOR",
                   "SGK1", "SLC16A3", "SLC2A1", "SLC2A4", "SLC2A9", "SLC5A1", "STK11", "TBC1D4", "TIGAR", "TP53", "TPI1", "TSC2"),
  "Gluconeogenesis" = c("BCKDHA", "BCKDHB", "FBP1", "FBP2", "G6CP2", "G6PC1", "G6PC3", "GLS", "GLS2", "GLUD1", "GLUD2", "GOT1",
                        "GOT2", "GPT", "ME1", "MMUT", "PC", "PCCA", "PCCB", "PCK1", "PCK2", "SIRT1", "SLC17A5", "SLC37A4"),
  "PPP" = c("G6PD", "H6PD", "PGD", "TALDO1", "TKT", "TKTL1", "TKTL2"),
  "OXPHOS_TCA" = c("ACO1", "ACO2", "CS", "FH", "IDH1", "IDH2", "IDH3", "MDH1", "MDH1B", "MDH2", "ME1", "ME2", "ME3", "MPC1",
                   "NDUFS2", "OGDH", "OXCT1", "OXCT2", "PDHA1", "PDHA2", "PDK1", "PDK2", "PDK3", "PDK4", "PDP1", "SCO2",
                   "SDHA", "SDHAF1", "SDHAF2", "SDHAF3", "SDHAF4", "SDHB", "SDHC", "SDHD", "SLC1A5", "SLC22A1", "SLC25A1", "UCP2", "PDPK1"),
  "One_Carbon_Metabolism" = c("BHMT", "DNMT1", "DNMT3A", "DNMT3B", "DNMT3L", "GSTM1", "MTHFD1", "MTHFD2", "MTHFR", "MTR",
                              "MTRR", "PHGDH", "PSPH", "SHMT1", "SHMT2", "TP53"),
  "Aminoacid_Metabolism" = c("AANAT", "ABAT", "ACAD8", "ACASB", "AGXT", "AGXT1", "AGXT2", "ALAD", "ALADS2", "ALAS2", "AMT",
                             "ARG1", "ARG2", "ASL", "ASNS", "ASNSD1", "ASS1", "BCAT1", "BCAT2", "BCKDHA", "BCKDHB", "BCKDK",
                             "BHMT", "CAD", "CBS", "CDO1", "CPS1", "CSAD", "CTH", "DBH", "DBT", "DDC", "DLD", "GAD1", "GAD2",
                             "GADL1", "GATD3", "GCAT", "GCDH", "GCH1", "GCLC", "GCLM", "GCSH", "GLDC", "GLS", "GLS2", "GLUL",
                             "GLYAT", "GNMT", "GOT1", "GOT1L1", "GOT2", "GPT", "GPT2", "GSR", "GSS", "HAL", "HDC", "IDO1",
                             "IDO2", "KYAT1", "LIPT1", "LIPT2", "MAT1A", "MAT2A", "MAT2B", "MMUT", "MPST", "MTHFR", "MTR",
                             "MTRR", "NAGS", "OAT", "ODC1", "OTC", "PAH", "PCCB", "PSAT1", "PTPRN", "SHMT1", "SHMT2", "SRR",
                             "TAT", "TDO2", "TGM1", "TGM2", "TH", "TIK2", "TPH1", "TPH2", "TST", "TYR", "UROC1")
)

# Convert to mapping frame
gene_long <- stack(metabolic_genes)
colnames(gene_long) <- c("Gene", "Pathway")
all_metabolic_genes <- unique(gene_long$Gene)

# Color scale definition (0 to 1)
col_fun <- colorRamp2(breaks = c(0, 0.5, 1), colors = c("#4575B4", "white", "#D73027"))

# ==============================================================================
# 4. Helper Function: Plot and Save Heatmap
# ==============================================================================
plot_metabolic_heatmap <- function(expr_matrix, pathway_name, cell_type_name, output_folder) {
  
  # Calculate dynamic height based on gene number
  n_genes <- nrow(expr_matrix)
  
  # Define base height plus dynamic component (approx 0.15 inch per gene)
  plot_height <- max(3, n_genes * 0.15) 
  
  file_name <- file.path(output_folder, paste0("Heatmap_", gsub(" ", "_", cell_type_name), "_", pathway_name, ".pdf"))
  
  pdf(file_name, width = 4, height = plot_height)
  
  ht <- Heatmap(
    expr_matrix,
    name = "Expression",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_side = "left",
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 8, fontface = "italic"),
    heatmap_legend_param = list(
      title = "LogNorm Expr",
      title_gp = gpar(fontsize = 10),
      grid_width = unit(4, "mm"),
      legend_height = unit(3.5, "cm")
    ),
    rect_gp = gpar(col = "grey90"),
    # Add borders to cells
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid.rect(x, y, w, h, gp = gpar(fill = NA, col = "grey90"))
    }
  )
  
  draw(ht, column_title = paste0(pathway_name, " - ", cell_type_name))
  dev.off()
}

# ==============================================================================
# 5. Loop Execution: Generate Heatmaps for CD4 and CD8 RTE
# ==============================================================================
message("Generating Metabolic Heatmaps...")

target_cell_types <- c("CD4 RTE", "CD8 RTE")

for (ctype in target_cell_types) {
  message(paste("Processing:", ctype))
  
  # 1. Subset Object
  # Ensure we only take relevant statuses
  cells_to_keep <- merged_all@meta.data %>%
    filter(cell_type == ctype, Status %in% c("Healthy", "RA", "PSA", "SPA")) %>%
    rownames()
  
  if(length(cells_to_keep) < 10) {
    warning(paste("Not enough cells for", ctype, "- skipping."))
    next
  }
  
  sub_obj <- subset(merged_all, cells = cells_to_keep)
  
  # 2. Extract Expression Data
  # Get data slot (LogNormalized)
  raw_data <- GetAssayData(sub_obj, assay = "RNA", slot = "data")
  
  # Intersect available genes
  valid_genes <- intersect(all_metabolic_genes, rownames(raw_data))
  expr_mat <- as.matrix(raw_data[valid_genes, ])
  
  # 3. Calculate Average Expression per Group
  group_info <- sub_obj@meta.data[colnames(expr_mat), "Status"]
  group_factor <- factor(group_info, levels = c("Healthy", "RA", "PSA", "SPA"))
  
  # Calculate row means per group
  expr_avg <- sapply(levels(group_factor), function(g) {
    cells_in_group <- names(group_info)[group_info == g]
    if(length(cells_in_group) == 0) return(rep(0, nrow(expr_mat)))
    if(length(cells_in_group) == 1) return(expr_mat[, cells_in_group])
    rowMeans(expr_mat[, cells_in_group])
  })
  
  colnames(expr_avg) <- c("Healthy", "RA", "PSA", "SPA")
  
  # 4. Loop through each Pathway and Plot
  for (path_name in names(metabolic_genes)) {
    path_genes <- metabolic_genes[[path_name]]
    
    # Filter genes present in our dataset
    genes_to_plot <- intersect(path_genes, rownames(expr_avg))
    
    if (length(genes_to_plot) == 0) {
      message(paste("Skipping", path_name, "- no matching genes found."))
      next
    }
    
    # Subset matrix for this pathway
    mat_subset <- expr_avg[genes_to_plot, , drop = FALSE]
    
    # Generate Heatmap
    plot_metabolic_heatmap(mat_subset, path_name, ctype, PLOT_DIR)
  }
}

message("Analysis Complete!")
message("1. GO/KEGG Enrichment tables and plots saved to: ", PLOT_DIR)
message("2. Metabolic Heatmaps saved to: ", PLOT_DIR)