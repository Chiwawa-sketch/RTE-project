# ==============================================================================
# Analysis Pipeline: scRNA-seq Reference, Deconvolution, and Metabolic Profiling
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(Matrix)
library(harmony)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(BayesPrism)
library(GSVA)
library(nnet)       # For multinomial logistic regression
library(survey)     # For weighted GLM
library(cobalt)     # For covariate balance check
library(readr)
library(data.table)

# Set global options
options(stringsAsFactors = FALSE)
set.seed(1234)

# Define working directory (Adjust as needed)
work_dir <- "./Analysis_Results"
if(!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
setwd(work_dir)

# ==============================================================================
# PART 1: Single-cell Reference Construction & RTE Definition
# ==============================================================================

# 1.1 Data Loading & Integration
# ------------------------------------------------------------------------------
base_path <- "./Raw_Data/scRNA_Seq" # Generic path
sample_folders <- list.dirs(path = base_path, full.names = FALSE, recursive = FALSE)
# Filter for relevant folders (Before/RA samples)
sample_folders <- sample_folders[grepl("^Before\\d$|^RA\\d$", sample_folders)]

seurat_list <- list()

for (sample in sample_folders) {
  data_path <- file.path(base_path, sample)
  
  # Check for 10X structure
  if(dir.exists(data_path)) {
    raw_data <- Read10X(data.dir = data_path)
    if (is.list(raw_data)) raw_data <- raw_data[["Gene Expression"]]
    
    seu <- CreateSeuratObject(counts = raw_data, project = sample, min.cells = 3, min.features = 200)
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
    
    # QC Filtering
    seu <- subset(seu, subset = percent.mt < 15 & nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 25000)
    
    # Normalization
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
    
    seurat_list[[sample]] <- seu
  }
}

# Merge and Integrate
combined <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list), project = "RA_Embedding")
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 30)

# Run Harmony
combined <- RunHarmony(combined, group.by.vars = "orig.ident")

# Clustering
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.8)

# 1.2 Annotation (Mapping Clusters to Cell Types)
# ------------------------------------------------------------------------------
# Mapping based on marker gene analysis
cluster_to_celltype <- c(
  "0" = "CD4 Naive", "9" = "CD8 Naive", "1" = "CD4 Tcm",
  "2" = "CD8 Tem", "15" = "CD8 Tem", 
  "4" = "NK Cells (CD56low)", "6" = "NK Cells (CD56low)", 
  "13" = "NK Cells (CD56low)", "14" = "NK Cells (CD56high)",
  "18" = "Hybrid", "21" = "pDC", "10" = "Non-classical Mono",
  "3" = "Classical Mono", "11" = "Classical Mono", "7" = "Classical Mono",
  "17" = "Classical Mono", "20" = "Classical Mono", "19" = "Platelets",
  "5" = "B Naive", "16" = "B Naive", "8" = "B Memory", "12" = "B Memory"
)

combined$cell_type <- plyr::mapvalues(
  as.character(combined$seurat_clusters),
  from = names(cluster_to_celltype),
  to = unname(cluster_to_celltype)
)

# 1.3 RTE Definition (SOX4+ Thresholding)
# ------------------------------------------------------------------------------
# Function to find peak density for thresholding
get_peak_threshold <- function(vec) {
  dens <- density(vec)
  max_idx <- which.max(dens$y)
  return(dens$x[max_idx])
}

meta <- combined@meta.data
rna_data <- GetAssayData(combined, assay = "RNA", slot = "data")
target_gene <- "SOX4"

# Iterate through Naive subsets to define RTE
for (ctype in c("CD4 Naive", "CD8 Naive")) {
  # Get expression vector for the specific cell type
  cells_of_interest <- rownames(meta)[meta$cell_type == ctype]
  expr_vec <- as.numeric(rna_data[target_gene, cells_of_interest])
  
  # Filter for non-zero expression to determine threshold
  expr_nonzero <- expr_vec[expr_vec > 0]
  
  if (length(expr_nonzero) > 0) {
    cutoff <- get_peak_threshold(expr_nonzero)
    
    # Identify RTE cells (Naive + SOX4 > Cutoff)
    rte_cells <- cells_of_interest[expr_vec > cutoff]
    
    # Update Metadata
    if (length(rte_cells) > 0) {
      new_label <- ifelse(ctype == "CD4 Naive", "CD4 RTE", "CD8 RTE")
      combined$cell_type[match(rte_cells, colnames(combined))] <- new_label
    }
  }
}

# Save References for Deconvolution
# CD8 Reference
cd8_ref <- subset(combined, subset = cell_type %in% c("CD8 Naive", "CD8 RTE", "CD8 Tem"))
saveRDS(cd8_ref, file.path(work_dir, "CD8_Reference_scRNA.rds"))

# CD4 Reference
cd4_ref <- subset(combined, subset = cell_type %in% c("CD4 Naive", "CD4 RTE", "CD4 Tcm"))
saveRDS(cd4_ref, file.path(work_dir, "CD4_Reference_scRNA.rds"))


# ==============================================================================
# PART 2: Bulk Deconvolution (BayesPrism)
# ==============================================================================

# Note: This block demonstrates CD8. Repeat similarly for CD4.
run_bayesprism_deconv <- function(sc_ref_path, bulk_matrix, output_prefix) {
  
  # Load Single Cell
  sce <- readRDS(sc_ref_path)
  sc_dat <- t(as.matrix(GetAssayData(sce, slot = "counts")))
  cell_types <- sce$cell_type
  
  # Intersect Genes
  common_genes <- intersect(colnames(sc_dat), rownames(bulk_matrix))
  sc_dat <- sc_dat[, common_genes]
  bulk_matrix <- bulk_matrix[common_genes, ]
  
  # Initialize BayesPrism
  bp <- new.prism(
    reference = sc_dat,
    mixture = t(bulk_matrix),
    input.type = "count.matrix",
    cell.type.labels = cell_types,
    cell.state.labels = cell_types, # Treating state same as type for this resolution
    key = NULL,
    outlier.cut = 0.01,
    outlier.fraction = 0.1
  )
  
  # Run Deconvolution
  bp_res <- run.prism(prism = bp, n.cores = 10) # Adjust cores as needed
  
  # Extract Proportions
  props <- get.fraction(bp = bp_res, which.theta = "final", state.or.type = "type")
  write.csv(props, file.path(work_dir, paste0(output_prefix, "_Proportions.csv")))
  
  # Extract RTE-specific Expression Profile (Gene x Sample)
  # This is crucial for the metabolic GSVA
  rte_type <- ifelse(grepl("CD8", output_prefix), "CD8 RTE", "CD4 RTE")
  
  # get.exp returns Sample x Gene, we transpose to Gene x Sample
  rte_exp_mat <- t(get.exp(bp = bp_res, state.or.type = "type", cell.name = rte_type))
  write.csv(rte_exp_mat, file.path(work_dir, paste0(output_prefix, "_RTE_Expression.csv")))
  
  return(list(proportions = props, rte_expression = rte_exp_mat))
}

# Load Bulk Data (Example Loading)
bulk_path <- "./Raw_Data/Bulk/GSE118829_Counts.txt"
bulk_df <- fread(bulk_path, data.table = FALSE)
rownames(bulk_df) <- bulk_df[,1]
bulk_mat <- as.matrix(bulk_df[,-1])

# Execute for CD8
cd8_results <- run_bayesprism_deconv(
  sc_ref_path = file.path(work_dir, "CD8_Reference_scRNA.rds"),
  bulk_matrix = bulk_mat,
  output_prefix = "CD8"
)

# ==============================================================================
# PART 3: Metabolic Pathway Analysis (GSVA)
# ==============================================================================

# 3.1 Define Custom Metabolic Gene Sets
# ------------------------------------------------------------------------------
metabolic_genes <- list(
  "Fatty Acid Metabolism" = c("ACAA1", "ACAA2", "ACACA", "ACACB", "ACADM", "ACADS", "ACADSB", "ACADVL", "ACAT1", "ACAT2",
                              "ACLY", "ACSL3", "ACSL4", "ACSL5", "ACSL6", "ACSM3", "ACSS1", "ACSS2", "ACSS3", "CD36", "CIDEC",
                              "CPT1A", "CPT1B", "CPT2", "CRCT2", "CROT", "DGAT1", "DGAT2", "ECHS1", "ELOVL1", "ELOVL3",
                              "ELOVL4", "ELOVL5", "ELOVL6", "ELOVL7", "FABP4", "FABP5", "FASN", "FFAR2", "FFAR3", "FOXO1",
                              "GPR109A", "GPR84", "HADH", "HADHA", "HADHB", "HNF4A", "ME1", "ME2", "ME3", "OXSM", "PC", "PCK1",
                              "PCK2", "PRKAA1", "PUN1", "SCD", "SCD5", "SIRT1", "SLC25A20", "SLC27A2", "SREBF1", "SREBF2"),
  
  "Glycolysis" = c("AKT1", "AKT2", "ALDOA", "ALDOB", "ALDOC", "BPGM", "EGFR", "GAPDH", "GPI", "HIF1A", "HK1", "HK2", "IGF1R",
                   "INSR", "LDHA", "LDHB", "MTOR", "MYC", "MYCN", "PFKFB3", "PFKL", "PFKM", "PFKP", "PGK1", "PGK2", "PIK3CD",
                   "PIK3IP1", "PIK3R1", "PIK3R2", "PIK3R3", "PKM1", "PKM2", "PRKAA1", "PRKAB1", "PRKAG1", "PTEN", "RPTOR",
                   "SGK1", "SLC16A3", "SLC2A1", "SLC2A4", "SLC2A9", "SLC5A1", "STK11", "TBC1D4", "TIGAR", "TP53", "TPI1", "TSC2"),
  
  "Glunconeogenesis" = c("BCKDHA", "BCKDHB", "FBP1", "FBP2", "G6CP2", "G6PC1", "G6PC3", "GLS", "GLS2", "GLUD1", "GLUD2", "GOT1",
                         "GOT2", "GPT", "ME", "MMUT", "PC", "PCCA", "PCCB", "PCK1", "PCK2", "SIRT1", "SLC17A5", "SLC37A4"),
  
  "PPP" = c("G6PD", "H6PD", "PGD", "TALDO1", "TKT", "TKTL1", "TKTL2"),
  
  "OXPHOS/TCA" = c("ACO1", "ACO2", "CS", "FH", "IDH1", "IDH2", "IDH3", "MDH1", "MDH1B", "MDH2", "ME1", "ME2", "ME3", "MPC1",
                   "NDUFS2", "OGDH", "OXCT1", "OXCT2", "PDHA1", "PDHA2", "PDK1", "PDK2", "PDK3", "PDK4", "PDP1", "SCO2",
                   "SDHA", "SDHAF1", "SDHAF2", "SDHAF3", "SDHAF4", "SDHB", "SDHC", "SDHD", "SLC1A5", "SLC22A1", "SLC25A1", "UCP2", "PDPK1"),
  
  "One-Carbon Metabolism" = c("BHMT", "DNMT1", "DNMT3A", "DNMT3B", "DNMT3L", "GSTM1", "MTHFD1", "MTHFD2", "MTHFR", "MTR",
                              "MTRR", "PHGDH", "PSPH", "SHMT1", "SHMT2", "TP53"),
  
  "Aminoacid Metabolism" = c("AANAT", "ABAT", "ACAD8", "ACASB", "AGXT", "AGXT1", "AGXT2", "ALAD", "ALADS2", "ALAS2", "AMT",
                             "ARG1", "ARG2", "ASL", "ASNS", "ASNSD1", "ASS1", "BCAT1", "BCAT2", "BCKDHA", "BCKDHB", "BCKDK",
                             "BHMT", "CAD", "CBS", "CDO1", "CPS1", "CSAD", "CTH", "DBH", "DBT", "DDC", "DLD", "GAD1", "GAD2",
                             "GADL1", "GATD3", "GCAT", "GCDH", "GCH1", "GCLC", "GCLM", "GCSH", "GLDC", "GLS", "GLS2", "GLUL",
                             "GLYAT", "GNMT", "GOT1", "GOT1L1", "GOT2", "GPT", "GPT2", "GSR", "GSS", "HAL", "HDC", "IDO1",
                             "IDO2", "KYAT1", "LIPT1", "LIPT2", "MAT1A", "MAT2A", "MAT2B", "MMUT", "MPST", "MTHFR", "MTR",
                             "MTRR", "NAGS", "OAT", "ODC1", "OTC", "PAH", "PCCB", "PSAT1", "PTPRN", "SHMT1", "SHMT2", "SRR",
                             "TAT", "TDO2", "TGM1", "TGM2", "TH", "TIK2", "TPH1", "TPH2", "TST", "TYR", "UROC1")
)

# 3.2 Run GSVA
# ------------------------------------------------------------------------------
# Input: RTE-specific expression matrix (Gene x Sample) from BayesPrism
rte_expr <- read.csv(file.path(work_dir, "CD8_RTE_Expression.csv"), row.names = 1)
rte_mat <- as.matrix(rte_expr)

# Log transformation if necessary (BayesPrism output is usually normalized counts)
# Adding pseudocount to avoid log(0)
rte_mat_log <- log2(rte_mat + 1)

# Set GSVA parameters
# method="gsva" or "ssgsea" depending on preference. 
# kcdf="Gaussian" for continuous data (log-transformed counts)
gsva_param <- gsvaParam(exprData = rte_mat_log, geneSets = metabolic_genes, kcdf = "Gaussian")
gsva_scores <- gsva(gsva_param, verbose = TRUE)

# Save Scores
write.csv(gsva_scores, file.path(work_dir, "CD8_RTE_GSVA_Metabolic_Scores.csv"))


# ==============================================================================
# PART 4: IPTW Adjustment and Statistical Analysis
# ==============================================================================

# 4.1 Load and Merge Clinical Data
# ------------------------------------------------------------------------------
# Load Metadata (Must contain: SampleID, Group, Age, Sex)
meta_df <- read.csv("./Raw_Data/Clinical_Metadata.csv") 

# Ensure Group Levels and Reference
meta_df$Group <- factor(meta_df$Group, levels = c("Healthy", "TN", "MTX", "TCZ", "IFX"))
meta_df$Group <- relevel(meta_df$Group, ref = "TN") # Treatment Naive as Reference

# Load Targets (Proportions and GSVA Scores)
props_df <- read.csv(file.path(work_dir, "CD8_Proportions.csv"), row.names = 1)
gsva_df <- read.csv(file.path(work_dir, "CD8_RTE_GSVA_Metabolic_Scores.csv"), row.names = 1)

# Transpose GSVA to have Samples as rows
gsva_t <- as.data.frame(t(gsva_df))
colnames(gsva_t) <- rownames(gsva_df)

# Merge all data
analysis_df <- merge(meta_df, props_df, by.x = "SampleID", by.y = "row.names")
analysis_df <- merge(analysis_df, gsva_t, by.x = "SampleID", by.y = "row.names")

# 4.2 Generalized IPTW Analysis Function
# ------------------------------------------------------------------------------
run_iptw_analysis <- function(data, target_vars, treatment_col = "Group", covariates = c("Age", "Sex")) {
  
  results_list <- list()
  
  # A. Propensity Score Estimation (Multinomial Logistic Regression)
  fmla_ps <- as.formula(paste(treatment_col, "~", paste(covariates, collapse = " + ")))
  ps_model <- multinom(fmla_ps, data = data, trace = FALSE)
  
  # Predict probabilities
  ps_probs <- fitted(ps_model)
  
  # Calculate Weights (1 / Probability of observed group)
  # Create index matrix to extract the prob corresponding to the actual group assignment
  actual_group_indices <- cbind(1:nrow(data), match(data[[treatment_col]], colnames(ps_probs)))
  prob_observed <- ps_probs[actual_group_indices]
  data$iptw_weight <- 1 / prob_observed
  
  # Optional: Stabilize weights if range is too large
  # marginal_probs <- table(data[[treatment_col]]) / nrow(data)
  # prob_marginal <- marginal_probs[data[[treatment_col]]]
  # data$iptw_weight <- prob_marginal / prob_observed
  
  # B. Iterate through each target variable (Cell types or Pathways)
  for (var in target_vars) {
    
    # Define Survey Design
    design <- svydesign(ids = ~1, data = data, weights = ~iptw_weight)
    
    # Weighted GLM (Gaussian family for continuous scores/proportions)
    # Using 'TN' as reference (set in factor levels previously)
    fmla_glm <- as.formula(paste("`", var, "` ~ ", treatment_col, sep = ""))
    fit <- svyglm(fmla_glm, design = design, family = gaussian())
    
    # Extract Statistics
    coef_summary <- summary(fit)$coefficients
    conf_intervals <- confint(fit)
    
    # Format output
    res_temp <- data.frame(
      Variable = var,
      Comparison = rownames(coef_summary),
      Estimate = coef_summary[, "Estimate"],
      StdError = coef_summary[, "Std. Error"],
      P_Value = coef_summary[, "Pr(>|t|)"],
      CI_Low = conf_intervals[, 1],
      CI_High = conf_intervals[, 2]
    )
    
    # Filter out Intercept
    res_temp <- res_temp[grepl(treatment_col, res_temp$Comparison), ]
    res_temp$Comparison <- gsub(treatment_col, "", res_temp$Comparison) # Clean name
    
    results_list[[var]] <- res_temp
  }
  
  return(do.call(rbind, results_list))
}

# 4.3 Execute Analysis
# ------------------------------------------------------------------------------

# A. Analyze Cell Proportions (specifically RTE)
target_cells <- c("CD8_RTE", "CD8_Naive") # Adjust based on exact column names from BayesPrism
res_proportions <- run_iptw_analysis(analysis_df, target_cells)
write.csv(res_proportions, file.path(work_dir, "IPTW_Results_CellProportions.csv"), row.names = FALSE)

# B. Analyze Metabolic Pathways (GSVA Scores)
target_pathways <- colnames(gsva_t)
res_metabolism <- run_iptw_analysis(analysis_df, target_pathways)
write.csv(res_metabolism, file.path(work_dir, "IPTW_Results_MetabolicScores.csv"), row.names = FALSE)

# 4.4 Visualization (Boxplots with P-values)
# ------------------------------------------------------------------------------
plot_results <- function(data, y_var, title) {
  p <- ggplot(data, aes_string(x = "Group", y = paste0("`", y_var, "`"), fill = "Group")) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
    scale_fill_brewer(palette = "Set2") +
    theme_classic() +
    labs(title = title, x = "Treatment Group", y = "Score / Proportion") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  return(p)
}

# Example Plot: CD8 RTE Proportion
pdf(file.path(work_dir, "Plot_CD8_RTE_Comparison.pdf"), width = 6, height = 5)
print(plot_results(analysis_df, "CD8_RTE", "CD8 RTE Proportion (IPTW Adjusted)"))
dev.off()

# Example Plot: Fatty Acid Metabolism
pdf(file.path(work_dir, "Plot_FattyAcid_Metabolism.pdf"), width = 6, height = 5)
print(plot_results(analysis_df, "Fatty Acid Metabolism", "Fatty Acid Metabolism Score (CD8 RTE)"))
dev.off()

print("Analysis Pipeline Completed.")