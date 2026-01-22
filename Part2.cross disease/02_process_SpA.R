# scripts/Part2_CrossDisease/02_process_SpA.R

library(Seurat)
library(dplyr)
library(Matrix)

# ==============================================================================
# 1. Configuration
# ==============================================================================
BASE_DIR <- "/path/to/ADT4DISEASES" # UPDATE THIS
CELLRANGER_DIR <- file.path(BASE_DIR, "cellranger_output")
OUTPUT_DIR <- file.path(BASE_DIR, "results/part2_cross_disease")
if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# The 5 specific Early Active patients defined in Step 01
spa_samples <- c("SPA1", "SPA2", "SPA3", "SPA4", "SPA5")

# ==============================================================================
# 2. Import & Quality Control
# ==============================================================================
message("Processing SpA Cohort (Early Active 1-5)...")

spa_list <- list()

for (sample_id in spa_samples) {
  # Path to Cell Ranger output
  matrix_dir <- file.path(CELLRANGER_DIR, sample_id, "outs", "filtered_feature_bc_matrix")
  
  if (!dir.exists(matrix_dir)) {
    warning(paste("Matrix not found for:", sample_id))
    next
  }
  
  counts <- Read10X(data.dir = matrix_dir)
  
  seu <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 3, min.features = 200)
  seu$sample <- sample_id
  seu$source <- "SpA_PRJNA1168183"
  seu$Status <- "SPA"
  
  # QC: < 10% Mitochondrial (Strict)
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mt < 10)
  
  # Pre-processing
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  
  spa_list[[sample_id]] <- seu
  message(paste("Processed:", sample_id))
}

# Merge and Save
spa_merged <- merge(x = spa_list[[1]], y = spa_list[-1], add.cell.ids = names(spa_list), project = "SpA_Cohort")
saveRDS(spa_merged, file = file.path(OUTPUT_DIR, "SpA_EarlyActive_processed.rds"))