R
# scripts/Part2_CrossDisease/03_integration.R
library(Seurat)
library(harmony)
library(dplyr)
# ==============================================================================
# 1. Configuration
# ==============================================================================
BASE_DIR <- "/path/to/ADT4DISEASES" 
OUTPUT_DIR <- file.path(BASE_DIR, "results/part2_cross_disease")
if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
# The specific individuals with the highest cell counts (as identified in Metadata)
# We hardcode them here to ensure reproducibility.
SELECTED_HC_PSA <- c("HC25", "HC26", "HC27", "PSA27", "PSA05", "PSA21")
# ==============================================================================
# 2. Load Datasets
# ==============================================================================
seurat_list <- list()
# --- A. Load SpA (Processed in Step 02) ---
spa_path <- file.path(OUTPUT_DIR, "SpA_EarlyActive_processed.rds")
if(file.exists(spa_path)) {
  seurat_list[["SPA"]] <- readRDS(spa_path)
} else { stop("SpA object not found.") }
# --- B. Load HC/PsA (With Selection & Metadata Fix) ---
hc_psa_path <- file.path(BASE_DIR, "filtered_subjects_subset_obj.rds")
if (file.exists(hc_psa_path)) {
  hc_psa <- readRDS(hc_psa_path)
  hc_psa$source <- "HC_PSA"
  # 1. REPORT SELECTION
  message("----------------------------------------------------")
  message("SELECTION REPORT (Top Cell Counts):")
  message(paste("Selected Subjects:", paste(SELECTED_HC_PSA, collapse=", ")))
  message("----------------------------------------------------")
  # 2. FILTER: Keep only selected subjects
  if ("Subject" %in% colnames(hc_psa@meta.data)) {
    hc_psa <- subset(hc_psa, subset = Subject %in% SELECTED_HC_PSA)
    # 3. METADATA FIX: Map 'Subject' to 'sample'
    # This ensures the sample column reads "HC25" instead of the original batch ID
    hc_psa$sample <- hc_psa$Subject
    message("Metadata fixed: 'sample' column updated with Subject IDs.")
  } else {
    stop("Error: 'Subject' column not found in HC/PsA object.")
  }
  # 4. QC
  hc_psa[["percent.mt"]] <- PercentageFeatureSet(hc_psa, pattern = "^MT-")
  seurat_list[["HC_PSA"]] <- subset(hc_psa, subset = percent.mt < 10 & nFeature_RNA > 200)
} else { stop("HC/PsA object not found.") }
# --- C. Load RA (E-MTAB-14192) ---
ra_samples <- c("RA1", "RA2", "RA3")
for (s in ra_samples) {
  path <- file.path(BASE_DIR, s)
  if (dir.exists(path)) {
    counts <- Read10X(path)
    seu <- CreateSeuratObject(counts, project = s)
    seu$sample <- s
    seu$source <- "RA_EMTAB14192"
    seu$Status <- "RA"
    
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
    seu <- subset(seu, subset = percent.mt < 10 & nFeature_RNA > 200)
    seu <- NormalizeData(seu, verbose = FALSE) %>% FindVariableFeatures(verbose = FALSE)
    seurat_list[[s]] <- seu
  }
}
# ==============================================================================
# 3. Merge Datasets
# ==============================================================================
message("Merging datasets (SpA + Selected HC/PsA + RA)...")
merged_all <- merge(x = seurat_list[[1]], y = seurat_list[-1], project = "CrossDisease")
# ==============================================================================
# 4. Integration & Processing
# ==============================================================================
message("Running Harmony Pipeline...")
DefaultAssay(merged_all) <- "RNA"
merged_all <- NormalizeData(merged_all, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE)
merged_all <- RunHarmony(merged_all, group.by.vars = "source", verbose = TRUE)
merged_all <- RunUMAP(merged_all, reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.8, verbose = FALSE)
# ==============================================================================
# 5. Save
# ==============================================================================
saveRDS(merged_all, file.path(OUTPUT_DIR, "merged_all_harmony.rds"))
message("Step 03 Complete! Integrated object saved.")