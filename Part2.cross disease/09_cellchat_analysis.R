library(Seurat)
library(CellChat)
library(patchwork)
library(dplyr)
library(ggplot2)

# ==============================================================================
# Configuration
# ==============================================================================
# Automatically determine project root relative to this script
args <- commandArgs(trailingOnly=FALSE)
file_arg <- grep("--file=", args, value=TRUE)
if (length(file_arg) > 0) {
  script_dir <- dirname(sub("--file=", "", file_arg))
} else {
  script_dir <- getwd()
}

base_dir <- normalizePath(file.path(script_dir, "..", ".."))
input_dir <- file.path(base_dir, "results", "part2_cross_disease")
output_dir <- file.path(input_dir, "CellChat_Analysis")

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# 1. Load Data
# ==============================================================================
# Load the final object from Step 05 (contains RTE annotations)
target_file <- file.path(input_dir, "merged_all_final_with_RTE.rds")

if(!file.exists(target_file)) {
  # Fallback to absolute path check if relative fails (for safety)
  if(file.exists("/mnt/data/home/tyu-a1ouoisq/ADT4DISEASES/results/part2_cross_disease/merged_all_final_with_RTE.rds")){
     target_file <- "/mnt/data/home/tyu-a1ouoisq/ADT4DISEASES/results/part2_cross_disease/merged_all_final_with_RTE.rds"
  } else {
     stop("Input file merged_all_final_with_RTE.rds not found. Please run Step 05 first.")
  }
}

message("Loading annotated Seurat object...")
seurat_obj <- readRDS(target_file)
Idents(seurat_obj) <- "cell_type"

# ==============================================================================
# 2. Define CellChat Analysis Function
# ==============================================================================
run_cellchat_custom <- function(object_subset, group_name, out_path) {
  
  message(paste0("\n>>> Running CellChat for: ", group_name))
  
  # A. Create CellChat Object
  data.input <- GetAssayData(object_subset, assay = "RNA", slot = "data")
  meta <- object_subset@meta.data
  
  # Ensure cell names match
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type")
  
  # B. Set Database
  CellChatDB <- CellChatDB.human 
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
  cellchat@DB <- CellChatDB.use
  
  # C. Preprocessing
  cellchat <- subsetData(cellchat) 
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # D. Inference
  # min.cells=10 ensures robustness
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # E. Save Results
  saveRDS(cellchat, file.path(out_path, paste0("cellchat_", group_name, ".rds")))
  
  # F. Export Tables (Supplementary Data)
  # This corresponds to Supplementary Data 2-5
  df_net <- subsetCommunication(cellchat)
  write.csv(df_net, file.path(out_path, paste0("SuppData_Interactions_", group_name, ".csv")))
  
  return(cellchat)
}

# ==============================================================================
# 3. Main Loop: Process Each Disease Status
# ==============================================================================
disease_groups <- c("Healthy", "RA", "PSA", "SPA")
cellchat_list <- list()

for (group in disease_groups) {
  if (group %in% unique(seurat_obj$Status)) {
    sub_obj <- subset(seurat_obj, subset = Status == group)
    
    # Run Pipeline
    tryCatch({
      cc_res <- run_cellchat_custom(sub_obj, group, output_dir)
      cellchat_list[[group]] <- cc_res
    }, error = function(e) {
      message(paste("Error processing", group, ":", e$message))
    })
  }
}

# ==============================================================================
# 4. Visualization (Circle Plots)
# ==============================================================================
message("\nGenerating Circle Plots...")
plot_dir <- file.path(output_dir, "CirclePlots")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

# Focus on RTEs as requested
target_groups <- c("CD4 RTE", "CD8 RTE")

for (group in names(cellchat_list)) {
  cc <- cellchat_list[[group]]
  groupSize <- as.numeric(table(cc@idents))
  
  # 1. Aggregated View
  pdf(file.path(plot_dir, paste0("Circle_Aggregated_", group, ".pdf")), width = 8, height = 8)
  netVisual_circle(cc@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                   title.name = paste(group, "- Total Interactions"))
  dev.off()
  
  # 2. Specific RTE Targets
  available_cells <- levels(cc@idents)
  
  for (target in target_groups) {
    if (target %in% available_cells) {
      # Incoming signals to RTE
      pdf(file.path(plot_dir, paste0("Circle_Target_", gsub(" ", "_", target), "_", group, ".pdf")), width = 10, height = 10)
      netVisual_circle(cc@net$count, targets.use = target, vertex.weight = groupSize, weight.scale = T, label.edge = F, 
                       title.name = paste(group, "- Incoming to", target))
      dev.off()
      
      # Outgoing signals from RTE
      pdf(file.path(plot_dir, paste0("Circle_Source_", gsub(" ", "_", target), "_", group, ".pdf")), width = 10, height = 10)
      netVisual_circle(cc@net$count, sources.use = target, vertex.weight = groupSize, weight.scale = T, label.edge = F, 
                       title.name = paste(group, "- Outgoing from", target))
      dev.off()
    }
  }
}

message("Step 09 Complete! Interaction tables and plots saved.")
