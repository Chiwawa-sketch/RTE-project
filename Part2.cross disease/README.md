# Part 2: Cross-Disease Transcriptional and Metabolic Profiling of RTEs

This module focuses on the comparative analysis of Recent Thymic Emigrants (RTEs) across three autoimmune arthritis conditions (Rheumatoid Arthritis, Psoriatic Arthritis, Spondyloarthritis) and Healthy controls. The analysis pipeline integrates multi-cohort single-cell data to uncover distinct immune signatures, metabolic programming, and intercellular communication networks.

## Data Sources

The analysis integrates single-cell RNA sequencing (scRNA-seq) data from the following cohorts:

- SpA (Spondyloarthritis): PRJNA1168183 (Selected 5 Early Active patients).
- RA (Rheumatoid Arthritis): E-MTAB-14192 (Treatment-naive).
- PsA (Psoriatic Arthritis) & Healthy: E-MTAB-9492.
  (Specific high-quality donors HC25, HC26, HC27, PSA27, PSA05, PSA21 were selected based on cell count depth to ensure balanced integration).

## Directory Structure

The analysis code is organized sequentially in scripts/Part2_CrossDisease/.

scripts/Part2_CrossDisease/
├── 01_upstream_cellranger.sh   # [Bash] Upstream processing of raw FastQ files
├── 02_process_SpA.R            # [R] QC and preprocessing of SpA cohort
├── 03_integration.R            # [R] Cohort merging, metadata harmonization, and Harmony integration
├── 04_annotation_plots.R       # [R] Cell type annotation and global UMAP visualization
├── 05_rte_analysis.R           # [R] RTE definition (SOX4+), thresholding, and composition analysis
├── 06_functional_analysis.R    # [R] GO/KEGG enrichment and Metabolic Pathway Heatmaps
├── 07_mebocost_analysis.py     # [Python] Metabolic cell-cell communication inference
├── 08_sankey_visualization.R   # [R] Visualization of metabolic communication (Sankey diagrams)
└── 09_cellchat_analysis.R      # [R] Ligand-Receptor interaction analysis (Supplementary Data 2-5)

## Analysis Workflow

### 1. Upstream Processing & QC
- Script: 01_upstream_cellranger.sh
  - Action: Maps reads to the GRCh38 reference using Cell Ranger.
- Script: 02_process_SpA.R
  - Action: Performs strict QC (<10% mitochondrial reads) on SpA samples and filters for "Early Active" patients.

### 2. Integration & Annotation
- Script: 03_integration.R
  - Action: Merges RA, SpA, and selected HC/PsA datasets. Applies Harmony to correct for batch effects.
- Script: 04_annotation_plots.R
  - Action: Annotates major cell lineages and generates Fig 3A (Global UMAP).

### 3. RTE Definition & Characterization
- Script: 05_rte_analysis.R
  - Action: Defines RTEs based on SOX4 expression (Valley method) and visualizes composition (Fig 3C, 3D).

### 4. Functional & Metabolic Profiling
- Script: 06_functional_analysis.R
  - Action: GO/KEGG enrichment and Metabolic Pathway Heatmaps.
- Script: 07_mebocost_analysis.py (Python)
  - Action: Infers metabolic cell-cell communication.
- Script: 08_sankey_visualization.R
  - Action: Visualizes metabolic flows via Sankey diagrams (Fig 3G, 3H).

### 5. Ligand-Receptor Analysis (CellChat)
- Script: 09_cellchat_analysis.R
  - Action:
    - Infers Ligand-Receptor interactions for each disease group.
    - Exports interaction tables (Supplementary Data 2-5).
    - Generates Circle Plots visualizing aggregated network strength and specific RTE targets.

## Key Outputs

- Fig 3A/B: Global UMAP and Proportions.
- Fig 3C/D: RTE Identification and Composition.
- Fig 3G/H: Metabolic Sankey Diagrams.
- Supplementary Data 2-5: CSV tables of significant Ligand-Receptor interactions (Step 09).
- Circle Plots: Visualization of intercellular communication networks (Step 09).

## Dependencies

R Packages:
- Seurat, Harmony, dplyr, ggplot2, patchwork
- CellChat, NMF, Circlize, ComplexHeatmap
- clusterProfiler, org.Hs.eg.db, ggalluvial

Python Packages:
- scanpy, mebocost, pandas

## Usage Notes

- Ensure the BASE_DIR path in all scripts matches your local environment.
- Script 09 requires the annotated object from Step 05 (merged_all_final_with_RTE.rds).
