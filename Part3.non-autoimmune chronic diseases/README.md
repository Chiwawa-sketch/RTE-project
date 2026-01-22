# Part 3: Comparison of RTE Remodeling in Other Chronic Diseases

This module extends the analysis to three comparator diseases (COPD, T2D, MG) to determine if the coordinated RTE remodeling observed in RA is specific to its pathology or a general feature of chronic inflammation.

## Directory Structure

scripts/Part3_DiseaseComparison/
├── 01_COPD_analysis.R       # [R] Analysis of COPD (GSE249584)
├── 02_MEBOCOST_COPD.py      # [Python] Metabolic communication in COPD
├── 03_T2D_analysis.R        # [R] Analysis of Type 2 Diabetes (GSE255566)
├── 04_MEBOCOST_T2D.py       # [Python] Metabolic communication in T2D
├── 05_MG_analysis.R         # [R] Analysis of Myasthenia Gravis (GSE227835)
└── 06_MEBOCOST_MG.py        # [Python] Metabolic communication in MG

## Workflow & Methodology

### 1. COPD Analysis (GSE249584)
- Script: 01_COPD_analysis.R
  - Processing: Integration of 15 samples (Control vs COPD).
  - Annotation: Manual annotation using canonical markers (e.g., Classical Mono, Naive B).
  - RTE Definition: SOX4 thresholds (CD4 > 1.25, CD8 > 1.18).
  - Functional: CellChat interaction analysis and irGSEA scoring (Radar plots).
  - Output: `COPD_scored.h5ad` for metabolic analysis.

- Script: 02_MEBOCOST_COPD.py
  - Analysis: Infers metabolic sender-receiver events for HC and COPD groups.

### 2. T2D Analysis (GSE255566)
- Script: 03_T2D_analysis.R
  - Processing: Integration of 6 samples (Con vs Mod).
  - Annotation: Automated Module Score annotation based on marker signatures.
  - RTE Definition: SOX4 thresholds (CD4 > 1.09, CD8 > 0.94).
  - Functional: CellChat and irGSEA scoring.
  - Output: `T2D_scored.h5ad`.

- Script: 04_MEBOCOST_T2D.py
  - Analysis: Metabolic communication inference for T2D vs HC.

### 3. MG Analysis (GSE227835)
- Script: 05_MG_analysis.R
  - Processing: Parallel loading of .txt.gz files, RPCA integration.
  - Annotation: Detailed manual mapping (including Regulatory T, Thymic-like T).
  - RTE Definition: SOX4 thresholds (CD4 > 0.98, CD8 > 0.93).
  - Functional: CellChat and irGSEA scoring.
  - Output: `GSE227835_scored.h5ad`.

- Script: 06_MEBOCOST_MG.py
  - Analysis: Metabolic communication inference for MG vs HC.

## Dependencies

- R: Seurat, CellChat, irGSEA, data.table, future, sceasy
- Python: scanpy, mebocost, pandas

## Key Outputs

- UMAP Plots: Visualizing cell type distributions across conditions.
- Marker Plots: Feature/Violin/Dot plots validating annotations (including Treg markers).
- Radar Plots: Comparing RTE functional states (Quiescence, Cytotoxicity, etc.).
- CellChat Plots: Circle plots showing differential interaction strength.
- MEBOCOST Tables: CSV files of predicted metabolic communication events.
