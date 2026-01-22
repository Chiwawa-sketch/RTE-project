import os
import scanpy as sc
import pandas as pd
import numpy as np
from mebocost import mebocost
import warnings
warnings.filterwarnings("ignore")

# ==============================================================================
# Configuration
# ==============================================================================
BASE_DIR = "/mnt/data/home/tyu-a1ouoisq/GSE227835_RAW"
DATA_FILE = os.path.join(BASE_DIR, "GSE227835_scored.h5ad")
OUTPUT_DIR = os.path.join(BASE_DIR, "MEBOCOST_Results")

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# ==============================================================================
# Helper Function
# ==============================================================================
def run_mebocost(adata_subset, group_name):
    print(f"\n>>> Running MEBOCOST for: {group_name}")
    
    # Run MEBOCOST
    comm_obj = mebocost.create_obj(
        adata = adata_subset,
        group_col = ['celltype'],
        met_est = 'mebocost',
        species = 'human'
    )
    
    # Standard parameters
    comm_obj.infer_comm(
        n_shuffle = 1000,
        seed = 12345,
        cutoff_exp = 0.05,
        cutoff_met = 0.05
    )
    
    # Save Results
    if comm_obj.comm_res is not None:
        res = comm_obj.comm_res.copy()
        res['Group'] = group_name
        
        # Filter for RTE interactions
        mask = (res['Sender'].str.contains('RTE')) | (res['Receiver'].str.contains('RTE'))
        res_filtered = res[mask].copy()
        
        out_csv = os.path.join(OUTPUT_DIR, f"MEBOCOST_{group_name}_RTE_Events.csv")
        res_filtered.to_csv(out_csv, index=False)
        print(f"Saved interactions to: {out_csv}")
        
    return comm_obj

# ==============================================================================
# Main Execution
# ==============================================================================
print(f"Loading data: {DATA_FILE}")
if not os.path.exists(DATA_FILE):
    print("Error: h5ad file not found. Please run Script 05 first.")
    exit(1)

adata = sc.read_h5ad(DATA_FILE)

# Ensure metadata is clean
if 'celltype' in adata.obs.columns:
    adata.obs['celltype'] = adata.obs['celltype'].astype(str)

# Split by Disease Group (HC vs MG)
hc_adata = adata[adata.obs['disease'] == 'HC'].copy()
mg_adata = adata[adata.obs['disease'] == 'MG'].copy()

# Run Analysis
if hc_adata.n_obs > 0:
    run_mebocost(hc_adata, "HC")
if mg_adata.n_obs > 0:
    run_mebocost(mg_adata, "MG")

print("\nStep 06 Complete! MG MEBOCOST analysis finished.")
