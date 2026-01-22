import os
import scanpy as sc
import pandas as pd
from mebocost import mebocost
import warnings

warnings.filterwarnings("ignore")

# ==============================================================================
# Configuration
# ==============================================================================
# Automatically determine project root relative to this script
# Assumes script is in: [Project_Root]/scripts/Part2_CrossDisease/
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(BASE_DIR, "results", "part2_cross_disease")
OUTPUT_DIR = os.path.join(DATA_DIR, "MEBOCOST")

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# ==============================================================================
# Data Loading
# ==============================================================================
h5ad_file = os.path.join(DATA_DIR, "merged_all.h5ad")

if not os.path.exists(h5ad_file):
    # Fallback to check absolute path if relative fails
    if os.path.exists("/mnt/data/home/tyu-a1ouoisq/ADT4DISEASES/results/part2_cross_disease/merged_all.h5ad"):
        h5ad_file = "/mnt/data/home/tyu-a1ouoisq/ADT4DISEASES/results/part2_cross_disease/merged_all.h5ad"
    else:
        raise FileNotFoundError(f"Input file not found: {h5ad_file}")

print(f"Loading data: {os.path.basename(h5ad_file)}")
adata = sc.read_h5ad(h5ad_file)

if adata.raw is not None:
    adata = adata.raw.to_adata()

# ==============================================================================
# Analysis Loop
# ==============================================================================
target_groups = ['Healthy', 'RA', 'PSA', 'SPA']

for group in target_groups:
    print(f"\nProcessing: {group}")
    
    if group not in adata.obs['Status'].unique():
        continue
        
    subset = adata[adata.obs['Status'] == group].copy()
    
    if subset.n_obs < 10:
        continue

    try:
        comm = mebocost.create_obj(
            adata = subset,
            group_col = 'cell_type',
            met_est = 'mebocost',
            species = 'human'
        )
        
        comm.infer_comm(
            n_shuffle = 1000, 
            seed = 12345, 
            cutoff_exp = 0.05, 
            cutoff_met = 0.05
        )

        if comm.comm_res is not None:
            res = comm.comm_res.copy()
            
            # Filter for RTE-associated interactions
            cols = res.columns
            if 'Sender' in cols and 'Receiver' in cols:
                mask = (res['Sender'].str.contains('RTE')) | (res['Receiver'].str.contains('RTE'))
                res_filtered = res[mask].copy()
                res_filtered['Status'] = group
                
                out_path = os.path.join(OUTPUT_DIR, f"MEBOCOST_Events_{group}.csv")
                res_filtered.to_csv(out_path, index=False)
                print(f"Output saved: {os.path.basename(out_path)}")

    except Exception as e:
        print(f"Skipping {group} due to error: {str(e)}")

print("\nAnalysis completed successfully.")
