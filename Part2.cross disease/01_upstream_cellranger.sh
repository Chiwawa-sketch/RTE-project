#!/bin/bash
# scripts/Part2_CrossDisease/01_upstream_cellranger.sh

# ==============================================================================
# Cell Ranger Processing for PRJNA1168183 (SpA Cohort)
# Target: Early Active Patients 1-5
# ==============================================================================

# --- 1. CONFIGURATION (PLEASE UPDATE THESE PATHS ON SERVER) ---

# Path to Cell Ranger Transcriptome Reference
# Example: /path/to/refdata-gex-GRCh38-2020-A
TRANSCRIPTOME="/path/to/your/refdata-gex-GRCh38-2020-A"

# Path containing the raw FastQ files
# Example: /path/to/raw_data/PRJNA1168183
FASTQ_DIR="/path/to/your/raw_data/PRJNA1168183"

# Output Directory
# Example: /path/to/ADT4DISEASES/cellranger_output
OUTPUT_DIR="/path/to/your/ADT4DISEASES/cellranger_output"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# --- 2. SAMPLE DEFINITION ---

# Define the 5 specific Early Active patients.
# KEY (Left) = The folder name for the output (e.g., SPA1)
# VALUE (Right) = The exact sample name prefix in your FastQ files (e.g., SRR1234567)
# Note: You must replace 'SRR_XXXX' with your actual FastQ prefixes.

declare -A SAMPLES
SAMPLES=(
    ["SPA1"]="SRR_ACCESSION_P1_REPLACE_ME"
    ["SPA2"]="SRR_ACCESSION_P2_REPLACE_ME"
    ["SPA3"]="SRR_ACCESSION_P3_REPLACE_ME"
    ["SPA4"]="SRR_ACCESSION_P4_REPLACE_ME"
    ["SPA5"]="SRR_ACCESSION_P5_REPLACE_ME"
)

# --- 3. EXECUTION LOOP ---

echo "Starting Cell Ranger Pipeline..."
echo "Output Directory: $OUTPUT_DIR"

for ID in "${!SAMPLES[@]}"; do
    SAMPLE_PREFIX="${SAMPLES[$ID]}"
    
    echo "------------------------------------------------"
    echo "Processing Sample: $ID"
    echo "FastQ Prefix:      $SAMPLE_PREFIX"
    echo "------------------------------------------------"
    
    # Navigate to output directory
    cd "$OUTPUT_DIR" || exit
    
    # Run Cell Ranger Count
    # Adjust --localcores and --localmem based on your server capacity
    cellranger count --id="$ID" \
                     --transcriptome="$TRANSCRIPTOME" \
                     --fastqs="$FASTQ_DIR" \
                     --sample="$SAMPLE_PREFIX" \
                     --localcores=16 \
                     --localmem=64
    
    if [ $? -eq 0 ]; then
        echo "SUCCESS: Sample $ID finished."
    else
        echo "ERROR: Sample $ID failed."
    fi
done

echo "------------------------------------------------"
echo "All jobs completed."
