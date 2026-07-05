#!/bin/bash

# =============================================================
# Pipeline: Antimicrobial resistance gene identification
# Tool: RGI (Resistance Gene Identifier) - bwt module
# Aligner: KMA
# Database: CARD (local)
# Parameters: --include_other_models --include_wildcard -n 30
# =============================================================

# --- Paths (adjust according to the environment) ---
WORK_DIR="/path/to/working/directory"   # directory where outputs will be saved
DATA_DIR="/path/to/data/directory"      # directory containing reads per sample

# --- Samples to process ---
# Add or remove samples as needed
samples=(
    "AB1097"
    "AB1235"
    # "AB1236"
    # ...
)

# --- Pipeline ---
cd "$WORK_DIR"

for sample_name in "${samples[@]}"; do
    sample_dir="${DATA_DIR}/${sample_name}"

    if [ -d "$sample_dir" ]; then
        R1=$(find "$sample_dir" -type f -name "*_1.fq.gz" | head -n 1)
        R2=$(find "$sample_dir" -type f -name "*_2.fq.gz" | head -n 1)

        if [[ -f "$R1" && -f "$R2" ]]; then
            echo "Processing sample: $sample_name"
            output_dir="${WORK_DIR}/output_${sample_name}"
            mkdir -p "$output_dir"

            rgi bwt \
                --read_one "$R1" \
                --read_two "$R2" \
                --output_file "${output_dir}/${sample_name}_rgi" \
                --aligner kma \
                --local \
                --include_other_models \
                --include_wildcard \
                -n 30
        else
            echo "No read files found for sample: $sample_name"
        fi
    else
        echo "Sample folder not found: $sample_name"
    fi
done

echo "Pipeline completed."
