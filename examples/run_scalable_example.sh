#!/bin/bash

# Example script for running the scalable workflow
# Adjust parameters based on your system and dataset

set -euo pipefail

# Configuration
INPUT_DIR="assemblies"
OUTPUT_DIR="results_scalable"
PROFILE="local_workstation"  # or "dgx_station" for high-performance systems

# Basic scalable mode for local workstation (12 cores, 64GB RAM)
echo "Running scalable mode for large dataset analysis..."

nextflow run main.nf \
  --input ${INPUT_DIR} \
  --outdir ${OUTPUT_DIR} \
  --scalable_mode true \
  --workflow_mode cluster \
  --mash_threshold 0.03 \
  --max_cluster_size 100 \
  --run_gubbins true \
  --gubbins_iterations 3 \
  --gubbins_tree_builder hybrid \
  --build_usher_mat true \
  --profile docker,${PROFILE} \
  -resume

echo "Analysis complete! Results in ${OUTPUT_DIR}"

# Optional: Generate summary report
if [ -f "${OUTPUT_DIR}/cluster_summary.txt" ]; then
    echo ""
    echo "Clustering Summary:"
    cat "${OUTPUT_DIR}/cluster_summary.txt"
fi

# Example for DGX Station (high-performance system)
# Uncomment and modify as needed:
#
# nextflow run main.nf \
#   --input ${INPUT_DIR} \
#   --outdir ${OUTPUT_DIR}_dgx \
#   --scalable_mode true \
#   --workflow_mode cluster \
#   --mash_threshold 0.025 \
#   --max_cluster_size 150 \
#   --run_gubbins true \
#   --gubbins_iterations 3 \
#   --build_usher_mat true \
#   --profile docker,dgx_station \
#   -resume

# Example for incremental updates (future implementation):
#
# # Add new samples to existing analysis
# nextflow run main.nf \
#   --input new_assemblies/ \
#   --workflow_mode place \
#   --existing_mat ${OUTPUT_DIR}/global.pb \
#   --outdir ${OUTPUT_DIR}_updated \
#   --profile docker,${PROFILE} \
#   -resume