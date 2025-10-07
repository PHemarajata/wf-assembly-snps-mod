#!/bin/bash

# Example script for running standalone tree grafting
# This demonstrates how to use graft_trees.py after a scalable workflow run

set -euo pipefail

echo "🌳 Standalone Tree Grafting Example"
echo "=================================="

# Check if required files exist
if [ ! -d "results" ]; then
    echo "❌ Error: 'results' directory not found"
    echo "Please run the scalable workflow first to generate the required input files"
    exit 1
fi

# Define expected input files
BACKBONE="results/backbone.treefile"
CLUSTERS_PATTERN="results/Clusters/**/cluster_*.final.treefile"
REPRESENTATIVES="results/cluster_representatives.tsv"

# Define output files
OUTDIR="results/tree_grafting"
OUTPUT_TREE="${OUTDIR}/final_grafted.treefile"
REPORT="${OUTDIR}/grafting_report.txt"
LOG="${OUTDIR}/grafting_log.txt"

# Create output directory
mkdir -p "${OUTDIR}"

# Check if backbone tree exists
if [ ! -f "$BACKBONE" ]; then
    echo "❌ Error: Backbone tree not found at $BACKBONE"
    echo "Make sure the scalable workflow completed successfully"
    exit 1
fi

# Check if cluster trees exist
CLUSTER_COUNT=$(find results/Clusters -name "cluster_*.final.treefile" 2>/dev/null | wc -l)
if [ "$CLUSTER_COUNT" -eq 0 ]; then
    echo "❌ Error: No cluster trees found matching pattern: $CLUSTERS_PATTERN"
    echo "Make sure the scalable workflow completed successfully"
    exit 1
fi

echo "📊 Found $CLUSTER_COUNT cluster trees"
echo "📁 Backbone tree: $BACKBONE"
echo "📋 Representatives file: $REPRESENTATIVES"
echo

# Check if Biopython is installed
if ! python -c "import Bio.Phylo" 2>/dev/null; then
    echo "❌ Error: Biopython is required but not installed"
    echo "Install with: pip install biopython"
    exit 1
fi

# Run the tree grafting script
echo "🚀 Running tree grafting..."
echo

# First, do a dry run to check what would be done
echo "Step 1: Dry run (checking inputs and planning grafting)"
./graft_trees.py \
    --backbone "$BACKBONE" \
    --clusters "$CLUSTERS_PATTERN" \
    --reps "$REPRESENTATIVES" \
    --out-tree "$OUTPUT_TREE" \
    --report "$REPORT" \
    --log "$LOG" \
    --rename-conflicts \
    --dry-run

echo
echo "Step 2: Actual tree grafting"

# Run the actual grafting
./graft_trees.py \
    --backbone "$BACKBONE" \
    --clusters "$CLUSTERS_PATTERN" \
    --reps "$REPRESENTATIVES" \
    --out-tree "$OUTPUT_TREE" \
    --report "$REPORT" \
    --log "$LOG" \
    --rename-conflicts \
    --parent-edge-mode keep

# Check results
if [ -f "$OUTPUT_TREE" ]; then
    echo
    echo "✅ Tree grafting completed successfully!"
    echo
    echo "📄 Results:"
    echo "  - Final tree: $OUTPUT_TREE"
    echo "  - Report: $REPORT"
    echo "  - Log: $LOG"
    echo
    
    # Get tree statistics
    echo "📊 Final tree statistics:"
    python -c "
import sys
try:
    from Bio import Phylo
    tree = Phylo.read('$OUTPUT_TREE', 'newick')
    print(f'  - Total tips: {len(tree.get_terminals())}')
    print(f'  - Tree depth: {tree.depth():.4f}')
    print(f'  - Total clades: {tree.count_terminals() + len(list(tree.get_nonterminals()))}')
except Exception as e:
    print(f'  - Error reading tree: {e}')
"
    
    echo
    echo "📋 Grafting summary:"
    cat "$REPORT"
    
else
    echo
    echo "❌ Tree grafting failed!"
    echo "Check the log file for details: $LOG"
    exit 1
fi

echo
echo "🎉 Tree grafting example completed!"
echo
echo "Next steps:"
echo "  1. Visualize the tree: Use FigTree, iTOL, or ggtree"
echo "  2. Validate the tree: Check tip labels and tree structure"
echo "  3. Compare results: Compare with the original backbone tree"