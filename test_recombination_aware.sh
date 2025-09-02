#!/bin/bash

# Test script for recombination-aware workflow
# This script validates the new workflow implementation

set -euo pipefail

echo "🧬 Testing Recombination-Aware Assembly SNPs Workflow"
echo "======================================================"

# Check if required directories exist
if [ ! -d "test_input" ]; then
    echo "❌ Error: test_input directory not found"
    echo "Please create test_input/ with some assembly files (.fa, .fasta, .fna)"
    exit 1
fi

# Count input files
input_count=$(find test_input -name "*.fa" -o -name "*.fasta" -o -name "*.fna" | wc -l)
echo "📁 Found $input_count assembly files in test_input/"

if [ "$input_count" -lt 3 ]; then
    echo "⚠️  Warning: Need at least 3 assemblies for meaningful phylogenetic analysis"
fi

# Check Nextflow installation
if ! command -v nextflow &> /dev/null; then
    echo "❌ Error: Nextflow not found. Please install Nextflow first."
    exit 1
fi

echo "✅ Nextflow found: $(nextflow -version | head -n1)"

# Check Docker installation
if ! command -v docker &> /dev/null; then
    echo "❌ Error: Docker not found. Please install Docker first."
    exit 1
fi

echo "✅ Docker found: $(docker --version)"

# Test workflow syntax
echo ""
echo "🔍 Testing workflow syntax..."
if nextflow run main.nf --help > /dev/null 2>&1; then
    echo "✅ Workflow syntax is valid"
else
    echo "❌ Workflow syntax error detected"
    exit 1
fi

# Run dry-run test
echo ""
echo "🧪 Running dry-run test..."
nextflow run main.nf \
    -c test_recombination_aware.config \
    --recombination_aware_mode \
    --input test_input/ \
    --outdir test_recombination_results \
    -preview \
    -dump-channels || {
    echo "❌ Dry-run failed"
    exit 1
}

echo "✅ Dry-run completed successfully"

# Offer to run full test
echo ""
echo "🚀 Ready to run full test?"
echo "This will execute the complete recombination-aware workflow."
echo ""
read -p "Run full test? (y/N): " -n 1 -r
echo

if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "🏃 Running full recombination-aware workflow test..."
    
    # Clean previous results
    if [ -d "test_recombination_results" ]; then
        echo "🧹 Cleaning previous results..."
        rm -rf test_recombination_results
    fi
    
    # Run the workflow
    nextflow run main.nf \
        -c test_recombination_aware.config \
        --recombination_aware_mode \
        --input test_input/ \
        --outdir test_recombination_results \
        -with-report test_recombination_results/pipeline_info/execution_report.html \
        -with-timeline test_recombination_results/pipeline_info/execution_timeline.html \
        -with-dag test_recombination_results/pipeline_info/pipeline_dag.html
    
    if [ $? -eq 0 ]; then
        echo ""
        echo "🎉 SUCCESS! Recombination-aware workflow completed successfully!"
        echo ""
        echo "📊 Results summary:"
        echo "- Output directory: test_recombination_results/"
        echo "- Final grafted tree: test_recombination_results/Final_Results/global_grafted.treefile"
        echo "- Execution report: test_recombination_results/pipeline_info/execution_report.html"
        echo ""
        echo "🔍 Key output files:"
        find test_recombination_results -name "*.treefile" -o -name "*.gff" -o -name "*.tsv" | head -10
    else
        echo "❌ Workflow execution failed"
        exit 1
    fi
else
    echo "ℹ️  Skipping full test. Use the following command to run manually:"
    echo ""
    echo "nextflow run main.nf \\"
    echo "  -c test_recombination_aware.config \\"
    echo "  --recombination_aware_mode \\"
    echo "  --input test_input/ \\"
    echo "  --outdir test_recombination_results"
fi

echo ""
echo "📚 For more information, see:"
echo "- RECOMBINATION_AWARE_WORKFLOW.md"
echo "- test_recombination_aware.config"
echo ""
echo "✨ Recombination-aware workflow testing complete!"