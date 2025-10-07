#!/usr/bin/env bash

# Simple wrapper script for running wf-assembly-snps-mod
# This provides an easy interface for common usage patterns

set -euo pipefail

SCRIPT_NAME="$(basename "$0")"

# Default parameters
PROFILE="docker"
MODE="standard"
INPUT=""
OUTDIR="results"
EXTRA_ARGS=""

usage() {
    cat << EOF
Usage: $SCRIPT_NAME [OPTIONS] --input INPUT_PATH

A wrapper script for running wf-assembly-snps-mod with common configurations.

REQUIRED:
    --input PATH            Path to input directory or samplesheet

OPTIONS:
    --outdir PATH          Output directory (default: results)
    --profile PROFILE      Execution profile (default: docker)
                          Options: docker, singularity, conda, local_workstation
    --mode MODE           Analysis mode (default: standard)
                          Options: standard, scalable, recombination
    --help                Show this help message

EXAMPLES:
    # Standard analysis with Docker
    $SCRIPT_NAME --input assemblies/ --outdir my_results

    # Scalable mode for large datasets  
    $SCRIPT_NAME --input assemblies/ --mode scalable --profile local_workstation

    # Recombination-aware analysis
    $SCRIPT_NAME --input assemblies/ --mode recombination

    # Custom parameters (passed through)
    $SCRIPT_NAME --input assemblies/ --mode scalable -- --mash_threshold 0.025 --max_cluster_size 80

For more advanced usage, run nextflow directly:
    nextflow run PHemarajata/wf-assembly-snps-mod --help

EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input)
            INPUT="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --profile)
            PROFILE="$2"
            shift 2
            ;;
        --mode)
            MODE="$2"
            shift 2
            ;;
        --help)
            usage
            exit 0
            ;;
        --)
            # Everything after -- gets passed to nextflow
            shift
            EXTRA_ARGS="$*"
            break
            ;;
        *)
            echo "Error: Unknown option $1"
            usage
            exit 1
            ;;
    esac
done

# Check required parameters
if [[ -z "$INPUT" ]]; then
    echo "Error: --input is required"
    usage
    exit 1
fi

# Build mode-specific parameters
MODE_PARAMS=""
case $MODE in
    standard)
        # Default mode, no additional parameters needed
        ;;
    scalable)
        MODE_PARAMS="--scalable_mode true"
        ;;
    recombination)
        MODE_PARAMS="--recombination_aware_mode true"
        ;;
    *)
        echo "Error: Unknown mode '$MODE'. Use: standard, scalable, or recombination"
        exit 1
        ;;
esac

# Construct and run the nextflow command
NEXTFLOW_CMD="nextflow run PHemarajata/wf-assembly-snps-mod"
NEXTFLOW_CMD+=" -profile $PROFILE"
NEXTFLOW_CMD+=" --input $INPUT"
NEXTFLOW_CMD+=" --outdir $OUTDIR"

if [[ -n "$MODE_PARAMS" ]]; then
    NEXTFLOW_CMD+=" $MODE_PARAMS"
fi

if [[ -n "$EXTRA_ARGS" ]]; then
    NEXTFLOW_CMD+=" $EXTRA_ARGS"
fi

echo "=========================================="
echo "🧬 wf-assembly-snps-mod Workflow Launcher"
echo "=========================================="
echo
echo "Mode:     $MODE"
echo "Profile:  $PROFILE"
echo "Input:    $INPUT"
echo "Output:   $OUTDIR"
echo
echo "Running command:"
echo "$NEXTFLOW_CMD"
echo
echo "=========================================="

# Execute the command
eval "$NEXTFLOW_CMD"