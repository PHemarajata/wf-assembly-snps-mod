# HPC-Specific Execution Scripts

This repository includes several HPC-specific execution scripts that were developed for institutional use. These scripts are provided for reference and may be adapted for other HPC environments.

## Available Scripts

### `_run_snp_identification.uge-nextflow`
- **Purpose**: Generic SNP identification wrapper for UGE (Univa Grid Engine)
- **Environment**: CDC institutional HPC clusters (Aspen/Rosalind)
- **Features**:
  - Automatic profile selection based on node number
  - Support for reference-based and reference-free analysis
  - Integrated error handling and logging

### `run_Parsnp_GENOMES.uge-nextflow`
- **Purpose**: Comprehensive genome panel analysis using Parsnp
- **Features**:
  - Full workflow execution with job management
  - Automatic file validation and filtering
  - Resource optimization for different dataset sizes
  - Comprehensive logging and error reporting

### `run_Parsnp_REFERENCE_vs_GENOMES.uge-nextflow`
- **Purpose**: Reference-based analysis comparing reference against genome panel
- **Features**:
  - Reference genome validation
  - Comparative analysis workflow
  - Optimized for reference-guided SNP calling

## Configuration Requirements

These scripts expect specific environment variables and infrastructure:

### Required Environment Variables
- `USER`: Username for scratch directory paths
- `HOSTNAME`: For automatic cluster detection
- `SNP_PACKAGE`: SNP calling method (default: parsnp)
- `GENOMES`: Input genome directory
- `OUT`: Output directory
- `REFERENCE`: Reference genome (for reference-based analysis)

### Required Modules
- `nextflow`: Nextflow execution environment

### Expected Infrastructure
- Scratch workspace at `/scicomp/scratch/${USER}/work`
- Email notification system configured for `${USER}@cdc.gov`
- UGE scheduler with proper resource allocation

## Adaptation for Other HPC Systems

To adapt these scripts for your HPC environment:

1. **Modify scheduler directives**: Update `#$` lines for your scheduler (SLURM, PBS, etc.)
2. **Update module loading**: Change `module load nextflow` to your environment setup
3. **Adjust paths**: Update scratch and work directory paths
4. **Configure profiles**: Modify profile selection logic in `nextflow.config`
5. **Update notifications**: Change email settings and notification methods

### Example SLURM Adaptation

```bash
#!/bin/bash
#SBATCH --job-name=wf-assembly-snps
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

module load nextflow

nextflow run PHemarajata/wf-assembly-snps-mod \
  -profile singularity \
  --input ${INPUT} \
  --outdir ${OUTPUT} \
  -work-dir ${SLURM_TMPDIR}/work
```

## Generic Usage

For general use without HPC-specific requirements, use the provided wrapper script:

```bash
./run_workflow.sh --input assemblies/ --mode scalable --profile docker
```

Or run Nextflow directly:

```bash
nextflow run PHemarajata/wf-assembly-snps-mod \
  -profile docker \
  --input assemblies/ \
  --outdir results \
  --scalable_mode true
```

## Support

These HPC scripts are provided as-is for reference. For support with general workflow usage, please refer to the main [README.md](README.md) and [documentation](docs/).

For HPC-specific adaptations, consult your local HPC support team and system documentation.