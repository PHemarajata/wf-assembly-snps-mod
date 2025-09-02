# Error Fixes Summary

## Issues Fixed

### 1. EXTRACT_CORE_SNPS Module Error ✅

**Problem:**
```
ERROR ~ Error executing process > 'ASSEMBLY_SNPS_SCALABLE:INTEGRATE_RESULTS:EXTRACT_CORE_SNPS (cluster_cluster_13)'
Caused by: No such variable: Extract -- Check script './workflows/../subworkflows/local/../../modules/local/extract_core_snps/main.nf' at line: 21
```

**Root Cause:**
The docstring `"""Extract core SNP positions from cluster alignment"""` inside the Python function was being interpreted by Nextflow as variable substitution because it was inside the script block.

**Fix Applied:**
- Replaced the triple-quoted docstring with a regular comment: `# Extract core SNP positions from cluster alignment`
- This prevents Nextflow from trying to interpret the docstring as variable substitution

### 2. Gubbins Empty Output Files ✅

**Problem:**
Gubbins processes were completing but producing zero-size output files:
- `cluster_18.filtered_polymorphic_sites.fasta` (0 bytes)
- `cluster_18.recombination_predictions.gff` (0 bytes)
- `cluster_18.node_labelled.final_tree.tre` (0 bytes)

**Root Cause:**
- Gubbins was failing silently due to insufficient variable sites or other issues
- The original module created empty files when Gubbins failed but didn't provide diagnostic information

**Fix Applied:**
- Added comprehensive diagnostic logging to show:
  - Number of sequences in alignment
  - Alignment length
  - Number of variable sites found
  - Gubbins exit codes
  - Output file sizes
- Added Python script to count variable sites before running Gubbins
- Added verbose output from Gubbins
- Better error handling and reporting when Gubbins fails

## Enhanced Diagnostics

### EXTRACT_CORE_SNPS Module
- Fixed variable substitution error
- Maintains all original functionality
- Now compatible with Nextflow's script parsing

### GUBBINS_CLUSTER Module
- Added detailed logging for debugging
- Reports sequence counts, alignment length, and variable sites
- Shows Gubbins exit codes and output file sizes
- Provides warnings when files are empty or missing
- Added `--verbose` flag to Gubbins for better error reporting

## Expected Behavior After Fixes

### EXTRACT_CORE_SNPS
- Should now run without the "No such variable: Extract" error
- Will successfully extract core SNPs from cluster alignments
- Produces proper output files for downstream integration

### GUBBINS_CLUSTER
- Will provide detailed diagnostic information in the log
- Will show why Gubbins is producing empty files (e.g., insufficient variable sites)
- Will continue to create empty files when Gubbins fails (to maintain pipeline flow)
- But now you'll know WHY the files are empty from the diagnostic output

## Testing the Fixes

To verify the fixes work:

1. **Check EXTRACT_CORE_SNPS**: The error should no longer occur and the process should complete successfully
2. **Check GUBBINS diagnostics**: Look at the process logs to see detailed information about why files might be empty

Example diagnostic output you should now see:
```
Starting Gubbins analysis for cluster cluster_18
Number of sequences in alignment: 5
Alignment length: 1234 characters
Variable sites found: 2
WARNING: Only 2 variable sites found, less than minimum required (5)
Gubbins may not produce meaningful results
```

This will help you understand whether clusters have sufficient variation for meaningful Gubbins analysis.

## Files Modified

1. `modules/local/extract_core_snps/main.nf` - Fixed docstring variable substitution error
2. `modules/local/gubbins_cluster/main.nf` - Added comprehensive diagnostics and error reporting

Both fixes maintain the original functionality while resolving the errors and providing better debugging information.