# Pipeline Fixes Summary

## Issues Fixed

### 1. Graft Tree Failure - Sample Name Mismatch

**Problem**: The graft tree process failed because representative sample names in the backbone tree didn't match the representative IDs in the mapping file.

**Root Cause**: 
- Backbone tree contained original sequence headers like `'IP-0030-7_S10_L001-SPAdes_1.fa'`
- Representative mapping file showed simplified IDs like `IP-0030`
- GRAFT_SUBTREES couldn't match these different naming conventions

**Fixes Applied**:

1. **Modified `modules/local/collect_representatives/main.nf`**:
   - Standardized sequence headers in the combined FASTA file
   - Now uses representative ID as the sequence header instead of original complex names
   - This ensures backbone tree uses consistent, simplified naming

2. **Enhanced `modules/local/graft_subtrees/main.nf`**:
   - Added robust representative matching logic
   - First tries exact match, then falls back to partial matching
   - Provides detailed logging of available labels when matching fails
   - Better error reporting for debugging

### 2. Resume Optimization for KEEP_INVARIANT_ATCG

**Problem**: The KEEP_INVARIANT_ATCG step was re-running unnecessarily when resuming the pipeline.

**Fixes Applied**:

1. **Enhanced `modules/local/keep_invariant_atcg/main.nf`**:
   - Added `cache 'lenient'` directive for better caching
   - Added `storeDir` to persist outputs in a dedicated cache directory
   - Implemented checksum-based optimization:
     - Creates MD5 checksum of input alignment
     - Skips processing if output exists with same input checksum
     - Stores checksum for future resume operations

2. **Added configuration parameters in `conf/params.config`**:
   - `enable_resume_optimization = true`: Control resume optimization
   - `work_cache_dir`: Configurable cache directory location
   - `recombination_aware_mode`: Enable recombination-aware workflow
   - Additional backbone tree and grafting parameters

## Configuration Changes

### New Parameters Added:

```groovy
// Resume optimization parameters
enable_resume_optimization        = true
work_cache_dir                    = "${params.outdir}/work_cache"

// Recombination-aware workflow parameters
recombination_aware_mode          = false
alignment_method                  = "snippy"
backbone_method                   = "parsnp"
backbone_fasttree_opts            = "-nt -gtr"
backbone_container                = null
rescale_grafted_branches          = false
```

## Expected Results

1. **Graft Tree Success**: 
   - Representative samples should now be found in backbone tree
   - Tree grafting should complete successfully
   - Better error messages if issues persist

2. **Faster Resume**:
   - KEEP_INVARIANT_ATCG step will skip processing if input hasn't changed
   - Cached outputs will be reused when resuming pipeline
   - Overall pipeline resume time should be significantly reduced

## Usage

To use the recombination-aware workflow with these fixes:

```bash
nextflow run main.nf \
  --input /path/to/assemblies \
  --outdir results \
  --recombination_aware_mode true \
  -resume
```

The `-resume` flag will now be much more effective, especially for the KEEP_INVARIANT_ATCG step.

## Testing Recommendations

1. Test the graft tree functionality with your existing data
2. Verify that resume works correctly by interrupting and restarting the pipeline
3. Check that the final grafted tree contains all expected samples
4. Monitor the grafting log files for any remaining issues

## Files Modified

- `modules/local/collect_representatives/main.nf`
- `modules/local/graft_subtrees/main.nf`
- `modules/local/keep_invariant_atcg/main.nf`
- `conf/params.config`