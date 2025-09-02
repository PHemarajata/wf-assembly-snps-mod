# Fix Summary: SKA_BUILD Input Tuple Mismatch Error

## Problem
The pipeline was failing with the error:
```
Input tuple does not match tuple declaration in process 'SKA_BUILD' -- offending value: [cluster_7, [[sample1, path1], [sample2, path2], ...]]
Path value cannot be null
```

## Root Cause
The clustering subworkflow was outputting data in the wrong format for the SKA_BUILD process:

**Expected by SKA_BUILD:**
```
tuple val(cluster_id), val(sample_ids), path(assemblies)
```

**Actually provided by clustering:**
```
tuple(cluster_id, [[sample_id, assembly], ...])
```

## Solution
Fixed the data transformation in the clustering subworkflow to output the correct tuple format:

### Changes Made:

1. **Fixed tuple format in clustering.nf:**
   - Changed from: `tuple(cluster_id, sample_assembly_pairs)`
   - Changed to: `tuple(cluster_id, sample_ids, assemblies)`

2. **Updated singleton handling:**
   - Fixed merge_singletons logic to properly separate sample_ids and assemblies
   - Ensured consistent tuple format across all code paths

3. **Improved join operation:**
   - Added explicit `by: 0` parameter to join operation for clarity

4. **Updated comments:**
   - Fixed channel format documentation throughout the codebase

### Files Modified:
- `subworkflows/local/clustering.nf`
- `subworkflows/local/clustered_snp_tree.nf`
- All duplicate files in nested modules directories

## Expected Result
The SKA_BUILD process should now receive properly formatted input tuples and the pipeline should run successfully in scalable mode.

## Testing
To test the fix, run the pipeline with scalable mode enabled:
```bash
nextflow run main.nf -profile docker --scalable_mode true --input <your_input> --outdir results
```