# Testing the Clustering Fix

## Problem
The workflow was completing successfully but not running the downstream phylogenetic processes (SKA_BUILD, SKA_ALIGN, IQTREE_FAST, GUBBINS_CLUSTER).

## Root Cause
The issue was that the clustering algorithm was producing only singleton clusters (each sample in its own cluster), and there was no handling for this case. Since phylogenetic analysis requires at least 3 samples, these singleton clusters were being passed to the CLUSTERED_SNP_TREE subworkflow but no processes were running.

## Fixes Applied

### 1. Fixed Channel Operations in CLUSTERING Subworkflow
- **File**: `subworkflows/local/clustering.nf`
- **Issue**: The `combine` operation was using incorrect indices for joining channels
- **Fix**: Changed from `combine(ch_assemblies, by: 1)` to proper `join` operation:
  ```groovy
  .map { cluster_id, sample_id -> tuple(sample_id, cluster_id) }
  .join(ch_assemblies)
  .map { sample_id, cluster_id, assembly -> tuple(cluster_id, sample_id, assembly) }
  ```

### 2. Added Singleton Cluster Filtering
- **File**: `subworkflows/local/clustering.nf`
- **Issue**: Singleton clusters (1 sample) cannot be used for phylogenetic analysis
- **Fix**: Added filtering to only process clusters with >1 sample:
  ```groovy
  .branch { cluster_id, sample_ids, assemblies ->
      multi_sample: sample_ids.size() > 1
      singleton: sample_ids.size() == 1
  }
  ```

### 3. Added New Parameter: `--merge_singletons`
- **File**: `conf/params.config`
- **Purpose**: Allow users to merge all singleton clusters into one large cluster
- **Default**: `false`
- **Usage**: `--merge_singletons` to enable

### 4. Added Informative Logging
- **Purpose**: Help users understand why processes aren't running
- **Features**:
  - Logs when singleton clusters are skipped
  - Provides suggestions for parameter adjustment
  - Shows cluster count summary

### 5. Updated Documentation
- **File**: `README.md`
- **Added**: Section explaining clustering parameters and troubleshooting

## Testing the Fix

To test if the fix works:

1. **For datasets with mostly singletons**:
   ```bash
   nextflow run . --scalable_mode --mash_threshold 0.05 --merge_singletons
   ```

2. **For datasets that should cluster better**:
   ```bash
   nextflow run . --scalable_mode --mash_threshold 0.1
   ```

3. **Check logs for clustering information**:
   - Look for "Found X clusters for phylogenetic analysis"
   - If 0 clusters, check suggestions in log messages

## Expected Behavior After Fix

1. **If multi-sample clusters exist**: Downstream processes (SKA_BUILD, etc.) will run
2. **If only singletons exist**: 
   - Without `--merge_singletons`: Clear warning messages, no downstream processes
   - With `--merge_singletons`: One large cluster created, downstream processes run
3. **Channel operations**: Should work correctly without UnboundLocalError

## Parameters for Troubleshooting

- `--mash_threshold`: Lower = more strict clustering, higher = more permissive
- `--merge_singletons`: Combine all singletons into one cluster
- `--max_cluster_size`: Split large clusters (shouldn't affect the original issue)