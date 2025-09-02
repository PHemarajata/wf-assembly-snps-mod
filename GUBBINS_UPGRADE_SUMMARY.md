# Gubbins Upgrade and Hybrid Tree Builder Implementation

## Summary of Changes

This document summarizes the changes made to upgrade Gubbins and implement proper hybrid tree building functionality.

## Container Update

**From:** `snads/gubbins@sha256:391a980312096f96d976f4be668d4dea7dda13115db004a50e49762accc0ec62` (Gubbins 3.1.4)
**To:** `quay.io/biocontainers/gubbins:3.3.5--py39pl5321he4a0461_0` (Gubbins 3.3.5)

## New Parameters Added

### conf/params.config
```groovy
// Gubbins optimization parameters
run_gubbins                       = true
gubbins_iterations                = 3
gubbins_use_hybrid                = true
gubbins_first_tree_builder        = "rapidnj"
gubbins_tree_builder              = "iqtree"
gubbins_min_snps                  = 5
```

## Hybrid Tree Builder Implementation

The hybrid approach now properly uses two tree builders as intended by Gubbins:

1. **First Tree Builder (rapidnj)**: Creates a fast neighbor-joining tree for initial topology
2. **Tree Builder (iqtree)**: Performs maximum likelihood optimization for final accuracy

### Command Structure
```bash
# When gubbins_use_hybrid = true
run_gubbins.py \
    --starting-tree input.tree \
    --prefix output_prefix \
    --first-tree-builder rapidnj \
    --tree-builder iqtree \
    --iterations 3 \
    --min-snps 5 \
    --threads 16 \
    alignment.fasta

# When gubbins_use_hybrid = false
run_gubbins.py \
    --starting-tree input.tree \
    --prefix output_prefix \
    --tree-builder iqtree \
    --iterations 3 \
    --min-snps 5 \
    --threads 16 \
    alignment.fasta
```

## Files Modified

1. **modules/local/gubbins_cluster/main.nf**
   - Updated container to Gubbins 3.3.5
   - Implemented hybrid tree builder logic
   - Added conditional command building

2. **modules/local/recombination_gubbins/main.nf**
   - Updated container to Gubbins 3.3.5
   - Implemented hybrid tree builder logic in shell script

3. **conf/params.config**
   - Added new hybrid tree builder parameters
   - Updated schema_ignore_params list

4. **docs/scalable-mode.md**
   - Updated parameter documentation
   - Updated optimization strategies description

5. **FIXES_APPLIED.md**
   - Added comprehensive documentation of the upgrade

## Benefits

- **Resolves KeyError**: Eliminates the `KeyError: 'hybrid'` error from Gubbins 3.1.4
- **Improved Accuracy**: IQ-TREE provides more accurate maximum likelihood inference
- **Speed Optimization**: RapidNJ provides fast initial tree, IQ-TREE refines for accuracy
- **Flexibility**: Users can disable hybrid mode or change tree builders
- **Modern Features**: Access to latest Gubbins features and bug fixes

## Usage Examples

### Enable Hybrid Mode (Default)
```bash
nextflow run main.nf \
  --scalable_mode true \
  --run_gubbins true \
  --gubbins_use_hybrid true \
  --gubbins_first_tree_builder rapidnj \
  --gubbins_tree_builder iqtree
```

### Use Single Tree Builder
```bash
nextflow run main.nf \
  --scalable_mode true \
  --run_gubbins true \
  --gubbins_use_hybrid false \
  --gubbins_tree_builder iqtree
```

### Alternative Tree Builders
```bash
# Use FastTree for speed
nextflow run main.nf \
  --scalable_mode true \
  --run_gubbins true \
  --gubbins_use_hybrid false \
  --gubbins_tree_builder fasttree

# Use RAxML for accuracy
nextflow run main.nf \
  --scalable_mode true \
  --run_gubbins true \
  --gubbins_use_hybrid false \
  --gubbins_tree_builder raxml
```

## Testing Recommendations

1. Test with small dataset to verify container pulls correctly
2. Verify hybrid mode produces expected output files
3. Compare results between hybrid and single tree builder modes
4. Test with different tree builder combinations
5. Validate performance improvements with larger datasets

## Troubleshooting

If you encounter issues:

1. **Container pull fails**: Verify internet connection and Docker/Singularity setup
2. **Tree builder not found**: Check Gubbins 3.3.5 supports your chosen tree builder
3. **Memory issues**: Reduce cluster sizes or disable hybrid mode for resource-constrained systems
4. **Performance issues**: Consider using FastTree instead of IQ-TREE for very large datasets