# SKA_BUILD Process Fix

## Problem
The SKA_BUILD process was failing with a Rust panic error because it couldn't find the input FASTA files. The error occurred in work directory `b96c91c2e9bc908d939f6d212ba795`.

## Root Cause Analysis

### What was happening:
1. Nextflow stages input files using symbolic links in the work directory
2. The `path(assemblies)` parameter creates these staged files with their original names
3. The SKA_BUILD script was creating a TSV input file using `${assembly}` which gave the full path object
4. SKA was trying to read files by their original names, but the staged files might have different names or broken symlinks

### Evidence from the failed work directory:
- **Input TSV content**: Referenced files like `ERS013364.fasta`, `IP-0087-7_S10_L001-SPAdes.fasta`
- **Staged files**: Symbolic links pointing to absolute paths that don't exist in the current environment
- **SKA error**: Rust panic indicating file access issues

## Solution

### The Fix:
Changed the SKA_BUILD process to use `${assembly.name}` instead of `${assembly}` when creating the input TSV file.

**Before:**
```groovy
def input_content = [sample_ids, assemblies].transpose().collect{ sample_id, assembly -> 
    "${sample_id}\t${assembly}" 
}.join('\n')
```

**After:**
```groovy
def input_content = [sample_ids, assemblies].transpose().collect{ sample_id, assembly -> 
    "${sample_id}\t${assembly.name}" 
}.join('\n')
```

### Additional improvements:
1. **Added file existence verification**: The script now checks which FASTA files are actually present in the work directory
2. **Better debugging**: Added logging to help troubleshoot file staging issues
3. **Updated both copies**: Fixed both the main module and the duplicate in the modules directory

## Files Modified:
- `modules/local/ska_build/main.nf`
- `modules/modules/local/ska_build/main.nf`

## Expected Result:
The SKA_BUILD process should now correctly reference the staged files and SKA should be able to access them without errors.

## Testing:
To test the fix, run the pipeline with scalable mode:
```bash
nextflow run main.nf -profile docker --scalable_mode true --input <your_input> --outdir results
```

The process should now complete successfully without the Rust panic error.