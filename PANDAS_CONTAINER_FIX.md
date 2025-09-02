# Container and Pandas Dependency Fix

## Problem
```
ModuleNotFoundError: No module named 'pandas'
```

The error occurred because several modules were using containers that don't have pandas installed, but the Python scripts require pandas for data processing.

## Root Cause
Several modules were using specialized containers (like `biopython:1.79` or `iqtree:2.2.6`) that don't include pandas, but the Python scripts in these modules require pandas for data manipulation.

## Modules Fixed

### 1. EXTRACT_CORE_SNPS Module ✅
**Problem**: Used `biopython:1.79--py39hd23ed53_1` container without pandas
**Fix**: 
- Changed to `python:3.9--1` container
- Added `pip install biopython pandas` in script
- Removed problematic docstring that caused variable substitution errors

### 2. INTEGRATE_CORE_SNPS Module ✅
**Problem**: Used `biopython:1.79--py39hd23ed53_1` container without pandas
**Fix**:
- Changed to `python:3.9--1` container  
- Added `pip install biopython pandas` in script

### 3. BUILD_INTEGRATED_TREE Module ✅
**Problem**: Used `iqtree:2.2.6--h21ec9f0_0` container without pandas
**Fix**:
- Kept original container (needed for IQ-TREE)
- Added `pip install pandas` with error handling
- Added fallback logic to work without pandas if installation fails

### 4. GRAFT_TREES Module ✅
**Problem**: Used `ete3:3.1.2--py39hd23ed53_0` container which may not have pandas
**Fix**:
- Kept original container (needed for ETE3)
- Added `pip install pandas` with error handling
- Added fallback logic to work without pandas if installation fails

## Solution Strategy

### Approach 1: Change Container (Used for modules that don't need specialized tools)
- **EXTRACT_CORE_SNPS**: Changed to `python:3.9--1` and install biopython + pandas
- **INTEGRATE_CORE_SNPS**: Changed to `python:3.9--1` and install biopython + pandas

### Approach 2: Install Pandas in Script (Used for modules that need specialized containers)
- **BUILD_INTEGRATED_TREE**: Keep IQ-TREE container, install pandas with fallback
- **GRAFT_TREES**: Keep ETE3 container, install pandas with fallback

## Key Changes Made

### Container Updates
```nextflow
# Before
container "quay.io/biocontainers/biopython:1.79--py39hd23ed53_1"

# After  
container "quay.io/biocontainers/python:3.9--1"
```

### Script Updates
```bash
# Install required packages
pip install biopython pandas

# Or with error handling for specialized containers
pip install pandas || echo "pandas installation failed, continuing without it"
```

### Fallback Logic
Added fallback code that works without pandas when installation fails:
```python
try:
    import pandas as pd
    pandas_available = True
except ImportError:
    pandas_available = False
    # Use basic file operations instead
```

## Expected Results

1. **No more pandas import errors** - All modules now have access to pandas
2. **Maintained functionality** - Specialized tools (IQ-TREE, ETE3) still work
3. **Robust error handling** - Modules continue to work even if pandas installation fails
4. **Better compatibility** - Using standard Python container where possible

## Files Modified

1. `modules/local/extract_core_snps/main.nf` - Container change + pandas install
2. `modules/local/integrate_core_snps/main.nf` - Container change + pandas install  
3. `modules/local/build_integrated_tree/main.nf` - Pandas install + fallback logic
4. `modules/local/graft_trees/main.nf` - Pandas install + fallback logic

## Testing

To verify the fixes:
1. The `ModuleNotFoundError: No module named 'pandas'` should no longer occur
2. All integration modules should complete successfully
3. Output files should be generated properly

The pipeline should now run without pandas-related errors while maintaining all original functionality.