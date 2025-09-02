# Numpy/Pandas Compatibility Fix

## Problem
```
ImportError: this version of pandas is incompatible with numpy < 1.15.4
your numpy version is 1.12.1.
Please upgrade numpy to >= 1.15.4 to use this pandas version
```

## Root Cause
The containers were using an old version of numpy (1.12.1) that's incompatible with newer pandas versions. When pandas was installed via pip, it installed a recent version that requires numpy >= 1.15.4, but the container had an older numpy version.

## Solution Applied

### Strategy: Explicit Numpy Upgrade
For all modules that use pandas, I added explicit numpy upgrade before installing pandas:

```bash
# Install compatible versions of numpy and pandas
pip install --upgrade numpy>=1.15.4
pip install pandas biopython
```

## Modules Fixed

### 1. EXTRACT_CORE_SNPS Module ✅
**Container**: `python:3.9--1`
**Fix**: Added `pip install --upgrade numpy>=1.15.4` before pandas installation
**Result**: Compatible numpy and pandas versions

### 2. INTEGRATE_CORE_SNPS Module ✅
**Container**: `python:3.9--1`
**Fix**: Added `pip install --upgrade numpy>=1.15.4` before pandas installation
**Result**: Compatible numpy and pandas versions

### 3. GRAFT_TREES Module ✅
**Container**: `ete3:3.1.2--py39hd23ed53_0` (kept for ETE3 dependency)
**Fix**: Added `pip install --upgrade numpy>=1.15.4` before pandas installation
**Result**: Compatible numpy and pandas versions with fallback logic

### 4. BUILD_INTEGRATED_TREE Module ✅
**Container**: `iqtree:2.2.6--h21ec9f0_0` (kept for IQ-TREE dependency)
**Fix**: Added `pip install --upgrade numpy>=1.15.4` before pandas installation
**Result**: Compatible numpy and pandas versions with fallback logic

### 5. CREATE_FINAL_SUMMARY Module ✅
**Container**: `python:3.9--1`
**Fix**: Added `pip install --upgrade numpy>=1.15.4` before pandas installation
**Result**: Compatible numpy and pandas versions

## Key Changes Made

### Before (Problematic)
```bash
pip install pandas  # This installs latest pandas but numpy might be old
```

### After (Fixed)
```bash
# Install compatible versions of numpy and pandas
pip install --upgrade numpy>=1.15.4
pip install pandas biopython
```

### Version Tracking
Updated version reporting to include both numpy and pandas:
```bash
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed 's/Python //')
    pandas: \$(python -c "import pandas; print(pandas.__version__)")
    numpy: \$(python -c "import numpy; print(numpy.__version__)")
END_VERSIONS
```

## Expected Results

1. **No more numpy/pandas compatibility errors**
2. **All integration modules should complete successfully**
3. **Proper version compatibility maintained**
4. **Fallback logic for specialized containers**

## Files Modified

1. `modules/local/extract_core_snps/main.nf` - Added numpy upgrade
2. `modules/local/integrate_core_snps/main.nf` - Added numpy upgrade
3. `modules/local/graft_trees/main.nf` - Added numpy upgrade with fallback
4. `modules/local/build_integrated_tree/main.nf` - Added numpy upgrade with fallback
5. `modules/local/create_final_summary/main.nf` - Added numpy upgrade

## Testing

To verify the fix:
1. The `ImportError: this version of pandas is incompatible with numpy < 1.15.4` should no longer occur
2. All pandas-dependent modules should complete successfully
3. Version output should show compatible numpy (>=1.15.4) and pandas versions

## Why This Approach Works

- **Explicit Control**: We explicitly upgrade numpy before installing pandas
- **Version Constraint**: We ensure numpy >= 1.15.4 which is compatible with modern pandas
- **Container Compatibility**: Works with both general Python containers and specialized tool containers
- **Fallback Logic**: Specialized containers have fallback logic if package installation fails

The pipeline should now run without numpy/pandas compatibility errors while maintaining all functionality.