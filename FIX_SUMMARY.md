# SKA_BUILD Channel Format Error Fix

## Problem
The pipeline was failing with the error:
```
Input tuple does not match tuple declaration in process `ASSEMBLY_SNPS_SCALABLE:CLUSTERED_SNP_TREE:SKA_BUILD` -- offending value: [cluster_3, [[sample_id, path], [sample_id, path], ...]]
```

And "Path value cannot be null" errors.

## Root Cause
The issue was a mismatch between the channel format expected by SKA_BUILD and what was being provided:

**SKA_BUILD expected:**
```nextflow
tuple val(cluster_id), val(sample_ids), path(assemblies)
```

**But was receiving:**
```nextflow
[cluster_id, [[sample_id, assembly_path], [sample_id, assembly_path], ...]]
```

The clustering workflow was creating nested arrays instead of separate lists of sample_ids and assemblies.

## Solution

### 1. Fixed clustering.nf channel transformation
**Before:**
```nextflow
ch_clustered_assemblies = ch_clustered_final
    .map { cluster_id, sample_ids, assemblies ->
        def sample_assembly_pairs = [sample_ids, assemblies].transpose()
        tuple(cluster_id, sample_assembly_pairs)
    }
```

**After:**
```nextflow
ch_clustered_assemblies = ch_clustered_final
    .map { cluster_id, sample_ids, assemblies ->
        tuple(cluster_id, sample_ids, assemblies)
    }
```

### 2. Fixed singleton merging logic
Updated the singleton merging to maintain the correct format with separate sample_ids and assemblies lists.

### 3. Simplified clustered_snp_tree.nf
Removed the unnecessary channel transformation since the input is now already in the correct format.

## Files Modified
1. `modules/subworkflows/local/clustering.nf`
2. `modules/subworkflows/local/clustered_snp_tree.nf`

## Result
The channel now passes data in the format that SKA_BUILD expects:
- `cluster_id` as a value
- `sample_ids` as a list of sample identifiers  
- `assemblies` as a list of file paths

This should resolve both the "Input tuple does not match" error and the "Path value cannot be null" errors.