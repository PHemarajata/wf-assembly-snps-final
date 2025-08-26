# SKA_BUILD Error Fix Summary

## Problem
The SKA_BUILD process was failing with error 101 because it was trying to access assembly files using absolute paths from other Nextflow work directories. These paths were not accessible within the container running the SKA_BUILD process.

## Root Cause
The issue was in how the clustering workflow was passing data to SKA_BUILD:
1. Assembly file paths were being passed as strings instead of actual file objects
2. SKA_BUILD was trying to access files using absolute paths that don't exist in its container
3. The input format was `tuple val(cluster_id), val(sample_list)` where sample_list contained string paths

## Solution
### 1. Updated SKA_BUILD module input format
Changed from:
```nextflow
input:
tuple val(cluster_id), val(sample_list)
```

To:
```nextflow
input:
tuple val(cluster_id), val(sample_ids), path(assemblies)
```

### 2. Updated input file creation
Changed the script to use local file names instead of absolute paths:
```nextflow
def input_content = [sample_ids, assemblies].transpose().collect{ sample_id, assembly -> 
    "${sample_id}\t${assembly}" 
}.join('\n')
```

### 3. Updated clustered_snp_tree workflow
Added transformation to convert the input channel format:
```nextflow
ch_ska_input = ch_clustered_assemblies
    .map { cluster_id, sample_assembly_pairs ->
        def sample_ids = sample_assembly_pairs.collect { it[0] }
        def assemblies = sample_assembly_pairs.collect { it[1] }
        tuple(cluster_id, sample_ids, assemblies)
    }
```

### 4. Fixed clustering workflow
- Added explicit `by: 0` parameter to join operation
- Fixed singleton merging to preserve file objects

## Files Modified
1. `modules/local/ska_build/main.nf`
2. `modules/modules/local/ska_build/main.nf` (duplicate)
3. `modules/subworkflows/local/clustered_snp_tree.nf`
4. `modules/subworkflows/local/clustering.nf`

## Result
Now SKA_BUILD will receive actual assembly files as input parameters, allowing it to access them locally within its container, rather than trying to access files via absolute paths from other work directories.