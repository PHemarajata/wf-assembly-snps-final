#!/bin/bash -euo pipefail
# Check if alignment has at least 3 sequences (minimum for tree building)
seq_count=$(grep -c "^>" cluster_3.aln.fa)

if [ $seq_count -lt 3 ]; then
    echo "WARNING: Alignment has only $seq_count sequences. IQ-TREE requires at least 3 sequences for tree building."
    echo "Skipping tree construction for cluster cluster_3"

    # Create empty output files to satisfy pipeline expectations
    touch cluster_3.treefile
    touch cluster_3.iqtree

    # Create versions file for small clusters
    echo '"ASSEMBLY_SNPS_SCALABLE:CLUSTERED_SNP_TREE:IQTREE_FAST":' > versions.yml
    echo '    iqtree: '$(iqtree2 --version 2>&1 | head -n1 | sed 's/^/    /') >> versions.yml

    exit 0
fi

# Run IQ-TREE with fast mode
iqtree2 \
    -s cluster_3.aln.fa \
    -st DNA \
    -m GTR+ASC \
    --fast \
    -nt AUTO \
    --prefix cluster_3 \


cat <<-END_VERSIONS > versions.yml
"ASSEMBLY_SNPS_SCALABLE:CLUSTERED_SNP_TREE:IQTREE_FAST":
    iqtree: $(iqtree2 --version 2>&1 | head -n1 | sed 's/^/    /')
END_VERSIONS
