#!/bin/bash -euo pipefail
# Check if alignment has at least 3 sequences (minimum for phylogenetic analysis)
seq_count=$(grep -c "^>" cluster_13.aln.fa)

if [ $seq_count -lt 3 ]; then
    echo "WARNING: Alignment has only $seq_count sequences. Gubbins requires at least 3 sequences for phylogenetic analysis."
    echo "Skipping Gubbins analysis for cluster cluster_13"

    # Create empty output files to satisfy pipeline expectations
    touch cluster_13.filtered_polymorphic_sites.fasta
    touch cluster_13.recombination_predictions.gff
    touch cluster_13.node_labelled.final_tree.tre

    # Create versions file for small clusters
    echo '"ASSEMBLY_SNPS_SCALABLE:CLUSTERED_SNP_TREE:GUBBINS_CLUSTER":' > versions.yml
    echo '    gubbins: '$(run_gubbins.py --version | sed 's/^/    /') >> versions.yml

    exit 0
fi

# Run Gubbins with optimized settings for cluster
run_gubbins.py \
    --starting-tree cluster_13.treefile \
    --prefix cluster_13 \
    --tree-builder hybrid \
    --iterations 3 \
    --min-snps 5 \
    --threads 16 \
     \
    cluster_13.aln.fa

cat <<-END_VERSIONS > versions.yml
"ASSEMBLY_SNPS_SCALABLE:CLUSTERED_SNP_TREE:GUBBINS_CLUSTER":
    gubbins: $(run_gubbins.py --version | sed 's/^/    /')
END_VERSIONS
