process GUBBINS_CLUSTER {
    tag "cluster_${cluster_id}"
    label 'process_high'
    container "snads/gubbins@sha256:391a980312096f96d976f4be668d4dea7dda13115db004a50e49762accc0ec62"

    input:
    tuple val(cluster_id), path(alignment), path(starting_tree)

    output:
    tuple val(cluster_id), path("${cluster_id}.filtered_polymorphic_sites.fasta"), emit: filtered_alignment
    tuple val(cluster_id), path("${cluster_id}.recombination_predictions.gff"), emit: recombination_gff
    tuple val(cluster_id), path("${cluster_id}.node_labelled.final_tree.tre"), emit: final_tree
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def iterations = params.gubbins_iterations ?: 3
    def tree_builder = params.gubbins_tree_builder ?: 'hybrid'
    def min_snps = params.gubbins_min_snps ?: 5
    """
    # Check if alignment has at least 3 sequences (minimum for phylogenetic analysis)
    seq_count=\$(grep -c "^>" $alignment)
    
    if [ \$seq_count -lt 3 ]; then
        echo "WARNING: Alignment has only \$seq_count sequences. Gubbins requires at least 3 sequences for phylogenetic analysis."
        echo "Skipping Gubbins analysis for cluster ${cluster_id}"
        
        # Create empty output files to satisfy pipeline expectations
        touch ${cluster_id}.filtered_polymorphic_sites.fasta
        touch ${cluster_id}.recombination_predictions.gff
        touch ${cluster_id}.node_labelled.final_tree.tre
        
        # Create versions file for small clusters
        echo '"${task.process}":' > versions.yml
        echo '    gubbins: '\$(run_gubbins.py --version | sed 's/^/    /') >> versions.yml
        
        exit 0
    fi

    # Run Gubbins with optimized settings for cluster
    run_gubbins.py \\
        --starting-tree $starting_tree \\
        --prefix ${cluster_id} \\
        --tree-builder $tree_builder \\
        --iterations $iterations \\
        --min-snps $min_snps \\
        --threads ${task.cpus} \\
        $args \\
        $alignment

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gubbins: \$(run_gubbins.py --version | sed 's/^/    /')
    END_VERSIONS
    """
}