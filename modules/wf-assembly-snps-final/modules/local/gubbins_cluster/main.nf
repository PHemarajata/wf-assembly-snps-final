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