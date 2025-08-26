process IQTREE_FAST {
    tag "cluster_${cluster_id}"
    label 'process_medium'
    container "quay.io/biocontainers/iqtree:2.2.6--h21ec9f0_0"

    input:
    tuple val(cluster_id), path(alignment)

    output:
    tuple val(cluster_id), path("${cluster_id}.treefile"), emit: tree
    tuple val(cluster_id), path("${cluster_id}.iqtree"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def model = params.iqtree_model ?: 'GTR+ASC'
    """
    # Run IQ-TREE with fast mode
    iqtree2 \\
        -s $alignment \\
        -st DNA \\
        -m $model \\
        --fast \\
        -nt AUTO \\
        --prefix ${cluster_id} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(iqtree2 --version 2>&1 | head -n1 | sed 's/^/    /')
    END_VERSIONS
    """
}