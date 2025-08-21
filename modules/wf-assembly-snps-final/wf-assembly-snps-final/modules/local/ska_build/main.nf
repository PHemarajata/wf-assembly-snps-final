process SKA_BUILD {
    tag "cluster_${cluster_id}"
    label 'process_medium'
    container "quay.io/biocontainers/ska2:0.3.7--h7d875b9_0"

    input:
    tuple val(cluster_id), val(sample_list)

    output:
    tuple val(cluster_id), path("${cluster_id}.skf"), emit: ska_file
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Create input file for SKA
    cat > ${cluster_id}_input.tsv <<EOF
${sample_list.collect{ "${it[0]}\t${it[1]}" }.join('\n')}
EOF

    # Build SKA file
    ska build \\
        -o ${cluster_id} \\
        $args \\
        ${cluster_id}_input.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska: \$(ska --version 2>&1 | head -n1 | sed 's/^/    /')
    END_VERSIONS
    """
}