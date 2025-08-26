process SKA_BUILD {
    tag "cluster_${cluster_id}"
    label 'process_medium'
    container "quay.io/biocontainers/ska2:0.3.7--h4349ce8_2"

    input:
    tuple val(cluster_id), val(sample_list)

    output:
    tuple val(cluster_id), path("${cluster_id}.skf"), emit: ska_file
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_content = sample_list.collect{ "${it[0]}\t${it[1]}" }.join('\n')
    """
    # Create input file for SKA
    cat > ${cluster_id}_input.tsv <<'EOFSKA'
${input_content}
EOFSKA

    # Verify input file was created correctly
    echo "Input file contents:"
    cat ${cluster_id}_input.tsv
    
    # Build SKA file
    ska build \\
        -o ${cluster_id} \\
        ${args} \\
        ${cluster_id}_input.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska: \$(ska --version 2>&1 | head -n1 | sed 's/^/    /')
    END_VERSIONS
    """
}