process MARKDUPLICATE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::picard=3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.0.0--hdfd78af_1':
        'biocontainers/picard:3.0.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path(".txt") , emit: metrics
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    picard MarkDuplicates \
        $args \
        I=$bam \
        O=${prefix}.bam \
        M=marked_dup_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    gatk: \$(gatk --version 2>1&)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    gatk: \$(gatk --version 2>1&)
    END_VERSIONS
    """
}
