process MARKDUPLICATE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gatk4-spark=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4-spark:4.4.0.0--hdfd78af_0':
        'biocontainers/gatk4-spark:4.4.0.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.txt"), optional: true, emit: metrics
    tuple val(meta), path("*.bai"), optional: true, emit: bai
    tuple val(meta), path("*.sbi"), optional: true, emit: sbi
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gatk MarkDuplicatesSpark \
        $args \
        -I $bam \
        -O ${prefix}.bam \
        -M ${meta.id}_marked_dup_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
