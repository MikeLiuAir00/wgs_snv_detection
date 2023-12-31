process ADDREADGROUPS {
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
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    picard AddOrReplaceReadGroups \
        I=$bam \
        O=${prefix}.bam \
        RGSM=$meta.id \
        RGID=$meta.rgid \
        RGLB=$meta.rglb \
        RGPU=$meta.rgpu \
        RGPL=$meta.rgpl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard AddOrReplaceReadGroups --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
    // picard: \$(picard AddOrReplaceReadGroups --version 2>1&)

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard AddOrReplaceReadGroups --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
    // not output in stdout picard: \$(picard AddOrReplaceReadGroups --version 2>1&)
}
