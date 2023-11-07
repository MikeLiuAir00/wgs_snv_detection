process BQSR {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
        'biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(ref)
    tuple val(meta3), path(dbvcf) // only support vcf file path now.

    output:
    path path(*.table)            , emit: recal_tbl
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gatk BaseRecalibrator \
        $args \
        -I $bam \
        -R $ref \
        --known-sites $dbvcf \
        -O ${prefix}.recal_data.table
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    gatk: \$(gatk --version 2>1&)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.recal_data.table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    gatk: \$(gatk --version 2>1&)
    END_VERSIONS
    """
}
