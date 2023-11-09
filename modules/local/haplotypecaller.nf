process HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
        'biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta),    path(bam)
    tuple val(meta1),   path(bai)
    tuple val(meta2),    path(ref)
    val(ifgvcf)

    output:
    tuple val(meta), path("*.vcf.gz"), optional: true, emit: vcf
    tuple val(meta), path("*.g.vcf.gz"), optional: true, emit: gvcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = ifgvcf ? '.g.vcf.gz': '.vcf.gz'
    def gvcf = ifgvcf ? '-ERC GVCF': ''
    """
    gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=2" \
        HaplotypeCaller  \
            -R $ref \
            -I $bam \
            -O ${prefix}${suffix} \
            $gvcf \
            $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>1&)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>1&)
    END_VERSIONS
    """
}
