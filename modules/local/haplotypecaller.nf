process HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_low'
    label 'process_long'

    conda "bioconda::gatk4=4.2.6.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.2.6.1--py36hdfd78af_1' }"

    input:
    tuple val(meta) ,   path(bam)
    tuple val(meta1),   path(bai)
    tuple val(meta2),   path(ref)
    tuple val(meta3),   path(fai)
    tuple val(meta4),   path(dict)
    path(dbsnp)
    path(dbsnp_tbi)
    path(intervals)
    val(ifgvcf)


    output:
    tuple val(meta), path("*.vcf.gz")   , optional: true, emit: vcf
    //tuple val(meta), path("*.g.vcf.gz") , optional: true, emit: gvcf
    tuple val(meta), path("*.tbi")      , optional:true, emit: tbi
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = ifgvcf ? '.g.vcf.gz': '.vcf.gz'
    def gvcf = ifgvcf ? '-ERC GVCF': ''
    def dbsnp_command = dbsnp ? "--dbsnp $dbsnp" : ""
    def interval_command = intervals ? "--intervals $intervals" : ""
    """
    gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=2" \
        HaplotypeCaller  \
            -R $ref \
            -I $bam \
            -O ${prefix}${suffix} \
            --tmp-dir . \
            $gvcf \
            $dbsnp_command \
            $interval_command \
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
    touch ${prefix}${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>1&)
    END_VERSIONS
    """
}
