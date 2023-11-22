process MKDBIMPORTLIST {
    tag "$meta.id"
    label 'process_single'
    debug true

    input:
    tuple val(meta), val(tbi), val(vcf)

    output:
    tuple val(meta), path("*.tsv"), optional: true, emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/bash
    printf "\$(basename $vcf)\\t$vcf" >> vcf_list.tsv
    printf "\$(basename $tbi)\\t$tbi" >> tbi_list.tsv
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
