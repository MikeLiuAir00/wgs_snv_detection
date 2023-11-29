include { GATK4_BASERECALIBRATOR as BQSR    } from '../../modules/nf-core/gatk4/baserecalibrator/main.nf'
include { GATK4_APPLYBQSR as APPLYBQSR           } from '../../modules/nf-core/gatk4/applybqsr/main.nf'

workflow BQSR_RUN {

    take:
    ch_bam
    ch_bai
    ch_fasta
    ch_fai
    ch_dict
    ch_vcf_db
    ch_vcf_db_tbi

    main:
    // data wrangle
    // ch_bam
    ch_bam
    | join (ch_bai)
    | map { meta, bam, bai ->
        [meta, bam, bai, []]
    }
    | set {ch_input}

    // ch_vcf_db
    BQSR(
        ch_input,
        ch_fasta.collect(),
        ch_fai.collect(),
        ch_dict.collect(),
        ch_vcf_db.collect(),
        ch_vcf_db_tbi.collect()
    )
    | set { ch_bqsr }

    ch_bam
    | join (ch_bai)
    | join (
        ch_bqsr.table
    )
    | map { meta, bam, bai, table ->
        [meta, bam, bai, table, []]
    }
    | set { ch_bam_in }

    // ch_bam_in.view()

    APPLYBQSR(
        ch_bam_in,
        ch_fasta.collect(),
        ch_fai.collect(),
        ch_dict.collect()
    )
    | set { ch_calibrated }
    ch_versions = Channel.empty()


    emit:
    calibrated      = ch_calibrated.bam

    versions = ch_versions                     // channel: [ versions.yml ]
}

