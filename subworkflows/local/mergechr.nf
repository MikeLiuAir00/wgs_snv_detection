include { GATK4_MERGEVCFS as MERGESNP   } from '../../modules/nf-core/gatk4/mergevcfs/main.nf'
include { GATK4_MERGEVCFS as MERGEINDEL   } from '../../modules/nf-core/gatk4/mergevcfs/main.nf'

workflow MERGECHR {

    take:
    ch_snp
    ch_indel
    ch_dict

    main:
    MERGESNP(ch_snp, ch_dict)
    | set { ch_merged_snp }

    MERGEINDEL(ch_indel, ch_dict)
    | set { ch_merged_indel }

    ch_versions = Channel.empty()

    emit:
    snp_vcf         = ch_merged_snp.vcf
    snp_tbi         = ch_merged_snp.tbi
    indel_vcf       = ch_merged_indel.vcf
    indel_tbi       = ch_merged_indel.tbi

    versions = ch_versions                     // channel: [ versions.yml ]
}

