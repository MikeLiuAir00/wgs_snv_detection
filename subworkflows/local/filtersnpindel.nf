include { GATK4_SELECTVARIANTS as SELECT_SNP} from '../../modules/nf-core/gatk4/selectvariants/main.nf'
include { GATK4_SELECTVARIANTS as SELECT_INDEL } from '../../modules/nf-core/gatk4/selectvariants/main.nf'
include { GATK4_SELECTVARIANTS as DROPFILTERED } '../../modules/nf-core/gatk4/selectvariants/main.nf'
include { GATK4_VARIANTFILTRATION as FILTER_SNP } from '../../modules/nf-core/gatk4/variantfiltration/main.nf'
include { GATK4_VARIANTFILTRATION as FILTER_INDEL } from '../../modules/nf-core/gatk4/variantfiltration/main.nf'

// Filter low qualit SNP and INDEL from raw calling result, and generate
// summary stat table for inspection and visualization
workflow FILTERSNPINDEL {

    take:
    ch_vcf      // vcf files
    ch_fasta    // reference genome
    ch_fai      // reference index
    ch_dict     // referecne dict

    main:
    // seperate snp and indel
    SELECT_SNP(ch_vcf, [], [])
    SELECT_INDEL(ch_vcf, [], [])

    // filter the vcf by cretaria
    // SNP
    SELECT_SNP.out.vcf
    | join (SELECT_SNP.out.tbi)
    | set { ch_snp }
    // INDEL
    SELECT_INDEL.out.vcf
    | join (SELECT_INDEL.out.tbi)
    | set { ch_indel }

    // apply hard filtering
    // FILTER_SNP(ch_snp, ch_fasta, ch_fai, ch_dict)
    // FILTER_INDEL(ch_indel, ch_fasta, ch_fai, ch_dict)

    // DROPFILTERED
    ch_versions = Channel.empty()



    emit:
    filtered_snp_vcf    = FILTER_SNP.out.vcf
    filtered_snp_tbi    = FILTER_SNP.out.tbi
    filtered_indel_vcf  = FILTER_INDEL.out.vcf
    filtered_indel_tbi  = FILTER_INDEL.out.tbi
    versions = ch_versions                     // channel: [ versions.yml ]
}

