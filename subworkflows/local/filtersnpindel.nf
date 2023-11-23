include { GATK4_SELECTVARIANTS as SELECT_SNP} from '../../modules/nf-core/gatk4/selectvariants/main.nf'
include { GATK4_SELECTVARIANTS as SELECT_INDEL } from '../../modules/nf-core/gatk4/selectvariants/main.nf'
include { GATK4_SELECTVARIANTS as DROPFILTEREDSNP } from '../../modules/nf-core/gatk4/selectvariants/main.nf'
include { GATK4_SELECTVARIANTS as DROPFILTEREDINDEL } from '../../modules/nf-core/gatk4/selectvariants/main.nf'
include { GATK4_VARIANTFILTRATION as FILTER_SNP } from '../../modules/nf-core/gatk4/variantfiltration/main.nf'
include { GATK4_VARIANTFILTRATION as FILTER_INDEL } from '../../modules/nf-core/gatk4/variantfiltration/main.nf'
include { VARIANTTOTABLE as SNPTOTABLE                       } from '../../modules/local/varianttotable.nf'
include { VARIANTTOTABLE as INDELTOTABLE         } from '../../modules/local/varianttotable.nf'

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
    SELECT_SNP(ch_vcf, [])
    SELECT_INDEL(ch_vcf, [])

    // filter the vcf by cretaria
    // SNP
    SELECT_SNP.out.vcf
    | join (SELECT_SNP.out.tbi)
    | map { meta, vcf, tbi ->
        meta = [id: meta.id, type: "snp"]
        [meta, vcf, tbi]
    }
    | set { ch_snp }
    // INDEL
    SELECT_INDEL.out.vcf
    | join (SELECT_INDEL.out.tbi)
    | map { meta, vcf, tbi->
        meta = [id: meta.id, type: "indel"]
        [meta, vcf, tbi]
    }
    | set { ch_indel }

    // hard filtering marking
    FILTER_SNP(ch_snp, ch_fasta.collect(), ch_fai.collect(), ch_dict.collect())
    FILTER_INDEL(ch_indel, ch_fasta.collect(), ch_fai.collect(), ch_dict.collect())

    // join tbi with vcf
    FILTER_SNP.out.vcf
    | join (FILTER_SNP.out.tbi)
    | set { ch_filter_snp }
    FILTER_INDEL.out.vcf
    | join (FILTER_INDEL.out.tbi)
    | set { ch_filter_indel }

    // generate stat table of snp and indel
    SNPTOTABLE(ch_filter_snp)
    INDELTOTABLE(ch_filter_indel)


    // drop entries with FILTERED mark
    DROPFILTEREDSNP(ch_filter_snp, [])
    | set { ch_filtered_snp }

    DROPFILTEREDINDEL(ch_filter_indel, [])
    | set { ch_filtered_indel }

    ch_versions = Channel.empty()



    emit:
    snp_vcf    = ch_filtered_snp.vcf
    snp_tbi    = ch_filtered_snp.tbi
    indel_vcf  = ch_filtered_indel.vcf
    indel_tbi  = ch_filtered_indel.tbi
    versions = ch_versions                     // channel: [ versions.yml ]
}

