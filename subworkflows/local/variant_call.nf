include { HAPLOTYPECALLER             } from '../../modules/local/haplotypecaller.nf'
include { GATK4_GENOMICSDBIMPORT      } from '../../modules/nf-core/gatk4/genomicsdbimport/main.nf'
include { GATK4_GENOTYPEGVCFS         } from '../../modules/nf-core/gatk4/genotypegvcfs/main.nf'
workflow VARIANT_CALL {

    take:
    ch_bam
    ch_bai
    ch_ref
    ch_fai
    ch_dict

    main:
    HAPLOTYPECALLER(
        ch_bam,
        ch_bai,
        ch_ref.collect(),
        ch_fai.collect(),
        ch_dict.collect(),
        [],
        [],
        [],
        true
    )
    | set { ch_callset }

    // wrangle data
    ch_callset.vcf
    | map {meta, vcf ->
            meta = [id:meta.id]
            [meta, vcf]
        }
    | collectFile(sort: true) { meta, vcf ->
        ["vcf.sample_map", meta.id + "\t" + vcf + '\n' ]
    }
    | set { vcf_in }

    vcf_in
    | map { f_vcf ->
        meta = [id: "sample list"]
        [meta, f_vcf]
    }
    | set { vcf_in }

    Channel.fromList((1..8).collect {"chr${it}"})
    | set { ch_interval }
    // create db on all chrs
    GATK4_GENOMICSDBIMPORT(vcf_in.collect(), ch_interval, [], [], false, false, true)
    | set { ch_merged }

    ch_merged.genomicsdb
    | map { meta, db ->
        [meta, db, [], [], []]
    }
    | set { ch_genomicsdb }
    GATK4_GENOTYPEGVCFS(
        ch_genomicsdb,
        ch_ref.collect(),
        ch_fai.collect(),
        ch_dict.collect(),
        [],
        []
    )
    | set { ch_genotype }
    ch_versions = Channel.empty()


    emit:
    vcf     = ch_genotype.vcf
    tbi     = ch_genotype.tbi

    versions = ch_versions                     // channel: [ versions.yml ]
}

