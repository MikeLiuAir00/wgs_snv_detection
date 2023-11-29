/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowSnvdetect.initialise(params, log)

def checkParamList = [
    params.input, params.fasta, params.annotation
]
for (param in checkParamList) {
    if (param) {
        file(param, checkIfExists: true)
    }
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { FASTP                       } from '../modules/nf-core/fastp/main'
include { SAMTOOLS_SORT               } from '../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_INDEX              } from '../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_FAIDX              } from '../modules/nf-core/samtools/faidx/main.nf'
include { SAMTOOLS_STATS              } from '../modules/nf-core/samtools/stats/main.nf'
include { SAMTOOLS_IDXSTATS           } from '../modules/nf-core/samtools/idxstats/main.nf'
include { MINIMAP2_INDEX              } from '../modules/nf-core/minimap2/index/main.nf'
include { MINIMAP2_ALIGN              } from '../modules/nf-core/minimap2/align/main.nf'

// GATK modules
include { GATK4_CREATESEQUENCEDICTIONARY } from '../modules/nf-core/gatk4/createsequencedictionary/main.nf'
include { ADDREADGROUPS                } from '../modules/local/addreadgroups.nf'
include { MARKDUPLICATE               } from '../modules/local/markduplicate.nf'
include { HAPLOTYPECALLER             } from '../modules/local/haplotypecaller.nf'
include { GATK4_GENOMICSDBIMPORT      } from '../modules/nf-core/gatk4/genomicsdbimport/main.nf'
include { GATK4_GENOTYPEGVCFS         } from '../modules/nf-core/gatk4/genotypegvcfs/main.nf'
include { GATK4_BASERECALIBRATOR as BQSR    } from '../modules/nf-core/gatk4/baserecalibrator/main.nf'
include { GATK4_APPLYBQSR as APPLYBQSR           } from '../modules/nf-core/gatk4/applybqsr/main.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

include { MKDBIMPORTLIST              } from '../modules/local/mkdbimportlist.nf'
include { DROPMINIMAPFIELD            } from '../modules/local/dropminimapfield.nf'


// subworkflow
include { FILTERSNPINDEL as filter_snp_indel     } from '../subworkflows/local/filtersnpindel.nf'
include { MERGECHR       as mergechr_bootstrap             } from '../subworkflows/local/mergechr.nf'
include { BQSR_RUN       as bqsr_run             } from '../subworkflows/local/bqsr_run.nf'
include { VARIANT_CALL   as post_calibration_call          } from '../subworkflows/local/variant_call.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SNVDETECT {

    // Read in samplesheet, validate and stage input files
    ch_versions = Channel.empty()
    input_ch = INPUT_CHECK(params.input)
    Channel.fromPath(params.fasta)
    | map{ file ->
        ref_id = file.name.tokenize('.')[0]
        meta = [id: ref_id]
        [meta, file]
    }
    | set { ref_ch }
    ref_ch.view()

    // Alignment START

    // Minimap2 index
    // make minimap2 index files from reference genome
    MINIMAP2_INDEX(ref_ch)
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

    // fastp quality filtering
    FASTP(input_ch.reads, [], false, false)
    ch_versions = ch_versions.mix(FASTP.out.versions)

    // minimap2 align
    // sorted bam output
    MINIMAP2_ALIGN(FASTP.out.reads, MINIMAP2_INDEX.out.index.collect(), true, false, false)
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    // samtools index
    // create index for alignment result
    SAMTOOLS_INDEX(MINIMAP2_ALIGN.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    // samtools stat
    SAMTOOLS_IDXSTATS(
        MINIMAP2_ALIGN.out.bam.join(SAMTOOLS_INDEX.out.bai)
        )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)

    SAMTOOLS_STATS(
        MINIMAP2_ALIGN.out.bam.join(SAMTOOLS_INDEX.out.bai),
        ref_ch.collect()
        )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    // Alignment Finished

    // GATK best practice workflow
    // join read group info with data channel by meta[sample_id]
    Channel.fromPath(params.readgroups)
    | splitCsv(header:true)
    | map {row ->
            rg = [row.RGSM, row.RGID, row.RGLB, row.RGPU, row.RGPL]
        }
    | set { readgroup_ch }

    MINIMAP2_ALIGN.out.bam
    | map { meta, bam -> tuple(meta.id, meta.single_end, bam) }
    | join ( readgroup_ch )
    | map { id, single_end, bam, rgid, rglb, rgpu, rgpl ->
            meta = [id:id, single_end:single_end, rgid:rgid, rglb:rglb, rgpu:rgpu, rgpl:rgpl]
            [meta, bam]
        }
    | set {ch_downstream }


    // Add read group info
    ADDREADGROUPS(ch_downstream)
    ch_versions = ch_versions.mix(ADDREADGROUPS.out.versions)
    // MarkDupliacte
    MARKDUPLICATE(ADDREADGROUPS.out.bam)
    ch_versions = ch_versions.mix(MARKDUPLICATE.out.versions)

    // Calibrate Base Quality Score
    // init-round --> HaplotypeCall --> filter --> get snpdb.vcf from filtered snv --> calibrate
    // initial round (obtain snv.db.vcf)
    // HplotypeCaller on uncalicrated data
    // HAPLOTYPECALLER(MARKDUPLICATE.out.bam)
    SAMTOOLS_FAIDX(ref_ch)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    // create reference dictionary
    GATK4_CREATESEQUENCEDICTIONARY(ref_ch)
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)

    // DROPMINIMAPFIELD(MARKDUPLICATE.out.bam)
    HAPLOTYPECALLER(
        MARKDUPLICATE.out.bam,
        MARKDUPLICATE.out.bai,
        ref_ch.collect(),
        SAMTOOLS_FAIDX.out.fai.collect(),
        GATK4_CREATESEQUENCEDICTIONARY.out.dict.collect(),
        [],
        [],
        [],
        true
        )

    // GenomicsDBImport
    // wrangle input
    HAPLOTYPECALLER.out.vcf
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
    | set {interval_ch}
    // create db on all chrs
    GATK4_GENOMICSDBIMPORT(vcf_in.collect(), interval_ch, [], [], false, false, true)

    // call snp indel on combined db
    GATK4_GENOMICSDBIMPORT.out.genomicsdb
    | map { meta, db ->
        [meta, db, [], [], []]
    }
    | set {genomicsdb_ch}
    GATK4_GENOTYPEGVCFS(
        genomicsdb_ch,
        ref_ch.collect(),
        SAMTOOLS_FAIDX.out.fai.collect(),
        GATK4_CREATESEQUENCEDICTIONARY.out.dict.collect(),
        [],
        []
    )

    GATK4_GENOTYPEGVCFS.out.vcf
    | join (GATK4_GENOTYPEGVCFS.out.tbi)
    | map{ meta, vcf, tbi ->
        chr = vcf.name =~ /(chr\d)_db/
        meta = [id: chr[0][1]]
        [meta, vcf, tbi]
    }
    | set { genotype_ch }

    // enter subworkflow
    // Filter low quality snv
    filter_snp_indel(genotype_ch, ref_ch, SAMTOOLS_FAIDX.out.fai, GATK4_CREATESEQUENCEDICTIONARY.out.dict)

    // wrangle data
    filter_snp_indel.out.snp_vcf
    | map { meta, vcf ->
        meta = [meta.type]
        [meta, vcf]
    }
    | groupTuple
    | map {meta, vcf ->
        meta = [id: meta[0]]
        [meta, vcf]
    }
    | set { ch_snp }

    filter_snp_indel.out.indel_vcf
    | map { meta, vcf ->
        meta = [meta.type]
        [meta, vcf]
    }
    | groupTuple
    | map {meta, vcf ->
        meta = [id: meta[0]]
        [meta, vcf]
    }
    | set { ch_indel }

    // combine all chrs
    mergechr_bootstrap(
        ch_snp,
        ch_indel,
        GATK4_CREATESEQUENCEDICTIONARY.out.dict
    )
    | set { merged }

    // wrangle data
    // vcf
    merged.snp_vcf
    | map { meta, vcf -> [[id: 'vcf_db'], vcf]}
    | concat (
        merged.indel_vcf.map({ meta, vcf -> [[id: 'vcf_db'], vcf]})
    )
    | groupTuple
    | set { ch_vcf_db }
    // tbi
    merged.snp_tbi
    | map { meta, tbi -> [[id: 'vcf_db'], tbi]}
    | concat (
        merged.indel_tbi.map({ meta, tbi -> [[id: 'vcf_db'], tbi]})
    )
    | groupTuple
    | set { ch_vcf_db_tbi }

    // basecall quality score calibration process /w bootstraped snp_db
    bqsr_run(
        MARKDUPLICATE.out.bam,
        MARKDUPLICATE.out.bai,
        ref_ch,
        SAMTOOLS_FAIDX.out.fai,
        GATK4_CREATESEQUENCEDICTIONARY.out.dict,
        ch_vcf_db,
        ch_vcf_db_tbi
    )
    | set { ch_bqsr_bam }

    post_calibration_call(
        ch_bqsr_bam.calibrated,
        MARKDUPLICATE.out.bai,
        ref_ch,
        SAMTOOLS_FAIDX.out.fai,
        GATK4_CREATESEQUENCEDICTIONARY.out.dict
    )
    | set { ch_unfiltered_genotype }
    ch_unfiltered_genotype.vcf.view()

    // HaplotypeCaller with calibrated bam

    // Filter low quality snv

    // combine vcf files


    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
