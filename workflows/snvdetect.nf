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
include { SAMTOOLS_IDXSTATS           } from '../modules/nf-core/samtools/idxstats/main.nf'
include { MINIMAP2_INDEX              } from '../modules/nf-core/minimap2/index/main.nf'
include { MINIMAP2_ALIGN              } from '../modules/nf-core/minimap2/align/main.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SNVDETECT {
    ch_versions = Channel.empty()
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    input_ch = INPUT_CHECK(params.input)
    Channel.fromPath(params.fasta)
    | map{ file ->
        ref_id = file.name.tokenize('.')[0]
        meta = [id: ref_id]
        [meta, file]
    }
    | set { ref_ch }
    ref_ch.view()
    //
    // MODULE: Run Test
    //
    // bwa index
    // make bwa index files from reference genome
    MINIMAP2_INDEX(ref_ch)
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

    // fastq quality filtering
    FASTP(input_ch.reads, [], false, false)
    ch_versions = ch_versions.mix(FASTP.out.versions)
    FASTP.out.reads.view()

    // minimap2 align
    // sorted bam output
    MINIMAP2_ALIGN(FASTP.out.reads, MINIMAP2_INDEX.out.index.collect(), true, false, false)
    MINIMAP2_ALIGN.out.bam.view()
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    // samtools index
    // create index for alignment result
    SAMTOOLS_INDEX(MINIMAP2_ALIGN.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    SAMTOOLS_INDEX.out.bai.view()

    // samtools stat
    SAMTOOLS_IDXSTATS(
        MINIMAP2_ALIGN.out.bam.join(SAMTOOLS_INDEX.out.bai)
        )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)
    SAMTOOLS_IDXSTATS.out.idxstats.view()

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
