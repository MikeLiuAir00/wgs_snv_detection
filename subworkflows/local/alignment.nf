// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SAMTOOLS_SORT      } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../modules/nf-core/samtools/index/main'
include { BWA_INDEX          } from '../../modules/nf-core/bwa/index/main.nf'
include { BWA_MEM            } from '../../modules/nf-core/bwa/index/main.nf'
workflow ALIGNMENT {

    take:
    read_pairs // channel: [ val(meta), [ sample_id, path ] ]
    reference // channel: [val(meta), path(ref_path) ]
    annotation // channel: [val(meta), path(annotation_path) ]

    main:
    // version log
    ch_versions = Channel.empty()

    // pipeline
    // build bwa index
    BWA_INDEX(reference)
    ch_versions = ch_versions.mix(BWA_INDEX.out.version.first())
    // bwa alignment sort with samtools
    BWA_MEM(read_pairs, BWA_INDEX.out, true)
    ch_versions = ch_versions.mix(BWA_MEM.out.version.first())


    emit:
    ref_idx  = BWA_INDEX.out.index
    bam      = BWA_MEM.out.bam           // channel: [ val(meta), [ bam ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

