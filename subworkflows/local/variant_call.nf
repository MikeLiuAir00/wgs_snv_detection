
workflow VARIANT_CALL {

    take:
    ch_bam
    ch_bai
    ch_ref
    ch_fai
    ch_dict
    main:

    ch_versions = Channel.empty()



    emit:

    versions = ch_versions                     // channel: [ versions.yml ]
}

