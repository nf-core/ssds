//
// Sort, index BAM file and run samtools stats, flagstat and idxstats
//
include { PICARD_SORTSAM } from '../../modules/nf-core/modules/picard/sortsam/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/modules/samtools/index/main'

workflow SORT_AND_INDEX_BAM {
    take:
    ch_bam // channel: [ val(meta), [ bam ] ]

    main:
    ch_versions = Channel.empty()
    
    PICARD_SORTSAM     ( ch_bam, 'coordinate' )
    SAMTOOLS_INDEX     ( PICARD_SORTSAM.out.bam )

    ch_versions = ch_versions.mix(SAMTOOLS_SORTSAM.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    bam      = PICARD_SORTSAM.out.bam // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai // channel: [ val(meta), [ bai ] ]
    versions = ch_versions
}
