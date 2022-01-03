//
// Sort, mark duplicates, index BAM file and get metrics
//
include { PICARD_SORTSAM }                from '../../modules/nf-core/modules/picard/sortsam/main'
include { PICARD_MARKDUPLICATES }         from '../../modules/nf-core/modules/picard/markduplicates/main'
include { SAMTOOLS_INDEX }                from '../../modules/nf-core/modules/samtools/index/main'

workflow MARKDUPLICATES_AND_INDEX {
    take:
    ch_bam // channel: [ val(meta), [ bam ] ]
    fasta
    
    main:
    ch_versions = Channel.empty()
    
    PICARD_SORTSAM     ( ch_bam, 'coordinate' )

    PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)   

    SAMTOOLS_INDEX (PICARD_MARKDUPLICATES.out.bam)\

    ch_versions = ch_versions.mix(PICARD_SORTSAM.out.versions.first())
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
    
    emit:
    bam       = PICARD_MARKDUPLICATES.out.bam              // channel: [ val(meta), [ bam ] ]
    bai       = SAMTOOLS_INDEX.out.bai                     // channel: [ val(meta), [ bai ] ]
    mdmetrics = PICARD_MARKDUPLICATES.out.metrics          // channel: [ val(meta), [ bam ] ]
    versions  = ch_versions          // path: *.version.txt

}
