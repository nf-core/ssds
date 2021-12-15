//
// Make bigwigs from SSDS fragment BED files
//
include { UCSC_BEDGRAPHTOBIGWIG }                   from '../../modules/nf-core/modules/ucsc/bedgraphtobigwig/main'
include { GENERATE_SSDS_COVERAGE }                  from '../../modules/local/modules/generate_ssds_coverage/main'

workflow MAKE_SSDS_BIGWIGS {
    take:
    ssds_beds // channel: [ val(meta), [ bam ] ]
    genome_index
    genome_windows
    
    main:
    ch_versions = Channel.empty()
        
    GENERATE_SSDS_COVERAGE(ssds_beds, genome_index, genome_windows)

    ch_coverage_bg = GENERATE_SSDS_COVERAGE.out.bedgraph.flatten().map {
                                                    def id   = it.name.replaceFirst('^(.+)\\.(FWD|REV|FR|TOT)\\..+$','$1.$2')
                                                    def name = it.name.replaceFirst('^(.+)\\.(FWD|REV|FR|TOT)\\..+$','$1')
                                                    def meta=[id:"${id}", name:"${name}", "single-end":"false"];
                                                    return [meta, it]}
                                           .groupTuple()
    
    UCSC_BEDGRAPHTOBIGWIG (ch_coverage_bg, genome_index)

    ch_versions = ch_versions.mix(GENERATE_SSDS_COVERAGE.out.versions.first())
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions.first())
    
    emit:
    bigwig   = UCSC_BEDGRAPHTOBIGWIG.out.bigwig    // channel: [ val(meta), [ bai ] ]
    versions = ch_versions
}
