//
// SSDS alignment 
//
include { MARKDUPLICATES_AND_INDEX as MD_AND_INDEX}     from '../../subworkflows/local/markduplicates_and_index_bam'
include { MARKDUPLICATES_AND_INDEX as MD_AND_INDEX_RAW} from '../../subworkflows/local/markduplicates_and_index_bam'

include { TRIMGALORE }                                  from '../../modules/nf-core/modules/trimgalore/main'
include { BWA_MEM }                                     from '../../modules/nf-core/modules/bwa/mem/main'
include { SAMTOOLS_VIEW as FILTER_SUPP_ALIGNMENTS}      from '../../modules/nf-core/modules/samtools/view/main'
include { PICARD_COLLECTMULTIPLEMETRICS }               from '../../modules/nf-core/modules/picard/collectmultiplemetrics/main'
include { PICARD_SORTSAM as PICARD_SORTSAM_QRY}         from '../../modules/nf-core/modules/picard/sortsam/main'
include { PARSESSDS }                                   from '../../modules/local/modules/parse_ssds_bam/main'
include { BEDTOOLS_SORT}                                from '../../modules/nf-core/modules/bedtools/sort/main'
include { SAMTOOLS_IDXSTATS }                           from '../../modules/nf-core/modules/samtools/idxstats/main'

workflow ALIGN_SS {
    take:
    fastq // channel: [ val(meta), [ bam ] ]
    bwa
    fasta
    
    main:
    
    ch_versions = Channel.empty()
    
    TRIMGALORE (fastq)
    BWA_MEM (TRIMGALORE.out.reads, bwa, false)
    FILTER_SUPP_ALIGNMENTS (BWA_MEM.out.bam, fasta)
    
    ch_unprocessed = FILTER_SUPP_ALIGNMENTS.out.bam.map {return [[id:"${it[0].id}.unprocessed_SSDS", single_end: it[0].single_end], it[1]]}
    MD_AND_INDEX_RAW(ch_unprocessed, fasta)
    
    PICARD_COLLECTMULTIPLEMETRICS(MD_AND_INDEX_RAW.out.bam, fasta)
    
    PICARD_SORTSAM_QRY (FILTER_SUPP_ALIGNMENTS.out.bam.map {return [[id:it[0].id.replaceAll(".unprocessed_SSDS",""), single_end: it[0].single_end], it[1]]}, 'queryname')
    PARSESSDS(PICARD_SORTSAM_QRY.out.bam)
    
    ch_ssds_bam_unsorted = PARSESSDS.out.bam.collect{it[1]}.flatten().map {
                                                 def id   = it.name.replaceFirst('^(.+)_(ssDNA|ssLow|type2|dsDNA|unclassified)\\..+$','$1.$2')
                                                 def name = it.name.replaceFirst('^(.+)_(ssDNA|ssLow|type2|dsDNA|unclassified)\\..+$','$1')
                                                 def type = it.name.replaceFirst('^(.+)_(ssDNA|ssLow|type2|dsDNA|unclassified)\\..+$','$2')
                                                 def meta = [id:"${id}", name:"${name}", type:"${type}", "single-end":"false"];
                                                 return [meta, it]}
                                            .groupTuple()
    
    ch_ssds_bed_unsorted = PARSESSDS.out.bed.collect{it[1]}.flatten().map {
                                                 def id   = it.name.replaceFirst('^(.+)_(ssDNA|ssLow|type2|dsDNA|unclassified)\\..+$','$1.$2')
                                                 def name = it.name.replaceFirst('^(.+)_(ssDNA|ssLow|type2|dsDNA|unclassified)\\..+$','$1')
                                                 def type = it.name.replaceFirst('^(.+)_(ssDNA|ssLow|type2|dsDNA|unclassified)\\..+$','$2')
                                                 def meta = [id:"${id}", name:"${name}", type:"${type}", "single-end":"false"];
                                                 return [meta, it]}
                                           .groupTuple()
    
    BEDTOOLS_SORT(ch_ssds_bed_unsorted, 'bed')
    
    MD_AND_INDEX (ch_ssds_bam_unsorted, fasta)

    ch_ssds_bam = MD_AND_INDEX.out.bam.join(MD_AND_INDEX.out.bai)    
    
    SAMTOOLS_IDXSTATS (ch_ssds_bam)
    
    ch_reports = TRIMGALORE.out.log.collect()
                 .mix(PICARD_COLLECTMULTIPLEMETRICS.out.metrics)
                 .mix(MD_AND_INDEX_RAW.out.mdmetrics)
                 .mix(PARSESSDS.out.report)
    
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())
    ch_versions = ch_versions.mix(FILTER_SUPP_ALIGNMENTS.out.versions.first())
    ch_versions = ch_versions.mix(MD_AND_INDEX_RAW.out.versions.first())
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())
    ch_versions = ch_versions.mix(PICARD_SORTSAM_QRY.out.versions.first())
    ch_versions = ch_versions.mix(PARSESSDS.out.versions.first())
    ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions.first())
    ch_versions = ch_versions.mix(MD_AND_INDEX.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.first())
    
    emit:
    bam      = ch_ssds_bam                // channel: [ val(meta), [ bam ] ]
    bed      = BEDTOOLS_SORT.out.sorted   // channel: [ val(meta), [ bai ] ]
    reports  = ch_reports
    versions = ch_versions           //    path: *.version.txt
}
