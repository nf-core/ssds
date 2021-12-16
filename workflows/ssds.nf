/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSsds.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check genome for ssds intervals
def spot_intervals_genome = params.spot_intervals_genome
if (!spot_intervals_genome) {
    spot_intervals_genome = params.genome ? params.genome : null
}

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

if (params.bwa)   {
    if (file(params.bwa).isFile()){ 
        ch_bwa = file(params.bwa).getParent() 
    } else {
        ch_bwa = file(params.bwa) 
    }
} 

ch_genome_index  = file("${params.fasta}.fai", checkIfExists: true) 
ch_makewin_input = [ [ id:"${params.genome}"],
                     file("${params.fasta}.fai", checkIfExists: true) ]

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    SSDS INTERVAL FILES
========================================================================================
*/
def ch_spot_intervals

if (params.spot_interval_beds["${spot_intervals_genome}"]){
    ch_spot_intervals = Channel.from(params.spot_interval_beds["${spot_intervals_genome}"])
                                   .map {[[id:it.id, 
                                           genome:spot_intervals_genome], 
                                           file(it.bed, checkIfExists: true)]}
}else{
    System.err.println("WARNING: No SSDS intervals found for SPoT analysis. SPoT analyses will not be performed.")
    ch_spot_intervals = Channel.empty()
}

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK }        from '../subworkflows/local/input_check'
include { ALIGN_SS }           from '../subworkflows/local/align_ss'
include { SORT_AND_INDEX_BAM } from '../subworkflows/local/sort_and_index_bam'
include { MAKE_SSDS_BIGWIGS }  from '../subworkflows/local/make_ssds_bigwigs'


/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { FASTQC}                       from '../modules/nf-core/modules/fastqc/main'
include { BWA_INDEX }                   from '../modules/nf-core/modules/bwa/index/main'
include { DEEPTOOLS_PLOTFINGERPRINT}    from '../modules/nf-core/modules/deeptools/plotfingerprint/main'
include { BEDTOOLS_MAKEWINDOWS }        from '../modules/nf-core/modules/bedtools/makewindows/main'
include { CALCULATESPOT }               from '../modules/local/modules/calculatespot/main.nf'
include { MULTIQC as MULTIQC_SSDS}      from '../modules/local/modules/multiqc_ssds/main'
/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SSDS {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    //
    // MODULE/SUBWORKFLOW: Index genome if required, then run SSDS alignment sub-workflow    
    //
    if (!params.bwa){
        BWA_INDEX (
            file(params.fasta)
        )
        ch_bwa = BWA_INDEX.out.index
        ch_versions = ch_versions.mix(BWA_INDEX.out.versions.ifEmpty(null))
    }
    
    ALIGN_SS(INPUT_CHECK.out.reads, ch_bwa, file(params.fasta))
    ch_versions = ch_versions.mix(ALIGN_SS.out.versions.ifEmpty(null))
    
    //
    // MODULE/SUBWORKFLOW: Run MAKE_SSDS_BIGWIGS subworkflow to generate SSDS coverage bigwigs; make genomic windows first
    //
    BEDTOOLS_MAKEWINDOWS (
        ch_makewin_input, false 
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions.ifEmpty(null))
    
    def ch_ssds_beds = ALIGN_SS.out.bed
                        .filter {file(it[1]).size() > 40000}
    
    MAKE_SSDS_BIGWIGS(ch_ssds_beds, ch_genome_index, BEDTOOLS_MAKEWINDOWS.out.tab)
    ch_versions = ch_versions.mix(MAKE_SSDS_BIGWIGS.out.versions.ifEmpty(null))
    
    //
    // MODULE: Run Calculatespot  
    // 
    ch_bams_for_spot = ALIGN_SS.out.bam.map { [["id":it[0].name, "single-end":it[0]."single-end"], it[1], it[2], it[1].name, it[0].type ] }
                                  .groupTuple(by:0)
    
    ch_intervals_for_spot = ch_spot_intervals.map{[[genome:it[0].genome],it[1]]}.groupTuple(by:0).collect()
    
    CALCULATESPOT (
        ch_bams_for_spot, ch_genome_index, ch_intervals_for_spot
    )
    ch_versions = ch_versions.mix(CALCULATESPOT.out.versions.ifEmpty(null))
    
    //
    // MODULE: Run Deeptools plotfingerprint
    //
    ch_fingerprint_bams = ALIGN_SS.out.bam.map { [["id":it[0].name, "single-end":it[0]."single-end"], it[1], it[2]] }
                                     .groupTuple(by:0)
    
    DEEPTOOLS_PLOTFINGERPRINT (
        ch_fingerprint_bams
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_PLOTFINGERPRINT.out.versions.ifEmpty(null))
    
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    
    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowSsds.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)
    
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN_SS.out.reports.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CALCULATESPOT.out.report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_PLOTFINGERPRINT.out.metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_PLOTFINGERPRINT.out.matrix.collect{it[1]}.ifEmpty([]))

    MULTIQC_SSDS (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC_SSDS.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC_SSDS.out.versions)
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
