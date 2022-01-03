process GENERATE_SSDS_COVERAGE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::pybedtools=0.8.2--py27h6a42192_1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/pybedtools:0.8.2--py27h6a42192_1' :
        'quay.io/biocontainers/pybedtools:0.8.2--py27h6a42192_1'}"

    input:
    tuple val(meta), path(bed)
    path genome_index
    tuple val(winmeta), path(windows_bed)
        
    output:
    path "*.bedgraph"  , emit: bedgraph
    path "versions.yml", emit: versions

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    """
    grep -w \\+ ${bed} >${prefix}.FWD.bed
    grep -w \\- ${bed} >${prefix}.REV.bed
    
    nF=`cat ${prefix}.FWD.bed |wc -l`
    nR=`cat ${prefix}.REV.bed |wc -l`
    
    make_SSDS_bedgraphs.py --fwd ${prefix}.FWD.bed \
                           --rev ${prefix}.REV.bed \
                           --g ${genome_index} \
                           --name ${prefix} \
                           --win ${windows_bed}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        make_SSDS_bedgraphs: \$(echo "0.8.2")
    END_VERSIONS
    """
}