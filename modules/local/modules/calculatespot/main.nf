process CALCULATESPOT {
    tag "${meta.id}"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::pybedtools=0.8.2--py27h6a42192_1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "quay.io/biocontainers/pybedtools:0.8.2--py27h6a42192_1" :
        "quay.io/biocontainers/pybedtools:0.8.2--py27h6a42192_1"}"

    input:
    tuple val(meta), path(reads_bam), path(reads_bai), val(bamids), val(names)
    path(index)
    tuple val(imeta), path(interval_beds)
    
    output:
    tuple val(meta), path('*report.txt'), emit: report
    path "versions.yml"                 , emit: versions

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    """    
    calculate_SPoT.py \\
        --reads_bam \"${bamids}\" \\
        --intervals all \\
        --name \"${names}\" \\
        --iname all \\
        --g ${index} \\
        --rand \\
        --o ${meta.id}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        calculate_SSDS_SPoT: \$(calculate_SPoT.py --vers)
    END_VERSIONS
    """
}
