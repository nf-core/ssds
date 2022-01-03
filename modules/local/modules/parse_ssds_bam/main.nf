process PARSESSDS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "mulled-v2-480c331443a1d7f4cb82aa41315ac8ea4c9c0b45" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-480c331443a1d7f4cb82aa41315ac8ea4c9c0b45:3e0fc1ebdf2007459f18c33c65d38d2b031b0052-0' :
        'quay.io/biocontainers/mulled-v2-480c331443a1d7f4cb82aa41315ac8ea4c9c0b45:3e0fc1ebdf2007459f18c33c65d38d2b031b0052-0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.bam')        , emit: bam
    tuple val(meta), path('*.bed')        , emit: bed
    tuple val(meta), path('*_report.txt') , emit: report
    path "versions.yml"                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix   = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    """
    parse_SSDS_BAM.py \\
        --bam $bam \\
        --name $prefix \\
        $args
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parse_SSDS_bam: \$(parse_SSDS_BAM.py --vers)
    END_VERSIONS
    """
}