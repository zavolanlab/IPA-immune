#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process TECTOOL {
    label "TECtool"
    container "docker://fgypas/tectool:0.4"

    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    
    input:
    tuple val(library), path(index)
    tuple val(library), path(bam)
    path annotation_gtf
    path polya_sites_bed
    path genome_fa
    
    output:
    tuple val(library), path('*'), emit: enriched_gtf

    script:
    """
    tectool \
        --annotation ${annotation_gtf} \
        --polyasites ${polya_sites_bed} \
        --bam ${bam} \
        --genome ${genome_fa} \
        --output_dir ${library} \
        &> ${library}_tectool.log
    """
}
