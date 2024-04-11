#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process TECTOOL {
    label "TECtool"
    container "docker://fgypas/tectool:0.4"
    
    tag { library }

    publishDir "${params.out_dir}", mode: 'copy', pattern: "*/*.tsv"
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    
    input:
    tuple val(library), path(index)
    tuple val(library), path(bam)
    path annotation_gtf
    path polya_sites_bed
    path genome_fa
    
    output:
    tuple val(library), path('*/*.gtf'), emit: enriched_gtf
    path '*/*.tsv', emit: tsv


    script:
    """
    tectool \
        --annotation ${annotation_gtf} \
        --polyasites ${polya_sites_bed} \
        --bam ${bam} \
        --genome ${genome_fa} \
        --output_dir ${library}_tectool \
        &> ${library}_tectool.log
    """
}

process TECTOOL_RENAME {
    
    tag { library }
    
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*_renamed.gtf"
    
    input:
    tuple val(library), path(gtf)
    
    output:
    tuple val(library), path('*_enriched.gtf'), emit: renamed_gtf


    script:
    """
    mv ${gtf} ${library}_enriched.gtf
    """
}
