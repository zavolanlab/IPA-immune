#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process BEDTOOLS_MERGE {

    label 'bedtools'
    
    tag { library }

    // publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*_merged.gtf"

    input:
    tuple val(library), file(bam)
    tuple val(library_1), file(gtf1)
    tuple val(library_2), file(gtf2)

    output:
    tuple val(library), path('*_merged.gtf'), emit: merged_gtf

    script:
    """
    cat ${gtf1} ${gtf2} | sort -k1,1 -k4,4n | uniq > ${library}_merged.gtf
    """
}