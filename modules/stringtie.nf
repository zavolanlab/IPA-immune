#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process STRINGTIE_QUANTIFY {

    label 'stringtie'
    
    tag { library }

    publishDir "${params.out_dir}", mode: 'copy', pattern: '*.gtf'

    input:
    tuple val(library), file(bam)
    tuple val(library), file(gtf)

    output:
    tuple val(library), path('*.gtf'), emit: stringie_gtf

    script:
    """
    stringtie ${bam} -e -G ${gtf} -o ${library}_stringtie_quantified.gtf
    """
}