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

process STRINGTIE_COUNT_MATRIX{
        
    label 'stringtie'
    
    tag { library }

    publishDir "${params.out_dir}", mode: 'copy', pattern: '*.csv'
    
    input:
    tuple val(library), file(gtf)

    output:
    tuple val(library), path('*.csv'), emit: STRINGTIE_COUNT_MATRIX

    script:
    """
    python ${project_dir}/modules/prepDE.py ${gtf} -t ${library}_transcript_count_matrix.csv
    """
}