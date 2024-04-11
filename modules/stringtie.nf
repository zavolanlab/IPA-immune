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
    tuple val(library), path('*.gtf'), emit: stringtie_gtf

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
    echo "${library} ${gtf}" > ${library}_sample_list.txt
    python ${projectDir}/modules/prepDE.py -i ${library}_sample_list.txt -t ${library}_transcript_count_matrix.csv
    echo "transcript_id,." > ${library}_novel_transcript_matrix.csv; grep "novel_" ${library}_transcript_count_matrix.csv >> ${library}_novel_transcript_matrix.csv
    """
}