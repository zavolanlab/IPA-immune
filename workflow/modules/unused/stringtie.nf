#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process STRINGTIE_QUANTIFY {

    label 'stringtie'
    
    tag { library }

    publishDir "${params.out_dir}", mode: 'copy', pattern: '*_stringtie_quantified.gtf'

    input:
    tuple val(library), file(bam)
    tuple val(library), file(gtf)

    output:
    tuple val(library), path('*_stringtie_quantified.gtf'), emit: stringtie_gtf

    script:
    """
    stringtie -o ${library}_stringtie_quantified.gtf -p ${params.threads_pe} -e -G ${gtf} ${bam}
    """
}

process STRINGTIE_COUNT_MATRIX{
        
    label 'stringtie'
    
    tag { library }

    publishDir "${params.out_dir}", mode: 'copy', pattern: '*_transcript_count_matrix.csv'
    
    input:
    tuple val(library), file(gtf)

    output:
    tuple val(library), path('*_transcript_count_matrix.csv'), emit: STRINGTIE_COUNT_MATRIX

    script:
    """
    echo "${library} ${gtf}" > ${library}_sample_list.txt
    python ${projectDir}/modules/prepDE.py -l 48 -i ${library}_sample_list.txt -t ${library}_transcript_count_matrix.csv
    echo "transcript_id,." > ${library}_novel_transcript_count_matrix.csv; grep "novel_" ${library}_transcript_count_matrix.csv >> ${library}_novel_transcript_count_matrix.csv
    """
}