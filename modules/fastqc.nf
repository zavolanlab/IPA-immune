#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process FASTQC {

    label 'fastqc'
    
    publishDir "${params.out_dir}", mode: 'copy', pattern: '*_fastqc.html'

    input:
    tuple val(library), file(reads)

    output:
    path '*_fastqc.html'

    script:
    """
    fastqc -t ${params.threads_se} ${reads}
    """
}