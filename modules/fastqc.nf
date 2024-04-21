#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process FASTQC {

    label 'fastqc'
    
    tag { library }
    
    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: '*_fastqc.html'

    input:
    tuple val(library), file(reads)

    output:
    path '*_fastqc.html'

    script:
    """
    fastqc -t ${params.threads_se} ${reads}
    """
}