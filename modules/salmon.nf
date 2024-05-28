#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SALMON_TRANSCRIPTOME {

    label 'salmon'
    
    tag { library }

    // publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: '*_transcriptome.fa'

    input:
    tuple val(library), file(gtf)
    path genome_fa

    output:
    tuple val(library), path('*_transcriptome.fa'), emit: transcriptome_fa

    script:
    """
    gffread \
        ${gtf} \
        -g ${genome_fa} \
        -w ${library}_transcriptome.fa
    """
}

process SALMON_INDEX {
        
    label 'salmon'
    
    tag { library }

    // publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: '*'
    
    input:
    tuple val(library), file(fasta)

    output:
    tuple val(library), path('*'), emit: salmon_index

    script:
    """
    salmon index \
        --transcripts ${fasta} \
        --index ${library}_transcripts_index \
        --keepDuplicates \
        --threads ${params.threads_pe}
    """
}

process SALMON_QUANTIFY {
        
    label 'salmon'
    
    tag { library }

    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: '*_quant.tsv'
    
    input:
    tuple val(library), path(index)
    tuple val(library_1), path(fastq_1)
    tuple val(library_2), path(fastq_2)
    
    output:
    tuple val(library), path('*_quant.tsv'), emit: salmon_counts

    script:
    """
    salmon quant \
        --index ${index} \
        --libType A \
        -1 ${fastq_1} -2 ${fastq_2} \
        --validateMappings \
        --seqBias \
        --output ${library}_transcript_quant \
        --threads ${params.threads_pe}
    mv ${library}_transcript_quant/quant.sf ${library}_quant.tsv
    """
}

process SALMON_FILTER_RESULT {
    
    tag { library }

    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: '*_filtered_quant.tsv'
    
    input:
    tuple val(library), file(gtf)
    path polya_sites_bed
    tuple val(library), file(tsv)

    output:
    tuple val(library), path('*_filtered_quant.tsv'), emit: salmon_filtered_counts

    script:
    """

    """
}