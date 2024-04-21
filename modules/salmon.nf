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
    gffread -w ${library}_transcriptome.fa -g ${genome_fa} ${gtf}
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
    salmon index -t ${fasta} -i ${library}_transcripts_index
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
    salmon quant -i ${index} -l A -1 ${fastq_1} -2 ${fastq_2} --validateMappings -o ${library}_transcript_quant
    mv ${library}_transcript_quant/quant.sf ${library}_quant.tsv
    """
}