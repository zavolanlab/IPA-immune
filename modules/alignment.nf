#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process STAR_INDEX_GENOME {

    label 'star'
    label 'indexing'
    
    // publishDir "${params.out_dir}/star_index_genome", mode: 'copy', pattern: 'starIndex'
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'

    input:
    path sequence

    output:
    path 'starIndex', emit: index
    path '*.log', emit: log

    script:
    """
    mkdir starIndex

    STAR --runThreadN ${params.threads_pe} \
        --runMode genomeGenerate \
        --genomeDir starIndex \
        --genomeFastaFiles ${sequence} \
        &> star_index_genome.log

    """
}

process STAR_ALIGN_PE {

    label "star"
    label "mapping"

    tag { library }

    publishDir "${params.out_dir}", mode: 'copy', pattern: "*.Aligned.sortedByCoord.out.bam"
    // publishDir "${params.out_dir}", mode: 'copy', pattern: "*.tab"
    // publishDir "${params.out_dir}", mode: 'copy', pattern: "*.Unmapped*"
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.out'

    input:
    tuple val(library), file(reads)
    path index

    output:
    tuple val(library), path('*.Aligned.sortedByCoord.out.bam'), emit: star_mapped_bam_tuple
    path '*.tab', emit: counts
    path '*.Unmapped*', emit: unmapped
    path '*.log', emit: log
    path '*.out', emit: out

    script:
    """
    STAR --runMode alignReads \
        --runThreadN ${params.threads_pe} \
        --genomeDir ${index} \
        --genomeLoad NoSharedMemory \
        --readFilesIn ${reads} \
        --limitOutSJcollapsed 5000000 \
        --outFileNamePrefix ${library}. \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM   SortedByCoordinate \
        --outSAMattributes All \
        --outBAMsortingThreadN 8 \
        --outSAMattrIHstart 0 \
        --outFilterType BySJout \
        --outFilterMultimapNmax 500000000 \
        --alignEndsType Local \
        --twopassMode None \
        &> ${library}_map_star.log
    """
}

process STAR_ALIGN_SE {

    label "star"
    label "mapping"
    
    tag { library } 

    publishDir "${params.out_dir}", mode: 'copy', pattern: "*.Aligned.sortedByCoord.out.bam"
    // publishDir "${params.out_dir}", mode: 'copy', pattern: "*.tab"
    // publishDir "${params.out_dir}", mode: 'copy', pattern: "*.Unmapped*"
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.out'

    input:
    tuple val(library), path(input_fastq)
    path index

    output:
    tuple val(library), path('*.Aligned.sortedByCoord.out.bam'), emit: star_mapped_bam
    path '*.tab', emit: counts
    path '*.Unmapped*', emit: unmapped
    path '*.log', emit: log
    path '*.out', emit: out

    script:
    """
    STAR --runMode alignReads \
        --runThreadN ${params.threads_se} \
        --genomeDir ${index} \
        --genomeLoad NoSharedMemory \
        --readFilesIn ${input_fastq} \
        --limitOutSJcollapsed 5000000 \
        --outFileNamePrefix ${library}. \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM   SortedByCoordinate \
        --outSAMattributes All \
        --outBAMsortingThreadN 8 \
        --outSAMattrIHstart 0 \
        --outFilterType BySJout \
        --outFilterMultimapNmax 500000000 \
        --alignEndsType Local \
        --twopassMode None \
        &> ${library}_map_star.log
    """
}
