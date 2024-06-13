#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process TECTOOL {
    label "TECtool"
    container "docker://fgypas/tectool:0.4"
    
    tag { library_split }

    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*/*.tsv"
    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*_enriched_annotation.gtf"
    publishDir "${params.log_dir}/${library}_logs", mode: 'copy', pattern: '*.log'
    
    input:
    tuple val(library), file(bam)
    tuple val(library_split), path(bam_split)
    path annotation_gtf
    path polya_sites_bed
    path genome_fa
    
    output:
    tuple val(library_split), path('*.gtf'), emit: enriched_gtf
    path '*/*.tsv', emit: tsv


    script:
    """
    samtools index -@ ${params.threads_se} -M ${bam_split}
    tectool \
        --annotation ${annotation_gtf} \
        --polyasites ${polya_sites_bed} \
        --bam ${bam_split} \
        --genome ${genome_fa} \
        --num_cores ${params.threads_se} \
        --output_dir ${library_split}_tectool &> ${library_split}_tectool.log
    mv ${library_split}_tectool/enriched_annotation.gtf ${library_split}_enriched_annotation.gtf  
    """
}

process TECTOOL_MERGE {

    label 'bedtools'
    
    tag { library }

    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*_merged.gtf"

    input:
    tuple val(library), file(bam)
    tuple val(library_1), path(gtf_files_1)
    tuple val(library_2), path(gtf_files_2)

    output:
    tuple val(library), path('*_merged.gtf'), emit: merged_gtf

    script:
    """
    echo -e "${gtf_files_1}\n${gtf_files_2}" > ${library}_tectool_annotation_files.tsv
    tectool_add_novel_transcripts_to_gtf_file \
        --list_of_gtf_files ${library}_tectool_annotation_files.tsv \
        --out-dir ${library}_tectool_merged_annotations
    mv ${library}_tectool_merged_annotations/merged_annotation.gtf ${library}_tectool_annotation_merged.gtf
    """
}