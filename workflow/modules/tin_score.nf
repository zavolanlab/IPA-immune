#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process TIN_GTF2BED {

    label 'gtf2bed'

    publishDir "${params.out_dir}", mode: 'copy', pattern: "*.bed"
    publishDir "${params.log_dir}", mode: 'copy', pattern: "*.log"

    input:
    path annotation_gtf

    output:
    path ('*.bed'), emit: transcripts_bed12

    script:
    """
    sed 's/transcript_type/transcript_biotype/g' ${annotation_gtf} > ${annotation_gtf}.biotype
    sort -k1,1 -k4,4n -k5,5nr ${annotation_gtf}.biotype > ${annotation_gtf}.sorted
    gtf2bed12 --gtf ${annotation_gtf}.sorted --bed12 full_transcripts_protein_coding.bed 2> ${annotation_gtf}_gtf2bed.log
    """
}

process CALCULATE_TIN_SCORES {

    label 'calculate_tin'

    tag { library }

    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.log_dir}/${library}_logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(library), path(bam)
    path transcripts_bed12

    output:
    path ('*.tsv'), emit: tin_scores_tsv

    script:
    """
    samtools index -@ ${params.threads_se} -M ${bam}
    calculate-tin.py -i ${bam} -r ${transcripts_bed12} --names ${library} -p ${params.threads_se} > ${library}_TIN_score.tsv 2> ${library}_tin_scores.log
    """
}