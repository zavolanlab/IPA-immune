#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SAMTOOLS_INDEX {

    label "samtools"
    
    tag { library }

    // publishDir "${params.out_dir}", mode: 'copy', pattern: "*.bai"

    input:
    tuple val(library), path(bam)

    output:
    tuple val(library), path('*.bai'), emit: index

    script:
    """
    samtools index -@ ${params.threads_se} -M ${bam}
    """
}

process SAMTOOLS_GET_UNIQUE_MAPPERS {

    label "samtools"
    
    tag { library }

    // publishDir "${params.out_dir}", mode: 'copy', pattern: "*.bam_filtered"

    input:
    tuple val(library), path(input_bam)

    output:
    tuple val(library), path('*.bam_filtered'), emit: filtered_bam_tuple

    script:
    """
    samtools view -@ ${params.threads_pe} -h -q 255 -u ${input_bam} | \
        samtools sort -@ ${params.threads_pe} -o output.bam_filtered
    """
}

process SAMTOOLS_GET_LOW_DUP_READS {

    label "samtools"

    tag { library }

    // publishDir "${params.out_dir}", mode: 'copy', pattern: "*.bam_low_dupl"

    input:
    tuple val(library), path(input_bam)

    output:
    tuple val(library), path('*.bam_low_dupl'), emit: bam_low_dupl_tupl

    script:
    """
    samtools collate -@ ${params.threads_pe} -O -u ${input_bam} | \
        samtools fixmate -@ ${params.threads_pe} -m -u - - | \
        samtools sort -@ ${params.threads_pe} -u - | \
        samtools markdup -@ ${params.threads_pe} --duplicate-count -t -S --include-fails - out.mkdupped_bam; \
    samtools sort -@ ${params.threads_pe} out.mkdupped_bam > out.sorted_mkdupped_bam; \
    samtools index -@ ${params.threads_pe} out.sorted_mkdupped_bam; \
    printf "1\\n2\\n3\\n4\\n5\\n6\\n7\\n8\\n9\\n10" > out.selected_dup_levels_file; \
    samtools view out.sorted_mkdupped_bam -@ ${params.threads_pe} -D dc:out.selected_dup_levels_file -u | samtools sort -@ ${params.threads_pe} - > out.deduplicated_bam_file_intermediate; \
    samtools view -@ ${params.threads_pe} out.deduplicated_bam_file_intermediate | awk -F"\\t" '{{print \$1}}' > out.selected_read_names_file; \
    samtools view -@ ${params.threads_pe} -D do:out.selected_read_names_file -u out.sorted_mkdupped_bam | samtools sort -@ ${params.threads_pe} - > out.low_duplicates_intermediate; \
    samtools merge -f -@ ${params.threads_pe} out.bam_low_dupl_nonsorted out.low_duplicates_intermediate out.deduplicated_bam_file_intermediate; \
    samtools sort -@ ${params.threads_pe} out.bam_low_dupl_nonsorted > out.bam_low_dupl
    """
}


process SAMTOOLS_BAM2FASTQ {

    label "samtools"
    
    tag { library }
    
    // publishDir "${params.out_dir}", mode: 'copy', pattern: "*_1.fastq"
    // publishDir "${params.out_dir}", mode: 'copy', pattern: "*_2.fastq"
    // publishDir "${params.out_dir}", mode: 'copy', pattern: "*_singleton.fastq"
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'

    input:
    tuple val(library), path(bam)

    output:  
    tuple val("${library}_1"), path("${library}_1.fastq"), emit: fastq1_tuple
    tuple val("${library}_2"), path("${library}_2.fastq"), emit: fastq2_tuple
    tuple val("${library}_singleton"), path("${library}_singleton.fastq"), emit: singleton_tuple

    script:
    """
    samtools fastq -@ ${params.threads_pe} -1 "${library}_1.fastq" -2 "${library}_2.fastq" -0 /dev/null -s "${library}_singleton.fastq" ${bam} &> ${library}_bam2fastq.log
    """
}

