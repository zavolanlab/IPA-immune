#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SAMTOOLS_INDEX {

    label "samtools"

	publishDir "${params.out_dir}", mode: 'copy', pattern: "*.bai"

    input:
    path bam

    output:
    path '*.bai', emit: bai

    script:
    """
    samtools index -@ ${params.threads} -M ${bam}
    """
}

process SAMTOOLS_GET_UNIQUE_MAPPERS {

    label "samtools"

	publishDir "${params.out_dir}", mode: 'copy', pattern: "*.bam_filtered"

    input:
    path input_bam

    output:
    path '*.bam_filtered', emit: filtered_bam

    script:
    """
    samtools view -@ ${params.threads} -h -q 255 -u ${input_bam} | \
        samtools sort -@ ${params.threads} -o output.bam_filtered
    """
}

process SAMTOOLS_GET_LOW_DUP_READS {

    label "samtools"

	publishDir "${params.out_dir}", mode: 'copy', pattern: "*.bam_low_dupl"

    input:
    path input_bam

    output:
    path '*.bam_low_dupl', emit: bam_low_dupl

    script:
    """
    samtools collate -@ ${params.threads} -O -u ${input_bam} | \
        samtools fixmate -@ ${params.threads} -m -u - - | \
        samtools sort -@ ${params.threads} -u - | \
        samtools markdup -@ ${params.threads} --duplicate-count -t -S --include-fails - out.mkdupped_bam; \
    samtools sort -@ ${params.threads} out.mkdupped_bam > out.sorted_mkdupped_bam; \
    samtools index -@ ${params.threads} out.sorted_mkdupped_bam; \
    printf "1\\n2\\n3\\n4\\n5\\n6\\n7\\n8\\n9\\n10" > out.selected_dup_levels_file; \
    samtools view out.sorted_mkdupped_bam -@ ${params.threads} -D dc:out.selected_dup_levels_file -u | samtools sort -@ ${params.threads} - > out.deduplicated_bam_file_intermediate; \
    samtools view -@ ${params.threads} out.deduplicated_bam_file_intermediate | awk -F"\\t" '{{print $1}}' > out.selected_read_names_file; \
    samtools view -@ ${params.threads} -D do:out.selected_read_names_file -u out.sorted_mkdupped_bam | samtools sort -@ ${params.threads} - > out.low_duplicates_intermediate; \
    samtools merge -f -@ ${params.threads} out.bam_low_dupl_nonsorted out.low_duplicates_intermediate out.deduplicated_bam_file_intermediate; \
    samtools sort -@ ${params.threads} out.bam_low_dupl_nonsorted > out.bam_low_dupl
    """
}


process SAMTOOLS_FASTQ {

    label "samtools"

	publishDir "${params.out_dir}", mode: 'copy', pattern: "*_1.fastq"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*_2.fastq"

    input:
    path bam

    output:
    path '*_1.fastq', emit: 1_fastq
    path '*_2.fastq', emit: 2_fastq

    script:
    """
    samtools fastq -@ ${params.threads} -1 paired_1.fastq -2 paired_2.fastq -0 /dev/null -s /dev/null
    """
}