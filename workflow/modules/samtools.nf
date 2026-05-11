#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SAMTOOLS_GET_UNIQUE_MAPPERS {
    label 'samtools'
    tag   { meta.id }
    // Use a dedicated cpus param if provided, otherwise mirror SE/PE threads
    cpus  { (params.samtools_cpus ?: (meta.single_end ? params.threads_se : params.threads_pe)) as int }

    publishDir "${params.out_dir}/${meta.id}_results", mode: 'copy', pattern: '*.unique_filtered.bam'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.unique_filtered.bam"), emit: filtered_bam_tuple

    script:
    def mapq = (params.unique_mapq ?: 255) as int
    """
    set -euo pipefail

    # Filter to unique mappers by MAPQ, keep header, stream uncompressed to speed sorting
    samtools view -@ ${task.cpus} -h -q ${mapq} -u ${bam} \
        | samtools sort -@ ${task.cpus} -o ${meta.id}.unique_filtered.bam
    """
}

process SAMTOOLS_GET_LOW_DUP_READS {
    label 'samtools'
    tag   { meta.id }
    cpus  { (params.samtools_cpus ?: (meta.single_end ? params.threads_se : params.threads_pe)) as int }

    publishDir "${params.out_dir}/${meta.id}_results", mode: 'copy', pattern: '*.low_dupl.bam'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.low_dupl.bam"), emit: bam_low_dupl_tuple

    script:
    def dupMax = (params.dup_max ?: 10) as int
    """
    set -euo pipefail

    # 1) Prepare & mark duplicates (annotate dc / do tags), keep everything
    samtools collate -@ ${task.cpus} -O -u ${bam} \
        | samtools fixmate -@ ${task.cpus} -m -u - - \
        | samtools sort    -@ ${task.cpus} -u - \
        | samtools markdup -@ ${task.cpus} --duplicate-count -t -S --include-fails - ${meta.id}.mkdup.bam

    samtools sort  -@ ${task.cpus} -o ${meta.id}.mkdup.sorted.bam ${meta.id}.mkdup.bam
    samtools index -@ ${task.cpus} ${meta.id}.mkdup.sorted.bam

    # 2) Allowed duplicate counts: 1..dupMax
    seq 1 ${dupMax} > ${meta.id}.dup_allowed.txt

    # 3) Representatives / uniques (dc in allowed set)
    samtools view -@ ${task.cpus} -D dc:${meta.id}.dup_allowed.txt -u ${meta.id}.mkdup.sorted.bam \
        | samtools sort -@ ${task.cpus} -o ${meta.id}.lowdup.reps.bam

    # 4) Collect representative read names (QNAME)
    samtools view -@ ${task.cpus} ${meta.id}.lowdup.reps.bam \
        | awk -F'\\t' '{print \$1}' > ${meta.id}.lowdup.repnames.txt

    # 5) Pull duplicates whose 'do' (dup-original) belongs to those reps
    samtools view -@ ${task.cpus} -D do:${meta.id}.lowdup.repnames.txt -u ${meta.id}.mkdup.sorted.bam \
        | samtools sort -@ ${task.cpus} -o ${meta.id}.lowdup.dups.bam

    # 6) Merge reps + their duplicates, and finalize
    samtools merge -f -@ ${task.cpus} ${meta.id}.lowdup.merged.unsorted.bam \\
        ${meta.id}.lowdup.dups.bam ${meta.id}.lowdup.reps.bam

    samtools sort -@ ${task.cpus} -o ${meta.id}.low_dupl.bam ${meta.id}.lowdup.merged.unsorted.bam
    """
}

// process SAMTOOLS_GET_UNIQUE_MAPPERS {

//     label "samtools"
    
//     tag { library }

//     // publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*.unique_filtered.bam"

//     input:
//     tuple val(library), path(input_bam)

//     output:
//     tuple val(library), path('*.unique_filtered.bam'), emit: filtered_bam_tuple

//     script:
//     """
//     samtools view -@ ${params.threads_pe} -h -q 255 -u ${input_bam} | \
//         samtools sort -@ ${params.threads_pe} -o ${library}.unique_filtered.bam
//     """
// }

// process SAMTOOLS_GET_LOW_DUP_READS {

//     label "low_dup"

//     tag { library }

//     publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*.low_dupl.bam"

//     input:
//     tuple val(library), path(input_bam)

//     output:
//     tuple val(library), path('*.low_dupl.bam'), emit: bam_low_dupl_tupl

//     script:
//     """
//     samtools collate -@ ${params.threads_pe} -O -u ${input_bam} | \
//         samtools fixmate -@ ${params.threads_pe} -m -u - - | \
//         samtools sort -@ ${params.threads_pe} -u - | \
//         samtools markdup -@ ${params.threads_pe} --duplicate-count -t -S --include-fails - out.mkdupped_bam; \
//     samtools sort -@ ${params.threads_pe} out.mkdupped_bam > out.sorted_mkdupped_bam; \
//     samtools index -@ ${params.threads_pe} out.sorted_mkdupped_bam; \
//     printf "1\\n2\\n3\\n4\\n5\\n6\\n7\\n8\\n9\\n10" > out.selected_dup_levels_file; \
//     samtools view out.sorted_mkdupped_bam -@ ${params.threads_pe} -D dc:out.selected_dup_levels_file -u | samtools sort -@ ${params.threads_pe} - > out.deduplicated_bam_file_intermediate; \
//     samtools view -@ ${params.threads_pe} out.deduplicated_bam_file_intermediate | awk -F"\\t" '{{print \$1}}' > out.selected_read_names_file; \
//     samtools view -@ ${params.threads_pe} -D do:out.selected_read_names_file -u out.sorted_mkdupped_bam | samtools sort -@ ${params.threads_pe} - > out.low_duplicates_intermediate; \
//     samtools merge -f -@ ${params.threads_pe} out.bam_low_dupl_nonsorted out.low_duplicates_intermediate out.deduplicated_bam_file_intermediate; \
//     samtools sort -@ ${params.threads_pe} out.bam_low_dupl_nonsorted > ${library}.low_dupl.bam
//     """
// }


process SAMTOOLS_BAM2FASTQ {

    label "samtools"
    
    tag { library }
    
    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*_1.fastq"
    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*_2.fastq"
    publishDir "${params.log_dir}/${library}_logs", mode: 'copy', pattern: '*.log'

    input:
    tuple val(library), path(bam)

    output:  
    tuple val("${library}_1"), path("${library}_1.fastq"), emit: fastq1_tuple
    tuple val("${library}_2"), path("${library}_2.fastq"), emit: fastq2_tuple
    path '*.log', emit: log

    script:
    """
    samtools sort -n -@ ${params.threads_pe} ${bam} -o ${library}.out.querysort.bam
    samtools fastq -@ ${params.threads_pe} -1 "${library}_1.fastq" -2 "${library}_2.fastq" -0 /dev/null -s /dev/null ${library}.out.querysort.bam &> ${library}_bam2fastq.log
    """
}

process SAMTOOLS_DEPTH {

    label "samtools"

    tag { library }

    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*.coverage.bed"

    input:
    tuple val(library), path(input_bam)

    output:
    tuple val(library), path('*.coverage.bed'), emit: depth_bed

    script:
    """
    samtools depth -@ ${params.threads_pe} -s ${input_bam} -o ${library}.coverage.bed
    """
}