#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process STAR_INDEX_GENOME {

    label 'star'
	label 'indexing'

    publishDir "${params.out_dir}/star_index_genome", mode: 'copy', pattern: 'starIndex'
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'

    input:
    path sequence

    output:
    path 'starIndex', emit: index
	path '*.log', emit: log

    script:
    """
    mkdir starIndex

    STAR --runThreadN ${params.threads} \
    	--runMode genomeGenerate \
        --genomeDir starIndex \
        --genomeFastaFiles ${sequence} \
		&> star_index_genome.log

    """
}

process STAR_ALIGN_PE {

    label "star"
	label "mapping"

    debug true

    publishDir "${params.out_dir}", mode: 'copy', pattern: "*.Aligned.sortedByCoord.out.bam"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*.tab"
	publishDir "${params.out_dir}", mode: 'copy', pattern: "*.Unmapped*"
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.out'

    input:
    path input_fastq_1
    path input_fastq_2
    path index

    output:
    path '*.Aligned.sortedByCoord.out.bam', emit: aligned
    path '*.tab', emit: counts
    path '*.Unmapped*', emit: unmapped
	path '*.log', emit: log
	path '*.out', emit: out

    script:
    """
    prefix=\$(echo "${input_fastq_1}" | sed 's/\\(\\.fastq\\.gz\\)*\$//') 

    STAR --runMode alignReads \
        --runThreadN ${params.threads} \
        --genomeDir ${index} \
        --genomeLoad NoSharedMemory \
        --readFilesIn ${input_fastq_1} ${input_fastq_2} \
        --limitOutSJcollapsed 5000000 \
        --outFileNamePrefix \${prefix}. \
        --readFilesCommand zcat \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM   SortedByCoordinate \
        --outSAMattributes All \
        --outBAMsortingThreadN 8 \
        --outSAMattrIHstart 0 \
        --outFilterType BySJout \
        --outFilterMultimapNmax 500000000 \
        --alignEndsType Local \
        --twopassMode None
        &> \${prefix}_map_star.log
    """
}

process STAR_ALIGN_SE {

    label "star"
	label "mapping"

    debug true

    publishDir "${params.out_dir}", mode: 'copy', pattern: "*.Aligned.sortedByCoord.out.bam"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*.tab"
	publishDir "${params.out_dir}", mode: 'copy', pattern: "*.Unmapped*"
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.out'

    input:
    path input_fastq
    path index

    output:
    path '*.Aligned.sortedByCoord.out.bam', emit: aligned
    path '*.tab', emit: counts
    path '*.Unmapped*', emit: unmapped
	path '*.log', emit: log
	path '*.out', emit: out

    script:
    """
    prefix=\$(echo "${input_fastq}") 

    STAR --runMode alignReads \
        --runThreadN ${params.threads} \
        --genomeDir ${index} \
        --genomeLoad NoSharedMemory \
        --readFilesIn ${input_fastq} \
        --limitOutSJcollapsed 5000000 \
        --outFileNamePrefix \${prefix}. \
        --readFilesCommand zcat \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM   SortedByCoordinate \
        --outSAMattributes All \
        --outBAMsortingThreadN 8 \
        --outSAMattrIHstart 0 \
        --outFilterType BySJout \
        --outFilterMultimapNmax 500000000 \
        --alignEndsType Local \
        --twopassMode None
        &> \${prefix}_map_star.log
    """
}