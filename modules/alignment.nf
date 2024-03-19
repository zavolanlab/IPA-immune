#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process MAP_STAR {

    label "star"
	label "mapping"

    debug true

    publishDir "${params.out_dir}", mode: 'copy', pattern: "*.Aligned.sortedByCoord.out.bam"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*.tab"
	publishDir "${params.out_dir}", mode: 'copy', pattern: "*.Unmapped*"
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'
    publishDir "${params.log_dir}", mode: 'copy', pattern: '*.out'

    input:
    each(path(reads))
    path index

    output:
    path '*.Aligned.sortedByCoord.out.bam', emit: aligned
    path '*.tab', emit: counts
    path '*.Unmapped*', emit: unmapped
	path '*.log', emit: log
	path '*.out', emit: out

    script:
    """
    for VAR in ${reads}
    do

        input=\$(basename \$VAR)
        prefix=\$(echo "\$input" | sed 's/\\(\\.fastq\\.gz\\)*\$//') 

        STAR --runMode alignReads \
            --runThreadN ${params.star_threads} \
            --genomeDir ${index} \
            --genomeLoad NoSharedMemory \
            --readFilesIn \$VAR \
            --limitOutSJcollapsed 5000000 \
            --outFileNamePrefix \$prefix. \
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

    done

    """

}

process INDEX_SAMTOOLS {

    label "samtools"

    debug true

	publishDir "${params.out_dir}", mode: 'copy', pattern: "*.bai"

    input:
    path bam

    output:
    path '*.bai', emit: bai

    script:
    """
    samtools index -M ${bam}
    """
}
