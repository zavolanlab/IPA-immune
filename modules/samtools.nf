process SAMTOOLS_INDEX {

    label "samtools"

    debug true

	publishDir "${params.out_dir}", mode: 'copy', pattern: "*.bai"

    input:
    path bam

    output:
    path '*.bai', emit: bai

    script:
    """
    samtools index -@ ${params.star_threads} -M ${bam}
    """
}


process SAMTOOLS_GET_UNIQUE_MAPPERS {

    label "samtools"

    debug true

	publishDir "${params.out_dir}", mode: 'copy', pattern: "*.bam"

    input:
    path bam

    output:
    path '*.bam', emit: filtered_bam

    script:
    """
    samtools view -@ ${params.star_threads} -h -q 255 -u ${bam} | \
        samtools sort -@ ${params.star_threads} -o {output.bam_filtered}
    """
}
