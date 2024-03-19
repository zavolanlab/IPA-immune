#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process TECTOOL {
    label "TECtool"
    container "docker://fgypas/tectool:0.4"

    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"

    input:

    path index_bai
    path sample_bam
    path annotation_gtf
    path polya_sites_bed
    path genome_fa

    script:
    """
    tectool \
        --annotation ${annotation_gtf} \
        --polyasites ${polya_sites_bed} \
        --bam ${sample_bam} \
        --genome ${genome_fa} \
        --output_dir ${params.out_dir} \
        --verbose
    """
}
