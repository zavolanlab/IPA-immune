#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process BEDTOOLS_INTERSECT {

    label 'bedtools'
    
    tag { library }

    input:
    tuple val(library), file(gtf1)
    tuple val(library), file(gtf2)

    output:
    tuple val(library), path('*.gtf'), emit: intersect_gtf

    script:
    """
    sort -k1,1 -k4,4n {gtf1} > enriched.sorted.1.gtf
    sort -k1,1 -k4,4n {gtf2} > enriched.sorted.2.gtf
    bedtools intersect -u -a enriched.sorted.1.gtf -b enriched.sorted.2.gtf > common_annotations.gtf
    """
}