#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process BEDTOOLS_MERGE {

    label 'bedtools'
    
    tag { library }

    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*_merged.gtf"

    input:
    tuple val(library), file(bam)
    tuple val(library_1), file(gtf1)
    tuple val(library_2), file(gtf2)

    output:
    tuple val(library), path('*_merged.gtf'), emit: merged_gtf

    script:
    """
    cat ${gtf1} ${gtf2} | sort -k1,1 -k4,4n | uniq > ${library}_merged.gtf
    """
}

process BEDTOOLS_FILTER_QUANT {

    label 'bedtools'
    
    tag { library }

    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*.csv"
    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*_overlaps_a.bed"
    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*_overlaps_a.gtf"

    input:
    tuple val(library), file(merged_gtf)
    path pas_bed
    tuple val(library), file(quant)

    output:
    tuple val(library), path('*.csv'), emit: filtered_quant

    script:
    """
    gtf2bed < ${merged_gtf} > ${library}_merged.bed
    bedtools intersect -a ${library}_merged.bed -b ${pas_bed} -wa > ${library}_overlaps_a.bed
    bedtools intersect -a ${library}_merged.bed -b ${pas_bed} -wb > ${library}_overlaps_b.bed
    awk '{print \$1"\t"\$7"\t"\$8"\t"\$2"\t"\$3"\t"\$5"\t"\$6"\t"\$9"\t"(substr(\$0, index(\$0,\$10)))}' ${library}_overlaps_a.bed > ${library}_overlaps_a.gtf

    python ${projectDir}/modules/filter_quant_results.py -g ${library}_overlaps_a.gtf -b ${library}_overlaps_b.bed -q ${quant} -o ${library}_filtered_quant.csv
    """
}