process TECTOOL {
    label "TECtool"
    container "docker://fgypas/tectool:0.4"

    publishDir params.outdir, mode:'copy'

    input:
    path sample_bam
    path annotation_gtf
    path polya_sites_bed
    path genome_fa

    output:
    path "results"

    script:
    """
    tectool \
        --bam $sample_bam \
        --annotation $annotation_gtf \
        --polyasites $polya_sites_bed \
        --genome $genome_fa \
        --output_dir results
    """
}

// process POST_RIBLAST {

//     publishDir params.outdir, mode:'copy'
    
//     input:
//     path riblast_output

//     output:
//     path "priming_sites.gtf", emit: gtf

//     script:
//     """
//     priming-site-predictor \
//         -r $riblast_output \
//         -o priming_sites.gtf \
//     """
// }