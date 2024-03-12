process TECTOOL {
    label "TECtool"
    container "docker://fgypas/tectool:0.4"

    beforeScript 'export SINGULARITY_BIND="$PWD/tests/test_data"'
    publishDir "${projectDir}", mode:'copy'

    input:

    path sample_bai
    path sample_bam
    path annotation_gtf
    path polya_sites_bed
    path genome_fa

    script:
    """
    echo $sample_bai;
    tectool \\
        --annotation $annotation_gtf \\
        --polyasites $polya_sites_bed \\
        --bam $sample_bam \\
        --genome $genome_fa \\
        --output_dir $params.outdir
    """
}
