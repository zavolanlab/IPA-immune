process INTRON_RETENTION{
        
    label 'intron_retention'
    
    tag { library }

    publishDir "${params.out_dir}", mode: 'copy', pattern: '*_intron_retention_matrix.csv'
    
    input:
    tuple val(library), path(bam)
    path annotation_gtf

    output:
    tuple val(library), path('*_intron_retention_matrix.csv'), emit: INTRON_RETENTION_MATRIX

    script:
    """
    python ${projectDir}/modules/intron_retention/IRworkflow.py -a ${annotation_gtf} -b ${bam} -o ${library}_intron_retention -f
    mv ${library}_intron_retention/*_result.csv ${library}_intron_retention_matrix.csv
    """
}