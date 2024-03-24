/* 
 * enables modules 
 */
nextflow.enable.dsl = 2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

log.info """\
 IPA-IMMUNE - N F   P I P E L I N E
 ===================================
 samples: ${params.input_fastq_1}, ${params.input_fastq_2}
 annotation GTF: ${params.annotation_gtf}
 polyA sites BED: ${params.polya_sites_bed}
 genome FASTA: ${params.genome_fa}
 outdir       : ${params.out_dir}
 """

// import modules
// include { STAR_INDEX_GENOME } from './modules/alignment.nf'
include { STAR_ALIGN_PE } from './modules/alignment.nf'
include { STAR_ALIGN_SE } from './modules/alignment.nf'

include { SAMTOOLS_INDEX } from './modules/samtools.nf'
include { SAMTOOLS_GET_UNIQUE_MAPPERS } from './modules/samtools.nf'
include { SAMTOOLS_GET_LOW_DUP_READS } from './modules/samtools.nf'
include { SAMTOOLS_FASTQ } from './modules/samtools.nf'

include { TECTOOL } from './modules/tectool.nf'


/* 
 * main script flow
 */
workflow {
    STAR_INDEX_GENOME(
        params.genome_fa
    )
    genome_index = STAR_INDEX_GENOME.out.index
    STAR_ALIGN_PE(
        params.input_fastq_1,
        params.input_fastq_2,
        genome_index
     )
    star_mapped_bam = STAR_ALIGN_PE.out.aligned
    SAMTOOLS_GET_UNIQUE_MAPPERS(
        star_mapped_bam
    )
    filtered_bam = SAMTOOLS_GET_UNIQUE_MAPPERS.out.filtered_bam
    SAMTOOLS_GET_LOW_DUP_READS(
        filtered_bam
    )
    bam_low_dupl = SAMTOOLS_GET_LOW_DUP_READS.out.bam_low_dupl
    SAMTOOLS_INDEX(
        bam_low_dupl
    )
    bam_low_dupl_index = SAMTOOLS_INDEX.out.bai
    SAMTOOLS_FASTQ(
        bam_low_dupl
    )
    1_fastq = SAMTOOLS_FASTQ.out.1_fastq
    2_fastq = SAMTOOLS_FASTQ.out.2_fastq
    // Combine these two into channels
    STAR_ALIGN_SE(
        1_fastq,
        params.star_genome_dir
    )
    1_star_mapped_bam = STAR_ALIGN_PE.out.aligned
    SAMTOOLS_INDEX(
        1_star_mapped_bam
    )
    1_star_mapped_bam_index = SAMTOOLS_INDEX.out.bai
        STAR_ALIGN_SE(
        2_fastq,
        params.star_genome_dir
    )
    2_star_mapped_bam = STAR_ALIGN_PE.out.aligned
    SAMTOOLS_INDEX(
        2_star_mapped_bam
    )
    2_star_mapped_bam_index = SAMTOOLS_INDEX.out.bai
    TECTOOL(
        1_star_mapped_bam_index,
        1_star_mapped_bam,
        params.annotation_gtf, 
        params.polya_sites_bed,
        params.genome_fa
    )
    TECTOOL(
        2_star_mapped_bam_index,
        2_star_mapped_bam,
        params.annotation_gtf, 
        params.polya_sites_bed,
        params.genome_fa
    )
}

/* 
 * completion handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}