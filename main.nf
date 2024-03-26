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
 samples: ${params.input_fastq}
 annotation GTF: ${params.annotation_gtf}
 polyA sites BED: ${params.polya_sites_bed}
 genome FASTA: ${params.genome_fa}
 outdir       : ${params.out_dir}
 """

// import modules
include { FASTQC } from './modules/fastqc.nf'
// include { STAR_INDEX_GENOME } from './modules/alignment.nf'
include { STAR_ALIGN_PE } from './modules/alignment.nf'
include { STAR_ALIGN_SE as ALIGN_FASTQ_1 } from './modules/alignment.nf'
include { STAR_ALIGN_SE as ALIGN_FASTQ_2 } from './modules/alignment.nf'
include { STAR_ALIGN_SE as ALIGN_FASTQ_3 } from './modules/alignment.nf'
include { SAMTOOLS_INDEX as INDEX_BAM_1 } from './modules/samtools.nf'
include { SAMTOOLS_INDEX as INDEX_BAM_2 } from './modules/samtools.nf'
include { SAMTOOLS_INDEX as INDEX_BAM_3 } from './modules/samtools.nf'
include { SAMTOOLS_GET_UNIQUE_MAPPERS } from './modules/samtools.nf'
include { SAMTOOLS_GET_LOW_DUP_READS } from './modules/samtools.nf'
include { SAMTOOLS_BAM2FASTQ } from './modules/samtools.nf'
include { TECTOOL as TECTOOL_1 } from './modules/tectool.nf'
include { TECTOOL as TECTOOL_2 } from './modules/tectool.nf'
include { TECTOOL as TECTOOL_3 } from './modules/tectool.nf'

input_fastq_ch = Channel.fromFilePairs(params.input_fastq)
genome_index_ch = channel.fromPath(params.genome_index)

/* 
 * main script flow
 */
workflow {
    // STAR_INDEX_GENOME(
    //     params.genome_fa
    // )
    // genome_index = STAR_INDEX_GENOME.out.index
    FASTQC(
        input_fastq_ch
    )
    STAR_ALIGN_PE(
        input_fastq_ch,
        genome_index_ch
    )
    star_mapped_bam_tuple = STAR_ALIGN_PE.out.star_mapped_bam_tuple
    SAMTOOLS_GET_UNIQUE_MAPPERS(
        star_mapped_bam_tuple
    )
    filtered_bam_tuple = SAMTOOLS_GET_UNIQUE_MAPPERS.out.filtered_bam_tuple
    SAMTOOLS_GET_LOW_DUP_READS(
        filtered_bam_tuple
    )
    bam_low_dupl_tupl = SAMTOOLS_GET_LOW_DUP_READS.out.bam_low_dupl_tupl
    
    // Convert BAM to 2 FASTQ file     
    SAMTOOLS_BAM2FASTQ(bam_low_dupl_tupl)
    fastq1_tuple = SAMTOOLS_BAM2FASTQ.out.fastq1_tuple
    fastq2_tuple = SAMTOOLS_BAM2FASTQ.out.fastq2_tuple
    singleton_tuple = SAMTOOLS_BAM2FASTQ.out.singleton_tuple
    
    // FASTQ1 from BAM
    ALIGN_FASTQ_1(
        fastq1_tuple,
        genome_index_ch
    )
    star_mapped_bam_1 = ALIGN_FASTQ_1.out.star_mapped_bam
    INDEX_BAM_1(
        star_mapped_bam_1
    )
    star_mapped_bam_index_1 = INDEX_BAM_1.out.index
    TECTOOL_1(
        star_mapped_bam_index_1,
        star_mapped_bam_1,
        params.annotation_gtf, 
        params.polya_sites_bed,
        params.genome_fa
    )
    // FASTQ2 from BAM
    ALIGN_FASTQ_2(
        fastq2_tuple,
        genome_index_ch
    )
    star_mapped_bam_2 = ALIGN_FASTQ_2.out.star_mapped_bam
    INDEX_BAM_2(
        star_mapped_bam_2
    )
    star_mapped_bam_index_2 = INDEX_BAM_2.out.index
    TECTOOL_2(
        star_mapped_bam_index_2,
        star_mapped_bam_2,
        params.annotation_gtf, 
        params.polya_sites_bed,
        params.genome_fa
    )
    // FASTQ_SINGLETON from BAM
    ALIGN_FASTQ_3(
        singleton_tuple,
        genome_index_ch
    )
    star_mapped_bam_singleton = ALIGN_FASTQ_3.out.star_mapped_bam
    INDEX_BAM_3(
        star_mapped_bam_singleton
    )
    star_mapped_bam_index_singleton = INDEX_BAM_3.out.index
    TECTOOL_3(
        star_mapped_bam_index_singleton,
        star_mapped_bam_singleton,
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
