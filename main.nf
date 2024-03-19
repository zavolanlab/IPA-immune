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
 samples: ${params.samples_fastq}
 STAR genome index: ${params.star_genome_dir}
 annotation GTF: ${params.tec_annotation_gtf}
 polyA sites BED: ${params.tec_polya_sites_bed}
 genome FASTA: ${params.tec_genome_fa}
 outdir       : ${params.out_dir}
 """

// import modules
include { MAP_STAR } from './modules/alignment.nf'
include { INDEX_SAMTOOLS } from './modules/alignment.nf'
include { TECTOOL } from './modules/tectool.nf'


// Define inputs
samples_fastq_ch = Channel.fromPath(params.samples_fastq, checkIfExists: true)
index_ch = Channel.fromPath (params.star_genome_dir, checkIfExists: true)
annotation_ch = Channel.fromPath(params.tec_annotation_gtf, checkIfExists: true)
polya_ch = Channel.value(params.tec_polya_sites_bed)
genome_ch = Channel.value(params.tec_genome_fa)

/* 
 * main script flow
 */
workflow {
    MAP_STAR(
        samples_fastq_ch,
        index_ch
     )
    star_mapped_bam = MAP_STAR.out.aligned
    INDEX_SAMTOOLS(
        star_mapped_bam
    )
    index_bai = INDEX_SAMTOOLS.out.bai
    TECTOOL(
        index_bai,
        star_mapped_bam,
        annotation_ch, 
        polya_ch,
        genome_ch
    )
    }

/* 
 * completion handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}