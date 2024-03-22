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
include { STAR_INDEX_GENOME } from './modules/alignment.nf'
include { STAR_ALIGN_PE } from './modules/alignment.nf'
include { SAMTOOLS_INDEX } from './modules/samtools.nf'

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
    SAMTOOLS_INDEX(
        star_mapped_bam
    )
    index_bai = SAMTOOLS_INDEX.out.bai
    // TECTOOL(
    //     index_bai,
    //     star_mapped_bam,
    //     annotation_ch, 
    //     polya_ch,
    //     genome_ch
    // )
    }

/* 
 * completion handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}