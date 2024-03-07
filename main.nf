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
 sample: ${params.sample_bam}
 annotation GTF: ${params.annotation_gtf}
 polyA sites BED: ${params.polya_sites_bed}
 genome FASTA: ${params.genome_fa}
 outdir       : ${params.outdir}
 """

// import modules
include { TECTOOL } from './modules/tectool.nf'


// Define inputs
sample_ch = Channel.fromPath( params.sample_bam )
annotation_ch = Channel.fromPath( params.annotation_gtf )
polya_ch = Channel.value( params.polya_sites_bed )
genome_ch = Channel.value( params.genome_fa )

/* 
 * main script flow
 */
workflow {
    TECTOOL( sample_ch, annotation_ch, polya_ch, genome_ch )
    }

/* 
 * completion handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}