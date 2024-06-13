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
 """

// import modules
include { FASTQC as FASTQC_FASTQ } from './modules/fastqc.nf'
include { FASTQC as FASTQC_BAM } from './modules/fastqc.nf'
include { STAR_ALIGN_PE } from './modules/alignment.nf'
include { STAR_ALIGN_SE as ALIGN_FASTQ_1 } from './modules/alignment.nf'
include { STAR_ALIGN_SE as ALIGN_FASTQ_2 } from './modules/alignment.nf'
include { SAMTOOLS_GET_UNIQUE_MAPPERS } from './modules/samtools.nf'
include { SAMTOOLS_GET_LOW_DUP_READS } from './modules/samtools.nf'
include { SAMTOOLS_BAM2FASTQ } from './modules/samtools.nf'
include { TECTOOL as TECTOOL_1 } from './modules/tectool.nf'
include { TECTOOL as TECTOOL_2 } from './modules/tectool.nf'
include { TECTOOL_MERGE } from './modules/tectool.nf'
include { BEDTOOLS_FILTER_QUANT } from './modules/bedtools.nf'
include { SALMON_TRANSCRIPTOME } from './modules/salmon.nf'
include { SALMON_INDEX } from './modules/salmon.nf'
include { SALMON_QUANTIFY } from './modules/salmon.nf'
include { INTRON_RETENTION } from './modules/intron_retention.nf'

genome_index_ch = Channel.fromPath(params.genome_index, checkIfExists: true).collect()
annotation_gtf_ch = Channel.fromPath(params.annotation_gtf, checkIfExists: true).collect()
polya_sites_bed_ch = Channel.fromPath(params.polya_sites_bed, checkIfExists: true).collect()
genome_fa_ch = Channel.fromPath(params.genome_fa, checkIfExists: true).collect()

// Subworkflow for preprocessing steps
workflow preprocessing {
    take:
        input_fastq

    main:
        FASTQC_FASTQ(input_fastq)
        STAR_ALIGN_PE(input_fastq, genome_index_ch)
        star_mapped_bam_tuple = STAR_ALIGN_PE.out.star_mapped_bam_tuple
        SAMTOOLS_GET_UNIQUE_MAPPERS(star_mapped_bam_tuple)
        filtered_bam_tuple = SAMTOOLS_GET_UNIQUE_MAPPERS.out.filtered_bam_tuple
        SAMTOOLS_GET_LOW_DUP_READS(filtered_bam_tuple)
        bam_low_dupl_tupl = SAMTOOLS_GET_LOW_DUP_READS.out.bam_low_dupl_tupl
        FASTQC_BAM(bam_low_dupl_tupl)

    emit:
        bam_low_dupl_tupl
}

// Subworkflow for TECtool analysis and downstream steps
workflow tectool_analysis {
    take:
        input_bam

    main:
        // Convert BAM to 2 FASTQ files
        SAMTOOLS_BAM2FASTQ(input_bam)
        fastq1_tuple = SAMTOOLS_BAM2FASTQ.out.fastq1_tuple
        fastq2_tuple = SAMTOOLS_BAM2FASTQ.out.fastq2_tuple

        // FASTQ1 from BAM
        ALIGN_FASTQ_1(fastq1_tuple, genome_index_ch)
        star_mapped_bam_1 = ALIGN_FASTQ_1.out.star_mapped_bam
        TECTOOL_1(star_mapped_bam_1, annotation_gtf_ch, polya_sites_bed_ch, genome_fa_ch)
        enriched_gtf_1 = TECTOOL_1.out.enriched_gtf

        // FASTQ2 from BAM
        ALIGN_FASTQ_2(fastq2_tuple, genome_index_ch)
        star_mapped_bam_2 = ALIGN_FASTQ_2.out.star_mapped_bam
        TECTOOL_2(star_mapped_bam_2, annotation_gtf_ch, polya_sites_bed_ch, genome_fa_ch)
        enriched_gtf_2 = TECTOOL_2.out.enriched_gtf

        // MERGE 2 enriched GTFs with TECtool script
        TECTOOL_MERGE(input_bam, enriched_gtf_1, enriched_gtf_2)
        merged_gtf = TECTOOL_MERGE.out.merged_gtf

        SALMON_TRANSCRIPTOME(merged_gtf, genome_fa_ch)
        transcriptome_fa = SALMON_TRANSCRIPTOME.out.transcriptome_fa
        SALMON_INDEX(transcriptome_fa)
        salmon_index = SALMON_INDEX.out.salmon_index
        SALMON_QUANTIFY(salmon_index, fastq1_tuple, fastq2_tuple)
        salmon_counts = SALMON_QUANTIFY.out.salmon_counts

    emit:
        salmon_counts
}

workflow intron_retention {
    take:
        input_bam

    main:
        INTRON_RETENTION(input_bam, annotation_gtf_ch)
        intron_retention = INTRON_RETENTION.out.INTRON_RETENTION_MATRIX

    emit:
        intron_retention
}

// Main workflow
workflow {
    if (params.run_mode == 'full') {
        input_fastq_ch = Channel.fromFilePairs(params.input_fastq, checkIfExists: true)
        input_fastq_ch.each {
            preprocessing(input_fastq_ch)
            tectool_analysis(preprocessing.out.bam_low_dupl_tupl)
            intron_retention(preprocessing.out.bam_low_dupl_tupl)
        }
    }
    if (params.run_mode == 'preprocessing') {
        input_fastq_ch = Channel.fromFilePairs(params.input_fastq, checkIfExists: true)
        input_fastq_ch.each {
            preprocessing(input_fastq_ch)
        }
    }
    if (params.run_mode == 'analysis') {
        input_bam_ch = Channel.fromPath(params.input_bam, checkIfExists: true).map { bam_path -> tuple(bam_path.baseName, bam_path) }
        input_bam_ch.each {
            tectool_analysis(input_bam_ch)
            intron_retention(input_bam_ch)
        }
    }
    if (params.run_mode == 'tectool') {
        input_bam_ch = Channel.fromPath(params.input_bam, checkIfExists: true).map { bam_path -> tuple(bam_path.baseName, bam_path) }
        input_bam_ch.each {
            tectool_analysis(input_bam_ch)
        }
    }
    if (params.run_mode == 'intron') {
        input_bam_ch = Channel.fromPath(params.input_bam, checkIfExists: true).map { bam_path -> tuple(bam_path.baseName, bam_path) }
        input_bam_ch.each {
            intron_retention(input_bam_ch)
        }
    }
}

/* 
 * completion handler
 */
workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}