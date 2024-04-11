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
 Samples: ${params.input_fastq}
 Genome index: ${params.genome_index}
 Annotation GTF: ${params.annotation_gtf}
 PolyA sites BED: ${params.polya_sites_bed}
 Genome FASTA: ${params.genome_fa}
 Output DIR: ${params.out_dir}
 """

// import modules
include { FASTQC as FASTQC_FASTQ } from './modules/fastqc.nf'
include { FASTQC as FASTQC_BAM } from './modules/fastqc.nf'
// include { STAR_INDEX_GENOME } from './modules/alignment.nf'
include { STAR_ALIGN_PE } from './modules/alignment.nf'
include { STAR_ALIGN_SE as ALIGN_FASTQ_1 } from './modules/alignment.nf'
include { STAR_ALIGN_SE as ALIGN_FASTQ_2 } from './modules/alignment.nf'
include { SAMTOOLS_INDEX as INDEX_BAM_1 } from './modules/samtools.nf'
include { SAMTOOLS_INDEX as INDEX_BAM_2 } from './modules/samtools.nf'
include { SAMTOOLS_GET_UNIQUE_MAPPERS } from './modules/samtools.nf'
include { SAMTOOLS_GET_LOW_DUP_READS } from './modules/samtools.nf'
include { SAMTOOLS_BAM2FASTQ } from './modules/samtools.nf'
include { TECTOOL as TECTOOL_1 } from './modules/tectool.nf'
include { TECTOOL_RENAME as TECTOOL_RENAME_1 } from './modules/tectool.nf'
include { TECTOOL_RENAME as TECTOOL_RENAME_2 } from './modules/tectool.nf'
include { TECTOOL as TECTOOL_2 } from './modules/tectool.nf'
include { BEDTOOLS_INTERSECT } from './modules/bedtools.nf'
include { STRINGTIE_QUANTIFY } from './modules/stringtie.nf'
include { STRINGTIE_COUNT_MATRIX } from './modules/stringtie.nf'

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
    FASTQC_FASTQ(
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
    FASTQC_BAM(
        bam_low_dupl_tupl
    )
    // Convert BAM to 2 FASTQ file     
    SAMTOOLS_BAM2FASTQ(bam_low_dupl_tupl)
    fastq1_tuple = SAMTOOLS_BAM2FASTQ.out.fastq1_tuple
    fastq2_tuple = SAMTOOLS_BAM2FASTQ.out.fastq2_tuple
    
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
    enriched_gtf_1 = TECTOOL_1.out.enriched_gtf
    TECTOOL_RENAME_1(
        enriched_gtf_1
    )
    renamed_gtf_1 = TECTOOL_RENAME_1.out.renamed_gtf
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
    enriched_gtf_2 = TECTOOL_2.out.enriched_gtf
    TECTOOL_RENAME_2(
        enriched_gtf_2
    )
    renamed_gtf_2 = TECTOOL_RENAME_2.out.renamed_gtf
    // MERGE 2 enriched GTFs and run STRINGTIE on original BAM
    BEDTOOLS_INTERSECT(
        renamed_gtf_1,
        renamed_gtf_2
    )
    intersect_gtf = BEDTOOLS_INTERSECT.out.intersect_gtf
    STRINGTIE_QUANTIFY(
        bam_low_dupl_tupl,
        intersect_gtf
    )
    stringtie_gtf = STRINGTIE_QUANTIFY.out.stringtie_gtf
    STRINGTIE_COUNT_MATRIX(stringtie_gtf)
}

/* 
 * completion handler
 */
workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}
