#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process FEATURECOUNTS_GENE_CDS {
  label 'featurecounts'
  tag { meta.id }

  // match your conventions
  publishDir "${params.out_dir}/${meta.id}_results", mode: 'copy', pattern: "${meta.id}.featureCounts.txt*"
  publishDir "${params.log_dir}/${meta.id}_logs",    mode: 'copy', pattern: '*.inputs.log'

  // allocate CPUs based on SE/PE (so -T matches scheduler)
  cpus { meta.single_end ? (params.threads_se as int) : (params.threads_pe as int) }

  input:
  tuple val(meta), path(bam)
  path gtf

  output:
  tuple val(meta), path("${meta.id}.featureCounts.txt"),           emit: counts
  tuple val(meta), path("${meta.id}.featureCounts.txt.summary"),   emit: summary
  path "*.log", emit: log

  script:
  // feature type: CDS or exon
  def feature_type = (params.fc_feature_type ?: 'CDS') as String
  // strand: 0 unstranded (default), 1 stranded, 2 reverse-stranded
  def strand = (params.fc_strand ?: 0) as int
  // paired-end flags when not single-end
  def peFlags = meta.single_end ? '' : '-p --countReadPairs'
  """
  set -euo pipefail

  {
    echo "Sample ID   : ${meta.id}"
    echo "BAM         : ${bam}"
    echo "GTF         : ${gtf}"
    echo "single_end  : ${meta.single_end}"
    echo "strand (-s) : ${strand}"
    echo "cpus (-T)   : ${task.cpus}"
  } > ${meta.id}.featureCounts.inputs.log

  featureCounts \\
    -F GTF \\
    -a ${gtf} \\
    -t ${feature_type} \\
    -g gene_id \\
    -s ${strand} \\
    ${peFlags} \\
    -T ${task.cpus} \\
    --verbose \\
    -o ${meta.id}.featureCounts.txt \\
    ${bam}
  """
}


/*
 * Gene-level quantification with RNA-SeQC (v2; conda package `rnaseqc`)
 * We run in strict mode (proper pairs, exonic only, edit distance ≤6, etc.).
 * We treat libraries as *unstranded* here (to mirror GTEx), regardless of featureCounts setting.
 *
 * Input : (meta, bam), gtf
 * Output: standardized gene counts / TPM (copied from RNA-SeQC output dir)
 */
process RNASEQC_COUNT {
  label 'rnaseqc'
  tag   { meta.id }

  // keep threads consistent with the rest of your pipeline
  cpus { meta.single_end ? (params.threads_se as int) : (params.threads_pe as int) }

  conda "${HOME}/miniconda3/envs/rnaseqc"

  // land standardized outputs with sample id in results; keep raw dir too
  publishDir "${params.out_dir}/${meta.id}_results", mode: 'copy', pattern: "${meta.id}.rnaseqc.dir"

  input:
  tuple val(meta), path(bam)
  path gtf

  output:
  tuple val(meta), path("${meta.id}.rnaseqc.gene_counts.gct"), emit: counts,  optional: true
  tuple val(meta), path("${meta.id}.rnaseqc.gene_tpm.gct"),    emit: tpm,     optional: true
  path "${meta.id}.rnaseqc.metrics.txt",                       emit: metrics, optional: true
  path "${meta.id}.rnaseqc.dir",                               emit: outdir

  script:
  // enforce unstranded here (as in GTEx); override-able via param if you want later
  def outdir = "${meta.id}.rnaseqc.dir"

  """
  set -euo pipefail

  # Run RNA-SeQC v2
  rnaseqc \\
    --sample ${meta.id} \\
    ${gtf} \\
    ${bam} \\
    ${outdir}

  """
}
