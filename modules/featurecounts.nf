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
    -t CDS \\
    -g gene_id \\
    -s ${strand} \\
    ${peFlags} \\
    -T ${task.cpus} \\
    --verbose \\
    -o ${meta.id}.featureCounts.txt \\
    ${bam}
  """
}
