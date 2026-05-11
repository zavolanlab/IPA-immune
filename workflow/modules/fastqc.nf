#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process FASTQC {
    tag { meta.id }
    cpus { params.threads_se }

    publishDir "${params.out_dir}/${meta.id}_results", mode: 'copy', pattern: '*_fastqc.html'
    publishDir "${params.log_dir}/${meta.id}_logs",    mode: 'copy', pattern: '*.inputs.log'

    input:
    tuple val(meta), path(reads)   // reads may be a single Path or a List<Path>

    output:
    tuple val(meta), path('*_fastqc.html'), emit: html
    tuple val(meta), path('*_fastqc.zip'),  emit: zip

    script:
    def files = (reads instanceof List) ? reads : [reads]
    def joined = files.collect { it.toString() }.join(' ')
    """
    set -euo pipefail

    # simple, robust debug log (no here-doc)
    printf 'Sample ID  : %s\n' "${meta.id}" > ${meta.id}.inputs.log
    printf 'Single-end : %s\n' "${meta.single_end}" >> ${meta.id}.inputs.log
    printf 'Reads (n=%s):\n' "\${#@}" >> ${meta.id}.inputs.log
    for f in ${joined}; do echo "\$f" >> ${meta.id}.inputs.log; done

    fastqc --threads ${task.cpus} ${joined} --outdir .
    """
}

process FASTQC_BAM {
    tag { meta.id }
    cpus { params.threads_se }

    publishDir "${params.out_dir}/${meta.id}_results", mode: 'copy', pattern: '*_fastqc.html'
    publishDir "${params.log_dir}/${meta.id}_logs",    mode: 'copy', pattern: '*.inputs.log'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*_fastqc.html'), emit: html
    tuple val(meta), path('*_fastqc.zip'),  emit: zip

    script:
    """
    set -euo pipefail

    {
        echo "Sample ID  : ${meta.id}"
        echo "BAM        : ${bam}"
    } > ${meta.id}.bam.inputs.log

    # FastQC can read BAM; if your version complains, switch to: samtools fastq "${bam}" | fastqc -o .
    fastqc --threads ${task.cpus} ${bam} --outdir .
    """
}

process MULTIQC_FASTQ {
  label 'multiqc'
  cpus { params.multiqc_cpus ?: 2 }
  tag { 'multiqc_fastq' }

  publishDir "${params.out_dir}/multiqc_fastq", mode: 'copy', pattern: '*.html'

  input:
  path(qc_files)   // a List<Path> from collect()

  output:
  path('multiqc_report.html'), emit: report
  path('multiqc_data'),        emit: data

  script:
  """
  multiqc --force --outdir . .
  """
}

process MULTIQC_BAM {
  label 'multiqc'
  cpus { params.multiqc_cpus ?: 2 }
  tag { 'multiqc_bam' }

  publishDir "${params.out_dir}/multiqc_bam", mode: 'copy', pattern: '*.html'

  input:
  path(qc_files)   // a List<Path> from collect()

  output:
  path('multiqc_report.html'), emit: report
  path('multiqc_data'),        emit: data

  script:
  """
  multiqc --force --outdir . .
  """
}
