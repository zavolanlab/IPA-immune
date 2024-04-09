# IPA-derived proteins expand the MHC-I-restricted immunopeptidome
A Nextflow pipeline for the multi-omics analysis of IPA-derived isoforms in cancer.

## Installation

### 1. Clone the repository

Go to the desired directory/folder on your file system, then clone/get the 
repository and move into the respective directory with:

```bash
git clone https://github.com/zavolanlab/IPA-immune
cd IPA-immune
```

### 2. Conda and Mamba installation

Workflow dependencies can be conveniently installed with the [Conda](https://docs.conda.io/projects/conda/en/stable/)
package manager. We recommend that you install [Miniconda](https://docs.anaconda.com/free/miniconda/miniconda-install/)
for your system (Linux). Be sure to select the Python 3 option. 

```bash
conda install -y mamba -n base -c conda-forge
```

### 3. Create environment

Install the remaining dependencies with:
```bash
mamba env create -f install/environment.yml
```

### 4. Activate environment

Activate the Conda environment with:

```bash
conda activate ipa-immune
```

## Workflow

The workflow makes use of [TECtool](https://github.com/balajtimate/TECtool).

Inputs:
1. Bulk RNA-Seq FASTQ files
2. GTF annotation file
3. Genome FASTA file
4. BED file with IPA sites


Outputs:
1. Table: raw number of reads supporting each PAS in the sample

## Running the workflow

To start the workflow, in your activated `ipa-immune` conda environment, run

```bash
nextflow main.nf -profile conda
```
> Currently, running the workflow is only supported with conda

For running on SLURM:
```bash
nextflow main.nf -profile slurm,conda
```

Unzip the test FASTQ files:
```bash
zcat /scicore/home/zavolan/GROUP/IPA/data/bulk_RNA_seq/test/ENCFF184CDV_ENCFF456OPJ.R1.fastq.gz > tests/ENCFF184CDV_ENCFF456OPJ.R1.fastq
zcat /scicore/home/zavolan/GROUP/IPA/data/bulk_RNA_seq/test/ENCFF184CDV_ENCFF456OPJ.R2.fastq.gz > tests/ENCFF184CDV_ENCFF456OPJ.R2.fastq
```

Running the workflow with setting the params for test files from sciCORE:
```bash
nextflow main.nf -profile slurm,conda \
    --input_fastq tests/ENCFF184CDV_*.R{1,2}.fastq \
    --genome_index /scicore/home/zavolan/GROUP/RBP_perturbational_networks/wf_runs/v1/output/STAR_indices/human/without_GTF/STAR_index \
    --annotation_gtf /scicore/home/zavolan/GROUP/Genomes/homo_sapiens/hg38_v42/gencode.v42.annotation.gtf \
    --polya_sites_bed /scicore/home/zavolan/GROUP/IPA/IPA_catalogue/SCINPAS_all_normal_q15Expr.bed \
    --genome_fa /scicore/home/zavolan/GROUP/Genomes/homo_sapiens/GRCh38.primary_assembly.genome.fa
```