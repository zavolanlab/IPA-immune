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