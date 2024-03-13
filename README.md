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

Update und activate the Conda environment with:

```bash
mamba env update -f install/environment-update.yml
conda activate ipa-immune
```

## Workflow

The workflow makes use of [TECtool](https://github.com/zavolanlab/TECtool).

Inputs:
1. Bulk RNA-Seq bam file
2. GTF annotation file
3. Genome FASTA file
4. BED file with IPA sites


Outputs:
1. Table: raw number of reads supporting each PAS in the sample

## TESTING

Download and uncompress testing data:

```bash
cd tests/
wget http://tectool.unibas.ch/data/test_data.tar.gz
tar xzvf test_data.tar.gz
```