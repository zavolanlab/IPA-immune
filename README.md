# IPA-immune
A Nextflow pipeline for the multi-omics analysis of IPA-derived isoforms in cancer.
> DISCLAIMER: The workflow is currently in development and still in an experimental stage

### Installation

#### 1. Clone the repository

Go to the desired directory/folder on your file system, then clone/get the 
repository and move into the respective directory with:

```bash
git clone https://github.com/zavolanlab/IPA-immune
cd IPA-immune
```

#### 2. Conda and Mamba installation

Workflow dependencies can be conveniently installed with the [Conda](https://docs.conda.io/projects/conda/en/stable/)
package manager. We recommend that you install [Miniconda](https://docs.anaconda.com/free/miniconda/miniconda-install/)
for your system (Linux). Be sure to select the Python 3 option. 

```bash
conda install -y mamba -n base -c conda-forge
```

#### 3. Create environment

Install the remaining dependencies with:
```bash
mamba env create -f install/environment.yml
```

#### 4. Activate environment

Activate the Conda environment with:

```bash
conda activate ipa-immune
```

### Workflow

The workflow makes use of 3 separate subworkflows:
1. Global quantification of intronic PAS usage (using [TECtool](https://github.com/balajtimate/TECtool))
2. Local quantification of intronic PAS
3. Local quantification of intron retention at splice sites

Inputs:
1. FASTQ bulk RNA-Seq files (paired) or BAM files
2. GTF annotation file 
3. FASTA genome file
4. BED file with IPA sites
5. STAR INDEX directory


Outputs:
1. Table: raw number of reads supporting each PAS in the sample

### Running the workflow
#### You have the choice of running the workflow in different configurations: 
(substitute one of the below options for the `<run_mode>` argument)

- `full`: default, to run the full workflow (this is computationally quite heavy and should be done in a cluster environment) - requires `input_fastq`
- `preprocessing`: to only run the preprocessing part of the workflow (alignment and filtering of low duplicate reads) - requires `input_fastq`
- `analysis`: to only run the postprocessing part of the workflow (quantification of IPA usage and intron retention) - requires `input_bam`
- `tectool`: to only run the IPA usage quantification, using TECtool - requires `input_bam`
- `intron`: to only run the intron retention quantification subworkflow - requires `input_bam`

In the case of `full` and `preprocessing` modes, the `input_fastq` is required, using a wildcard character, e.g.: `--input_fastq='test_data\*{1,2}.fastq'`

In the case of `analysis`, `tectool` and `intron` modes, the `input_bam` is required, using a wildcard character, e.g.: `--input_bam='test_data\*.bam'`


#### To start the workflow, in your activated `ipa-immune` conda environment, run

```bash
nextflow main.nf -profile conda 
    --genome_index <genome_index>
    --annotation_gtf <annotation_gtf>
    --polya_sites_bed <polya_sites_bed>
    --genome_fa <genome_fa>
    --out_dir <out_dir>
    --logs_dir <logs_dir>
    --run_mode <run_mode>
    --input_fastq <input_fastq> \ --input_bam <input_bam>
```
> Currently, running the workflow is only supported with conda

For running on SLURM:
```bash
nextflow main.nf -profile slurm,conda 
    --genome_index <genome_index>
    --annotation_gtf <annotation_gtf>
    --polya_sites_bed <polya_sites_bed>
    --genome_fa <genome_fa>
    --out_dir <out_dir>
    --logs_dir <logs_dir>
    --run_mode <run_mode>
    --input_fastq <input_fastq> \ --input_bam <input_bam>
```

### Testing

Download and uncompress testing data:

```bash
cd tests/
wget http://tectool.unibas.ch/data/test_data.tar.gz
tar xzvf test_data.tar.gz
```

## Subworkflow descriptions:
### Intron retention (by [@Zhihan Zhu](https://github.com/Jade0904))
IRworkflow.py is used to extract the count of reads that supports (1) intron retention events (2) splicing at spliced sites.

```
python IRworkflow.py --annotation <gtf file> --bam <bam file> --out <output folder>
```

**Possible options**:

`--filtered_input`, `-f`: If specify this option, the step of filtering multimappers will be skipped. By default, when not specifying this option, multimappers will be filtered. However, filtering is often recommended, and we offer this option in case the filtering process is done beforehand and you would like to save some computational efforts.

`--read_orientation`, `-r`: Options are `SR` or `SF`, in which "SR" stands for first read being opposite strand, and "SF" stands for first read being the same strand. Default is `SR`.

**Inputs** are:

(1) gtf file: annotation file in gtf format.

(2) bam file: output from STAR.

(3) output folder: for example "xx/xx/out".

**Outputs** are:

(1) An sj-like file extracting the start/end position of introns and their strand (similar to SJ.out.tab);

(2) All the spliced sites from both gtf file and sj file (merged.bed);

(3) The filtered bed file from bam file (filteredReads.bed);

(4) The filtered bed flle "groupby" the position of the reads, using the "name" column to save the count of the same reads (filteredNameAsCount.bed);

(5,6) The output of "bedtools intersect" (intersect.bed & intersect.log);

**(7) The result table containing the count of the reads of both events (result.csv). (This is the main result!)**

*Other files:*

sjFromSAM.awk: dependent awk script. Should be put under the same path as IRworkflow.py.
