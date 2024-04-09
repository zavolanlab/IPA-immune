IRworkflow.py is used to extract the count of reads that supports (1) intron retention events (2) splicing at spliced sites.

```
python IRworkflow.py --a <gtf file> --bam <bam file> --out <output folder>
```

**Possible options**:

`--f`: If specify this option, the step of filtering multimappers will be skipped. By default, when not specifying this option, multimappers will be filtered. However, filtering is often recommended, and we offer this option in case the filtering process is done beforehand and you would like to save some computational efforts.

`--read_orientation`: Options are `SR` or `SF`, in which "SR" stands for first read being opposite strand, and "SF" stands for first read being the same strand. Default is `SR`.

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

myenv1.yml: the environment file.

sjFromSAM.awk: dependent awk script. Should be put under the same path as IRworkflow.py.
