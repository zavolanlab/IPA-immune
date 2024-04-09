import os
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import subprocess

pd.options.mode.chained_assignment = None

parser = ArgumentParser(description="Intron retention quantification arguments",
                        formatter_class=RawTextHelpFormatter)


# gtf file input
parser.add_argument('--a', type = str, required = True,
                    help = "The file path in gtf format of the annotation file")

# bam file input
parser.add_argument('--bam', type = str, required = True,
                    help = "The file path in bam format output by STAR")

# output table of Spliced Sites
parser.add_argument('--out', type = str, required = True,
                    help = "The output path for the table of the spliced sites")

# whether to filter multimappers
parser.add_argument('--f', help = "Specify this option to skip filtering multimappers",
                    action = "store_true")

# specify libType
### "SF" stands for "first read same strand"
### "SR" stands for "first read negative strand"
### input reads should be paired
parser.add_argument('--read_orientation', choices = ['SF', 'SR'], default = 'SR',
                    help = "Specify whether the first read is on the same strand or the opposite strand")

options = parser.parse_args()



# extract file name from path
filename = os.path.basename(options.bam)
filename = filename.split('.')
filename = filename[0]



# set gtf path and column names
gtf_path = options.a
gtf_columns = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

# read gtf file
gtf_df = pd.read_csv(gtf_path, sep = '\t', comment = '#',
                     header = None, names = gtf_columns)

# extract exons
exons = gtf_df[gtf_df['feature'] == 'exon']

# add transcript ID column
exons['transcript_id'] = exons.apply(lambda x:x['attribute'].split('transcript_id "')[1].split('";')[0],1)

# add exon number column
exons['exon_number'] = exons.apply(lambda x:x['attribute'].split('exon_number ')[1].split(';')[0],1)
exons['exon_number'] = pd.to_numeric(exons['exon_number'])

# add two bool columns to identify if exons are the first/last one of a certain transcript
min_exon = exons.groupby('transcript_id')['exon_number'].transform('min')
max_exon = exons.groupby('transcript_id')['exon_number'].transform('max')
exons['is_min'] = exons['exon_number'] == min_exon
exons['is_max'] = exons['exon_number'] == max_exon

# forward and reverse strand
forward = exons[exons['strand'] == '+']
reverse = exons[exons['strand'] == '-']



### extract the spliced sites ###

# forward and extracting start position
# delete the first exon in every transcript because they are the starting point of genes
forward_A = forward[forward['is_min'] == False]
forward_A = forward_A[['chr', 'start', 'strand']] # forward strands start with 3'SS (Acceptor)
forward_A.rename(columns = {'start': 'end'}, inplace = True)
forward_A['start'] = forward_A['end'] - 1
forward_A['score'] = 0
forward_A['name'] = forward_A['chr'] + '_' + forward_A['end'].astype(str) + '_A_K' # known
forward_A = forward_A[['chr', 'start', 'end', 'name', 'score', 'strand']]
forward_A = forward_A.drop_duplicates()

# forward and extracting end position
# delete the last exon in every transcript because they are the ending point of genes
forward_D = forward[forward['is_max'] == False]
forward_D = forward_D[['chr', 'end', 'strand']] # forward strands end with 5'SS (Donor)
forward_D['start'] = forward_D['end'] - 1
forward_D['score'] = 0
forward_D['name'] = forward_D['chr'] + '_' + forward_D['end'].astype(str) + '_D_K' # known
forward_D = forward_D[['chr', 'start', 'end', 'name', 'score', 'strand']]
forward_D = forward_D.drop_duplicates()

# reverse and extracting start position
# delete the last exon in every transcript because they are the ending point of genes
reverse_D = reverse[reverse['is_max'] == False]
reverse_D = reverse_D[['chr', 'start', 'strand']] # reverse strands start with D
reverse_D.rename(columns = {'start': 'end'}, inplace = True)
reverse_D['start'] = reverse_D['end'] - 1
reverse_D['score'] = 0
reverse_D['name'] = reverse_D['chr'] + '_' + reverse_D['end'].astype(str) + '_D_K' # known
reverse_D = reverse_D[['chr', 'start', 'end', 'name', 'score', 'strand']]
reverse_D = reverse_D.drop_duplicates()

# reverse and extracting end position
# delete the first exon in every transcript because they are the starting point of genes
reverse_A = reverse[reverse['is_min'] == False]
reverse_A = reverse_A[['chr', 'end', 'strand']] # reverse strands end with A
reverse_A['start'] = reverse_A['end'] - 1
reverse_A['score'] = 0
reverse_A['name'] = reverse_A['chr'] + '_' + reverse_A['end'].astype(str) + '_A_K' # known
reverse_A = reverse_A[['chr', 'start', 'end', 'name', 'score', 'strand']]
reverse_A = reverse_A.drop_duplicates()


# the result dataframe
gtf2bed = pd.concat([forward_A, forward_D, reverse_A, reverse_D])
gtf2bed = gtf2bed.sort_values(by = ['chr', 'start', 'name'])



# set SJ-like file path
sj_basename = filename + ".sj"
sj_path = os.path.join(options.out, sj_basename)

# this script is originally from
# "https://github.com/alexdobin/STAR/blob/master/extras/scripts/sjFromSAMcollapseUandM.awk"
awk_script = "sjFromSAM.awk"

# generate sj-like file from bam file
command = f"samtools view -h {options.bam} | awk -f {awk_script} | sort -V > {sj_path}"
subprocess.run(command, shell = True)

# set column names
sj_columns = ['chr', 'start', 'end', 'readReverse', 'mateReverse', 'firstInPair', 'secondInPair']

# read SJ-like file
sj_df = pd.read_csv(sj_path, sep = '\t', names = sj_columns)


# convert strand from numbers to '+' and '-'
def convert_strand_SR(row):
    readReverse = row['readReverse']
    mateReverse = row['mateReverse']
    firstInPair = row['firstInPair']
    secondInPair = row['secondInPair']
    forward1 = (readReverse==0) & (mateReverse==1) & (firstInPair==0) & (secondInPair==1)
    forward2 = (readReverse==1) & (mateReverse==0) & (firstInPair==1) & (secondInPair==0)
    reverse1 = (readReverse==0) & (mateReverse==1) & (firstInPair==1) & (secondInPair==0)
    reverse2 = (readReverse==1) & (mateReverse==0) & (firstInPair==0) & (secondInPair==1)
    if (forward1 | forward2):
        return '+'
    elif (reverse1 | reverse2):
        return '-'
    else:
        return 0

def convert_strand_SF(row):
    readReverse = row['readReverse']
    mateReverse = row['mateReverse']
    firstInPair = row['firstInPair']
    secondInPair = row['secondInPair']
    forward1 = (readReverse==0) & (mateReverse==1) & (firstInPair==1) & (secondInPair==0)
    forward2 = (readReverse==1) & (mateReverse==0) & (firstInPair==0) & (secondInPair==1)
    reverse1 = (readReverse==0) & (mateReverse==1) & (firstInPair==0) & (secondInPair==1)
    reverse2 = (readReverse==1) & (mateReverse==0) & (firstInPair==1) & (secondInPair==0)
    if (forward1 | forward2):
        return '+'
    elif (reverse1 | reverse2):
        return '-'
    else:
        return 0

if options.read_orientation=='SR':
    sj_df['strand'] = sj_df.apply(convert_strand_SR, axis=1)
elif options.read_orientation=='SF':
    sj_df['strand'] = sj_df.apply(convert_strand_SF, axis=1)


sj_df = sj_df[['chr', 'start', 'end', 'strand']] # keep only useful columns
sj_df = sj_df.drop_duplicates()

# mark undefined strands
duplicated = sj_df.duplicated(subset = ['chr', 'start', 'end'], keep = False)
sj_df.loc[duplicated, 'strand'] = 0

# discard the undefined strands
sj_df = sj_df[sj_df['strand'] != 0]

# forward and reverse strand
sj_forward = sj_df[sj_df['strand'] == '+']
sj_reverse = sj_df[sj_df['strand'] == '-']



### extract the spliced sites ###

# SJ.out.tab use the start and end of *intron*

# forward and extracting start position
sj_forward_D = sj_forward[['chr', 'start', 'strand']] # forward introns start with 5'SS (Donor)
sj_forward_D.rename(columns = {'start': 'end'}, inplace = True)
sj_forward_D['end'] = sj_forward_D['end'] - 1 # from the first position of intron to the last position of exon
sj_forward_D['start'] = sj_forward_D['end'] - 1
sj_forward_D['score'] = 0
sj_forward_D['name'] = sj_forward_D['chr'] + '_' + sj_forward_D['end'].astype(str) + '_D_U' # unknown
sj_forward_D = sj_forward_D[['chr', 'start', 'end', 'name', 'score', 'strand']]
sj_forward_D = sj_forward_D.drop_duplicates()

# forward and extracting end position
sj_forward_A = sj_forward[['chr', 'end', 'strand']] # forward introns end with 3'SS (Acceptor)
sj_forward_A['end'] = sj_forward_A['end'] + 1 # from the last position of intron to the first position of exon
sj_forward_A['start'] = sj_forward_A['end'] - 1
sj_forward_A['score'] = 0
sj_forward_A['name'] = sj_forward_A['chr'] + '_' + sj_forward_A['end'].astype(str) + '_A_U' # unknown
sj_forward_A = sj_forward_A[['chr', 'start', 'end', 'name', 'score', 'strand']]
sj_forward_A = sj_forward_A.drop_duplicates()

# reverse and extracting start position
sj_reverse_A = sj_reverse[['chr', 'start', 'strand']] # reverse introns start with 3'SS (Accepter)
sj_reverse_A.rename(columns = {'start': 'end'}, inplace = True)
sj_reverse_A['end'] = sj_reverse_A['end'] - 1 # from the first position of intron to the last position of exon
sj_reverse_A['start'] = sj_reverse_A['end'] - 1
sj_reverse_A['score'] = 0
sj_reverse_A['name'] = sj_reverse_A['chr'] + '_' + sj_reverse_A['end'].astype(str) + '_A_U' # unknown
sj_reverse_A = sj_reverse_A[['chr', 'start', 'end', 'name', 'score', 'strand']]
sj_reverse_A = sj_reverse_A.drop_duplicates()

# reverse and extracting end position
sj_reverse_D = sj_reverse[['chr', 'end', 'strand']] # reverse introns end with 5'SS (Donor)
sj_reverse_D['end'] = sj_reverse_D['end'] + 1 # from the last position of intron to the first position of exon
sj_reverse_D['start'] = sj_reverse_D['end'] - 1
sj_reverse_D['score'] = 0
sj_reverse_D['name'] = sj_reverse_D['chr'] + '_' + sj_reverse_D['end'].astype(str) + '_D_U' # unknown
sj_reverse_D = sj_reverse_D[['chr', 'start', 'end', 'name', 'score', 'strand']]
sj_reverse_D = sj_reverse_D.drop_duplicates()


# the result dataframe
sj2bed = pd.concat([sj_forward_A, sj_forward_D, sj_reverse_A, sj_reverse_D])
sj2bed = sj2bed.sort_values(by = ['chr', 'start', 'name'])



# the spliced sites dataframes have the same columns:
# bed_columns = ['chr', 'start', 'end', 'name', 'score', 'strand']

merged_df = pd.concat([gtf2bed, sj2bed])
merged_df = merged_df.sort_values(by = ['chr', 'start', 'name'])


# define whether the splice sites are known
def is_unknown(name):
    if name.endswith('_K'):
        return 0
    elif name.endswith('_U'):
        return 1
    return 2

merged_df['is_unknown'] = merged_df['name'].apply(is_unknown)


# delete duplicated rows
merged_df.sort_values(by = ['chr', 'start', 'end', 'is_unknown'], ascending = True, inplace = True)
merged_df.drop_duplicates(subset = ['chr', 'start', 'end', 'strand'], keep = 'first', inplace = True)
merged_df.drop('is_unknown', axis = 1, inplace = True)

# save merged spliced sites
merged_basename = filename + "_merged.bed"
merged_path = os.path.join(options.out, merged_basename)
merged_df.to_csv(merged_path, sep = '\t', index = False, header = False)



# filter multimappers
# convert bam to bed using "bedtools bamtobed"
# retain only useful information

### if "--f" is specified, multimapper filtering will be skipped
### please be sure that your input bam file is already filtered before you specify "--f"
### this is only to save computational efforts

filtered_reads_basename = filename + "_filteredReads.bed"
filtered_reads_path = os.path.join(options.out, filtered_reads_basename)

if options.f:
    command = f"bedtools bamtobed -split -i {options.bam} | awk 'BEGIN{{OFS=\"\\t\"}} {{$4=\"-\"; $5=0; print}}' > {filtered_reads_path}"
else:
    command = f"samtools view -h -q 255 {options.bam} | samtools view -h -b | bedtools bamtobed -split -i | awk 'BEGIN{{OFS=\"\\t\"}} {{$4=\"-\"; $5=0; print}}' > {filtered_reads_path}" 

subprocess.run(command, shell = True)



# read filtered and spliced sites file
bed_columns = ['chr', 'start', 'end', 'name', 'score', 'strand']
filtered_df = pd.read_csv(filtered_reads_path, sep = '\t', names = bed_columns)
filtered_df.sort_values(by = ['chr', 'start', 'end'], ascending = True, inplace = True)
SS_df = pd.read_csv(merged_path, sep = '\t', names = bed_columns)
SS_df.sort_values(by = ['chr', 'start', 'end'], ascending = True, inplace = True)

# groupby bed by chr, start, end, strand to get the number of reads mapped to the same genomic location
filtered_grouped = filtered_df.groupby(['chr', 'start', 'end', 'strand'])
# add 'count' column, saving the duplicated counts
filtered_countDup = filtered_grouped.size().reset_index(name = 'count')
# since name is not important, save count information in "name" column
filtered_countDup.rename(columns = {'count': 'name'}, inplace = True)
# add a "score" column to keep the format
filtered_countDup['score'] = 0
# reorder the column
filtered_countDup = filtered_countDup[['chr', 'start', 'end', 'name', 'score', 'strand']]
# save the filtered count file
filtered_count_basename = filename + "_filteredNameAsCount.bed"
filtered_count_path = os.path.join(options.out, filtered_count_basename)
filtered_countDup.to_csv(filtered_count_path, sep = '\t', index = False, header = False)


# bedtools intersect command
intersect_basename = filename + "_intersect.bed"
intersect_path = os.path.join(options.out, intersect_basename)
intersect_log_basename = filename + "_intersect.log"
intersect_log_path = os.path.join(options.out, intersect_log_basename)
command = f"bedtools intersect -wa -wb -s -a {merged_path} -b {filtered_count_path} -sorted 1>{intersect_path} 2>{intersect_log_path}"
subprocess.run(command, shell = True)

# read only the columns needed
intersect_columns = ['chr', 'SS_start', 'SS_end', 'SS_name', 'SS_score', 'strand', 'r_chr', 'r_start', 'r_end', 'r_num', 'r_score', 'r_strand']
col_interest = ['chr', 'SS_start', 'SS_end', 'SS_name', 'r_start', 'r_end', 'r_num', 'strand']
intersect_df = pd.read_csv(intersect_path, sep = '\t', header = None, names = intersect_columns, usecols = col_interest)

# define whether the spliced sites are known
intersect_df['is_unknown'] = intersect_df['SS_name'].apply(is_unknown)
intersect_df = intersect_df[['chr', 'SS_start', 'SS_end', 'is_unknown', 'strand', 'r_start', 'r_end', 'r_num']]


### intron retention ###

# find intron retention event
ir_df = intersect_df[(intersect_df['SS_start']-intersect_df['r_start']>=8)&(intersect_df['r_end']-intersect_df['SS_end']>=8)]
# count reads that support intron retention event
ir_df['count'] = ir_df.groupby(['chr', 'SS_start', 'SS_end', 'strand'])['r_num'].transform('sum')


### spliced reads ###

# reads spliced at SS
SS_reads = intersect_df[(intersect_df['SS_start']==intersect_df['r_start'])|(intersect_df['SS_end']==intersect_df['r_end'])]
# count reads that support spliced event
SS_reads['count'] = SS_reads.groupby(['chr', 'SS_start', 'SS_end', 'strand'])['r_num'].transform('sum')



# intron retention
ir_table = ir_df.drop_duplicates(subset = ['chr', 'SS_start', 'SS_end', 'is_unknown', 'strand', 'count'])
ir_table = ir_table[['chr', 'SS_start', 'SS_end', 'is_unknown', 'strand', 'count']]
ir_table.rename(columns = {'count': 'ir_count'}, inplace = True)

# reads that support spliced sites
SS_reads_table = SS_reads.drop_duplicates(subset = ['chr', 'SS_start', 'SS_end', 'is_unknown', 'strand', 'count'])
SS_reads_table = SS_reads_table[['chr', 'SS_start', 'SS_end', 'is_unknown', 'strand', 'count']]
SS_reads_table.rename(columns = {'count': 'SS_reads_count'}, inplace = True)

# merge two events
result = pd.merge(ir_table, SS_reads_table, how = "outer", on = ['chr', 'SS_start', 'SS_end', 'is_unknown', 'strand'])
result = result.fillna(0)
result[['ir_count', 'SS_reads_count']] = result[['ir_count', 'SS_reads_count']].astype(int)
# keep only useful columns
result = result[['chr', 'SS_end', 'is_unknown', 'strand', 'ir_count', 'SS_reads_count']]
result.rename(columns = {'SS_end': 'SS'}, inplace = True)

# save the result table
result_basename = filename + "_result.csv"
result_path = os.path.join(options.out, result_basename)
result.to_csv(result_path, index = False)