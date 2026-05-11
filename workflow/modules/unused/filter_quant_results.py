import argparse
import gtfparse
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter exons based on specified files.")
    parser.add_argument("-g", "--gtf_file", help="Path to the GTF file.")
    parser.add_argument("-b", "--bed_file", help="Path to the BED file.")
    parser.add_argument("-q", "--quant_file", help="Path to the quant file.")
    parser.add_argument("-o", "--output_file", help="Path to the output CSV file.")
    return parser.parse_args()

def parse_bed_with_pandas(bed_file):
    bed_df = pd.read_csv(
        bed_file, sep='\t', index_col=False, header=None, names=[
            'seqname', 'start', 'end', 'gene_id', 'score',
            'strand', 'source', 'feature', 'frame', 'attributes',
            'pas_chr', 'pas_start', 'pas_end',
            'pas_name', 'pas_score', 'pas_strand'
        ]
    )
    return bed_df[['pas_chr', 'pas_start', 'pas_end', 'pas_name', 'pas_score', 'pas_strand']]

def filter_exons(df):
    positive_strand = df[df['strand'] == '+']
    positive_strand_filtered = positive_strand[(positive_strand['strand'] == positive_strand['pas_strand']) & 
                                               (positive_strand['end'] >= positive_strand['pas_start']) & 
                                               (positive_strand['end'] <= positive_strand['pas_end'])]
    negative_strand = df[df['strand'] == '-']
    negative_strand_filtered = negative_strand[(negative_strand['strand'] == negative_strand['pas_strand']) & 
                                               (negative_strand['start'] >= negative_strand['pas_start']) & 
                                               (negative_strand['start'] <= negative_strand['pas_end'])]
    filtered_df = pd.concat([positive_strand_filtered, negative_strand_filtered])
    return filtered_df

def find_terminal_exons(transcripts_df):
    start_terminal_exons = transcripts_df[(transcripts_df['feature'] == 'exon') & (transcripts_df['start'].isin(transcripts_df[transcripts_df['feature'] == 'transcript']['start']))]
    end_terminal_exons = transcripts_df[(transcripts_df['feature'] == 'exon') & (transcripts_df['end'].isin(transcripts_df[transcripts_df['feature'] == 'transcript']['end']))]
    terminal_exons = pd.concat([start_terminal_exons, end_terminal_exons]).drop_duplicates()
    return terminal_exons

def main(gtf_file, bed_file, quant_file, output_file):
    quant_df = pd.read_csv(quant_file, sep='\t')
    quant_df.columns = ['transcript_id', 'length', 'effective_length', 'tpm', 'counts'] 
    
    gtf_df = gtfparse.read_gtf(gtf_file)
    pas_sites_df = parse_bed_with_pandas(bed_file)
    bed_df = pd.read_csv(bed_file, sep='\t', index_col=False, header=None)
    gtf_df.index = pas_sites_df.index
    merged_df = gtf_df.join(pas_sites_df)
    
    transcripts_df = pd.merge(merged_df, quant_df, on='transcript_id', how='left')
    transcripts_df = transcripts_df.dropna(subset=['counts']).query('counts != 0')
    transcripts_df = transcripts_df.sort_values(by=['seqname', 'start'])
    
    terminal_exons = find_terminal_exons(transcripts_df)
    
    filtered_df = filter_exons(terminal_exons)
    
    # Save the filtered dataframe to a CSV file
    filtered_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    args = parse_arguments()
    main(args.gtf_file, args.bed_file, args.quant_file, args.output_file)
