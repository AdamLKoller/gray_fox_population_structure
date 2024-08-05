import pandas as pd
import numpy as np
import argparse


def parser_arguments():
    par = argparse.ArgumentParser()
    parser = par.add_argument_group("required arguments")
    parser.add_argument(
        "-g", "--geno", help="path to geno call file from ANGSD", required=True
    )
    parser.add_argument(
        "-b", "--bam", help="path to the file containing list of bams used as ANGSD input", required=True
    )
    parser.add_argument(
        "-r", "--relate", help="path to the ngsrelate output file", required=True
    )
    parser.add_argument(
        "-c", "--cutoff", help="Relatedness cutoff", required=True
    )
    parser.add_argument(
        "-s", "--summary", help="Summary table keeping track of number of SNPs and number of individuals", required=True
    )
    parser.add_argument(
        "-o", "--output", help="path to write output file", required=True
    )
    args = par.parse_args()

    return args


def count_missing(row):
    return (row == 'NN').sum()

def get_highest_missing_id(row):
    if row['missing_percent_a'] > row['missing_percent_b']:
        return row['Sample_ID_a']
    else:
        return row['Sample_ID_b']
    
def apply_rcr_filter(args):
    """
    Filters out close relatives (rab > 0.45) and lower quality replicates (more missing sites)
    """
    
    bam = pd.read_table(args.bam, header=None)
    bam.columns = ['Sample_ID']
    bam['Sample_ID'] = bam['Sample_ID'].apply(lambda x: x.split('/')[-1].split('.')[0])
    
    

    geno_df = pd.read_table(args.geno, header=None)
    geno_df = geno_df.drop([0, 1, len(geno_df.columns)-1], axis=1).T
    geno_df.index = range(len(geno_df))
    geno_df['missing_count'] = geno_df.apply(lambda row: count_missing(row), axis=1)
    geno_df['missing_percent'] = geno_df['missing_count'] / (len(geno_df.columns) - 1)
    
    ID_to_perc = bam.merge(geno_df, left_index=True, right_index=True)[['Sample_ID','missing_percent']]
    
    rel_df = pd.read_table(args.relate)
    rel_df['ida'] = rel_df['ida'].apply(lambda x: x.rstrip('.bam'))
    rel_df['idb'] = rel_df['idb'].apply(lambda x: x.rstrip('.bam'))
    rel_df = rel_df.merge(ID_to_perc.add_suffix('_a'), left_on='ida', right_on='Sample_ID_a')
    rel_df = rel_df.merge(ID_to_perc.add_suffix('_b'), left_on='idb', right_on='Sample_ID_b')
    rel_df = rel_df[['Sample_ID_a','Sample_ID_b','missing_percent_a','missing_percent_b','rab']]
    
    close_relatives = rel_df.loc[rel_df.rab >= float(args.cutoff)]
    close_relatives['to_exclude'] = close_relatives.apply(get_highest_missing_id, axis = 1)
    close_relatives.to_csv(args.summary, sep='\t')
    
    for idx, row in close_relatives.iterrows():
        print(f"From dyad {row['Sample_ID_a']} and {row['Sample_ID_b']}, {row['to_exclude']} will be excluded as it has the lowest number of sites")
    
    to_include = list(set(bam['Sample_ID']) - set(close_relatives['to_exclude'])) 
    
    out_file = open(args.output, 'w')
    for s_id in to_include:
        out_file.write(f'data/merged_bams/{s_id}.bam\n')
        

def main():
    args = parser_arguments()
    return apply_rcr_filter(args)


if __name__ == "__main__":
    main()
