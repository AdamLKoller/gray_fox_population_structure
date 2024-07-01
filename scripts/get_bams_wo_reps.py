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
        "-m", "--meta", help="path to the meta data file", required=True
    )
    parser.add_argument(
        "-o", "--output", help="path to write output file", required=True
    )
    args = par.parse_args()

    return args


def count_missing(row):
    return (row == 'NN').sum()

def get_highest_missing_id(row):
    if row['missing_count'] > row['missing_count_a']:
        return row['Sample_ID']
    else:
        return row['Sample_ID_a']
    


def get_bams_wo_replicates(args):
    """
    Excludes replicates with the lowest number of sites and prepares to rerun ANGSD.
    Also excludes samples with > 50% missing data.
    """
    
    meta = pd.read_csv(args.meta)
    meta = meta.sort_values(by='Sample_ID')
    meta.index = range(len(meta))

    geno_df = pd.read_table(args.geno, header=None)
    geno_df = geno_df.drop([0, 1, len(geno_df.columns)-1], axis=1).T
    geno_df.index = range(len(geno_df))
    n_SNPs = len(geno_df.columns)

    geno_df['missing_count'] = geno_df.apply(lambda row: count_missing(row), axis=1)
    geno_df = geno_df[['missing_count']]
    geno_df = meta.merge(geno_df, left_index=True, right_index=True)

    reps = geno_df.merge(geno_df.add_suffix('_a'), left_on='replicate', right_on='Sample_ID_a')
    reps['worse_sample'] = reps.apply(lambda row: get_highest_missing_id(row), axis=1)

    print(f"Filtering {len(meta)} samples")
    meta_filtered = meta[~meta['Sample_ID'].isin(reps['worse_sample'])]
    print(f"After removing lower quality replicates, {len(meta_filtered)} samples remain")
    meta_filtered = meta_filtered[~meta_filtered['Sample_ID'].isin(geno_df.loc[geno_df.missing_count > (n_SNPs / 2)]['Sample_ID'])]
    print(f"After removing low quality samples, {len(meta_filtered)} samples remain")
    
    out_file = open(args.output, 'w')
    for s_id in list(meta_filtered['Sample_ID']):
        out_file.write(f'data/merged_bams/{s_id}.bam\n')

    

def main():
    args = parser_arguments()
    return get_bams_wo_replicates(args)


if __name__ == "__main__":
    main()
