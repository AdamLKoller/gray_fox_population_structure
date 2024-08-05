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
        "-f", "--filter", help="percentage of SNPs required to keep individual", required=True
    )
    
    parser.add_argument(
        "-o", "--output", help="path to write output file", required=True
    )
    args = par.parse_args()

    return args


def count_missing(row):
    return (row == 'NN').sum()


def apply_imiss_filter(args):
    """
    Keeps samples with less than provided threshold percentage of missing data
    Ex. if -f is set to 0.8, then individuals missing data at more than 80% of SNPs are removed, 
    and individuals missing data at less than 80% of sites are retained.
    """
    
    bam = pd.read_table(args.bam, header=None)
    bam.columns = ['Sample_ID']
    bam['Sample_ID'] = bam['Sample_ID'].apply(lambda x: x.split('/')[-1].split('.')[0])

    geno_df = pd.read_table(args.geno, header=None)
    geno_df = geno_df.drop([0, 1, len(geno_df.columns)-1], axis=1).T
    geno_df.index = range(len(geno_df))
    
    geno_df['missing_count'] = geno_df.apply(lambda row: count_missing(row), axis=1)
    geno_df['missing_percent'] = geno_df['missing_count'] / (len(geno_df.columns) - 1)
    geno_df_w_ID = bam.merge(geno_df, left_index=True, right_index=True)
    geno_df_filtered = geno_df_w_ID.loc[geno_df_w_ID.missing_percent < float(args.filter)] 
    
    print(f"From the original {len(geno_df)} samples, {len(geno_df) - len(geno_df_filtered)} samples were removed as they were missing data at more than {float(args.filter) * 100}% of sites")
    
    out_file = open(args.output, 'w')
    for s_id in list(geno_df_filtered['Sample_ID']):
        out_file.write(f'data/merged_bams/{s_id}.bam\n')
        
    print("Filter threshold:", args.filter)
    print("Number of SNPs:", len(geno_df.columns)-2)
    print("Number of individuals prior to filtering:", len(geno_df))
    print("Number of individuals post filtering:", len(geno_df_filtered))
        


def main():
    args = parser_arguments()
    return apply_imiss_filter(args)


if __name__ == "__main__":
    main()
