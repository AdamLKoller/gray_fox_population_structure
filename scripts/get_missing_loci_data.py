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
        "-c", "--column", help="Name of output column", required = True
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
    return (row == 'NN').sum() / len(row)
    
def get_missing_loci_data(args):
    """
    Creates a file with the sample ID and the number of missing loci
    """
    
    meta = pd.read_csv(args.meta)
    meta = meta.sort_values(by='Sample_ID')
    meta.index = range(len(meta))

    geno_df = pd.read_table(args.geno, header=None)
    geno_df = geno_df.drop([0, 1, len(geno_df.columns)-1], axis=1).T
    geno_df.index = range(len(geno_df))
    n_SNPs = len(geno_df.columns)

    geno_df[args.column] = geno_df.apply(lambda row: count_missing(row), axis=1)
    geno_df = geno_df[[args.column]]
    geno_df = meta.merge(geno_df, left_index=True, right_index=True)

    geno_df[['Sample_ID', args.column]].to_csv(args.output)
    

def main():
    args = parser_arguments()
    return get_missing_loci_data(args)


if __name__ == "__main__":
    main()
