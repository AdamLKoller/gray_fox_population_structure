import pandas as pd
from collections import Counter
import numpy as np
import argparse


def parser_arguments():
    par = argparse.ArgumentParser()
    parser = par.add_argument_group("required arguments")
    parser.add_argument(
        "-g", "--geno", help="path to input geno file from ANGSD", required=True
    )
    parser.add_argument(
        "-o", "--output", help="path to write output file", required=True
    )
  
    args = par.parse_args()

    return args


def bases_to_number(bases: str, major_allele: str) -> int:
    if bases == "NN":
        return 9
    else:
        return 2 - bases.count(major_allele)


def geno_to_lfmm(args):
    """
    Transforms ANGSD genotype data to a matrix with 0s, 1s, and 2s
    based on recessive vs heterozygous vs dominant alleles
    missing data is encoded with a 9
    """
    df = pd.read_table(args.geno, header=None)
    df = df.drop([0, 1, len(df.columns)-1], axis=1).T
    df.index = range(len(df))
    
    major_alleles = []
    for column_name in df.columns:
        column = df[column_name].to_list()
        bases = "".join(column)
        allele_frequencies = dict(Counter(bases))
        major_allele = max(zip(allele_frequencies.values(), allele_frequencies.keys()))[
            1
        ]
        major_alleles.append(major_allele)

    for column_name, index in zip(df.columns, range(len(df.columns))):
        df[column_name] = df[column_name].apply(
            bases_to_number, major_allele=major_alleles[index]
        )

    df.to_csv(args.output, sep=' ', index = False, header=False)

    return df


def main():
    args = parser_arguments()
    return geno_to_lfmm(args)


if __name__ == "__main__":
    main()
