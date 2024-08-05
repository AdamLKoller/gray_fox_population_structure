import pandas as pd
import argparse
import numpy as np


def parser_arguments():
    par = argparse.ArgumentParser()
    parser = par.add_argument_group("required arguments")
    parser.add_argument(
        "-g", "--geno", help="path to input geno file from ANGSD", required=True
    )
    parser.add_argument(
        "-m", "--meta", help="path to the meta data file", required=True
    )
    parser.add_argument(
        "-o", "--output", help="path to write output file", required=True
    )
    parser.add_argument(
        "-b", "--bam", help="path to the file containing list of bams used as ANGSD input", required=True
    )
    args = par.parse_args()

    return args


def prepare_fst_input(args):
    """
    Transforms ANGSD genotype data to input file for hierfstat to calculate pairwise Fst
    """
    df = pd.read_table(args.geno, header=None)
    df = df.drop([0, 1, len(df.columns)-1], axis=1).T
    df = df.replace({'AA':11,'AG':12,'AC':13,'AT':14,
                     'CA':31,'CG':32,'CC':33,'CT':34,
                     'GA':21,'GG':22,'GC':23,'GT':24,
                     'TA':41,'TG':42,'TC':43,'TT':44,'NN':np.nan})
    df.index = range(len(df))
    
    bam = pd.read_table(args.bam, header=None)
    bam.columns = ['Sample_ID']
    bam['Sample_ID'] = bam['Sample_ID'].apply(lambda x: x.split('/')[-1].split('.')[0])
    
    df_meta = pd.read_csv(args.meta)
    bam = pd.merge(bam, df_meta, left_on='Sample_ID', right_on='Sample_ID')
    bam = bam[["POP"]]
    
    df = pd.merge(bam, df, left_index=True, right_index=True)
   
    df.to_csv(args.output)
    
    return df


def main():
    args = parser_arguments()
    return prepare_fst_input(args)


if __name__ == "__main__":
    main()
