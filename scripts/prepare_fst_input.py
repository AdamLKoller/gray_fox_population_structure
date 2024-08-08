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
        "-q", "--qmatrix", help="path to the q matrix file", required=True
    )
    parser.add_argument(
        "-o", "--output", help="path to write output file", required=True
    )
    args = par.parse_args()

    return args


def prepare_fst_input(args):
    """
    Transforms ANGSD genotype data to input file for hierfstat to calculate pairwise Fst
    """
    df = pd.read_table(args.geno, header=None)
    df = df.drop([0, 1, len(df.columns)-1], axis=1).T
    df.index = range(len(df))
   
    df = df.replace({'AA':11,'AG':12,'AC':13,'AT':14,
                     'CA':31,'CG':32,'CC':33,'CT':34,
                     'GA':21,'GG':22,'GC':23,'GT':24,
                     'TA':41,'TG':42,'TC':43,'TT':44,'NN':np.nan})

    q_matrix = pd.read_csv(args.qmatrix, index_col=0)
    q_matrix = q_matrix.apply(pd.to_numeric)
    q_matrix['max_prob'] = q_matrix.iloc[:, 1:].max(axis=1)
    q_matrix['cluster_assignment'] = q_matrix.iloc[:, 1:-1].idxmax(axis=1)
    q_matrix['cluster_assignment'] = q_matrix['cluster_assignment'].apply(lambda x: int(x.replace('prop_','')))
    
    q_matrix = q_matrix[['cluster_assignment']]
    
    df = pd.merge(q_matrix, df, left_index=True, right_index=True)
    
    df.to_csv(args.output)
    
    return df


def main():
    args = parser_arguments()
    return prepare_fst_input(args)


if __name__ == "__main__":
    main()
