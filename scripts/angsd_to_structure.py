import pandas as pd
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


def get_left_base(geno):
    return geno[0]
    
def get_right_base(geno):
    return geno[1]
    

def angsd_to_structure(args):
    """
    Transforms ANGSD genotype data to input file for STRUCTURE
    """
    df = pd.read_table(args.geno, header=None)
    df = df.drop([0, 1, len(df.columns)-1], axis=1).T
    
    df_l = df.applymap(get_left_base)
    df_r = df.applymap(get_right_base)
    df_cat = pd.concat([df_l, df_r], ignore_index=False)
    df_cat = df_cat.sort_index()
    
    df_cat = df_cat.replace({'A':1,'T':2,'C':3,'G':4,'N':-1})
    df_cat.to_csv(args.output,sep='\t',index=False,header=True)
    
    return df_cat


def main():
    args = parser_arguments()
    return angsd_to_structure(args)


if __name__ == "__main__":
    main()
