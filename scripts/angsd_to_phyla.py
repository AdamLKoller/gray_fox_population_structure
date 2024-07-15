import pandas as pd
import argparse
from collections import Counter


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
    args = par.parse_args()

    return args


def bases_to_number(bases: str, major_allele: str) -> int:
    if bases == "NN":
        return '-'
    else:
        return "1" if bases != major_allele*2 else "0"
    

def angsd_to_phyla(args):
    """
    Transforms ANGSD genotype data to input file for iqtree
    """
    df = pd.read_table(args.geno, header=None)
    df = df.drop([0, 1, len(df.columns)-1], axis=1).T
    df.index = range(len(df))
    df_meta = pd.read_csv(args.meta)
    df_meta = df_meta[df_meta.to_exclude == False]
    df_meta = df_meta[["Sample_ID"]]
    df = pd.merge(df, df_meta, left_index=True, right_index=True)
    df.index = df.Sample_ID
    df = df.drop(
        [
            "Sample_ID"
            
        ],
        axis=1,
    )

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
    
    phyla_out = open(args.output, 'w')
    
    phyla_out.write(f'{df.shape[0]}\t{df.shape[1]}\n')
    
    for index, row in df.iterrows():
        row = ''.join(row.tolist())
        phyla_out.write(f'{index}\t{row}\n')
    

def main():
    args = parser_arguments()
    return angsd_to_phyla(args)


if __name__ == "__main__":
    main()
