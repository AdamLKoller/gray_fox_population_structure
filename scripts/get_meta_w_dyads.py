import pandas as pd
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
        "-r", "--relate", help="path to the NgsRelate output file", required=True
    )
    parser.add_argument(
        "-o", "--output", help="path to write output file", required=True
    )
    args = par.parse_args()

    return args


def count_missing(row):
    return (row == 'NN').sum()

def get_highest_missing_id(row):
    if row['missing_count_a'] > row['missing_count_b']:
        return row['Sample_ID_a']
    else:
        return row['Sample_ID_b']
    
def to_exclude_or_not_to_exclude(row, exclusion_list):
    return row['Sample_ID'] in exclusion_list



def get_meta_w_dyads(args):
    """
    Uses output from NgsRelate to add a new column to the meta file `to_exclude`
    which excludes samples from population genetics analyses which are a member of 
    highly related dyads (defined by >0.45 relatedness). The member with the lowest
    number of genomic sites is removed and the member with the highest is kept.
    """
    
    rel_df = pd.read_table(args.relate)
    meta = pd.read_csv(args.meta)
    geno_df = pd.read_table(args.geno, header=None)

    rel_df = rel_df.merge(meta.add_suffix('_a'), left_on='a', right_index=True )
    rel_df = rel_df.merge(meta.add_suffix('_b'), left_on='b', right_index=True )

    geno_df = geno_df.drop([0, 1, len(geno_df.columns)-1], axis=1).T
    geno_df.index = range(len(geno_df))

    geno_df['missing_count'] = geno_df.apply(lambda row: count_missing(row), axis=1)
    geno_df = geno_df[['missing_count']]

    geno_df = meta.merge(geno_df, left_index=True, right_index=True)[['Sample_ID', 'missing_count']]

    filtered_df = rel_df.loc[rel_df.rab>0.45][['Sample_ID_a','Sample_ID_b','rab']]
    filtered_df = filtered_df.merge(geno_df.add_suffix('_a'), left_on='Sample_ID_a',right_on='Sample_ID_a')
    filtered_df = filtered_df.merge(geno_df.add_suffix('_b'), left_on='Sample_ID_b',right_on='Sample_ID_b')
    filtered_df
    filtered_df['to_exclude'] = filtered_df.apply(get_highest_missing_id, axis = 1)

    meta['to_exclude'] = meta.apply(lambda row: to_exclude_or_not_to_exclude(row, list(filtered_df['to_exclude'])), axis = 1)
    
    meta.to_csv(args.output)
    
    return meta

def main():
    args = parser_arguments()
    return get_meta_w_dyads(args)


if __name__ == "__main__":
    main()
