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
        "-b", "--bam", help="path to file containing list of bams used by ANGSD", required=True
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
    bam = pd.read_csv(args.bam, header=None)

    meta = meta.sort_values(by='Sample_ID')
    meta.index = range(len(meta))

    bam[0] = bam[0].apply(lambda x: x.rstrip('.bam').split('/')[-1])
    meta = meta.merge(bam, left_on='Sample_ID', right_on=0)
    meta.index = range(len(meta))

    rel_df['ida'] = rel_df['ida'].apply(lambda x: x.rstrip('.bam'))
    rel_df['idb'] = rel_df['idb'].apply(lambda x: x.rstrip('.bam'))
    rel_df = rel_df.merge(meta.add_suffix('_a'), left_on='ida', right_on='Sample_ID_a')
    rel_df = rel_df.merge(meta.add_suffix('_b'), left_on='idb', right_on='Sample_ID_b')


    geno_df = geno_df.drop([0, 1, len(geno_df.columns)-1], axis=1).T
    geno_df.index = range(len(geno_df))
    geno_df['missing_count'] = geno_df.apply(lambda row: count_missing(row), axis=1)
    geno_df = geno_df[['missing_count']]
    geno_df = meta.merge(geno_df, left_index=True, right_index=True)


    print(f"Filtering {len(meta)} samples")
    filtered_df = rel_df.loc[rel_df.rab>0.45][['Sample_ID_a','Sample_ID_b','rab']]
    print(f"{len(filtered_df)} dyads detected (rab>0.45)")
    filtered_df = filtered_df.merge(geno_df.add_suffix('_a'), left_on='Sample_ID_a',right_on='Sample_ID_a')
    filtered_df = filtered_df.merge(geno_df.add_suffix('_b'), left_on='Sample_ID_b',right_on='Sample_ID_b')
    filtered_df['to_exclude'] = filtered_df.apply(get_highest_missing_id, axis = 1)

    for idx, row in filtered_df.iterrows():
        print(f"From dyad {row['Sample_ID_a']} and {row['Sample_ID_b']}, {row['to_exclude']} will be excluded as it has the lowest number of sites")


    meta['to_exclude'] = meta.apply(lambda row: to_exclude_or_not_to_exclude(row, list(filtered_df['to_exclude'])), axis = 1)
    
    meta.to_csv(args.output)

    return meta

def main():
    args = parser_arguments()
    return get_meta_w_dyads(args)


if __name__ == "__main__":
    main()
