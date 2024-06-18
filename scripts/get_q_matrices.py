import pandas as pd
import argparse
import re

def parser_arguments():
    par = argparse.ArgumentParser()
    parser = par.add_argument_group("required arguments")
    parser.add_argument(
        "-i", "--input", help="path to STRUCTURE output file", required=True
    )

    parser.add_argument(
        "-o", "--output", help="path to write output file", required=True
    )
    parser.add_argument(
        "-k", "--pops", help="number of populations (k)", required=True
    )
    args = par.parse_args()

    return args




def get_q_matrices(args):
    """
    Extracts q matrices from structure _f output files
    """
    
    
    structure_in_file = open(args.input)
    line = structure_in_file.readline()
    while not line.startswith('Inferred ancestry of individuals:'):
        line = structure_in_file.readline()
    
    line = structure_in_file.readline()
    line = structure_in_file.readline()
    columns = ['prop_missing']
    for i in range(1,int(args.pops)+1):
        columns.append(f'prop_{i}')
        
    df = pd.DataFrame(columns=columns)
    
    while not line.startswith('\n'):
        row = re.split(r'\s+', line.strip().replace(':', ''))
        row[1] = float(row[1].replace('(','').replace(')','')) / 100
        row = [float(x) for x in row]
        df.loc[int(row[0])-1]=row[1:] 
        line = structure_in_file.readline()
            
    df.to_csv(args.output, index = True, header=True)
    
    return df


def main():
    args = parser_arguments()
    return get_q_matrices(args)


if __name__ == "__main__":
    main()
