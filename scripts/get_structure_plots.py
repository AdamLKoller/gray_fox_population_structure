import pandas as pd
import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

def parser_arguments():
    par = argparse.ArgumentParser()
    parser = par.add_argument_group("required arguments")
    parser.add_argument(
        "-i", "--input", help="basename to q matrices", required=True
    )
    
    parser.add_argument(
        "-o", "--output", help="path to write output figure", required=True
    )
    parser.add_argument(
        "-k", "--maxpops", help="display figures for k = 2 to maxpops (upper limit)", required=True
    )
    args = par.parse_args()

    return args



def get_structure_plots(args):
    """
    Produces structure plots
    """
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(16, 18))

    colors = ['orangered', 'mediumturquoise', 'magenta', 'gold']

    for i, axis in enumerate([ax1, ax2, ax3]):

        k = i + 2

        # Preparing data
        df = pd.read_csv(f'{args.input}{k}', index_col=0)
        df = df.apply(pd.to_numeric)
        df['max_prob'] = df.iloc[:, 1:].max(axis=1)
        df['cluster_assignment'] = df.iloc[:, 1:-1].idxmax(axis=1)
        df = df.sort_values(by=['cluster_assignment', 'max_prob'], ascending=[True, False])

        # Plotting data
        df.iloc[:, 1:-2].plot(kind='bar', stacked=True, ax=axis, width=1.0, color = colors)
        axis.set_xticks([])
        axis.legend().set_visible(False)
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)
        axis.spines['bottom'].set_visible(False)
        axis.spines['left'].set_visible(False)
        axis.set_ylabel("q", fontname='Arial', fontsize=16, fontweight='bold', fontstyle='italic')
        axis.tick_params(axis='both', width=2)

        axis.set_title(f'K = {k}', fontsize=18, fontweight='bold', fontstyle='italic')

    plt.tight_layout()
    
    # Save to PNG
    plt.savefig(args.output)
    

def main():
    args = parser_arguments()
    return get_structure_plots(args)


if __name__ == "__main__":
    main()
