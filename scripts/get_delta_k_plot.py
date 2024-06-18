import pandas as pd
import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

def parser_arguments():
    par = argparse.ArgumentParser()
    parser = par.add_argument_group("required arguments")
    parser.add_argument(
        "-i", "--input", help="path to summary file from STRUCTURE harvester", required=True
    )
    
    parser.add_argument(
        "-o", "--output", help="path to write output figure", required=True
    )
    
    args = par.parse_args()

    return args


def comma_format(x, pos):
    return '{:,.0f}'.format(x)

def get_delta_k_plot(args):
    """
    Produces delta k figure
    """
    
    structure_summary_file = open(args.input)
    line = structure_summary_file.readline()
    while not line.startswith('# K	Reps	mean est. LnP(Data)	stdev est. LnP(Data)'):
        line = structure_summary_file.readline()
        
    line = structure_summary_file.readline()
    
    df = pd.DataFrame(columns=['K', 'mean est. LnP(Data)', 'stdev est. LnP(Data)'])
    
    i = 0
    while not line.startswith('\n'):
        k,reps,mean,stdev = line.strip().split('\t')
        k,reps,mean,stdev = int(k), int(reps), float(mean), float(stdev)
        df.loc[i]= [k, mean,stdev]
        i += 1
        line = structure_summary_file.readline()
    
    
    df['delta_k'] = df['mean est. LnP(Data)'].diff()
    
    # Create a figure with two subplots side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # First plot
    ax1.errorbar(df['K'], df['mean est. LnP(Data)'], yerr=df['stdev est. LnP(Data)'], fmt='o',
                 color='blue', ecolor='lightgray', capsize=5, capthick=2)
    ax1.set_xlabel('K', fontname='Arial', fontweight='bold', fontstyle='italic', fontsize=16)
    ax1.set_ylabel(r'$\mathbf{Mean\ } \mathit{\mathbf{L(K)}}$', fontname='Arial', fontsize=16)
    ax1.yaxis.set_major_formatter(FuncFormatter(comma_format))
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_linewidth(2)
    ax1.spines['left'].set_linewidth(2)
    ax1.tick_params(axis='both', width=2)
    ax1.text(0.05, 0.95, '(a)', transform=ax1.transAxes, fontsize=16, fontweight='bold')

    # Second plot
    ax2.plot(df['K'], df['delta_k'], marker='o', color='blue', linestyle='-', linewidth=2)
    ax2.set_xlabel('K', fontname='Arial', fontweight='bold', fontstyle='italic', fontsize=16)
    ax2.set_ylabel(r'$\mathbf{\Delta\ K}$', fontname='Arial', fontsize=16, fontweight='bold')
    ax2.yaxis.set_major_formatter(FuncFormatter(comma_format))
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_linewidth(2)
    ax2.spines['left'].set_linewidth(2)
    ax2.tick_params(axis='both', width=2)
    ax2.axhline(y=0, color='black', linestyle='--', linewidth=0.5)
    ax2.text(0.05, 0.95, '(b)', transform=ax2.transAxes, fontsize=16, fontweight='bold')
    #ax2.set_ylim(top=100000)

    # Adjust layout
    plt.tight_layout()
    
    # Save to PNG
    plt.savefig(args.output)
    

def main():
    args = parser_arguments()
    return get_delta_k_plot(args)


if __name__ == "__main__":
    main()
