import pandas as pd
import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import matplotlib

def parser_arguments():
    par = argparse.ArgumentParser()
    parser = par.add_argument_group("required arguments")
    parser.add_argument(
        "-i", "--input", help="path to evanno.txt file from STRUCTURE harvester", required=True
    )
    
    parser.add_argument(
        "-o", "--output", help="path to write output figure", required=True
    )
    
    args = par.parse_args()

    return args


def comma_format(x, pos):
    return '{:,.0f}'.format(x)

def LegendVertical(Ax, Rotation=90, XPad=0, YPad=0, **LegendArgs):
    if Rotation not in (90,270):
        raise NotImplementedError('Rotation must be 90 or 270.')

    # Extra spacing between labels is needed to fit the rotated labels;
    # and since the frame will not adjust to the rotated labels, it is
    # disabled by default
    DefaultLoc = 'center right' if Rotation==90 else 'center right'
    ArgsDefaults = dict(loc=DefaultLoc, labelspacing=4, frameon=False)
    Args = {**ArgsDefaults, **LegendArgs}

    Handles, Labels = Ax.get_legend_handles_labels()
    if Rotation==90:
        # Reverse entries
        Handles, Labels = (reversed(_) for _ in (Handles, Labels))
    AxLeg = Ax.legend(Handles, Labels, **Args)

    LegTexts = AxLeg.get_texts()
    LegHandles = AxLeg.legend_handles

    for L,Leg in enumerate(LegHandles):
        if type(Leg) == matplotlib.patches.Rectangle:
            BBounds = np.ravel(Leg.get_bbox())
            BBounds[2:] = BBounds[2:][::-1]
            Leg.set_bounds(BBounds)

            LegPos = (
                # Ideally,
                #    `(BBounds[0]+(BBounds[2]/2)) - AxLeg.handletextpad`
                # should be at the horizontal center of the legend patch,
                # but for some reason it is not. Therefore the user will
                # need to specify some padding.
                (BBounds[0]+(BBounds[2]/2)) - AxLeg.handletextpad + XPad,

                # Similarly, `(BBounds[1]+BBounds[3])` should be at the vertical
                # top of the legend patch, but it is not.
                (BBounds[1]+BBounds[3])+YPad
            
            )
            print('hello')

        elif type(Leg) == matplotlib.lines.Line2D:
            LegXY = Leg.get_xydata()[:,::-1]
            Leg.set_data(*(LegXY[:,_] for _ in (0,1)))

            LegPos = (
                LegXY[0,0] - AxLeg.handletextpad + XPad,
                max(LegXY[:,1]) + YPad
            )
            
            print(LegPos)
            

        elif type(Leg) == matplotlib.collections.PathCollection or type(Leg) == matplotlib.collections.LineCollection:
            LegPos = (
                Leg.get_offsets()[0][0] + XPad,
                Leg.get_offsets()[0][1] + YPad,
            )
            print(Leg.get_offsets()[0][0])
            print(Leg.get_offsets()[0][1])
        else:
            print(type(Leg))
            print(Leg.get_offsets()[0][0])
            raise NotImplementedError('Legends should contain Rectangle, Line2D or PathCollection.')

        PText = LegTexts[L]
        PText.set_verticalalignment('bottom')
        PText.set_rotation(Rotation)
        PText.set_x(LegPos[0])
        PText.set_y(LegPos[1])
  

    return(None)





def get_delta_k_plot(args):
    """
    Produces delta k figure
    """
    
    structure_summary_file = open(args.input)
    line = structure_summary_file.readline()
    while not line.startswith("# K	Reps	Mean LnP(K)	Stdev LnP(K)	Ln'(K)	|Ln''(K)|	Delta K"):
        line = structure_summary_file.readline()
        
    line = structure_summary_file.readline()
    
    df = pd.DataFrame(columns=['K', 'Reps', 'Mean LnP(K)','Stdev LnP(K)',"Ln'(K)", "|Ln''(K)|", 'delta_k'])
    
    i = 0
    while not line.strip() == '':
        k,reps,mean,stdev,first_order,second_order,delta_k = [float(x) if x != 'NA' else pd.NA for x in line.strip().split('\t') ]
        df.loc[i]= [k,reps,mean,stdev,first_order,second_order,delta_k]
        i += 1
        line = structure_summary_file.readline()
    
    
    print(df)
    
    # Create a figure and a primary axis
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot the first series
    ax1.errorbar(df['K'], df['Mean LnP(K)'], yerr=df['Stdev LnP(K)'], fmt='o',
                 color='black', ecolor='black', capsize=5, capthick=2, linestyle='-',linewidth=2, label = r'$\mathbf{L(K)}$')
    ax1.set_xlabel('K', fontname='Arial', fontweight='bold', fontstyle='italic', fontsize=16)
    
       
    #ax1.set_ylabel(r'$\mathbf{Mean\ } \mathit{\mathbf{L(K)}}$', fontname='Arial', fontsize=16)
    
    ax1.yaxis.set_major_formatter(FuncFormatter(comma_format))
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_linewidth(2)
    ax1.spines['left'].set_linewidth(2)
    ax1.tick_params(axis='both', width=2)
    LegendVertical(ax1, 90, XPad=-50,YPad=-70, bbox_to_anchor=[-0.02, 0.67], prop = { "size": 16})
    #ax1.text(0.05, 0.95, '(a)', transform=ax1.transAxes, fontsize=16, fontweight='bold')

    # Create a secondary y-axis
    ax2 = ax1.twinx()
    ax2.plot(df['K'], df['delta_k'], marker='o', markerfacecolor='none', color='gray', linestyle='-', linewidth=2, label=r'$\mathbf{\Delta\ K}$')
    #ax2.set_ylabel(r'$\mathbf{\Delta\ K}$', fontname='Arial', fontsize=16, fontweight='bold')
    ax2.yaxis.set_major_formatter(FuncFormatter(comma_format))
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_linewidth(2)
    ax2.spines['bottom'].set_linewidth(2)
    ax2.spines['left'].set_visible(False)
    ax2.tick_params(axis='both', width=2)
    ax2.axhline(y=0, color='black', linestyle='--', linewidth=0.5)
    
    LegendVertical(ax2, 90, XPad=-65,YPad=-90, bbox_to_anchor=[1.2, 0.6], prop = { "size": 16})

    
    # Save to PNG
    plt.savefig(args.output)

    

def main():
    args = parser_arguments()
    return get_delta_k_plot(args)


if __name__ == "__main__":
    main()
