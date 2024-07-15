import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from matplotlib.patches import Ellipse
import numpy as np
import argparse
import matplotlib.gridspec as gridspec
import geopandas as gpd

def parser_arguments():
    par = argparse.ArgumentParser()
    parser = par.add_argument_group("required arguments")
    parser.add_argument(
        "-g", "--geno", help="path to genotype matrix file", required=True
    )
    parser.add_argument(
        "-q", "--qmatrix", help="path to q matrix file", required=True
    )
    parser.add_argument(
        "-o", "--output", help="path to write output file", required=True
    )
    parser.add_argument(
        "-k", "--numpops", help="number of populations assumed", required=True
    )
    parser.add_argument(
        "-m", "--meta", help="path to metadata file", required=True
    )
    
    args = par.parse_args()

    return args

def plot_ellipse(ax, data, color, label):
    """
    Plot an ellipse representing the 95% confidence interval of the data points.
    """
    mean = np.mean(data, axis=0)
    cov = np.cov(data.T)
    lambda_, v = np.linalg.eig(cov)
    lambda_ = np.sqrt(lambda_)
    
    # Scaling factor for the 95% confidence interval
    scaling_factor = np.sqrt(5.991)  # Chi-squared value for 2 degrees of freedom and 95% confidence
    
    # Calculate the angle of the ellipse
    angle = np.degrees(np.arctan2(*v[:, 0][::-1]))
    
    # Width and height are 2 standard deviations, scaled to the 95% confidence interval
    width, height = 2 * lambda_ * scaling_factor
    
    ellipse = Ellipse(xy=mean,
                      width=width, height=height,
                      angle=angle, color=color, linestyle='--', linewidth=1, fill=False)
    
    ax.add_patch(ellipse)
    ax.scatter(data[:, 0], data[:, 1], s=4, color=color)

def get_pca_plot(ax, args):
    """
    Produces PCA plot colored for a given q matrix
    """
    X = pd.read_csv(args.geno, index_col=0)
    
    reduction_obj = PCA(n_components=2, random_state=0)
    projection = reduction_obj.fit_transform(X)
    projection = pd.DataFrame(projection, columns=['PCA1', 'PCA2'])
    
    q_matrix = pd.read_csv(args.qmatrix, index_col=0)
    q_matrix = q_matrix.apply(pd.to_numeric)
    q_matrix['max_prob'] = q_matrix.iloc[:, 1:].max(axis=1)
    q_matrix['cluster_assignment'] = q_matrix.iloc[:, 1:-1].idxmax(axis=1)
    
    merged_df = pd.merge(projection, q_matrix, left_index=True, right_index=True)
    
    meta = pd.read_csv(args.meta)
    meta = meta[meta.to_exclude == False] 
    meta.index = range(len(meta))
    merged_df = pd.merge(merged_df, meta, left_index=True, right_index=True)
    
    k = int(args.numpops)
    
    
    colors = ['gold', 'mediumturquoise', 'orangered','magenta',
              'orange', 'red', 'aqua', 'hotpink', 'lime', 'blue', 'khaki','sienna']
    for i, pop in enumerate(merged_df.POP.unique()):
        cluster = merged_df.loc[merged_df.POP == pop]
        if len(cluster) > 3:
            ax.scatter(x='PCA1', y='PCA2', color=colors[i-1], s=4, data=cluster, label=pop)
            #plot_ellipse(ax, cluster[['PCA1', 'PCA2']].values, colors[i], f'{pop}')
    
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.axvline(x=0, color='black', linewidth=0.5)
    
    print(f"Axes 1 explains {reduction_obj.explained_variance_[0]}% of the variance")
    print(f"Axes 2 explains {reduction_obj.explained_variance_[1]}% of the variance")
    
    ax.set_xlabel(f'Axis 1 ({round(reduction_obj.explained_variance_[0],2)}%)', fontname='Arial', fontweight='bold', fontsize=10)
    ax.set_ylabel(f'Axis 2 ({round(reduction_obj.explained_variance_[1],2)}%)', fontname='Arial', fontweight='bold', fontsize=10)
    
    legend = ax.legend(title='Cluster')
    plt.setp(legend.get_title(), fontsize=12, fontweight='bold')

def get_structure_plot(ax, args):
    """
    Produces structure bar plot
    """
    # Preparing data
    df = pd.read_csv(args.qmatrix, index_col=0)
    df = df.apply(pd.to_numeric)
    df['max_prob'] = df.iloc[:, 1:].max(axis=1)
    df['cluster_assignment'] = df.iloc[:, 1:-1].idxmax(axis=1)
    
    meta = pd.read_csv(args.meta)
    meta = meta[meta.to_exclude == False] 
    meta.index = range(len(meta))
    df = pd.merge(df, meta, left_index=True, right_index=True)
    
    df = df.sort_values(by=['POP','cluster_assignment', 'max_prob'], ascending=[True, True, False])

    # Plotting data
    colors = ['gold', 'mediumturquoise', 'orangered','magenta',
              'orange', 'red', 'aqua', 'hotpink', 'lime', 'blue']
    
    df.plot(kind='bar', y=[x for x in df.columns if x.startswith('prop_') and x != 'prop_missing'], stacked=True, ax=ax, width=1.0, color=colors)
    
    # Getting x ticks
    xticks = []
    xticklabels = []
    boundaries = []
    current_pop = df['POP'].iloc[0]
    last_pop_idx = 0

    for idx, pop in enumerate(df['POP']):
        if pop != current_pop:
            xticks.append(last_pop_idx)
            boundaries.append(last_pop_idx + (idx - last_pop_idx) // 2)
            xticklabels.append(current_pop)
            current_pop = pop
            last_pop_idx = idx

    # Add the last population
    xticks.append(last_pop_idx)
    boundaries.append(last_pop_idx + (len(df) - last_pop_idx) // 2)
    xticklabels.append(current_pop)
    
    ax.set_xticks(xticks)
    ax.set_xticklabels([''] * len(xticks))  # Set empty labels for boundary ticks

    # Set labels at the center of each population
    for boundary, label in zip(boundaries, xticklabels):
        ax.text(boundary, -0.01, label, ha='center', va='top', rotation=90, fontname='Arial', fontsize=10, transform=ax.get_xaxis_transform())

    ax.legend().set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_ylabel("q", fontname='Arial', fontsize=10, fontweight='bold', fontstyle='italic')
    ax.tick_params(axis='both', width=2)

def get_pie_map(ax, args):
    """
    Produces geographic map with each individual as a pie chart showing their proportions of ancestry.
    """
    k = int(args.numpops)
    usa_shapefile = gpd.read_file(r'./data/shapefiles/s_08mr23.shp')
    subspecies_shapefile = gpd.read_file(r'./data/shapefiles/subspecies_22780.shp')

    colors = ['gold', 'mediumturquoise', 'orangered','magenta',
              'orange', 'red', 'aqua', 'hotpink', 'lime', 'blue']
    
    meta = pd.read_csv(args.meta)
    meta = meta[meta.to_exclude == False] 
    meta.index = range(len(meta))
    
    qmatrix = pd.read_csv(args.qmatrix, index_col=0)
    qmatrix = pd.merge(qmatrix, meta, left_index=True, right_index=True)
    
    usa_shapefile.plot(ax=ax, facecolor='White', edgecolor='gray', linewidth=0.5)
    subspecies_shapefile.plot(ax=ax, facecolor='White', edgecolor='black', linewidth=2, linestyle='--', alpha=0.5)
    
    ax.set_ylim([30, 48.5])
    ax.set_xlim([-98, -78])

    ax.set_xlabel('Longitude', fontweight='bold', fontsize=10)
    ax.set_ylabel('Latitude', fontweight='bold', fontsize=10)

    ax.text(-97.5, 40, 'U. c. ocythous', fontsize=10, fontstyle='italic')
    ax.text(-87.5, 32.5, 'U. c. floridanus', fontsize=10, fontstyle='italic')
    ax.text(-87.5, 35.0, 'U. c. cinereoargenteus', fontsize=10, fontstyle='italic')

    for idx, row in qmatrix.iterrows():
        proportions = row[1:k+1]
        ax.pie(proportions, colors=colors, radius=0.25, 
               center=(row['longitude'], row['latitude']), wedgeprops={'clip_on': True}, frame=True)

def get_combined_plots(args):
    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1])
                 
    # Add subplots to the specified positions in the grid
    ax1 = fig.add_subplot(gs[0, 0])  # First row, first column
    ax2 = fig.add_subplot(gs[0, 1])  # First row, second column
    ax3 = fig.add_subplot(gs[1, 0]) # Second row, first column

    # Create the individual plots directly on the combined figure's axes
    get_pca_plot(ax1, args)
    get_pie_map(ax2, args)
    get_structure_plot(ax3, args)
                     
    plt.tight_layout()
                     
    plt.savefig(args.output)
    plt.show()

def main():
    args = parser_arguments()
    get_combined_plots(args)

if __name__ == "__main__":
    main()
