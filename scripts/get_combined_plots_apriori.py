import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from matplotlib.patches import Ellipse
import numpy as np
import argparse
import matplotlib.gridspec as gridspec
import geopandas as gpd
import seaborn as sns
from shapely.geometry import Point, Polygon

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
    q_matrix['Cluster'] = q_matrix['cluster_assignment'].apply(lambda x: x.split('_')[-1])
    
    merged_df = pd.merge(projection, q_matrix, left_index=True, right_index=True)
    
    meta = pd.read_csv(args.meta)
    meta = meta.rename(columns={'POP': "Sampling Group"})
    meta = meta[meta.to_exclude == False] 
    meta.index = range(len(meta))
    merged_df = pd.merge(merged_df, meta, left_index=True, right_index=True)
    
    k = int(args.numpops)
    
    colors = ['gold', 'mediumturquoise', 'orangered', 'magenta',
              'orange', 'red', 'aqua', 'hotpink', 'lime', 'blue', 'khaki', 'sienna']

    markers = ['o','d','v','^','<','>','s','P','*','D', 'X', 'p']

    
    for i, pop in enumerate(merged_df["Sampling Group"].unique()):
        data = merged_df.loc[merged_df["Sampling Group"] == pop]
        ax.scatter(x='PCA1',y='PCA2', data=data, color = colors[i], label=pop, s=4)
        plot_ellipse(ax, data[['PCA1', 'PCA2']].values, colors[i], pop)
    
    #scatter = sns.scatterplot(data=merged_df, x='PCA1', y='PCA2', hue='Cluster',
    #                          style='Sampling Group', palette=colors[:k], ax=ax,
    #                         markers=markers)
    
    
    
    #for cluster_num in range(1,k+1):
    #    cluster = merged_df.loc[q_matrix.Cluster == str(cluster_num)]
    #    if len(cluster) > 3:
    #        plot_ellipse(ax, cluster[['PCA1', 'PCA2']].values, colors[cluster_num-1], '')
    
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.axvline(x=0, color='black', linewidth=0.5)
    
    print(f"Axes 1 explains {reduction_obj.explained_variance_[0]}% of the variance")
    print(f"Axes 2 explains {reduction_obj.explained_variance_[1]}% of the variance")
    
    ax.set_xlabel(f'Axis 1 ({round(reduction_obj.explained_variance_[0], 2)}%)', fontname='Arial', fontweight='bold', fontsize=10)
    ax.set_ylabel(f'Axis 2 ({round(reduction_obj.explained_variance_[1], 2)}%)', fontname='Arial', fontweight='bold', fontsize=10)
    
    # Remove the original legend
    #ax.get_legend().remove()
    
    ax.set_xlim(-5,14)
    
    legend = ax.legend(title='Cluster')
    plt.setp(legend.get_title(), fontsize=12, fontweight='bold')
    
    # Update legend
#     handles, labels = ax.get_legend_handles_labels()
    
#     # Separate handles and labels for clusters and sampling groups
#     cluster_handles = handles[1:k+1]
#     cluster_labels = labels[1:k+1]
#     sampling_group_handles = handles[k+2:]
#     sampling_group_labels = labels[k+2:]
    
#     # Create custom legend with two columns
#     from matplotlib.legend import Legend

#     cluster_legend = Legend(ax, cluster_handles, cluster_labels, title='Cluster', loc='upper right', bbox_to_anchor=(0.6, 1), fontsize=10, frameon=False)
#     sampling_group_legend = Legend(ax, sampling_group_handles, sampling_group_labels, title='Sampling Group', loc='upper right', bbox_to_anchor=(1, 1), fontsize=10, frameon=False)
    
#     ax.add_artist(cluster_legend)
#     ax.add_artist(sampling_group_legend)
    
#     plt.setp(cluster_legend.get_title(), fontsize=12, fontweight='bold')
#     plt.setp(sampling_group_legend.get_title(), fontsize=12, fontweight='bold')
    
#     for legend in [cluster_legend, sampling_group_legend]:
#         legend.get_frame().set_facecolor('none')
#         legend.get_frame().set_edgecolor('none')
        
    

    

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
    
    df['POP'] = pd.Categorical(df['POP'], ['MN_WI_ND','IA_NE','MO','KS_OK','AR','LA','IL_IN','MS','KY_TN', 'OH_WV','SC_NC','MI'])
    
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
    #apriori_pops = gpd.read_file(r'./data/shapefiles//XYRed2017POP.shp')
    

    colors = ['gold', 'mediumturquoise', 'orangered','magenta',
              'orange', 'red', 'aqua', 'hotpink', 'lime', 'blue']
    
    meta = pd.read_csv(args.meta)
    meta = meta[meta.to_exclude == False] 
    meta.index = range(len(meta))

    qmatrix = pd.read_csv(args.qmatrix, index_col=0)
    qmatrix = qmatrix.apply(pd.to_numeric) #
    qmatrix['max_prob'] = qmatrix.iloc[:, 1:].max(axis=1) #
    qmatrix = pd.merge(qmatrix, meta, left_index=True, right_index=True)
    
    usa_shapefile.plot(ax=ax, facecolor='White', edgecolor='gray', linewidth=0.5)
    subspecies_shapefile.plot(ax=ax, facecolor='White', edgecolor='black', linewidth=2, linestyle='--', alpha=0.5)
    #apriori_pops.plot(ax=ax, facecolor='gray',edgecolor='black', linewidth=2, 
    #                      linestyle = '--', alpha=0.1)
    
    # Plotting apriori pops
    for population, group in qmatrix.groupby('POP'):
        points = [Point(lon, lat) for lon, lat in zip(group['longitude'], group['latitude'])]
        poly = gpd.GeoSeries(points).union_all().convex_hull
        x, y = poly.exterior.xy
        ax.fill(x, y, edgecolor='black', facecolor='gray', linewidth=1.5, alpha = 0.1)
    
    ax.set_ylim([30, 48.5])
    ax.set_xlim([-98, -77.5])

    ax.set_xlabel('Longitude', fontweight='bold', fontsize=10)
    ax.set_ylabel('Latitude', fontweight='bold', fontsize=10)

    ax.text(-97.5, 40, 'U. c. ocythous', fontsize=10, fontstyle='italic')
    ax.text(-87.5, 32.5, 'U. c. floridanus', fontsize=10, fontstyle='italic')
    ax.text(-87.5, 35.0, 'U. c. cinereoargenteus', fontsize=10, fontstyle='italic')

    
    for idx, row in qmatrix.iterrows():
        
        if row['max_prob'] < 0.75: #
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
