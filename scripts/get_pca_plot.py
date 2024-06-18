import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from matplotlib.patches import Ellipse
import numpy as np
import argparse

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
        "-k", "--numpops", help="number of populations assumed", required = True
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



def get_pca_plot(args):
    """
    Produces PCA plot colored for a given q matrix

    """
    X = pd.read_csv(args.geno, index_col=0)
    reduction_obj = PCA(n_components = 2, random_state = 0)
    projection = reduction_obj.fit_transform(X)
    projection = pd.DataFrame(projection, columns=['PCA1','PCA2'])
    
    q_matrix = pd.read_csv(args.qmatrix, index_col=0)
    q_matrix = q_matrix.apply(pd.to_numeric)
    q_matrix['max_prob'] = q_matrix.iloc[:, 1:].max(axis=1)
    q_matrix['cluster_assignment'] = q_matrix.iloc[:, 1:-1].idxmax(axis=1)
    
    merged_df = pd.merge(projection, q_matrix, left_index=True, right_index=True)
    
    k = int(args.numpops)
    
    colors = ['orangered', 'mediumturquoise', 'magenta', 'gold', 
             'orange','red','aqua','hotpink','lime', 'blue']
    
    plt.figure(figsize=(10, 6))
    ax = plt.gca()
    
    for i in range(1,k+1):
        cluster = merged_df.loc[merged_df.cluster_assignment == f'prop_{i}']
        if len(cluster) > 3:
            plt.scatter(x='PCA1',y='PCA2', color = colors[i-1], s = 4, data=cluster, 
                       label=str(i))
            plot_ellipse(ax, cluster[['PCA1', 'PCA2']].values, colors[i-1], f'Cluster {i}')

    ax = plt.gca()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axhline(y=0, color='black', linewidth = 0.5)
    ax.axvline(x=0, color='black', linewidth = 0.5)
    ax.set_xlabel(f'Axis 1 ({round(reduction_obj.explained_variance_[0],2)}%)', fontname='Arial', fontweight='bold', fontsize=10)
    ax.set_ylabel(f'Axis 2 ({round(reduction_obj.explained_variance_[1],2)}%)', fontname='Arial', fontweight='bold', fontsize=10)
    
    legend = plt.legend(title='Cluster')
    plt.setp(legend.get_title(), fontsize=12, fontweight='bold')
    
    plt.savefig(args.output)
    
    
    

def main():
    args = parser_arguments()
    return get_pca_plot(args)


if __name__ == "__main__":
    main()
