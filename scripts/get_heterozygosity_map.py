import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
from pykrige.ok import OrdinaryKriging
import argparse

def parser_arguments():
    par = argparse.ArgumentParser()
    parser = par.add_argument_group("required arguments")
    parser.add_argument(
        "-g", "--geno", help="path to geno file from ANGSD", required=True
    )
    parser.add_argument(
        "-m", "--meta", help="path to q matrix file", required=True
    )
    parser.add_argument(
        "-o", "--output", help="path to write output file", required=True
    )

    args = par.parse_args()

    return args

def is_heterozygous(genotype):
    if genotype == 'NN':
        return np.nan
    return int(genotype[0] != genotype[1])


def get_heterozygosity_map(args):
    """
    Uses Kriging to interpolate heterozygosity rates over the map

    """
    df = pd.read_table(args.geno, header=None)
    df = df.drop([0, 1, len(df.columns)-1], axis=1).T
    df.index = range(len(df))
    df_hetero = df.applymap(is_heterozygous)
    df_hetero['Heterozygosity rate'] = df_hetero.mean(axis=1)

    meta = pd.read_csv(args.meta)
    meta = meta[meta.to_exclude == False] 


    df_hetero = pd.merge(df_hetero, meta, left_index=True, right_index=True )


    # Create a grid of points over the map area
    grid_lon = np.linspace(df_hetero['longitude'].min(), df_hetero['longitude'].max(), 500)
    grid_lat = np.linspace(df_hetero['latitude'].min(), df_hetero['latitude'].max(), 500)
    grid_lon, grid_lat = np.meshgrid(grid_lon, grid_lat)

    # Perform Kriging interpolation
    OK = OrdinaryKriging(
        df_hetero['longitude'], df_hetero['latitude'], df_hetero['Heterozygosity rate'],
        variogram_model='linear', verbose=False, enable_plotting=False
    )
    z, ss = OK.execute('grid', grid_lon[0], grid_lat[:, 0])

    # Load shapefiles
    usa_shapefile = gpd.read_file(r'./data/shapefiles/s_08mr23.shp')
    subspecies_shapefile = gpd.read_file(r'./data/shapefiles/subspecies_22780.shp')

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 7))

    # Plot shapefiles
    usa_shapefile.plot(ax=ax, facecolor='White', edgecolor='gray', linewidth=0.5)
    subspecies_shapefile.plot(ax=ax, facecolor='None', edgecolor='black', linewidth=2, linestyle='--', alpha=0.5)

    # Plot the interpolated surface
    c = ax.contourf(grid_lon, grid_lat, z, levels=20, cmap='jet', alpha=0.6)

    # Add color bar
    cbar = plt.colorbar(c, ax=ax)
    cbar.set_label('Heterozygosity rate', fontweight='bold', fontsize=12)
    
    # Plot points
    ax.scatter(df_hetero['longitude'], df_hetero['latitude'], color='black', s=5)

    # Set axis limits and labels
    ax.set_ylim([30, 48.5])
    ax.set_xlim([-98, -78])
    ax.set_xlabel('Longitude', fontweight='bold', fontsize=12)
    ax.set_ylabel('Latitude', fontweight='bold', fontsize=12)

    # Add text annotations
    ax.text(-97.5, 40, 'U. c. ocythous', fontsize=12, fontstyle='italic')
    ax.text(-87.5, 32.5, 'U. c. floridanus', fontsize=12, fontstyle='italic')
    ax.text(-87.5, 35.0, 'U. c. cineroargentus', fontsize=12, fontstyle='italic')
    
    plt.savefig(args.output)
    
    
    

def main():
    args = parser_arguments()
    return get_heterozygosity_map(args)


if __name__ == "__main__":
    main()
