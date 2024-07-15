import pandas as pd
import argparse
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point, Polygon, box

def parser_arguments():
    par = argparse.ArgumentParser()
    parser = par.add_argument_group("required arguments")
    parser.add_argument(
        "-q", "--qmatrix", help="path to q matrix", required=True
    )
    
    parser.add_argument(
        "-m", "--meta", help="path to meta file", required=True
    )
    
    parser.add_argument(
        "-o", "--output", help="path to write output figure", required=True
    )
    parser.add_argument(
        "-k", "--numpops", help="display figures for k assumed populations", required=True
    )
    args = par.parse_args()

    return args



def get_pie_map(args):
    """
    Produces geographic map with each individual as a pie chart showing their proportions of ancestry.
    """
    
    k = int(args.numpops)
    usa_shapefile = gpd.read_file(r'./data/shapefiles/s_08mr23.shp')
    subspecies_shapefile = gpd.read_file(r'./data/shapefiles/subspecies_22780.shp')
    
    fig, ax = plt.subplots(figsize=(7,7))

    colors = ['orangered', 'mediumturquoise', 'magenta', 'gold', 
             'orange','red','aqua','hotpink','lime', 'blue'] 
    
    meta = pd.read_csv(args.meta)
    meta = meta[meta.to_exclude == False] 
    meta.index = range(len(meta))
    
    qmatrix = pd.read_csv(args.qmatrix, index_col=0)
    qmatrix = pd.merge(qmatrix, meta, left_index=True, right_index=True )
    
    usa_shapefile.plot(ax=ax, facecolor='White',edgecolor='gray', linewidth=0.5)
    subspecies_shapefile.plot(ax=ax, facecolor='White',edgecolor='black', linewidth=2, 
                              linestyle = '--', alpha=0.5)
    
    ax.set_ylim([30, 48.5])
    ax.set_xlim([-98, -78])

    ax.set_xlabel('Longitude', fontweight = 'bold', fontsize = 12)
    ax.set_ylabel('Latitude', fontweight = 'bold', fontsize = 12)

    ax.text(-97.5, 40, 'U. c. ocythous',fontsize=12, fontstyle = 'italic')
    ax.text(-87.5, 32.5, 'U. c. floridanus',fontsize=12, fontstyle = 'italic')
    ax.text(-87.5, 35.0, 'U. c. cinereoargenteus',fontsize=12, fontstyle = 'italic')

    for idx, row in qmatrix.iterrows():
        proportions = row[1:k+1]
        ax.pie(proportions, colors=colors, radius = 0.25, 
               center = (row['longitude'], row['latitude']), wedgeprops={'clip_on':True}, frame=True)
      
    # Save to PNG
    plt.savefig(args.output)
    

def main():
    args = parser_arguments()
    return get_pie_map(args)


if __name__ == "__main__":
    main()
