import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import argparse

def parser_arguments():
    par = argparse.ArgumentParser()
    parser = par.add_argument_group("required arguments")
    parser.add_argument(
        "-m", "--meta", help="path to the meta data file", required=True
    )
    parser.add_argument(
        "-o", "--output", help="path to write output file", required=True
    )
    args = par.parse_args()

    return args


def add_apriori_pops(args):
    """
    Adds column 'POP' to meta data which contains the apriori sampling group.
    Apriori groups were defined as shapefile polygons. Points within polygons 
    are assigned to that group. Some points are outside of predefined polygons and
    are handled seperately.
    """
    
    apriori_pops = gpd.read_file(r'data/shapefiles/XYRed2017POP.shp') 
    meta = pd.read_csv(args.meta)
    geometry = [Point(xy) for xy in zip(meta['longitude'], meta['latitude'])]
    meta_geo = gpd.GeoDataFrame(meta, geometry=geometry)

    apriori_pops.set_crs(epsg=4326, inplace=True, allow_override=True)
    meta_geo.set_crs(epsg=4326, inplace=True)

    apriori_pops['geometry'] = apriori_pops['geometry'].buffer(0.0001)

    meta_with_pop = gpd.sjoin(meta_geo, apriori_pops[['POP', 'geometry']], how="left", op='within')
    meta_with_pop = meta_with_pop.drop(['geometry','index_right'], axis=1)

    meta_with_pop.loc[(meta_with_pop.latitude.round(3) == 30.863)
                                     & (meta_with_pop.longitude.round(3) == -91.799), 'POP'] = 'LA'

    meta_with_pop.loc[(meta_with_pop.latitude == 47.578635)
                                     & (meta_with_pop.longitude == -92.514570), 'POP'] = 'MN_WI_ND'

    meta_with_pop.loc[(meta_with_pop.latitude == 43.153671)
                                     & (meta_with_pop.longitude == -90.084763), 'POP'] = 'MN_WI_ND'

    meta_with_pop.loc[(meta_with_pop.latitude == 36.291422)
                                     & (meta_with_pop.longitude == -91.855491), 'POP'] = 'MO'

    meta_with_pop.to_csv(args.output)
    
    return meta_with_pop

def main():
    args = parser_arguments()
    return add_apriori_pops(args)


if __name__ == "__main__":
    main()
