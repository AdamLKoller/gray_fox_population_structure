import random
import pandas as pd
from geopy.distance import distance
from geopy.point import Point
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


def jitter_coordinates(row, angle, distance_miles=20):
    origin = Point(row['latitude'], row['longitude'])
    new_point = distance(miles=distance_miles).destination(point=origin, bearing=angle)
    return new_point.latitude, new_point.longitude

def jitter_locations(args):
    """
    Jitters samples with duplicate coordinates to help in visualizing samples.
    """
    
    meta = pd.read_csv(args.meta)
    duplicated_coords = meta[meta.duplicated(subset=['latitude', 'longitude'], keep=False)]
    for idx, row in duplicated_coords.iterrows():
            new_lat, new_long = jitter_coordinates(row, random.randint(0,360))
            meta.loc[meta.Sample_ID == row.Sample_ID, 'latitude'] = new_lat
            meta.loc[meta.Sample_ID == row.Sample_ID, 'longitude'] = new_long
            
    meta.to_csv(args.output)

    
def main():
    args = parser_arguments()
    return jitter_locations(args)


if __name__ == "__main__":
    main()
