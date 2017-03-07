#!/usr/bin/env python

"""
Identify the closest variant of a given type to each cluster.
Useful for filtering things that are too close to fiducial marks for example.

Warnings: The distance matrix computation is catastrophically memory intensive
for a whole chip of data. Limit core use to prevent overloading the server.

Inputs:
   CPannot file
   variant_ID of interest

Outputs:
   pickled 'nearest variant' file

Ben Ober-Reynolds
"""



import os
import sys
import argparse
import pandas as pd
import numpy as np
import gc
from scipy import spatial
from joblib import Parallel, delayed

def main():
    ################ Parse input parameters ################

    #set up command line argument parser
    parser = argparse.ArgumentParser(description='script for identifying the nearest cluster of a specific variant, for each other variant.')
    group = parser.add_argument_group('required arguments:')
    group.add_argument('-ca', '--cpannot_file', required=True,
                        help='The CPannot file (cluster_IDs and associated variant_IDs).')
    group.add_argument('-v', '--variant_ID_of_interest', required=True,
                        help='The variant ID for which you want to find distances to (probably the fiducial mark ID)')
    group = parser.add_argument_group('optional arguments')
    group.add_argument('-od','--output_directory', default="",
                        help='output directory for output file (default is current directory)')
    group.add_argument('-of','--output_filename', default="",
                        help='prefex for the calculated distance file (default is CPannot filename with variant_ID appended)')
    group.add_argument('-n','--num_cores', type=int, default=1,
                        help='number of cores to use (uses lots of memory!)')

    # Some fixed parameters:
    VOI_warning_threshold = 1000
    variant_col_header = 'variant_ID'
    index_col_header = 'index'
    tile_col_header = 'tile'
    location_col_header = 'location'

    # print help if no arguments provided
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    #parse command line arguments
    args = parser.parse_args()

    # Fiducial mark variant_ID:
    VOI_ID = args.variant_ID_of_interest

    # If no output directory given, use current directory
    output_dir = args.output_directory
    if output_dir == "":
        output_dir = './'
    if not os.path.isdir(output_dir):
        print "Error: invalid output directory selection. Exiting..."
        sys.exit()

    # Set output filename
    output_filename = args.output_filename
    if output_filename == "":
        output_filename = os.path.basename(args.cpannot_file).split('.')[0] + "_{}_distances.pkl".format(VOI_ID)

    output_filename = output_dir + '/' + output_filename
    print output_filename

    # Number of cores to use for distance matrix calculation
    numCores = args.num_cores


    # Read in data:
    print "Reading in CPannot file {}...".format(args.cpannot_file)
    annot_df = pd.read_pickle(args.cpannot_file)
    num_variants = len(annot_df[annot_df[variant_col_header] == VOI_ID].index)
    if num_variants < VOI_warning_threshold:
        print "Only found {} clusters with {} {}!\
         Is this the correct ID? Continue?".format(num_fiducials, variant_col_header, VOI_ID)
        answer = raw_input('[y/n]')
        if answer.lower() != 'y':
            sys.exit()
    
    print "Data loaded successfully. Preparing for variant proximity filtering..."
    # free up memory
    del annot_df
    gc.collect()

    min_dist_df = compute_min_distances(annot_df, VOI_ID, variant_col_header, 
            index_col_header, tile_col_header, location_col_header, numCores)
        

    min_dist_df.to_pickle(output_dir)

    
    
def extractLocation(cluster_id):
    # Extract the location (x, y) from a cluster ID
    split_list = cluster_id.split(':')
    x = int(split_list[5])
    y = int(split_list[6])
    return [x, y]

def extractTile(cluster_id):
    # Extract the tile number from a cluster ID
    return int(cluster_id.split(':')[4][2:])


def preallocate_dict(keys):
    # Pre-allocate a dictionary with the keys required.
    d = {}
    for key in keys:
        d[key] = []
    return d


def compute_min_distances(annot_df, VOI_ID, variant_col_header, index_col_header, 
    tile_col_header, location_col_header, numCores):
    """
    A very hairy function for calculating the real distance for each cluster to
    its closest fiducial mark. Use of pandas slicing allows for significant speedup of 
    some sections at the cost of programatic elegance...
    """
    # separate the VOI and non-VOI clusters:
    no_VOI_df = annot_df.loc[annot_df[variant_col_header] != VOI_ID]

    # Filter by proximity to VOI clusters:
    VOI_indicees = pd.Series(annot_df.loc[annot_df[variant_col_header] == VOI_ID].index.values)
    non_VOI_indicees = pd.Series(no_VOI_df.index.values)

    # Get the locations and tiles for non-fiducial clusters and fiducial clusters
    non_VOI_locations = non_VOI_indicees.apply(func=extractLocation)
    non_VOI_tiles = non_VOI_indicees.apply(func=extractTile)

    VOI_locations = VOI_indicees.apply(func=extractLocation)
    VOI_tiles = VOI_indicees.apply(func=extractTile)


    # Construct some data frames with the previously obtained tile and location information
    # these intermediate structures will speed up dictionary creation in the next step.
    non_VOI_loc_df = pd.DataFrame({index_col_header: non_VOI_indicees, tile_col_header: non_VOI_tiles, 
        location_col_header: non_VOI_locations})
    non_VOI_loc_df = non_VOI_loc_df[[index_col_header, tile_col_header, location_col_header]]
    non_VOI_loc_df = non_VOI_loc_df.set_index(index_col_header)

    VOI_loc_df = pd.DataFrame({index_col_header: VOI_indicees, tile_col_header: VOI_tiles, 
        location_col_header: VOI_locations})
    VOI_loc_df = VOI_loc_df[[index_col_header, tile_col_header, location_col_header]]
    VOI_loc_df = VOI_loc_df.set_index(index_col_header)

    # Construct dictionaries keyed by tile:
    # Originally thought using dictionaries would reduce memory usage during parallel processing... they don't though.
    # distance matrix computation
    all_tiles = sorted(list(set(VOI_tiles)))
    VOI_loc_dict = preallocate_dict(all_tiles)
    non_VOI_index_dict = preallocate_dict(all_tiles)
    non_VOI_loc_dict = preallocate_dict(all_tiles)

    # Fill in the dictionaries (the series maintain their original sorting)
    for tile in all_tiles:
        non_VOI_index_dict[tile] = non_VOI_loc_df[non_VOI_loc_df[tile_col_header] == tile].index.values.tolist()
        non_VOI_loc_dict[tile] = non_VOI_loc_df[non_VOI_loc_df[tile_col_header] == tile][location_col_header].tolist()
        VOI_loc_dict[tile] = VOI_loc_df[VOI_loc_df[tile_col_header] == tile][location_col_header].tolist()

    # Free up as much memory as possible...
    del VOI_indicees
    del non_VOI_indicees
    del VOI_locations
    del non_VOI_locations
    del VOI_tiles
    del non_VOI_tiles
    del VOI_loc_df
    del non_VOI_loc_df
    print "garbage collected again..."
    gc.collect()

    # Use the dictionaries for parallel computation of distance matrices:
    if numCores > 1:
        print "Calculating distances in parallel with {} cores...".format(numCores)
        frames = (Parallel(n_jobs=numCores, verbose = 10)(delayed(get_VOI_proximities)(\
            non_VOI_index_dict[tile], 
            non_VOI_loc_dict[tile],
            VOI_loc_dict[tile], tile) for tile in all_tiles))
    else:
        print "Calculating distances on a single core..."
        frames = [get_VOI_proximities(
            non_VOI_index_dict[tile], 
            non_VOI_loc_dict[tile],
            VOI_loc_dict[tile], tile) for tile in all_tiles]

    return pd.concat(frames)


def get_VOI_proximities(non_VOI_clusterIDs, non_VOI_loc_list, VOI_loc_list, tile):
    # Construct a data frame of distances to closest fiducial mark for 
    # each cluster
    print "calculating distances for tile {}".format(tile)
    dm = spatial.distance_matrix(non_VOI_loc_list, VOI_loc_list)
    min_distances = [min(x) for x in dm]
    dist_df = pd.DataFrame({'cluster_ID': non_VOI_clusterIDs, 
        'min_dist': min_distances}).set_index('cluster_ID')
    return dist_df



if __name__ == '__main__':
    main()