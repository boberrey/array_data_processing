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
    group.add_argument('-v', '--variant_IDs_of_interest', required=True,
                        help='A file containing variant IDs for which you want to calculate min distance to (min distance to closest of any included variant_ID)')
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
    VOI_IDs = args.variant_IDs_of_interest

    # variants to calculate distance to:
    VOI_IDs = []
    with open(args.variant_IDs_of_interest, 'r') as f:
        for line in f:
            VOI_IDs.append(line.strip())
        

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
        output_filename = os.path.basename(args.cpannot_file).split('.')[0] + "_VOI_distances.pkl"

    output_filename = output_dir.strip('/') + '/' + output_filename

    # Number of cores to use for distance matrix calculation
    numCores = args.num_cores

    # Read in data:
    print "Reading in CPannot file {}...".format(args.cpannot_file)
    annot_df = pd.read_pickle(args.cpannot_file)
    num_variants = len(annot_df[annot_df[variant_col_header].isin(VOI_IDs)].index)
    if num_variants < VOI_warning_threshold:
        print "Only found {} clusters with in variant ID list!".format(num_variants)
        print "Is this correct? Continue with distance calculation?"
        answer = raw_input('[y/n]')
        if answer.lower() != 'y':
            sys.exit()
    
    print "Data loaded successfully. Preparing for variant proximity filtering..."

    # Calculate distances:
    min_dist_df = compute_min_distances(annot_df, VOI_IDs, variant_col_header, 
            index_col_header, tile_col_header, location_col_header, numCores)
        
    # Save output as pickle
    print "Saving output..."
    min_dist_df.to_pickle(output_filename)
    print "Done."


def extract_x(cluster_id):
    # Extract x position
    split_list = cluster_id.split(':')
    return int(split_list[5])

def extract_y(cluster_id):
    # Extract y position
    split_list = cluster_id.split(':')
    return int(split_list[6])

def extractTile(cluster_id):
    # Extract the tile number from a cluster ID
    return int(cluster_id.split(':')[4][2:])


def preallocate_dict(keys):
    # Pre-allocate a dictionary with the keys required.
    d = {}
    for key in keys:
        d[key] = []
    return d


def compute_min_distances(annot_df, VOI_IDs, variant_col_header, index_col_header, 
    tile_col_header, location_col_header, numCores):
    """
    A very hairy function for calculating the real distance for each cluster to
    its closest fiducial mark. Use of pandas slicing allows for significant speedup of 
    some sections at the cost of programatic elegance...
    """
    # separate the VOI and non-VOI clusters:
    no_VOI_df = annot_df.loc[~annot_df[variant_col_header].isin(VOI_IDs)]

    # Filter by proximity to VOI clusters:
    VOI_indicees = pd.Series(annot_df.loc[annot_df[variant_col_header].isin(VOI_IDs)].index.values)
    non_VOI_indicees = pd.Series(no_VOI_df.index.values)

    non_VOI_x = non_VOI_indicees.apply(func=extract_x)
    non_VOI_y = non_VOI_indicees.apply(func=extract_y)
    non_VOI_tiles = non_VOI_indicees.apply(func=extractTile)

    VOI_x = VOI_indicees.apply(func=extract_x)
    VOI_y = VOI_indicees.apply(func=extract_y)
    VOI_tiles = VOI_indicees.apply(func=extractTile)



    # Construct some data frames with the previously obtained tile and location information
    # these intermediate structures will speed up dictionary creation in the next step.

    non_VOI_loc_df = pd.DataFrame({index_col_header: non_VOI_indicees, tile_col_header: non_VOI_tiles, 
        'x': non_VOI_x, 'y': non_VOI_y})
    non_VOI_loc_df = non_VOI_loc_df[[index_col_header, tile_col_header, 'x', 'y']]
    non_VOI_loc_df = non_VOI_loc_df.set_index(index_col_header)

    VOI_loc_df = pd.DataFrame({index_col_header: VOI_indicees, tile_col_header: VOI_tiles, 
        'x': VOI_x, 'y': VOI_y})
    VOI_loc_df = VOI_loc_df[[index_col_header, tile_col_header, 'x', 'y']]
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
        non_VOI_index_dict[tile] = non_VOI_loc_df[non_VOI_loc_df[tile_col_header] == tile].index.values
        non_VOI_loc_dict[tile] = non_VOI_loc_df[non_VOI_loc_df[tile_col_header] == tile][['x', 'y']].values
        VOI_loc_dict[tile] = VOI_loc_df[VOI_loc_df[tile_col_header] == tile][['x', 'y']].values

    # Free up as much memory as possible...
    del VOI_indicees
    del non_VOI_indicees
    del VOI_x
    del VOI_y
    del non_VOI_x
    del non_VOI_y
    del VOI_tiles
    del non_VOI_tiles
    del VOI_loc_df
    del non_VOI_loc_df
    del annot_df
    print "garbage collected..."
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