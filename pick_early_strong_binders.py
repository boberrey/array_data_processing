#!/usr/bin/env python

"""
Isolate clusters that appear to have strong binding at early time points.
This is the first selection criteria, based only on signal percentiles.
Selection is also based on distance from fiducial marks.

Warnings: The distance matrix computation is catastrophically memory intensive
for a whole chip of data. Limit core use to prevent overloading the server.

Inputs:
   CPsignal file
   CPannot file

Outputs:
   filtered fastq files

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
    parser = argparse.ArgumentParser(description='script for isolating specific clusters from fastq files')
    group = parser.add_argument_group('required arguments:')
    group.add_argument('-cs', '--cpseries_file', required=True,
                        help='The CPseries file containing cluster signal info.')
    group.add_argument('-ca', '--cpannot_file', required=True,
                        help='The CPannot file.')
    group.add_argument('-fv', '--fiducial_variant_ID', required=True,
                        help='The variant ID associated with fiducial marks')
    group = parser.add_argument_group('optional arguments')
    group.add_argument('-pd','--previous_distance_file', default="",
                        help='A previously calculated min-distance file (will save a lot of time)')
    group.add_argument('-sd','--save_distance_file', default="yes",
                        help='If not providing a previously calculated distance file, do you want to save the one that will be calculated? [y/n], defualt y')
    group.add_argument('-pc','--percentile_cutoff', type=float, default=97.5,
                        help='The percentile cutoff for selection (default is 97.5)')
    group.add_argument('-id','--image_indicees', default='1,2,3',
                        help='Which image indicees to use. Default is the first 3 ("1,2,3").')
    group.add_argument('-set','--set_operation', default='union',
                        help='Whether to use the union of points passing percentile cutoff,\
                         or the intersection. (union/intersect), default is union.')
    group.add_argument('-od','--output_directory', default="",
                        help='output directory for output file (default is current directory)')
    group.add_argument('-op','--output_prefix', default="strong_binders",
                        help='prefex for the identified strong binders')
    group.add_argument('-n','--num_cores', type=int, default=1,
                        help='number of cores to use (should be same as number of tiles)')

    # Some fixed parameters:
    fiducial_warning_threshold = 1000
    variant_col_header = 'variant_ID'
    index_col_header = 'index'
    tile_col_header = 'tile'
    location_col_header = 'location'
    distance_filename = 'tst_dist.pkl'

    # print help if no arguments provided
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    #parse command line arguments
    args = parser.parse_args()

    # If no output directory given, use input directory
    output_dir = args.output_directory
    if output_dir == "":
        output_dir = './'
    if not os.path.isdir(output_dir):
        print "Error: invalid output directory selection. Exiting..."
        sys.exit()
    # Set output prefix
    output_prefix = args.output_prefix

    # Fiducial mark variant_ID:
    fid_var_ID = args.fiducial_variant_ID

    # Number of cores to use for distance matrix calculation
    numCores = args.num_cores

    # Check if distance file provided:
    # the distance df flag will be False otherwise
    
    if args.previous_distance_file != "":
        min_dist_df = pd.read_pickle(args.previous_distance_file)
    else:
        min_dist_df = False
        print "No previous distance file provided. Distances will be calculated."

    # Read in data:
    print "Reading in CPannot file {}...".format(args.cpannot_file)
    annot_df = pd.read_pickle(args.cpannot_file)
    num_fiducials = len(annot_df[annot_df[variant_col_header] == fid_var_ID].index)
    if num_fiducials < fiducial_warning_threshold:
        print "Only found {} fiducial marks with {} {}!\
         Is this the correct ID? Continue?".format(num_fiducials, variant_col_header, fid_var_ID)
        answer = raw_input('[y/n]')
        if answer.lower() != 'y':
            sys.exit()
    print "Reading in CPseries file {}...".format(args.cpseries_file)
    series_df = pd.read_pickle(args.cpseries_file)
    
    print "Data loaded successfully. Formatting for fiducial proximity filtering..."
    # merge to incorporate variant_IDs, then free up memory
    merged_df = annot_df.merge(series_df, how='inner', left_index=True, right_index=True)
    del annot_df
    del series_df
    print "garbage collected..."
    gc.collect()

    if not min_dist_df:
        # This takes a very long time if you're using a whole chip of data.
        min_dist_df = compute_min_distances(merged_df, fid_var_ID, variant_col_header, 
            index_col_header, tile_col_header, location_col_header, numCores)

    print(min_dist_df.head())
    min_dist_df.to_pickle(distance_filename)




    
    
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


def compute_min_distances(merged_df, fid_var_ID, variant_col_header, index_col_header, 
    tile_col_header, location_col_header, numCores):
    """
    A very hairy function for calculating the real distance for each cluster to
    its closest fiducial mark. Use of pandas slicing allows for significant speedup of 
    some sections at the cost of programatic elegance...
    """
    # separate the fiducial and non-fiducial clusters:
    no_fiducial_df = merged_df.loc[merged_df[variant_col_header] != fid_var_ID]

    # Filter by proximity to fiducial marks:
    fiducial_indicees = pd.Series(merged_df.loc[merged_df[variant_col_header] == fid_var_ID].index.values)
    non_fiducial_indicees = pd.Series(no_fiducial_df.index.values)

    # Get the locations and tiles for non-fiducial clusters and fiducial clusters
    non_fiducial_locations = non_fiducial_indicees.apply(func=extractLocation)
    non_fiducial_tiles = non_fiducial_indicees.apply(func=extractTile)

    fiducial_locations = fiducial_indicees.apply(func=extractLocation)
    fiducial_tiles = fiducial_indicees.apply(func=extractTile)


    # Construct some data frames with the previously obtained tile and location information
    # these intermediate structures will speed up dictionary creation in the next step.
    non_fid_loc_df = pd.DataFrame({index_col_header: non_fiducial_indicees, tile_col_header: non_fiducial_tiles, 
        location_col_header: non_fiducial_locations})
    non_fid_loc_df = non_fid_loc_df[[index_col_header, tile_col_header, location_col_header]]
    non_fid_loc_df = non_fid_loc_df.set_index(index_col_header)

    fid_loc_df = pd.DataFrame({index_col_header: fiducial_indicees, tile_col_header: fiducial_tiles, 
        location_col_header: fiducial_locations})
    fid_loc_df = fid_loc_df[[index_col_header, tile_col_header, location_col_header]]
    fid_loc_df = fid_loc_df.set_index(index_col_header)

    # Construct dictionaries keyed by tile:
    # Reformatting as dictionaries significantly reduces memory overhead during parallel
    # distance matrix computation
    all_tiles = sorted(list(set(fiducial_tiles)))
    fid_loc_dict = preallocate_dict(all_tiles)
    non_fid_index_dict = preallocate_dict(all_tiles)
    non_fid_loc_dict = preallocate_dict(all_tiles)

    # Fill in the dictionaries (the series maintain their original sorting)
    for tile in all_tiles:
        non_fid_index_dict[tile] = non_fid_loc_df[non_fid_loc_df[tile_col_header] == tile].index.values.tolist()
        non_fid_loc_dict[tile] = non_fid_loc_df[non_fid_loc_df[tile_col_header] == tile][location_col_header].tolist()
        fid_loc_dict[tile] = fid_loc_df[fid_loc_df[tile_col_header] == tile][location_col_header].tolist()

    # Free up as much memory as possible...
    del fiducial_indicees
    del non_fiducial_indicees
    del fiducial_locations
    del non_fiducial_locations
    del fiducial_tiles
    del non_fiducial_tiles
    del fid_loc_df
    del non_fid_loc_df
    print "garbage collected again..."
    gc.collect()

    # Use the dictionaries for parallel computation of distance matrices:
    if numCores > 1:
        print "Calculating distances in parallel with {} cores...".format(numCores)
        frames = (Parallel(n_jobs=numCores, verbose = 10)(delayed(getFiducialMarkProximities)(\
            non_fid_index_dict[tile], 
            non_fid_loc_dict[tile],
            fid_loc_dict[tile], tile) for tile in all_tiles))
    else:
        print "Calculating distances on a single core..."
        frames = [getFiducialMarkProximities(
            non_fid_index_dict[tile], 
            non_fid_loc_dict[tile],
            fid_loc_dict[tile], tile) for tile in all_tiles]

    return pd.concat(frames)


def getFiducialMarkProximities(non_fid_cluster_IDs, non_fid_loc_list, fid_loc_list, tile):
    # Construct a data frame of distances to closest fiducial mark for 
    # each cluster
    print "calculating distances for tile {}".format(tile)
    dm = spatial.distance_matrix(non_fid_loc_list, fid_loc_list)
    min_distances = [min(x) for x in dm]
    dist_df = pd.DataFrame({'cluster_ID': non_fid_cluster_IDs, 
        'min_dist': min_distances}).set_index('cluster_ID')
    return dist_df



if __name__ == '__main__':
    main()