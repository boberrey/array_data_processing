#!/usr/bin/env python

"""
Isolate clusters that appear to have strong binding at early time points.
This is the first selection criteria, based only on signal value ranks.
Clusters can be filtered out based on variant_ID, proximity to a specific
variant, and quantification success.


Inputs:
   CPsignal file
   CPannot file

Outputs:
   CPannot file of suspected strong binders

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
    group.add_argument('-n','--n_cutoff', type=int, required=True,
                        help='cutoff for selection size.')
    
    group = parser.add_argument_group('optional arguments')
    group.add_argument('-ev', '--excluded_variants', default="",
                        help='A file containing variant IDs to be exluded from selection')
    group.add_argument('-pd','--previous_distance_file', default="",
                        help='A previously calculated min-distance file')
    # gaussian sd ranges from 0.95 to 1.55 pixels. Conversion from pixels to 'sequencer units' is 3.7.
    # In practice, optical abberations can complicate this cutoff. For example, you get more overlap of 
    # fiducial marks and 'strong binders' at edges of images. This means you may want to use a much more cautious 
    # cutoff than you initially thought.
    group.add_argument('-md','--min_distance', type=int, default=50,
                        help='The minimum distance in the provided distance (if provided) that will be considered. Default = 50')
    group.add_argument('-id','--image_indicees', default='1,2,3',
                        help='Which image indicees to use. Default is the first 3 excluding the baseline ("1,2,3"). Note that the first quantified image will be "0"')
    group.add_argument('-od','--output_directory', default="",
                        help='output directory for output file (default is current directory)')
    group.add_argument('-bf','--binder_filename', default="",
                        help='output filename for suspected binders')
    group.add_argument('-gnb','--get_non_binders', default="y",
                        help='should a selection of non-binders be picked? [y/n] default "y". Selection size will be the same as the binder selection.')
    group.add_argument('-nbf','--non_binder_filename', default="",
                        help='output filename for non_binder sample')

    # Some fixed parameters:
    variant_col_header = 'variant_ID'
    index_col_header = 'index'

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

    # Set output filenames
    binder_filename = args.binder_filename
    if binder_filename == "":
        binder_filename = os.path.basename(args.cpseries_file).split('.')[0] + "_early_binders.pkl"
    binder_filename = output_dir.strip('/') + '/' + binder_filename
    non_binder_filename = args.non_binder_filename
    if non_binder_filename == "":
        non_binder_filename = os.path.basename(args.cpseries_file).split('.')[0] + "_random_clusters.pkl"
    non_binder_filename = output_dir.strip('/') + '/' + non_binder_filename

    # variants to be excluded from selection:
    excluded_variants = []
    if args.excluded_variants != "":
        with open(args.excluded_variants, 'r') as f:
            for line in f:
                excluded_variants.append(line.strip())

    # Images to be used for selection:
    image_indicees = [int(x) for x in args.image_indicees.split(',')]

    # Check if distance file provided:
    # the distance df flag will be False otherwise
    if args.previous_distance_file != "":
        min_dist_df = pd.read_pickle(args.previous_distance_file)
    else:
        min_dist_df = pd.DataFrame()
        print "No previous distance file provided."

    # Read in data:
    print "Reading in CPannot file {}...".format(args.cpannot_file)
    annot_df = pd.read_pickle(args.cpannot_file)
    
    print "Reading in CPseries file {}...".format(args.cpseries_file)
    series_df = pd.read_pickle(args.cpseries_file)
    
    print "Data loaded successfully."
    # merge to incorporate variant_IDs, then free up memory
    merged_df = annot_df.merge(series_df, how='inner', left_index=True, right_index=True)
    del series_df
    gc.collect()


    # Filter by minimum distance to previously flagged variants:
    if not min_dist_df.empty:
        merged_df = filter_by_min_distance(merged_df, min_dist_df, args.min_distance)

    # Filter by excluded variants:
    if excluded_variants:
        merged_df = filter_by_excluded_variants(merged_df, excluded_variants, variant_col_header)

    # Determine ranks for indicated images, then sort by average rank
    merged_df = aggregate_and_sort_data(merged_df, image_indicees)

    # Take the top n clusters as the suspected strong binders
    strong_binder_indicees = list(merged_df.head(args.n_cutoff).index.values)
    strong_binder_annot_df = annot_df.loc[strong_binder_indicees,:]

    # save the strong binders
    print "Saving suspected strong binder file: {}".format(binder_filename)
    strong_binder_annot_df.to_pickle(binder_filename)

    # Do we also want a random selection of binders?
    if args.get_non_binders == 'y':
        random_indicees = np.random.choice(list(merged_df.index.values), 
            size=len(strong_binder_indicees), replace=False)
        random_annot_df = annot_df.loc[random_indicees,:]
        print "Saving random cluster file: {}".format(non_binder_filename)
        random_annot_df.to_pickle(non_binder_filename)

    print "Done."






def filter_by_min_distance(merged_df, min_dist_df, min_distance):
    """
    Filter data by proximity to previously flagged clusters (probably fiducial marks).
    """
    pre_len = len(merged_df.index)
    merged_df = merged_df.merge(min_dist_df, how='inner', left_index=True, right_index=True)
    merged_df = merged_df[merged_df['min_dist'] > min_distance]
    post_len = len(merged_df.index)
    percent = round(100*float(pre_len - post_len)/pre_len , 3)
    print "Distance cutoff of {} removed {} of {} clusters ({}%).".format(min_distance, 
        pre_len - post_len, pre_len, percent)
    del min_dist_df
    gc.collect()
    return merged_df


def filter_by_excluded_variants(merged_df, excluded_variants, variant_col_header):
    """
    Filter data by variants that should be excluded
    """
    pre_len = len(merged_df.index)
    merged_df = merged_df[~merged_df[variant_col_header].isin(excluded_variants)]
    post_len = len(merged_df.index)
    percent = round(100*float(pre_len - post_len)/pre_len , 3)
    print "Excluded variant filtration removed {} of {} clusters ({}%).".format(pre_len - post_len, 
        pre_len, percent)
    return merged_df


def aggregate_and_sort_data(merged_df, image_indicees):
    """
    Take the average of the signal ranks for the indicated image indicees
    """
    pre_len = len(merged_df.index)
    merged_df = merged_df[image_indicees].dropna()
    post_len = len(merged_df.index)
    percent = round(100*float(pre_len - post_len)/pre_len , 3)
    print "Unquantified value filtration removed {} of {} clusters ({}%).".format(pre_len - post_len, 
        pre_len, percent)
    merged_df = merged_df.rank(axis=0)
    merged_df = pd.DataFrame(merged_df.mean(1).sort_values(ascending=False))
    return merged_df



    
    




if __name__ == '__main__':
    main()