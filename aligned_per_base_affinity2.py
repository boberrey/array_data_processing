#!/usr/bin/env python

"""
Calculate the per-base median affinity from an an array experiment on a fragmented,
aligned library.

Requires bedtools

Note: Python 3

Inputs:
   CPfitted output file from Sarah's 'fitSingleClusters.py' script
   insert bed file of aligned library from Ben's align_fastqs_array_lib.py'

Outputs:
   Per-base median affinities 

Ben Ober-Reynolds
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import subprocess
import time
from joblib import Parallel, delayed


def main():  
    # set up command line argument parser
    parser = argparse.ArgumentParser(description='script for calculating the \
        median per-base affinity of an aligned array library')
    group = parser.add_argument_group('required arguments:')
    group.add_argument('-cf', '--CPfitted_file', required=True,
        help='CPfitted file containing single-cluster fit data')
    group.add_argument('-bf', '--bed_file', required=True,
        help='bed file containing cluster "inserts"')
    group.add_argument('-g', '--genome_file', required=True,
        help='The genome file for experiment organism')
    group = parser.add_argument_group('optional arguments')
    group.add_argument('--sense', action='store_true',
        help='use flag if sense orientation should be maintained.')
    group.add_argument('-od', '--output_directory',
        help='directory to output per-base measurements. Default is bed_file directory')
    group.add_argument('-os', '--output_suffix', type=str, default='per_base_affinity',
        help='output prefix for per-base files (default is "per_base_affinity")')
    group.add_argument('-mc', '--median_cutoff', type=int, default=5,
        help='fewest number of clusters allowed for calculating medians (default = 5)')
    group.add_argument('-n', '--max_cores', type=int, default=20,
        help='maximum number of cores to process data on. Will try to be number of \
        chromosomes/regions detected in bed file')

    # print help if no arguments provided
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    # parse command line arguments
    args = parser.parse_args()
    max_cores = args.max_cores

    # Parameters:
    rsq_cutoff = 0.7
    outfile_ext = '.tsv'
    

    # If no output directory given, use bed file directory
    output_dir = args.output_directory
    if not output_dir:
        output_dir = os.path.dirname(args.bed_file)
    if output_dir == '':
        output_dir = '.'
    if not os.path.isdir(output_dir):
        print("Error: invalid output directory selection. Exiting...")
        sys.exit()

    # Set output prefix name
    output_suffix = args.output_suffix
    
    # Read in data
    print("Reading in data...")
    bed_df = pd.read_csv(args.bed_file, sep='\t', header=None)
    bed_df.columns = ['region', 'start', 'stop', 'cluster_ID', 'score', 'strand']

    # Get all regions (useful downstream)
    all_regions = set(bed_df.region)

    fit_df = pd.read_pickle(args.CPfitted_file)
    fit_df = fit_df.apply(pd.to_numeric)

    # Filter fit file by relevant columns:
    print("Filtering data...")
    cols_to_keep = ['dG', 'rsq']
    filtered_fit_df = fit_df[cols_to_keep]
    
    # Merge data frames:
    merged_df = bed_df.merge(filtered_fit_df, how='left', 
        left_on='cluster_ID', right_index=True)

    # What to do with poor-quality fits? Force them to affinity threshold
    # or remove them entirely? Clusters that look like higher affinity interactions
    # get better rsq in general, so removing poor fits may bias medians toward
    # higher affinity values... Really what we want is to have a heuristic
    # that forces clusters to the affinity threshold during the fitting.
    # I know from data exploration that imposing an rsq cutoff removes a huge
    # fraction of the data (e.g. more than 2/3 clusters)

    # For now, I'll write both options

    # Filter option:
    
    filt_merg_df = merged_df[merged_df.rsq > rsq_cutoff]
    cleaned_data = filt_merg_df.dropna(axis=0, subset=['dG'])
    

    """
    # Force option:
    dG_threshold = max(merged_df.dG)
    mask = (merged_df.rsq < rsq_cutoff) | (np.isnan(merged_df.dG))
    cleaned_data = merged_df
    cleaned_data.loc[mask, 'dG'] = dG_threshold
    cleaned_data.loc[mask, 'rsq'] = 0
    """
    # If indicated, maintain sense going forward
    print("Writing temporary files...")
    data_frames_dict = {}
    if args.sense:
        cleaned_data_plus = cleaned_data[cleaned_data.strand == '+']
        cleaned_data_minus = cleaned_data[cleaned_data.strand == '-']
        data_frames_dict['pos'] = cleaned_data_plus
        data_frames_dict['neg'] = cleaned_data_minus
    else:
        data_frames_dict['all'] = cleaned_data

    # Write all relevant temporary files:
    temp_data_filenames = {}
    for kind, frame in data_frames_dict.items():
        temp_data_filenames[kind] = output_dir + '/temp_data_' + kind + '.bed'
        frame.to_csv(temp_data_filenames[kind], sep='\t', header=None, index=None)

    # use bedtools to get coverage files for each saved file:
    print("Calculating coverage...")
    coverage_data_dict = {}
    for kind, filename in temp_data_filenames.items():
        command_list = ['bedtools', 'genomecov', '-i', filename, '-g', args.genome_file, '-d']
        print("Running command: \n\t {}".format(' '.join(command_list)))
        stdout_collect = subprocess.check_output(command_list).decode("utf-8")
        # Process stdout
        coverage_data_dict[kind] = clean_stdout_result(stdout_collect)
        # Remove temporary files
        subprocess.call(['rm', temp_data_filenames[kind]])

    # Split all data by region/chromosome for parallel processing
    # WARNING: The memory costs of this operation may be large. Be careful...
    region_split_data_dict = split_df_dict_by_regions(data_frames_dict, all_regions)
    region_split_cov_dict = split_df_dict_by_regions(coverage_data_dict, all_regions)
    
    ### Parallel Processing ###

    # Compute per-base median affinity in parallel
    sub_regions = list(region_split_cov_dict.keys())
    numCores = len(sub_regions)
    if numCores > args.max_cores:
        numCores = args.max_cores

    median_series = []
    if numCores > 1:
        print("Calculating per-base medians on {} cores...".format(numCores))
        median_series = (Parallel(n_jobs=numCores, verbose=10)\
            (delayed(get_base_affinities)(region_split_cov_dict[reg], 
                region_split_data_dict[reg], args.median_cutoff, reg) for reg in sub_regions))
    else:
        print("Calculating per-base medians on a single core...")
        median_series = [get_base_affinities(region_split_cov_dict[reg], 
                region_split_data_dict[reg], args.median_cutoff, reg) for reg in sub_regions]

    # Re-split median results by strand, if necessary
    master_series = {}
    for kind in data_frames_dict.keys():
        master_series[kind] = []

    for i, region in enumerate(sub_regions):
        for kind in master_series.keys():
            if kind in region:
                master_series[kind].append(median_series[i])

    # Concat all strand-relevant series together (recombine regions)
    for kind, s_list in master_series.items():
        master_series[kind] = pd.concat(s_list).rename("med_dG")

    # Finally, join per-base medians with coverage info and save output
    for kind, cov_df in coverage_data_dict.items():
        to_save = cov_df.join(master_series[kind])
        outfile = output_dir + '/' + os.path.basename(args.bed_file)
        outfile = os.path.splitext(outfile)[0] + '_' + args.output_suffix + kind + outfile_ext
        to_save.to_csv(outfile, sep='\t', header=None, index=None)



"""

    # Get clusters from all CPseq files
    cluster_sets = []
    if numCores > 1:
        print("Pulling cluster IDs from {} CPseq files on {} cores...".format(
            len(CPseq_list), numCores))
        cluster_sets = (Parallel(n_jobs=numCores, verbose=10)\
            (delayed(get_clusters_to_keep)(
                CPseq_file, filter_list) for CPseq_file in CPseq_list))
    else:
        cluster_sets = [get_clusters_to_keep(
            CPseq_file, filter_list) for CPseq_file in CPseq_list]
"""







def clean_stdout_result(stdout_result):
    """
    Cleanup the byte string resulting from a subprocess.check_output call
    from bedtools genomecov.
    Inputs:
        stdout_result (str) - the byte string collected from subprocess call
    Outputs:
        sanitized_df (DataFrame) - data frame formatted from byte string
    """
    cov_lists = [x.split('\t') for x in stdout_result.split('\n')[:-1]]
    cov_df = pd.DataFrame(cov_lists)
    cov_df.columns = ['region', 'base_pos', 'coverage']
    cov_df[['base_pos', 'coverage']] = cov_df[['base_pos', 'coverage']].apply(pd.to_numeric)
    return cov_df


def split_df_dict_by_regions(df_dict, all_regions):
    """
    Split a dataframe dict into sub-dataframes based on region.
    Inputs:
        df_dict (dict) - the data frame dict
    Outputs:
        region_dict (dict) - further subdivided region df dict
    """
    region_dict = {}
    for kind, df in df_dict.items():
        for region in all_regions:
            hard_copy = df[df.region == region].copy()
            region_dict[kind + '_' + region] = hard_copy
    return region_dict


def get_base_affinities(cov_df, data_df, median_cutoff, region):
    """
    Wrapper function around per-base median calculation.
    Inputs:
        cov_df (DataFrame) - the region-focused coverage dataframe
        data_df (DataFrame) - the region-focused data dataframe
        median_cutoff (int) - the minimum coverage per base for calculating
            median
    Outputs:
        medians (Series) - the median affinity per base
    """
    start_time = time.time()
    medians = cov_df.apply(lambda x: find_base_median(x, data_df, median_cutoff), axis=1)
    print("Finished region {} in {} minutes".format(region, 
        round((time.time() - start_time)/60.0, 2)))
    return medians


def find_base_median(row, data_df, median_cutoff):
    """
    Calculating the median affinity in this way is pretty slow and memory intensive. 
    Is there a better way to do this? (This is at least an order of magnitude
    faster than the most naive per-base slicing method)
    Inputs:
        row (Series) - the current row
        data_df (DataFrame) - the data dataframe
        median_cutoff (int) - the minimum coverage per base for calculating
            median
    Outputs:
        median (float) - the median affinity for this base
    """
    if row.coverage < median_cutoff:
        return float('NaN')
    # Start positions are sorted. Assume chromosomes are split.
    region = row.region
    base = row.base_pos
    # slice start by first insert that stops after base
    start_idx = data_df.stop.values.searchsorted(base, side='left')
    # slice end by first insert that starts after base
    end_idx = data_df.start.values.searchsorted(base, side='left')
    narrow_search = data_df[start_idx:end_idx]
    # Still need the refined search since fragments are of random length
    mini_df = narrow_search[(narrow_search.start < base) & (narrow_search.stop > base)]
    return mini_df.dG.median()



if __name__ == '__main__':
    main()
