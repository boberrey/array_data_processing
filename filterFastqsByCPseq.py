#!/usr/bin/env python

"""
Filter a set of fastq files by specific clusters from a set of CPseq files.
This script was written for use in analyzing array libraries that require
alignment to a reference genome as part of their data processing. 

Note: Python 3

Inputs:
   directory of CPseq files from which to pick clusters
   directory of fastq files to filter

Outputs:
   filtered fastq files

Ben Ober-Reynolds
"""

import os
import sys
import argparse
from joblib import Parallel, delayed


def main():  
    # set up command line argument parser
    parser = argparse.ArgumentParser(description='script for isolating specific \
        clusters from fastq files, based on a set of CPseq files')
    group = parser.add_argument_group('required arguments:')
    group.add_argument('-sd', '--CPseq_dir', required=True,
        help='directory containing CPseq files')
    group.add_argument('-fd', '--fastq_directory', required=True,
        help='directory containing fastq files')
    group = parser.add_argument_group('optional arguments')
    group.add_argument('-f', '--filters', type=str, nargs='+',
        help='Which filter(s) from the CPseq files to extract. Separate multiple filters by a space.')
    group.add_argument('-od', '--output_directory',
        help='output directory for filtered fastq files (default is original \
            fastq_directory)')
    group.add_argument('-op', '--output_prefix', type=str, default='CPseq_filtered',
        help='output prefix for filtered fastq files (default is "CPseq_filtered")')
    group.add_argument('-n', '--num_cores', type=int, default=1,
        help='number of cores to use (should be same as number of fastq \
            files)')

    # print help if no arguments provided
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    # parse command line arguments
    args = parser.parse_args()
    numCores = args.num_cores

    # Pre-defined variables, constants, and settings
    fastq_extension = 'fastq'
    CPseq_extension = 'CPseq'
    default_prefix = 'CPseq_filter'

    # Check input directories
    CPseq_dir = args.CPseq_dir
    if not os.path.isdir(CPseq_dir):
        print("Error: invalid CPseq directory selection. Exiting...")
        sys.exit()
    fastq_dir = args.fastq_directory
    if not os.path.isdir(fastq_dir):
        print("Error: invalid fastq directory selection. Exiting...")
        sys.exit()

    # If no output directory given, use input directory
    output_dir = args.output_directory
    if not output_dir:
        output_dir = fastq_dir
    if not os.path.isdir(output_dir):
        print("Error: invalid output directory selection. Exiting...")
        sys.exit()

    # Set output prefix name
    output_prefix = args.output_prefix
    
    # Gather CPseq files:
    print("Finding CPseq files in directory {}".format(CPseq_dir))
    CPseq_list = find_files_in_directory(CPseq_dir, 
        extensionList=[CPseq_extension])

    # Gather fastq files:
    print("Finding fastq files in directory {}".format(fastq_dir))
    fastq_list = find_files_in_directory(fastq_dir, 
        extensionList=[fastq_extension])

    # Pick filters to use
    filter_list = args.filters
    if not filter_list:
        print("No filters provided. Extracting all clusters from CPseq files.")
    else:
        print("Extracting clusters with filters: {}".format(filter_list))

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

    # Combine clusters into one set
    # (Here's a handy trick: use a '*' before a list/tuple to break it apart, so 
    # each member becomes a new argument)
    all_clusters = set.union(*cluster_sets)
    print("Identified {} clusters for filtering fastqs".format(len(all_clusters)))

    # Adjust number of cores down to the number of fastq files, if necessary
    if numCores > len(fastq_list):
        numCores = len(fastq_list)
    
    # loop thorugh fastq files in parallel or in sequence
    results = []
    if numCores > 1:
        print("Filtering fastq files on {} cores...".format(numCores))
        results = (Parallel(n_jobs=numCores, verbose=10)\
            (delayed(filter_fastq)(all_clusters, fastq_file, output_prefix, 
            output_dir, fastq_extension) for fastq_file in fastq_list))
    else:
        results = [filter_fastq(all_clusters, fastq_file, output_prefix, 
            output_dir, fastq_extension) for fastq_file in fastq_list]
    
    # Report results of filtering:
    for result in results:
        print("file {} has {} clusters, filtered down from {}".format(
            result[0], result[1], result[2]))


def find_files_in_directory(dirPath, extensionList=None, 
                            excludedExtensionList=None):
    """
    Locate files in a given directory path. Optionally, desired files are 
    identified as matching one of the extension types provided in 
    'extensionList'
    Input: 
        dirPath (str) - path to directory
        extensionList (list) - list of acceptable extensions
        excludedExtensionList (list) - list of unacceptable extensions
    Output: 
        fileList (list) - list of found files (with path)
    """
    def extension_match(filename, extensionList=None):
        # from CPlibs
        if extensionList is not None:
            for currExt in extensionList:
                if filename.lower().endswith(currExt.lower()):
                    return True
        return False

    dirList = os.listdir(dirPath)
    fileList = []
    for currFilename in dirList:
        if (extension_match(currFilename, extensionList) 
        and not extension_match(currFilename, excludedExtensionList)): 
            fileList.append(dirPath+currFilename)
    if len(dirList) == 0:
        print('\tNONE FOUND')
    else:
        for filename in fileList:
            print("found:\t\t{}".format(filename))
        return fileList


def get_clusters_to_keep(CPseq_file, filter_list):
    """
    Generate a set of clusters extracted from a provided filename
    Input: 
        CPseq_file (str) - CPseq file from which to find clusters
        filter_list (list) - list of acceptable filters for a cluster to have
            Note: if no filters were provided, this will just return all clusters.
    Output: 
        cluster_set (set) - set of clusters to keep
    """
    cluster_set = set()
    if not filter_list:
        with open(CPseq_file, 'r') as f:
            for line in f:
                cluster = line.split()[0]
                cluster_set.add(cluster)
    else:
        with open(CPseq_file, 'r') as f:
            for line in f:
                # Note that in the case of a cluster with no filter, this will 
                # assign the R1 sequence to clust_filt. This is not exactly what
                # we want, but it shouldn't affect the script's behavior.
                cluster, clust_filt = line.split()[:2]
                for filt in filter_list:
                    if filt in clust_filt:
                        cluster_set.add(cluster)
    return cluster_set



def filter_fastq(filter_set, fastq_filename, output_prefix, output_dir, 
    fastq_extension):
    """
    filter a fastq file by clusters that exist in the cluster set, then
    save the filtered file as a new file
    (Note: I tried this function using biopython tools first, but it was an 
        order of magnitude slower than writing my own fastq parser. From what
        I could find online, this is likely because the biopython 
        implementations of the parsing and writing functions include 
        significantly more error checking than required for standard
        4-line fastq's)
    Input: filter_set, fastq_filename, fastq_identifier, output_prefix
    Output: saved filtered file
    """
    # Example fastq format:
    # @M00653:218:000000000-AYC5G:1:1101:20964:1096 1:N:0:1
    # CNTATAATGATTCTTATTGACCAAAAAGCTGACAATTCACTTATTTTGCTTGACTATTTATTATACTTTCA
    # +
    # C#8BCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGG
    
    fastq_basename = os.path.splitext(os.path.basename(fastq_filename))[0]
    new_filename = output_dir + output_prefix + '_' + fastq_basename + '.' + \
        fastq_extension

    # get the total number of lines:
    total_lines = 0
    with open(fastq_filename, 'r') as f:
        total_lines = sum(1 for line in f)

    # Every four lines is a new cluster:
    num_clusters = int(total_lines/4)
    cluster_count = 0

    # Loop through file in chunks of four lines:
    with open(fastq_filename, 'r') as infile, open(new_filename, 'w') as outfile:
        for chunk in range(num_clusters):
            cluster = infile.readline()
            seq = infile.readline()
            spacer = infile.readline()
            qual_score = infile.readline()

            # The first character of the cluster_ID in the fastq ('@') is not
            # present in the cluster_ID of the filter set.
            if cluster[1:].split()[0] in filter_set:
                cluster_count += 1
                outfile.write(cluster)
                outfile.write(seq)
                outfile.write(spacer)
                outfile.write(qual_score)
    
    # return the new filename, the starting number of clusters,
    # and the number of clusters kept
    return [new_filename, cluster_count, num_clusters]


if __name__ == '__main__':
    main()
