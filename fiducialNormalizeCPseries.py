#!/usr/bin/env python

"""
Normalize a CPseries file by the median fiducial signals within the file

Note: Python 2

Inputs:
   CPseries file (pickle)
   CPannot file (pickle)

Outputs:
   Normalized CPseries file (pickle)

Ben Ober-Reynolds
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np


def main():  
    # set up command line argument parser
    parser = argparse.ArgumentParser(description='Script for normalizing a CPseries \
        file by the median fiducial signal at each point.')
    group = parser.add_argument_group('required arguments:')
    group.add_argument('-cs', '--CPseries', required=True,
        help='CPseries file')
    group.add_argument('-ca', '--CPannot', required=True,
        help='Corresponding CPannot file')
    group = parser.add_argument_group('optional arguments')
    group.add_argument('-f', '--fiducial_variant_ID', type=str, default='11111111',
        help='The variant ID Corresponding to fiducial clusters (default is "11111111")')
    group.add_argument('-o', '--output_file',
        help='Output filename. Default is the input filename with "normalized" pre-appended')

    # print help if no arguments provided
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    # parse command line arguments
    args = parser.parse_args()

    series_file = os.path.abspath(args.CPseries)
    annot_file = os.path.abspath(args.CPannot)

    # Read in provided files
    print "Reading in CPseries file: {}".format(series_file)
    series_df = pd.read_pickle(series_file)
    print "Reading in CPannot file: {}".format(annot_file)
    annot_df = pd.read_pickle(annot_file)
    merged_df = annot_df.merge(series_df, how='inner', left_index=True, right_index=True)

    # Get median fiducial signal at each point
    fid_df = merged_df[merged_df.iloc[:,0] == args.fiducial_variant_ID]
    fid_medians = fid_df.iloc[:,1:].median(axis=0, skipna=True)
    print "Fiducial median signals:"
    for i in fid_medians:
        print round(i, 3)

    fid_med_med = np.nanmedian(fid_medians)

    # Divide all series by fiducial median, then rescale to overall median fiducial signal
    print "Normalizing CPseries to fiducial medians..."
    normalized_df = series_df.apply(lambda x: (x/fid_medians)*fid_med_med, 1)

    out_file = os.path.dirname(series_file) + "/" + "normalized_" + os.path.basename(series_file)
    if args.output_file is not None:
        out_file = args.output_file

    # Write to pickle
    normalized_df.to_pickle(out_file)
    print "Done."


if __name__ == '__main__':
    main()
