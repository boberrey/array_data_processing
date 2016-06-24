#!/usr/bin/env python
""" Assign variant IDs to clusters from a CPseries file. 

 Variant IDs are assigned based on matching sequence. 
 The user may specify which bases of which sequence are to be used for assigning variants.

 Inputs:
   CPseries files

 Outputs:
   CPseries files with variants labeled
   File containing list of variant IDs and associated sequences

 Ben Ober-Reynolds"""

import sys
import os
import argparse
import cpfiletools
import random
import string
import subprocess
import pandas as pd
import numpy as np

import time


### MAIN ###


def main():
	################ Parse input parameters ################

	#set up command line argument parser
	parser = argparse.ArgumentParser(description='Script for assigning unique IDs to variants in CPseries files')
	group = parser.add_argument_group('required arguments')
	group.add_argument('-sd', '--series_directory', required=True,
	                    help='directory that holds the CPseries files that need variant IDs')
	group.add_argument('-sc', '--seq_column', required=True,
	                    help='which column in the CPseries file you want to use for assigning variants')

	group = parser.add_argument_group('optional arguments for processing data')
	group.add_argument('-st','--seq_start', default=0,
	                    help='start position within sequence for matching')
	group.add_argument('-ed','--seq_end', default=0,
	                    help='end position within sequence for matching (note: will use whole sequence if position "0" given')
	group.add_argument('-lb','--label', default="ID_ed",
	                    help='label attached to output files')
	group.add_argument('-od','--output_directory', default="",
	                    help='output directory for series files with labeled variants (default will use series_directory)')
	group.add_argument('-if','--ID_file', default="ID_file.txt",
	                    help='file name for the list of IDs and corresponding sequences')

	if not len(sys.argv) > 1:
	    parser.print_help()
	    sys.exit()


	#parse command line arguments
	args = parser.parse_args()

	# If no output directory given, use input directory
	if args.output_directory == "":
		args.output_directory = args.series_directory

	# This script will run through each of the provided CPseries files sequentially in order to 
	# ensure that each variant gets assigned only one variant ID.

	CPseriesFiles = cpfiletools.find_files_in_directory(args.series_directory, ['CPseries'])

	numLines = 0
	for seriesFile in CPseriesFiles:
		numLines += int(subprocess.check_output(("wc -l {} | ".format(
			os.path.join(args.series_directory, seriesFile))+" awk \'{print $1}\'"), shell=True).strip())

	start = time.time()
	randID_set = set()
	print "Generating random IDs..."
	while len(randID_set) < numLines:
		randID = ''.join([random.choice(string.ascii_uppercase + string.digits) for n in range(8)])	# 8^36 ~ 3.25e32 possible IDs
		randID_set.add(randID)
	print "ID generation: {0:.2f} seconds".format(time.time() - start)


	# This dictionary will contain all the variants assigned, keyed by sequence match
	# The entries in variant dict will be three-element lists, the first is the ID, the second is the filter
	# associated with that variant (if any), and the third is the number of times that variant has been seen
	variantDict = {}

	seq_col = int(args.seq_column)
	strt = int(args.seq_start)
	end = int(args.seq_end)
	fileNum = 1

	# Loop through each CPseries file to assign variants:
	for seriesFile in CPseriesFiles:
		print "Working on file: {}...{} of {}".format(seriesFile, fileNum, len(CPseriesFiles))
		labeled_filename = os.path.join(args.output_directory, ".".join(['_'.join([os.path.splitext(seriesFile)[0], args.label]), 'CPseries']))
		# Time each loop for now:
		start = time.time()
		# Read in CPseries file as pandas df
		series_df = pd.read_table(os.path.join(args.series_directory, seriesFile), header=None)
		
		# If no end range provided, use entire sequence length
		seq_length = len(series_df.iloc[0,seq_col])
		if end <= strt:
			end = seq_length

		# Fill in list of IDs to be used as new column
		IDs = []
		total_rows = len(series_df.index)

		# Iterate through entire CPseries file:
		for row in range(total_rows):
			seq = series_df.iloc[row, seq_col][strt:end]
			# If sub-sequence has already been seen, assign existing ID
			if seq in variantDict:
				IDs.append(variantDict[seq][0])
				variantDict[seq][2] += 1	# Count how many times a variant has been seen
			else:
				newID = randID_set.pop()
				IDs.append(newID)
				variantDict[seq] = [newID, series_df.iloc[row, 1] ,1]
			# Curtis' cool progress bar:
			cpfiletools.update_progress(row, total_rows)

		# Add in new ID column: (currently puts it next to the filter column)
		series_df.insert(loc=2, column="IDs", value=IDs)
		np.savetxt(labeled_filename, series_df.values, fmt='%s', delimiter='\t')
		print "finished file: {0:.2f} seconds".format(time.time() - start)
		fileNum += 1

	# Now write a file containing the key for all the assigned IDs:
	print "Creating ID file: {}...".format(args.ID_file)
	variant_df = pd.DataFrame(variantDict).transpose()
	seqs = list(variant_df.index)
	variant_df.insert(loc=0, column="sequence", value=seqs)
	sorted_df = variant_df.sort([2, "sequence"], ascending=[False, True])	# Sort by number of variants, then by sequence
	np.savetxt(os.path.join(args.output_directory, args.ID_file), sorted_df.values, fmt='%s', delimiter='\t')
	print "Done"



def printList(lst):
	for l in lst:
		print "\t{}".format(l)

if __name__ == '__main__':
    main()