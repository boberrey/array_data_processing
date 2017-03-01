#!/usr/bin/env python
""" 
Assign variant IDs to create a CPannot file for use in Sarah's pipeline

 Variant IDs are assigned based on matching sequence. 
 The user may specify which bases of which sequence are to be used for assigning variants.

 Inputs:
   CPseq files

 Outputs:
   CPannot file
   File containing list of variant IDs and associated sequences (ID file)

 Ben Ober-Reynolds, boberrey@stanford.edu
 20160816
 """

import sys
import os
import argparse
import cpfiletools
import random
import string
import subprocess
import pandas as pd
import numpy as np
import re

import time


### MAIN ###


def main():
	################ Parse input parameters ################

	#set up command line argument parser
	parser = argparse.ArgumentParser(description='Script for assigning unique IDs to variants in CPseq files')
	group = parser.add_argument_group('required arguments')
	group.add_argument('-sd', '--seq_directory', required=True,
	                    help='directory that holds the CPseq files that need variant IDs')
	group.add_argument('-sc', '--seq_column', required=True,
	                    help='which column in the CPseq file you want to use for assigning variants')

	group = parser.add_argument_group('optional arguments for processing data')
	group.add_argument('-pi','--previous_ID_file', default="",
	                    help='An ID file previously created for variants expected in the new CPseq files')
	group.add_argument('-fk','--filters_to_use', default="",
	                    help='Which filters should be kept. Separate by commas: filter1,filter2,filter3,etc. If you want to use clusters without a filter, include "blank" (filter1,filter2,blank). Default is to use all filters.')
	group.add_argument('-st','--seq_start', default=0,
	                    help='start position within sequence for matching. Will use beginning of sequence if none specified.')
	group.add_argument('-ed','--seq_end', default=0,
	                    help='end position within sequence for matching. Will use end of sequence if none specified.')
	group.add_argument('-lb','--label', default="ID_ed",
	                    help='label attached to output files. Default is "ID_ed"')
	group.add_argument('-od','--output_directory', default="",
	                    help='output directory for series files with labeled variants (default will use seq_directory)')
	group.add_argument('-if','--ID_file', default="ID_file.txt",
	                    help='file name for the list of IDs and corresponding sequences. Default is "ID_file.txt"')

	if not len(sys.argv) > 1:
	    parser.print_help()
	    sys.exit()


	#parse command line arguments
	args = parser.parse_args()

	# If no output directory given, use input directory
	if args.output_directory == "":
		args.output_directory = args.seq_directory
	output_directory = args.output_directory
	if not os.path.isdir(output_directory):
		print "Error: invalid output directory selection. Exiting..."
		sys.exit()

	# Crate a set of filters to be kept:
	filters = set(args.filters_to_use.split(','))
	if "blank" in filters:
		filters.remove("blank")
		filters.add("no_filter")	#I coerce the pandas dataframes to contain 'no_filter' instead of NaN's
	if filters.pop() == "":
		filters.add("all")
	print filters
	# This script will run through each of the provided CPseq files sequentially in order to 
	# ensure that each variant gets assigned only one variant ID.

	CPseqFiles = cpfiletools.find_files_in_directory(args.seq_directory, ['.CPseq'])

	numLines = 0
	for seqFile in CPseqFiles:
		numLines += int(subprocess.check_output(("wc -l {} | ".format(
			os.path.join(args.seq_directory, seqFile))+" awk \'{print $1}\'"), shell=True).strip())

	start = time.time()
	randID_set = set()
	print "Generating random IDs..."
	while len(randID_set) < numLines:
		randID = ''.join([random.choice(string.ascii_uppercase + string.digits) for n in range(8)])	# 36^8 ~ 2.8e12 possible IDs
		randID_set.add(randID)
	print "ID generation: {0:.2f} seconds".format(time.time() - start)


	# This dictionary will contain all the variants assigned, keyed by sequence match
	# The entries in variant dict will be three-element lists, the first is the ID, the second is the filter
	# associated with that variant (if any), and the third is the number of times that variant has been seen
	variantDict = {}

	# If a previous ID file was provided, it will pre-populate the variantDict.
	# Note: it is up to the user to ensure that seq_column, seq_start and seq_end match those used to 
	# create the previous ID file!
	if args.previous_ID_file != "":
		with open(args.previous_ID_file, 'r') as f:
			for line in f:
				seq, ID, filtr, n = line.split()
				variantDict[seq] = [ID, filtr, int(n)]

	
	fileNum = 1
	CPannot_df = pd.DataFrame()

	# Loop through each CPseq file to assign variants:
	for seqFile in CPseqFiles:
		print "Working on file: {}...{} of {}".format(seqFile, fileNum, len(CPseqFiles))
		# Time each loop for now:
		start = time.time()
		# Read in CPseq file as pandas df
		seq_df = pd.read_table(os.path.join(args.seq_directory, seqFile), header=None)
		seq_df = seq_df.fillna('no_filter')

		print "length pre-filter " + str(len(seq_df))
		# filter df by filters to keep (if any)
		if "all" not in filters:
			seq_df = seq_df[seq_df.iloc[:,1].isin(filters)]
		print "length post-filter " + str(len(seq_df))
		
		# set sequence selection parameters:
		seq_col = int(args.seq_column) - 1	# Allow for intuitive column selection (i.e. start at 1)
		if seq_col < 0 or seq_col > len(seq_df.columns):
			print "Error: invalid seq column selected. Out of range. Must be within {} and {}".format(1, len(seq_df.columns))
			sys.exit()
		
		# Test to ensure provided column contains sequence data:
		test_seq = seq_df.iloc[0,seq_col]
		if not re.match("^[a-zA-Z]+$", test_seq):
			print "Error: provided column does not contain sequence data, e.g. {}".format(test_seq)
			sys.exit()

		# Test to ensure start and end sequence positions are valid:
		seq_length = len(seq_df.iloc[0,seq_col])
		strt = int(args.seq_start)
		if strt < 0 or strt > seq_length - 1:
			print "Error: invalid start position selected. Must be positive and less than seq length"
			sys.exit()
		end = int(args.seq_end)
		if end < strt or end > seq_length:
			print "Error: invalid end position selected. Must be greater than start position and <= seq length"
			sys.exit()
			
		# If no end range provided, use entire sequence length
		if end == 0:
			end = seq_length

		# Fill in list of IDs to be used as new column
		clusterIDs = []
		IDs = []
		total_rows = len(seq_df.index)

		# Iterate through entire CPseq file:
		for row in range(total_rows):
			seq = seq_df.iloc[row, seq_col][strt:end]
			# If sub-sequence has already been seen, assign existing ID
			if seq in variantDict:
				IDs.append(variantDict[seq][0])
				variantDict[seq][2] += 1	# Count how many times a variant has been seen
			else:
				newID = randID_set.pop()
				IDs.append(newID)
				variantDict[seq] = [newID, seq_df.iloc[row, 1] ,1]
			clusterIDs.append(seq_df.iloc[row, 0])
			# Curtis' cool progress bar:
			cpfiletools.update_progress(row, total_rows)

		# Start making the CPannot file:
		if fileNum == 1:
			CPannot_df = pd.DataFrame({"cluster_ID":clusterIDs, "variant_ID":IDs})
		else:
			CPannot_df = pd.concat([CPannot_df, pd.DataFrame({"cluster_ID":clusterIDs, "variant_ID":IDs})])

		print "finished file: {0:.2f} seconds".format(time.time() - start)
		fileNum += 1

	# Save the CPannot file as a pickle
	CPannotFilename = "_".join(longestSubstring(CPseqFiles).split("_")[:-1])+".CPannot.pkl"
	print "Creating CPannot.pkl file: {}...".format(CPannotFilename)
	CPannot_df = CPannot_df.set_index("cluster_ID")
	
	CPannot_df.to_pickle(args.output_directory+CPannotFilename)


	
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

def longestSubstring(lst):
	# Return the longest substring shared by a list of strings
	# Note: 'longest substring' is a famous CS problem, this function
	# is simplified in that matches must begin at the beginning of each string
	# (and this is probably not the most elegant solution either...)
	substr = ""
	match = True
	while match:
		letter_to_match = lst[0][0]
		matches = []
		for index in range(len(lst)):
			if len(lst[index]) >= 1:
				letter_in_question = lst[index][0]
			else:
				match = False
				break
			if len(lst[index]) > 1:
				lst[index] = lst[index][1:]
			else:
				match = False
				break
			if letter_in_question == letter_to_match:
				matches.append(True)
			else:
				matches.append(False)
		if all(matches):
			substr = substr + letter_to_match
		else:
			match = False
	return substr


if __name__ == '__main__':
    main()
