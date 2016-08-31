#!/usr/bin/env python
""" Translate a CPannot file into a new file that matches variant IDs to sequences

This script will take a .CPannot file and a directory of CPseq files, then 
create a new file that relates the variant IDs in the .CPannot file to the 
actual sequences in the CPseq file.

Inputs:
  .CPannot.pkl file
  Directory of CPseq files

Outputs:
  CPseries files with variants labeled
  File containing list of variant IDs and associated sequences

Ben Ober-Reynolds"""

import sys
import os
import argparse
import cpfiletools
import string
import pandas as pd
import numpy as np
import re
import time
from fittinglibs import fileio
from joblib import Parallel, delayed



### MAIN ###


def main():
	################ Parse input parameters ################

	#set up command line argument parser
	parser = argparse.ArgumentParser(description='script for relating variantIDs from CPannot file to sequence')
	group = parser.add_argument_group('required arguments')
	group.add_argument('-a', '--annot_file', required=True,
						help='A .CPannot.pkl file')
	group.add_argument('-sd', '--seq_dir', required=True,
						help='A directory of .CPseq files')

	group = parser.add_argument_group('optional arguments for running script')
	group.add_argument('-l','--length', default="short",
						help='translate in "long" or "short" format: long will include every cluster in CPannot file, short will only show variantID and one sequence match. Default = short')
	group.add_argument('-sc','--seq_cols', default="3",
						help='Which sequence columns to output in CPtranslate file. May use multiple columns for long format (seqarate by commas). key: 3 = r1, 5 = r2, 7 = i7, 9 = i5')
	group.add_argument('-od','--output_dir', default=os.getcwd(),
						help='Output directory. default is current directory')
	group.add_argument('-n','--num_cores', default=1,
						help='How many cores to use for parallel processing')

	if not len(sys.argv) > 1:
		parser.print_help()
		sys.exit()

	##### parse command line arguments #####
	args = parser.parse_args()

	annot_file = args.annot_file
	seq_dir = args.seq_dir

	length = args.length
	if length != "short" and length != "long":
		print "Error: length must be either 'short' or 'long'. Exiting..."
		sys.exit()

	seq_cols = [0] + [int(n) - 1 for n in args.seq_cols.split(',')]
	seq_col_names = assign_names(seq_cols)

	output_dir = args.output_dir
	if not os.path.isdir(output_dir):
		print "Error: output directory is invalid. Exiting..."
		sys.exit()
	if output_dir[-1] != '/':
		output_dir = output_dir + '/'

	num_cores = int(args.num_cores)

	########################################


	# Read in CPannot file
	print "Reading in CPannot file..."
	start = time.time()
	annot_df = fileio.loadFile(annot_file)
	print "file loaded: {0:.2f} seconds\n".format(time.time() - start)

	# Read in CPseq files as a concatenated data frame
	print "Reading in CPseq files..."
	start = time.time()
	seq_files = cpfiletools.find_files_in_directory(seq_dir, ['.CPseq'])
	print "found CPseq files: "
	cpfiletools.printList(seq_files)
	seq_df = pd.DataFrame()
	for seq_file in seq_files:
		new_df = pd.read_csv(seq_dir+seq_file, sep='\t', index_col=0, usecols=seq_cols, header=None)
		seq_df = pd.concat([seq_df, new_df])
	seq_df.columns = seq_col_names
	print str(len(seq_files)) + " files loaded: {0:.2f} seconds\n".format(time.time() - start)
	
	# Merge the data frames 
	print "Merging data frames..."
	start = time.time()
	merged_df = annot_df.merge(seq_df, how='left', left_index=True, right_index=True)
	print "Merged: {0:.2f} seconds\n".format(time.time() - start)


	# Save long format CPtranslate if requested
	if length == "long":
		print "Saving long format CPtranslate.pkl..."
		start = time.time()
		filename = os.path.basename(annot_file).rstrip('.CPannot.pkl')+".long.CPtranslate.pkl"
		print "filename = "+filename
		merged_df.to_pickle(output_dir+filename)
		print "Saved: {0:.2f} seconds\n".format(time.time() - start)

	# Create the short format CPtranslate:
	if length == "short":
		print "Generating short format CPtranslate..."
		start = time.time()
		# Make a list of unique variant_IDs

		grouped_variants = merged_df.groupby('variant_ID')

		all_variants = (Parallel(n_jobs=num_cores, verbose=10)
						(delayed(fetch_sequences)(name, group, seq_col_names) for name, group in grouped_variants))

		short_df = pd.DataFrame(all_variants)
		short_df.columns = ['variant_ID', 'count']+seq_col_names
		print "short format generated: {0:.2f} seconds\n".format(time.time() - start)
		print short_df.head()

		print "Saving short format CPtranslate.pkl..."
		start = time.time()
		filename = os.path.basename(annot_file).rstrip('.CPannot.pkl')+".short.CPtranslate.pkl"
		short_df.to_pickle(output_dir+filename)
		print "Saved: {0:.2f} seconds\n".format(time.time() - start)



		
def fetch_sequences(variant_ID, variant_df, seq_col_names):
	"""
	Gather the consensus sequences from a given variant_IDs
	Inputs: 
	variant_df: Data frame of a single variant: DataFrame
	seq_col_names: column names of sequences: list
	Output:
	list containing variant_ID and all associated consensus sequences
	"""
	variant_info = [variant_ID, len(variant_df.index)]
	for seq_name in seq_col_names:
		seqs = variant_df[seq_name].tolist()
		consensus = dumb_consensus(seqs)
		variant_info.append(consensus)
	return variant_info



def dumb_consensus(seq_list, threshold=0.7):
	"""
	A very simple consensus sequence finder. Assumes that all sequences are same 
	length and already 'aligned'
	"""
	consensus = ""
	seq_len = len(seq_list[0])
	for index in range(seq_len):
		bases = []
		for seq in seq_list:
			bases.append(seq[index])
		A = float(bases.count('A'))
		G = float(bases.count('G'))
		C = float(bases.count('C'))
		T = float(bases.count('T'))
		N = float(bases.count('N'))
		total = A + G + C + T + N
		if A / total > threshold:
			consensus = consensus + "A"
		elif G / total > threshold:
			consensus = consensus + "G"
		elif C / total > threshold:
			consensus = consensus + "C"
		elif T / total > threshold:
			consensus = consensus + "T"
		else:
			consensus = consensus + "N"
	return consensus




def assign_names(lst):
	names = []
	count = 1
	for col in lst:
		if col == 0:
			pass
		elif col == 2:
			names.append('read_1')
		elif col == 4:
			names.append('read_2')
		elif col == 6:
			names.append('index_i7')
		elif col == 8:
			names.append('index_i5')
		else:
			names.append('unknown_{}'.format(count))
			count += 1
	return names





if __name__ == '__main__':
    main()