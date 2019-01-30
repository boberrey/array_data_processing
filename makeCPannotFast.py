#!/usr/bin/env python
""" 
Assign variant IDs to create a CPannot file for use in Sarah's pipeline

Variant IDs are assigned based on a provided key of variants

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
import string
import cpfiletools
import pandas as pd
import numpy as np
import time
from joblib import Parallel, delayed

import time


### Global Vars ###
read1_primers = [
"GGATCCAGGAACGTCTTCCATACAACCTCCTTACTACAT", # Stall sequence
"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  # TruSeqR1 sequence
]

# 20171103: We now have multiple R2 primers. I've modified this script to allow for
# multiple R2 primers when running in R1 mode, but maybe the smart thing is just to run
# in R2 mode from now on?
#read2_primer = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG"	# TruSeqR2
#read2_primer = "GATCGTGAGCCGGACCCAGCGTTGAGAAGAGGCAAAG"  # TtEndoR2

read2_primers = [
"AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG",  # TruSeqR2
"GATCGTGAGCCGGACCCAGCGTTGAGAAGAGGCAAAG",  # TtEndoR2
"CGGACGCGGGAAGACAGAATAGAAGTGTACGCCGTCG",  # let7-pShortR2
"CACAGGGAGCGAGGTATACCGAGAGCTGACAGCGTAG"   # let7-pLongR2
]

primer_overlap = 15
#max_seq_length = 70
max_seq_length = 50

clusterID_column = 0
r1_column = 2
r2_column = 4

# Trim bases (may get more annotations if you trim the first n bases of variant
# sequences and r1:
trim_length = 0

transtab = string.maketrans("ACGT", "TGCA")

### MAIN ###

def main():
	start = time.time()
	################ Parse input parameters ################

	#set up command line argument parser
	parser = argparse.ArgumentParser(description='Script for generating a \
		CPannot file based on previously designed variants')
	group = parser.add_argument_group('required arguments')
	group.add_argument('-sd', '--seq_directory', required=True,
		help='directory that holds the CPseq files that need variant IDs')
	group.add_argument('-vt', '--variant_table', required=True,
		help='A tab-delimited table containing the variant information \
		(first column sequence, second column variant ID)')
	group.add_argument('-r', '--read_num', type=int, required=True,
		help='which read to use for matching variants')
	group = parser.add_argument_group('optional arguments for processing data')
	group.add_argument('-od','--output_directory',
		help='output directory for series files with labeled \
		variants (default will use seq_directory)')
	group.add_argument('-p','--paired', action='store_true',
		help='flag to indicate if variants must be read in both directions (more conservative)')
	group.add_argument('-n','--num_cores', type=int, default=19,
                        help='number of cores to use')

	if not len(sys.argv) > 1:
	    parser.print_help()
	    sys.exit()

	#parse command line arguments
	args = parser.parse_args()
	numCores = args.num_cores

	# If no output directory given, use current directory
	if not args.output_directory:
		args.output_directory = "./"
	output_directory = args.output_directory
	if not os.path.isdir(output_directory):
		print "Error: invalid output directory selection. Exiting..."
		sys.exit()

	# Construct variant dict:
	print "Reading in variant dict: {}".format(args.variant_table)
	variant_dict = get_variant_dict(args.variant_table, args.read_num)

	# Find CPseqs in seq_directory:
	print "Finding CPseq files in directory: {}".format(args.seq_directory)
	CPseqFiles = cpfiletools.find_files_in_directory(args.seq_directory, ['.CPseq'])

	if numCores > 1:
		print "Annotating clusters in parallel on {} cores...".format(numCores)
		annotated_cluster_lists = (Parallel(n_jobs=numCores, verbose=10)\
			(delayed(annotate_clusters)(
				args.seq_directory + CPseq, variant_dict, args.read_num, args.paired) for CPseq in CPseqFiles))
	else:
		print "Annotating clusters on a single core"
		annotated_cluster_lists = [annotate_clusters(
			args.seq_directory + CPseq, variant_dict, args.read_num, args.paired) for CPseq in CPseqFiles]

	# Combine cluster lists:
	print "Formatting and saving CPannot file..."
	all_annotations = []
	map(all_annotations.extend, annotated_cluster_lists)
	CPannot_df = pd.DataFrame(all_annotations)
	try:
		CPannot_df.columns = ['cluster_ID', 'variant_ID']
	except:
		print "No variants annotated!"
		sys.exit()
	
	# Save the CPannot file as a pickle
	CPannotFilename = "_".join(longestSubstring(CPseqFiles).split("_")[:-1])+".CPannot.pkl"
	print "Creating CPannot.pkl file: {}...".format(CPannotFilename)
	CPannot_df = CPannot_df.set_index("cluster_ID")
	
	CPannot_df.to_pickle(output_directory+CPannotFilename)
	print "Done. {} minutes".format(round((time.time() - start)/60, 2))





def get_variant_dict(filename, read_num):
	"""
	Read in a variant table and extract the necessary information for 
	constructing the variant dict:
	Inputs:
		filename (str) - the filename for the variant dict
	Outputs:
		variant_dict (dict) - the variant dict, keyed by sequence, 
		with variant IDs as values
	"""
	variant_dict = {}
	with open(filename, 'r') as f:
		for line in f:
			split_line = line.split('\t')
			seq = split_line[0]
			variant_ID = split_line[1]
			if read_num == 1:
				if len(seq) > max_seq_length:
					seq = seq[:max_seq_length]
			else:
				seq = rev_comp(seq)
				if len(seq) > max_seq_length:
					seq = seq[:max_seq_length]

			if seq not in variant_dict:
				# Give priority to variant annotations that come first in list
				variant_dict[seq] = variant_ID
			else:
				variant_dict[seq] = variant_dict[seq] + ';' + variant_ID
	#print_dict(variant_dict)
	return variant_dict


def annotate_clusters(CPseq_filename, variant_dict, read_num, paired):
	"""
	Annotate cluster IDs with their appropriate variants
	Inputs:
		CPseq_filename (str) - the CPseq filename
		variant_dict (dict) - the variant dict
	Outputs:
		annotated_clusters (list) - list with annotated clusters
	"""
	annotated_clusters = []
	with open(CPseq_filename, 'r') as f:
		for line in f:
			split_line = line.split('\t')
			clusterID = split_line[clusterID_column]
			read1 = split_line[r1_column][trim_length:]
			read2 = split_line[r2_column][trim_length:]

			# Get the insert sequence from paired reads:
			if read_num == 1:
				#insert_seq = get_insert_seq(read1, read2, read2_primer)
				insert_seq = get_insert_seq(read1, read2, read2_primers, paired)
			else:
				#insert_seq = get_insert_seq(read2, read1, rev_comp(read1_primer))
				insert_seq = get_insert_seq(read2, read1, read1_primers, paired)

			# if insert seq not in variant dict, continue to next line
			if not insert_seq in variant_dict:
				continue

			# If still going, it means there is a match, so add that annotation
			annotated_clusters.append([clusterID, variant_dict[insert_seq]])
	return annotated_clusters

'''
def get_insert_seq(readA, readB, primer):
	"""
	Find the insert sequence of two paired reads. If no overlap is found, will
	return false
	Inputs:
		readA (str) - the read of focus
		readB (str) - the other read
	Outputs:
		insert_seq (str) - the insert sequence or the whole read if no insert
			found
	"""
	primer_seq = primer[:primer_overlap]
	insert_end = readA.find(primer_seq)
	# If the primer sequence isn't found, just return the whole read 1
	if insert_end < 0:
		return readA[:max_seq_length]
	insert = readA[:insert_end]
	revB = rev_comp(readB)
	# It seems like there is some difficulty in matching the full length thing...
	if revB.find(insert[1:-1]) > 0:
		return insert
	else:
		return readA[:max_seq_length]
'''


def get_insert_seq(readA, readB, primers, paired):
	"""
	Find the insert sequence of two paired reads. If no overlap is found, will
	return false
	Inputs:
		readA (str) - the read of focus
		readB (str) - the other read
	Outputs:
		insert_seq (str) - the insert sequence or the whole read if no insert
			found
	"""
	primer_seqs = [primer[:primer_overlap] for primer in primers]
	insert_end = -1
	for p in primer_seqs:
		insert_end = readA.find(p)
		if insert_end >= 0:
			break
	# If the primer sequence isn't found, just return the whole read 1
	if insert_end < 0:
		return readA[:max_seq_length]
	insert = readA[:insert_end]
	if paired:
		# If paired flag was indicated, check that the insert was read in both directions
		revB = rev_comp(readB)
		# It seems like there is some difficulty in matching the full length thing...
		if revB.find(insert[1:-1]) >= 0:
			return insert
		else:
			return readA[:max_seq_length]
	return insert


def rev_comp(seq):
	# Reverse complement a sequence
	return seq.translate(transtab)[::-1]

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


def print_dict(d):
	for key, val in d.items():
		print "{}:\t{}".format(key, val)

if __name__ == '__main__':
    main()
