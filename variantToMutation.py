#!/usr/bin/env python
""" 
Create a new file relating variant ID's to a mutation notation relative to a wild-type sequence. 
Currently, this will only look for point mutants, no indels.

Really this should be some alignment-based check, but I'm going to do it fast and dirty for now.

 Variant IDs are assigned based on matching sequence. 
 The user may specify which bases of which sequence are to be used for assigning variants.

 Inputs:
   File containing variants related to sequences (e.g. a CPtranslate file)
   A consensus or WT sequence to relate mutant to
   Optional: A start and end position in the sequence that should match the consensus sequence

 Outputs:
   File relating variant IDs to mutations from the consensus sequence

 Ben Ober-Reynolds, boberrey@stanford.edu
 20161027
 """

import sys
import os
import argparse
import pandas as pd
from Bio import pairwise2
from joblib import Parallel, delayed

import time


### MAIN ###

def main():
	################ Parse input parameters ################

	#set up command line argument parser
	parser = argparse.ArgumentParser(description='Script for translating sequences into mutation annotations from a consensus sequence')
	group = parser.add_argument_group('required arguments')
	group.add_argument('-ct', '--CP_translate_file', required=True,
	                    help='A pickle file that relates variants to sequences')
	group.add_argument('-cs', '--consensus_sequence', required=True,
	                    help='The sequence you want to relate mutants to')

	group = parser.add_argument_group('optional arguments for processing data')
	# May not need to specify a start and end position since we're doing this with alignments
	group.add_argument('-st','--seq_start', default=0, type=int,
	                    help='start position within sequence for matching. Will use beginning of sequence if none specified.')
	group.add_argument('-ed','--seq_end', default=0, type=int,
	                    help='end position within sequence for matching (inclusive). Will use end of sequence if none specified.')
	group.add_argument('-od','--output_directory', default="",
	                    help='output directory for series files with labeled variants (default will use current directory)')
	group.add_argument('-of','--output_file', default="",
	                    help='file name for the mutation annotation file. Default is input file name + ".CPmut"')
	group.add_argument('-gp','--gap_penalty', default=1, type=float, metavar="GP",
	                    help='The gap penalty to use when aligning mutant sequences to the consensus sequence. Default is 5.')
	group.add_argument('-ep','--extension_penalty', default=0.5, type=float, metavar="EP",
	                    help='The extension penalty to use when aligning mutant sequences to the consensus sequence. Default is 0.5.')
	group.add_argument('-cd','--counting_direction', default="f",
	                    help='The direction to count mutations from in the consensus sequence ["f" for forward or "r" for reverse]. Default is "f"')
	group.add_argument('-n', '--numCores', default=20, type=int, metavar="N",
                   help='number of cores. default = 20')

	if not len(sys.argv) > 1:
	    parser.print_help()
	    sys.exit()

	start = time.time()
	print "Loading data..."
	################################################
	### Parse, error check, and format arguments ###
	################################################
	args = parser.parse_args()
	translate_df = pd.read_pickle(args.CP_translate_file)
	# trim to include just the variant_ID and the sequence:
	translate_df = translate_df.iloc[:,[0,2]]
	# Trim sequences if given range:
	if args.seq_end != 0:
		print "Trimming sequences to index range {} to {}...".format(args.seq_start, args.seq_end)
		translate_df.iloc[:,1] = translate_df.iloc[:,1].apply(lambda x: x[args.seq_start:args.seq_end + 1])
	consensus_sequence = args.consensus_sequence
	output_dir = args.output_directory
	if output_dir == "":
		output_dir = os.getcwd()
	if not os.path.isdir(output_dir):
		print "given output directory is not a directory. Exiting..."
		sys.exit()
	output_file = args.output_file
	if output_file == "":
		output_file = os.path.basename(args.CP_translate_file)
		# need to split extension twice, since it's a pickle file...
		output_file = os.path.splitext(output_file)[0]
		output_file = os.path.splitext(output_file)[0]
	output_file = output_dir + "/" + output_file + ".CPmut.pkl"
	gap_penalty = -args.gap_penalty
	extension_penalty = -args.extension_penalty
	counting_direction = args.counting_direction
	if counting_direction != "f" and counting_direction != "r":
		print "Invalid counting direction, must be either 'f' or 'r'. Exiting..."
		sys.exit()
	numCores = args.numCores

	################################################

	# First, convert the data frame into a two-element list for easier parallelization

	# zip converts zip([1, 3, 5], [2, 4, 6]) = [(1,2), (3,4), (5,6)]
	variant_seq_list = zip(translate_df.iloc[:,0].tolist(), translate_df.iloc[:,1].tolist())

	# Annotate in parallel:
	print "{} variants to annotate...".format(len(variant_seq_list))
	print "generating annotations in parallel, cores = {}".format(numCores)
	annotations = (Parallel(n_jobs=numCores, verbose=10)
					(delayed(unpackAndAnnotate)
						(variant_seq_list[i], consensus_sequence, gap_penalty, extension_penalty, counting_direction)
							for i in range(len(variant_seq_list))))

	# Converting from a list of tuples into a dataframe is fast and easy:
	annotation_df = pd.DataFrame(annotations, columns=['variant_ID', 'mutations', 'num_mutations'])
	annotation_df.set_index('variant_ID', inplace=True)

	print "Annotated {} mutations in {} seconds.".format(len(annotation_df.index), time.time() - start)
	print "Saving data..."

	annotation_df.to_pickle(output_file)
	print "Done."








def unpackAndAnnotate(var_seq_tuple, consensus_sequence, gap_penalty, extension_penalty, counting_direction):
	# Wrapper function to unpack the tuples being fed into the parallelization.
	variant_ID, mutant_sequence = var_seq_tuple
	annotation, num_mutations = annotateMutation(consensus_sequence, mutant_sequence, gap_penalty, extension_penalty, counting_direction)
	return (variant_ID, annotation, num_mutations)



def annotateMutation(consensus_sequence, mutant_sequence, gap_penalty, extension_penalty, counting_direction):
	"""
	Create a mutation annotation for a mutant sequence to its consensus sequence.
	Inputs: 
		consensus_sequence: a string containing the consensus sequence in 5'-3' orientation
		mutant_sequence: a string containing the mutant sequence in 5'-3' orientation
	outputs:
		A two element list: first element is the mutation notations, second is the number of mutations

	"""
	alignment = pairwise2.align.localms(consensus_sequence, mutant_sequence, 1, 0, gap_penalty, extension_penalty, one_alignment_only=True)
	# This segment was running into issues when it incountered sequences that were so wrong it couldn't align at all.
	# My attempt at catching those errors:
	try:
		aligned_consensus, aligned_mutant = alignment[0][0:2]
	except IndexError:
		return ("NO_ALIGNMENT", float('nan'))
	else:
		trimmed_consensus, trimmed_mutant = trimEndGaps(aligned_consensus, aligned_mutant)
		return findMutations(trimmed_consensus, trimmed_mutant, counting_direction)
	"""
	aligned_consensus, aligned_mutant = alignment[0][0:2]
	trimmed_consensus, trimmed_mutant = trimEndGaps(aligned_consensus, aligned_mutant)
	return findMutations(trimmed_consensus, trimmed_mutant, counting_direction)
	"""



def trimEndGaps(aligned_consensus, aligned_mutant):
	"""
	Remove end gaps from both the aligned consensus sequence and the aligned mutant sequence
	Imputs: 
		aligned_consensus: the alignment string for the consensus_sequence
		aligned_mutant: the alignment string for the mutant sequences
	Outputs:
		A two element tuple: first element is the trimmed_consensus, second is the trimmed_mutant
	"""
	n_leading_gaps = 0
	n_trailing_gaps = 0
	while aligned_consensus[0] == "-":
		n_leading_gaps += 1
		aligned_consensus = aligned_consensus[1:]
	while aligned_consensus[-1] == "-":
		n_trailing_gaps += 1
		aligned_consensus = aligned_consensus[:-1]
	trimmed_consensus = aligned_consensus
	trimmed_mutant = aligned_mutant[n_leading_gaps:len(aligned_mutant)-n_trailing_gaps]
	return trimmed_consensus, trimmed_mutant


def findMutations(trimmed_consensus, trimmed_mutant, counting_direction):
	"""
	Right now it only identifies point mutants. Obviously that isn't ideal. Fix this soon.
	Also, for short alignments, need a way to prevent it from reporting spurious alignments. 
	What's a fair way to do this? Discard annotations where >1/2 of the consensus is mutant?
	Inputs: 
		trimmed_consensus: the trimmed consensus sequence
		trimmed_mutant: the trimmed mutant sequence
	Output:
		A two element tuple: first element is the mutation annotations, second is the number of mutations
	"""
	mutations = ""
	count = 0 
	if counting_direction == "r":
		trimmed_consensus = invertString(trimmed_consensus)
		trimmed_mutant = invertString(trimmed_mutant)
	for i in range(len(trimmed_consensus)):
		consensus = trimmed_consensus[i]
		mutant = trimmed_mutant[i]
		if mutant != consensus:
			count += 1
			# Currently will count first base as 'base 1'
			mutations = mutations + str(i+1) + mutant + ":"
	if count > len(trimmed_consensus)/2:
		return "UPPER_LIM", float('nan')
	else:
		# Trim off the last ':'
		return mutations[:-1], count


def invertString(string):
	# invert a string
	inverted = ""
	while string != "":
		l = string[0]
		inverted = l + inverted
		string = string[1:]
	return inverted




if __name__ == '__main__':
    main()