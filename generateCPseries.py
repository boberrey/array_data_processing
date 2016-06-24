#!/usr/bin/env python
""" Convert CPfluors into CPseries.

 This script requires you to already have quantified images with another pipeline.

 Inputs:
   Sequence data (.CPseq files)
   CPfluor files (.CPfluor files or directories)

 Outputs:
   CPseries files 

 Ben Ober-Reynolds, modified from a script by Sarah Denny """

import sys
import os
import argparse
import cpfiletools
from joblib import Parallel, delayed


### MAIN ###


def main():
	################ Parse input parameters ################

	#set up command line argument parser
	parser = argparse.ArgumentParser(description='Script for generating CPseries files from CPseq and CPfluor files')
	group = parser.add_argument_group('required arguments for processing data')
	group.add_argument('-fs', '--filtered_CPseqs', required=True,
	                    help='directory that holds the filtered sequence data (CPseq)')
	group.add_argument('-bs', '--bsCPfluors', required=True,
	                    help='directory containing binding series CPfluor files')

	group = parser.add_argument_group('optional arguments for processing data')
	group.add_argument('-od','--output_dir', default="CPseries",
	                    help='save output files to here. default = ./CPseries')
	group.add_argument('-ar','--allRNA', default='',
	                    help='directory containing allRNA CPfluor files')
	group = parser.add_argument_group('other settings')
	group.add_argument('-n','--num_cores', type=int, default=20,
	                    help='maximum number of cores to use. default=20')

	if not len(sys.argv) > 1:
	    parser.print_help()
	    sys.exit()


	#parse command line arguments
	args = parser.parse_args()
	numCores = args.num_cores

	# import CPseq filtered files split by tile
	print "Finding CPseq files in directory {}...".format(args.filtered_CPseqs)

	# Gather all of the CPseq files in the 'filtered_CPseqs' file directory
	CPseqFilenames = cpfiletools.find_files_in_directory(args.filtered_CPseqs, ['CPseq'])
	if len(CPseqFilenames) < 1:
		print "Error: No CPseq files found in directory: " + args.filtered_CPseqs
		sys.exit()

	print "Found CPseq files: "
	printList(CPseqFilenames)
	# Create a dictionary of the CPseq files keyed by tile
	CPseqDict = cpfiletools.make_tile_dict(CPseqFilenames, args.filtered_CPseqs)
	tileList = CPseqDict.keys()

	# Gather all of the CPfluor files for all RNA images, if provided
	allRNA_Dict = {}
	if args.allRNA != '':
		print "Finding allRNA CPfluor files in directory {}...".format(args.allRNA)
		allRNAfilenames = cpfiletools.find_files_in_directory(args.allRNA, ['CPfluor'])
		print "Found allRNA files: "
		printList(allRNAfilenames)
		if len(allRNAfilenames) < 1:
			print "Error: no CPfluor files found in directory: " + args.allRNA
		allRNA_Dict = cpfiletools.make_tile_dict(allRNAfilenames, args.allRNA)
	else:
		for tile in tileList:
			allRNA_Dict[tile] = ''

	# Gather all of the CPfluor files for creating the cluster binding series
	print "Finding binding series CPfluor files in directory {}...".format(args.bsCPfluors)
	bindingSeriesList = cpfiletools.find_files_in_directory(args.bsCPfluors, ['CPfluor'])
	print "Found CPfluor files: "
	printList(bindingSeriesList)
	bindingSeriesDict = cpfiletools.make_tile_dict_multiple(bindingSeriesList, args.bsCPfluors)



	# Make sure output directory is ready:
	outputDirectory = args.output_dir
	if os.path.isdir(outputDirectory):
		print "Output directory {} already exists".format(outputDirectory)
	else:
		outputDirectory = os.path.join(os.getcwd(), outputDirectory)
		print "Making output directory: {}".format(outputDirectory)
		os.makedirs(outputDirectory)

	# Make CPseries files
	
	CPseriesDict = {}
	for tile, fileName in CPseqDict.items():
		path, baseFile = os.path.split(fileName)
		CPseriesDict[tile] = os.path.join(outputDirectory, baseFile.split('.')[0]+'.CPseries')
	
	# Make CPseries files in parallel:
	print "Making CPseries files..."
	(Parallel(n_jobs=numCores, verbose = 10)
		(delayed(cpfiletools.generate_CPseries_files)
			(CPseqDict[tile], 
			allRNA_Dict[tile], 
			bindingSeriesDict[tile], 
			CPseriesDict[tile], 
			tile)
		for i, tile in enumerate(tileList)))
	print "Done"
	


def printList(lst):
	for l in lst:
		print "\t{}".format(l)



if __name__ == '__main__':
    main()