#!/usr/bin/env python

"""
Script for running MATLAB analysis of multiple tiles (all RNA images)

Inputs:


Outputs:


Ben Ober-Reynolds
"""



import os
import sys
import argparse
import re
import cpfiletools

# Parallelization
from joblib import Parallel, delayed
import multiprocessing



def main():
	parser = argparse.ArgumentParser(description="script for parallelizing the matlab 'AnalyseImage' script")
	group = parser.add_argument_group('required arguments for processing data')
	group.add_argument('-cd', '--CPseq_directory', required=True,
		help='directory that holds the .CPseq files')
	group.add_argument('-id', '--image_files_directory', required=True,
		help='directory that holds the image files')
	group.add_argument('-n', '--numCores', required=True,
		help='number of cores to use (recommended same as number of tiles')
	group = parser.add_argument_group('optional arguments for processing data')
	group.add_argument('-if', '--image_flag', default="",
		help='flag to designate which images to use')
	group.add_argument('-fs', '--filterSubsets', default="",
		help='filter subsets to use in data registration')
	group.add_argument('-od', '--output_directory', default="",
		help='output directory for CPfluor files (if not provided, uses image_files_directory)')


	# print help if no arguments provided
	if len(sys.argv) <= 1:
		parser.print_help()
		sys.exit()

	#parse command line arguments
	args = parser.parse_args()
	numCores = int(args.numCores)

	# Get the absolute paths for everything
	args.CPseq_directory = os.path.abspath(args.CPseq_directory)+"/"
	args.image_files_directory = os.path.abspath(args.image_files_directory)+"/"
	if args.output_directory == "":
		args.output_directory = args.image_files_directory
	else:
		args.output_directory = os.path.abspath(args.output_directory)+"/"

	# temporary paths to use while pipeline migration is incomplete
	tempPaths = [
	"/raid1/lab/ben/array_data/miR-21_oligo_array_expts/imagesForAnalysis/scripts/temp_MATLAB",
	"/raid1/lab/ben/array_data/miR-21_oligo_array_expts/imagesForAnalysis/scripts/temp_MATLAB/CPlibs",
	"/raid1/lab/ben/array_data/miR-21_oligo_array_expts/imagesForAnalysis/scripts/temp_MATLAB/CPscripts"
	]

	CPseqFiles = cpfiletools.find_files_in_directory(args.CPseq_directory, ['CPseq'])
	
	# Filter images to use according to provided flag
	allImages = cpfiletools.find_files_in_directory(args.image_files_directory, ['tif'])
	imageFiles = []
	if args.image_flag != "":
		for filename in allImages:
			if args.image_flag in filename:
				imageFiles.append(filename)
	else:
		imageFiles = allImages

	print "using image files:"
	print_list(sorted(imageFiles))

	# Create CPseq and imageFile dictionaries according to tile
	CPseqDict = {}
	for CPseqFile in CPseqFiles:
		tile = cpfiletools.get_tile_number_from_filename(CPseqFile)
		CPseqDict[tile] = args.CPseq_directory + CPseqFile
	imageDict = {}
	for imageFile in imageFiles:
		tile = cpfiletools.get_tile_number_from_filename(imageFile)
		imageDict[tile] = args.image_files_directory + imageFile
	allTiles = sorted(imageDict.keys())


	dataScaling = 'MiSeq_to_TIRFStation1'
	filterSubsets = "''"
	if args.filterSubsets != "":
		filterSubsets = "{{'{}'}}".format(args.filterSubsets)

	workingPath = args.output_directory

	print "CPseqDict = "
	print CPseqDict
	print "imageDict = "
	print imageDict

	# Parallelize running get_registration_offset from cpfiletools (which spawns matlab jobs)
	(Parallel(n_jobs=len(imageDict.keys()), verbose=10)
		(delayed(cpfiletools.get_registration_offset)(CPseqDict[tile], imageDict[tile],
			dataScaling, filterSubsets, tempPaths) for tile in sorted(imageDict.keys())))

	# Gather .roff files
	allRoffFiles = cpfiletools.find_files_in_directory(args.image_files_directory, ['roff'])
	roffDict = {}
	for roffFile in allRoffFiles:
		tile = cpfiletools.get_tile_number_from_filename(roffFile)
		roffDict[tile] = args.image_files_directory + roffFile
	print roffDict


	# Parallelize running analyse_image from cpfiletools (which spawns matlab jobs)
	(Parallel(n_jobs=len(roffDict.keys()), verbose=10)
		(delayed(cpfiletools.analyse_image)(CPseqDict[tile],imageDict[tile], 
			dataScaling, filterSubsets, roffDict[tile], tempPaths) for tile in sorted(roffDict.keys())))




def print_list(lst):
	for i in lst:
		print i



if __name__ == '__main__':
	main()

