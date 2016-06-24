#!/usr/bin/env python

"""
Script for running MATLAB analysis of multiple tiles

NOTE: This script needs a lot of work to get it more usable. See analyseTilesRNA.py
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



def main():
	# This script is going to be very ugly for now. 
	numCores = 18

	tempPaths = [
	"/raid1/lab/ben/array_data/miR-21_oligo_array_expts/imagesForAnalysis/scripts/temp_MATLAB",
	"/raid1/lab/ben/array_data/miR-21_oligo_array_expts/imagesForAnalysis/scripts/temp_MATLAB/CPlibs",
	"/raid1/lab/ben/array_data/miR-21_oligo_array_expts/imagesForAnalysis/scripts/temp_MATLAB/CPscripts"
	]

	CPseqFiles = [
	"/raid1/lab/ben/array_data/miR-21_oligo_array_expts/imagesForAnalysis/seqdata/AHGTK_ALL_tile{:03}_Bottom_filtered.CPseq".format(i) for i in range(1,19)
	]

	imageListFiles = [
	"/raid1/lab/ben/array_data/miR-21_oligo_array_expts/imagesForAnalysis/8nM_all_tile{:03}.ipf".format(i) for i in range(1,19)
	]

	darkImageIndex = 1
	registrationImageIndex = 10
	dataScaling = 'MiSeq_to_TIRFStation1'
	filterSubsets = "{'PERFECT_C'}"
	workingPath = "/raid1/lab/ben/array_data/miR-21_oligo_array_expts/imagesForAnalysis/CPfluors"

	print "tempPaths:"
	print_list(tempPaths)
	print "CPseqFiles:"
	print_list(CPseqFiles)
	print "imageListFiles:"
	print_list(imageListFiles)

	# Parallelize running analyse_series from cpfiletools (which spawns matlab jobs)
	(Parallel(n_jobs=numCores, verbose=10)
		(delayed(cpfiletools.analyse_series)(CPseqFiles[i],imageListFiles[i], 
		1, 11, dataScaling, filterSubsets, workingPath, tempPaths) for i in range(18)))
		



def print_list(lst):
		for i in lst:
			print i



if __name__ == '__main__':
    main()

