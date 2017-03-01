#!/usr/bin/env python
""" Make empty CPfluor files corresponding to images that
were not successfully quantified. These will be used in processData.py 
to ensure unquantified data points are placed in the proper positions in a image series.

This script depends on both the images and the fluor 
files having their original timestamps!

 Inputs:
   Directory of images with timestamps
   Directory of CPfluors with timestamps

 Outputs:
   Phony CPfluor files corresponding to images that were 
   not successfully quantified 

 Ben Ober-Reynolds, boberrey@stanford.edu
 20160806
 """

import sys
import os
import argparse
import cpfiletools


### MAIN ###


def main():
	################ Parse input parameters ################

	#set up command line argument parser
	parser = argparse.ArgumentParser(description='Script for generating phony CPfluors for unquantified images')
	group = parser.add_argument_group('required arguments:')
	group.add_argument('-id', '--image_dir', required=True,
	                    help='directory that holds the all the images on which quantification was attempted (successful or not)')
	group.add_argument('-fd', '--fluor_dir', required=True,
	                    help='directory containing CPfluor files that were generated')
	group.add_argument('-sd', '--seq_dir', required=True,
						help='directory that contains the CPseq files for this experiment.')

	group = parser.add_argument_group('optional arguments for processing data')
	group.add_argument('-od','--output_dir', default="fluor_dir",
	                    help='where the output files will be saved. Default is the fluor_dir.')
	group.add_argument('-fl','--flag', default="phony",
	                    help='optional flag to be inserted at the front of phony CPfluor file names.')


	if not len(sys.argv) > 1:
	    parser.print_help()
	    sys.exit()


	#parse command line arguments
	args = parser.parse_args()

	# check that output directory is valid:
	if args.output_dir == "fluor_dir":
		output_dir = args.fluor_dir
	else:
		if os.path.isdir(args.output_dir):
			output_dir = args.output_dir
		else:
			print "Error: output directory "+args.output_dir+" is not a directory. Exiting..."
			sys.exit()

	# import fluor files
	print "Finding fluor files in directory "+args.fluor_dir+" ..."
	fluorFilenames = cpfiletools.find_files_in_directory(args.fluor_dir, ['.CPfluor'])
	if len(fluorFilenames) < 1:
		print "Error: No fluor files found in directory: " + args.fluor_dir
		sys.exit()

	# import image files
	print "Finding image files in directory "+args.image_dir+" ..."
	imageFilenames = cpfiletools.find_files_in_directory(args.image_dir, ['.tif', '.tiff'])
	if len(imageFilenames) < 1:
		print "Error: No image files found in directory: " + args.image_dir
		sys.exit()

	# find the relevant CPseq files:
	print "Finding CPseq files in directory "+args.seq_dir+" ..."
	seqFilenames = cpfiletools.find_files_in_directory(args.seq_dir, ['.CPseq'])
	if len(seqFilenames) < 1:
		print "Error: No CPseq files found in directory: " + args.seq_dir
		sys.exit()

	# Make a set of timestamps from the fluor files
	# This script assumes that no two images will have the same timestamp

	fluorTimestamps = set()
	for filename in fluorFilenames:
		fluorTimestamps.add(getTimestamp(filename))

	# Now identify which images do not have corresponding CPfluor files:
	lonelyImageFiles = []
	for filename in imageFilenames:
		timestamp = getTimestamp(filename)
		if timestamp not in fluorTimestamps:
			lonelyImageFiles.append(filename)

	if len(lonelyImageFiles) < 1:
		print "No need for phony files. Exiting..."
		sys.exit()

	# Make a CPseq dict keyed by tile number:
	seq_dict = cpfiletools.make_tile_dict(seqFilenames, args.seq_dir)

	# Now make the new phony files
	for filename in lonelyImageFiles:
		root, ext = os.path.splitext(filename)
		newFluorName = args.flag + filename.strip(ext) + ".CPfluor"
		# find the CPseq file relevant to this image:
		tile = cpfiletools.get_tile_number_from_filename(filename)
		cpseq = seq_dict[tile]
		with open(output_dir+'/'+newFluorName, 'w') as outfile, open(cpseq, 'r') as infile:
			for line in infile:
				cluster_ID = line.split()[0]
				outfile.write(cluster_ID+':0:0.000000:0.000000:0.000000:0.000000\n')
		print "Generated phony file: "+newFluorName

		


def getTimestamp(filename):
	# Return the timestamp from the end of a filename as a string
	# (assumes that the timestamp is the last part of the filename
	# and is separated from the rest of the name by an '_')
	root, ext = os.path.splitext(filename)
	try:
		timestamp=filename.strip(ext).split('_')[-1]
	except ValueError:
		print "ERROR: no timestamp on file: {}".format(filename)
		sys.exit()
	return timestamp









if __name__ == '__main__':
    main()