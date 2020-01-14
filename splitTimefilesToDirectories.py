#!/usr/bin/env python
""" Split a folder of images into multiple directories 
for use in classic quantification pipeline

Currently will assume that you have as the same number
of images for each tile. 

 Inputs:
   Directory with image files(.tif/.tiff files)

 Outputs:
   Image files moved into new directories such that
   each new directory has only one image per tile

 Ben Ober-Reynolds, boberrey@stanford.edu 
 Started 20160802, last changed 20160830
 """

import sys
import os
import argparse
import cpfiletools

def main():
	################ Parse input parameters ################

	#set up command line argument parser
	parser = argparse.ArgumentParser(description='Script for splitting a directory of timed files (e.g. images or CPfluors) into multiple directories')
	group = parser.add_argument_group('required arguments for processing data')
	group.add_argument('-fd', '--file_directory', required=True,
	                    help='directory that holds files to be split')
	group.add_argument('-ext', '--extension', required=True,
	                    help='extension of files to be split')

	group = parser.add_argument_group('optional arguments for processing data')
	group.add_argument('-p','--prefix', default="set",
	                    help='prefix for new directories. default = set')
	group.add_argument('-od','--output_directory', default='file_directory',
	                    help='directory in which new directories will be made')
	group.add_argument('-a','--action', default='l',
	                    help='what to do with the images (m = move, l = symbolic link). Default is to link.')
	

	if not len(sys.argv) > 1:
	    parser.print_help()
	    sys.exit()

	#parse command line arguments
	args = parser.parse_args()
	if args.action != "m" and args.action != "l":
		print "Error: action must be either 'm' (move) or 'l' (link)!"
		sys.exit()

	# Gather the files in the provided image directory
	print "Finding files in directory {}...".format(args.file_directory)


	imageFiles = cpfiletools.find_files_in_directory(args.file_directory, [args.extension])
	if len(imageFiles) < 1:
		print "Error: no files found in directory {} matching extension {}. Exiting...".format(args.file_directory, args.extension)
		sys.exit()


	# Make a dictionary of all the image files keyed by tile number
	imageDirectory = os.path.abspath(args.file_directory)

	imageDict = cpfiletools.make_tile_dict_multiple(imageFiles, imageDirectory)
	tileList = imageDict.keys()

	numImagesPerTile = len(imageDict[tileList[0]])

	# now make new directories to hold split images:
	if args.output_directory == 'file_directory':
		outputPath = args.file_directory
	else:
		outputPath = args.output_directory
		if not os.path.exists(outputPath):
			print "Error: directory {} does not exist!".format(outputPath)

	newDirList = []
	for n in range(numImagesPerTile):
		dirname = outputPath+args.prefix+"{:02}".format(n+1)
		os.mkdir(dirname)
		newDirList.append(dirname)
		print "made directory: {}".format(dirname)


	# Now that directories are made, move images into those directories (or link)
	count = 0
	while count < numImagesPerTile:
		for tile in tileList:
			fullFileName = imageDict[tile].pop(0)
			prevPath, fileName = os.path.split(fullFileName)
			if args.action == "m":
				os.rename(fullFileName, newDirList[count]+"/"+fileName)
			if args.action == "l":
				os.symlink(fullFileName, newDirList[count]+"/"+fileName)
		count += 1

	print "Files split successfully"






def printDict(d):
	for key in d.keys():
		images = d[key]
		print "tile {}".format(key)
		for image in images:
			print "\t" + image 






if __name__ == '__main__':
    main()

