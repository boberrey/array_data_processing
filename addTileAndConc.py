#!/usr/bin/env python
""" 
Add the tile number and concentration to a CPseries file 
as columns at the end of the file.

 This script requires you to already have quantified images with another pipeline.

 Inputs:
   CPseries file

 Outputs:
   CPseries file with tile number and concentration

 Ben Ober-Reynolds, boberrey@stanford.edu
 20160808
 """

import sys
import os
import argparse
import subprocess
import re


### MAIN ###


def main():
	################ Parse input parameters ################

	#set up command line argument parser
	parser = argparse.ArgumentParser(description='adding tile number and concentration to a CPseries file')
	group = parser.add_argument_group('required arguments for processing data')
	group.add_argument('-f', '--CPseries_file', required=True,
	                    help='a CPseries file')

	group = parser.add_argument_group('optional arguments')
	group.add_argument('-c','--conc', default="",
	                    help='concentration to add to CPseries file')
	group.add_argument('-nf','--new_file_flag', default="tc",
	                    help='flag to add to new file')


	if not len(sys.argv) > 1:
	    parser.print_help()
	    sys.exit()


	#parse command line arguments
	args = parser.parse_args()

	filename = os.path.abspath(args.CPseries_file)
	if not os.path.isfile(filename):
		print "Error: file "+filename+" is not a valid file. Exiting..."
		sys.exit()

	tilenum = get_tile_number_from_filename(filename)
	conc = args.conc

	basename = os.path.basename(filename)
	path = filename.strip(basename)
	newfilename = path+args.new_file_flag+"_"+basename

	subprocess.call("awk \'{{print $0 \"\t{}\t{}\"}}\'".format(tilenum, conc)+"<{} > {}".format(filename, newfilename), shell=True)

	print "Made new file: "+newfilename





def get_tile_number_from_filename(inFilename):
    """
    (from cpfiletools)
    Extract the tile number from a provided filename based on the presence of
    'tile###'
    Input: filename (string)
    Output: three digit tile number (string)
    """
    # from CPlibs
    (path,filename) = os.path.split(inFilename) #split the file into parts
    (root,ext) = os.path.splitext(filename)
    matches = re.findall('tile[0-9]{1,3}',root.lower())
    tileNumber = ''
    if matches != []:
        tileNumber = '{:03}'.format(int(matches[-1][4:]))
    return tileNumber


if __name__ == '__main__':
    main()