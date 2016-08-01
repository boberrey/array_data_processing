#!/usr/bin/env python
""" Split a folder of images into multiple directories 
for use in classic quantification pipeline


 Inputs:
   Directory with image files(.CPseq files)
   CPfluor files (.CPfluor files or directories)

 Outputs:
   CPseries files 

 Ben Ober-Reynolds, boberrey@stanford.edu """

import sys
import os
import argparse
import cpfiletools

def main():
	################ Parse input parameters ################

	#set up command line argument parser
	parser = argparse.ArgumentParser(description='Script for splitting a directory of images into multiple directories')
	group = parser.add_argument_group('required arguments for processing data')
	group.add_argument('-id', '--image_directory', required=True,
	                    help='directory that holds images to be split (.tif)')

	group = parser.add_argument_group('optional arguments for processing data')
	group.add_argument('-p','--prefix', default="set",
	                    help='prefix for new directories. default = set')
	group.add_argument('-od','--output_directory', default='image_directory',
	                    help='directory in which new directories will be made')
	group.add_argument('-a','--action', default='l',
	                    help='what to do with the images (m = move, l = symbolic link). Default is to link.')
	

	if not len(sys.argv) > 1:
	    parser.print_help()
	    sys.exit()






if __name__ == '__main__':
    main()