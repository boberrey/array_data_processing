#!/usr/bin/env python

"""
Create input files for new image analysis software

Inputs:
   directory containing all images from array experiments
   (optional) which tiles to use

Outputs:
   image analysis input files 

Ben Ober-Reynolds
"""

# This version only creates an output file that contains tile_ID, path, and time (in that order)

# Need to think of a better way to handle time stamps. Also, it sounds like we're going to be 
# thinking more about how tile ID's work. Perhaps I need to add an ID field?

# I've editided this version to output a 'time table' for hacky analysis stuff


import os
import sys
import argparse
import re
import cpfiletools


def main():
    ################ Parse input parameters ################

    #set up command line argument parser
    parser = argparse.ArgumentParser(description='script for generating input files for image stack quantification')
    group = parser.add_argument_group('required arguments for processing data')
    group.add_argument('-id', '--input_directory', required=True,
                        help='directory that holds the image files of an array experiment')
    group = parser.add_argument_group('optional arguments for processing data')
    group.add_argument('-tl', '--tile_list', default="",
                        help='which tiles to use when generating input files (default is all)')
    group.add_argument('-od','--output_directory', default="",
                        help='save output files to here. default = input directory')
    group.add_argument('-op','--output_prefix', default="",
                        help='optional output file prefix')
    group.add_argument('-bf','--baseline_flag', default="",
                        help='flag denoting image files that contain baseline measurements')
    group.add_argument('-ef','--experiment_flag', default="",
                        help='flag denoting image files that contain experimental measurements')

    


    
    # print help if no arguments provided
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()
    
    #parse command line arguments
    args = parser.parse_args()
    
    if args.output_directory == "":
        args.output_directory = args.input_directory
    
    # is there an os-insensitive way to add the slash at the end of this? Does it matter?
    args.absPath = os.path.abspath(args.input_directory)+"/"
    
    # add underscore for better formatting
    if args.output_prefix != "":
        args.output_prefix = args.output_prefix + "_"
    
    ################ Make input files ################
    
    # Gather all image files in input directory and extract all tiles
    allFiles = cpfiletools.find_files_in_directory(args.input_directory, ['tif'])
    allTiles = set()
    for filename in allFiles:
        allTiles.add(cpfiletools.get_tile_number_from_filename(filename))
    
    # decide which tiles you want to use for making inputFiles
    if args.tile_list == "":
        tilesToUse = set(allTiles)
    else:
        tilesToUse = set(parse_tile_input(args.tile_list)) & allTiles
    
    # Tile dictionary for storing file data later on
    tileDict = {}
    for tile in tilesToUse:
        tileDict[tile] = []
    
    # Make a list of files (filtered by tile) that will be used to create input files
    filteredFiles = []
    for filename in allFiles:
        if cpfiletools.get_tile_number_from_filename(filename) in tilesToUse \
        and (args.baseline_flag in filename or args.experiment_flag in filename):
            filteredFiles.append(filename)
            print "will use:\t{}".format(filename)
    
    # Make separate lists for differently flagged files
    baselineFiles = []
    expFiles = []
    
    for filename in filteredFiles:
        if args.baseline_flag != "" and args.baseline_flag in filename:
            baselineFiles.append(filename)
        if args.experiment_flag != "" and args.experiment_flag in filename:
            expFiles.append(filename)
    
    if len(baselineFiles) < 1 and len(expFiles) < 1:
	print "ERROR: no tiles selected!"
	sys.exit()

    # Add all baseline files to the tile dictionary
    if len(baselineFiles) > 0:
	add_data_to_tile_dict(tileDict, args, baselineFiles, args.baseline_flag, 0)
    
    # Add all experimental files to the tile dictionary
    if len(expFiles) > 0:
	minTimeStamp = cpfiletools.parse_timestamp_from_filename(expFiles[0])
    	# assumes that the first experimental image timestamp (over all tiles)
    	# is the pseudo-zero timestamp (THIS NEEDS A MORE ELEGANT SOLUTION)
    	for filename in expFiles:
            if cpfiletools.parse_timestamp_from_filename(filename) < minTimeStamp:
                minTimeStamp = cpfiletools.parse_timestamp_from_filename(filename)
    	add_data_to_tile_dict(tileDict, args, expFiles, args.experiment_flag, minTimeStamp)
    
    # Time table for use in hacky analysis (20160201)
    timeTable = {}

    # sort output, add sequence number to experimental file entries, and print all files
    for tile in sorted(tilesToUse):
        tileDict[tile].sort()
        count = 1
        timeTable[tile] = []    # Fill in time table for each tile
        for filedata in tileDict[tile]:
            timeTable[tile].append(filedata.timestamp)
            if args.experiment_flag in filedata.ID:
                filedata.ID = filedata.ID+"_"+str(count)
                count += 1
        filename = args.output_prefix + "tile" + tile + ".ipf"
        with open(args.output_directory+filename, 'w') as f:
            header = "{}".format("time")
            f.write(header + "\n")
            for filedata in tileDict[tile]:
                f.write("{}\n".format(filedata))
            f.write("\n")
        print "Successfully made file: {}".format(filename)

    # Print out the time Table (20160201)
    with open("timeTable.txt", 'w') as f:
        tiles = sorted(timeTable.keys())
        for tile in tiles:
            f.write(tile)
            for time in timeTable[tile]:
                f.write("\t"+str(time))
            f.write("\n")
    print "successfully made file: timeTable.txt"

    
def parse_tile_input(tileList):
    """
    Parse the input of a ',' and '-' delimited tile list into an actual list of tiles
    (with leading zeros as needed)
    Input: tile list with runs of tiles split by '-' and groups/single tiles split by ','
    Output: list of correctly formatted tiles
    """
    commaSplit = tileList.split(',')
    newTileList = []
    formattedTileList = []
    for group in commaSplit:
        seq = group.split('-')
        if len(seq) > 1:
            tilesInSeq = range(int(seq[0]), int(seq[1]) + 1)
            newTileList += tilesInSeq
        else:
            newTileList.append(int(seq[0]))
    for tileNumber in newTileList:
        fomattedTileNum = '{:03}'.format(tileNumber)
        formattedTileList.append(fomattedTileNum)
    return formattedTileList


def add_data_to_tile_dict(tileDict, args, fileList, flag, minTimeStamp):
    """
    Add FileData class instances to the tile dictionary.
    Inputs: Tile Dictionary, args, file list, flag for file list, minTimeStamp
    Output: Tile Dictionary with newly added FileData
    """
    for filename in fileList:
        tile = cpfiletools.get_tile_number_from_filename(filename)
        ID = flag+"_tile"+tile
        path = args.absPath + filename
        if minTimeStamp == 0:
            timestamp = 0
        else:
            timestamp = cpfiletools.parse_timestamp_from_filename(filename)
            # Add one second to all timestamps to prevent having two zero timepoints on tile 1
            timestamp = cpfiletools.get_time_delta(timestamp, minTimeStamp) + 1
        tileDict[tile].append(FileData(ID, path, timestamp))


def print_list(list):
    for i in list:
        print i

class FileData:
    """
    class for storing information about a given file.
    """
    def __init__(self, ID, path, timestamp):
        self.ID = ID
        self.path = path
        self.timestamp = timestamp

    
    def __cmp__(self, other):
        # FileData comparisons are based on the timestamp
        if self.timestamp < other.timestamp:
            return -1
        if self.timestamp > other.timestamp:
            return 1
        if self.timestamp == other.timestamp:
            return 0
    
    
    def __str__(self):
        # Printing FileData is every field in tab-delimitted format
        return "{}\t{}\t{}".format(self.ID, self.path, self.timestamp).strip()
    
    def __hash__(self):
        # Need to define a hash method for custom classes if you want to
        # use them in hash structures (sets/dictionaries)
        return hash(ID+str(timestamp))
    



if __name__ == '__main__':
    main()
