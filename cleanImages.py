#!/usr/bin/env python

"""
Run the matlab cleanImageScript.m script

Note: Python 3

Inputs:
   directory of images to clean
   

Outputs:
   cleaned images 

Ben Ober-Reynolds
"""

import os
import sys
import argparse
import re
import uuid
import subprocess
import time
from collections import OrderedDict
from joblib import Parallel, delayed


##### Gloval vars #####
offset_scale_x = -3.7
offset_scale_y = 3.7


def main():  
    # set up command line argument parser
    parser = argparse.ArgumentParser(description='script for cleaning blobs from images')
    group = parser.add_argument_group('required arguments:')
    group.add_argument('-id', '--image_directory', required=True,
        help='directory containing images to clean')
    group.add_argument('-cd', '--cleaned_image_directory', required=True,
        help='directory that will contain cleaned images')
    group.add_argument('-gv','--global_vars_path', required=True, 
        help='path to the directory in which the "GlobalVars.m" parameter file \
        for the run can be found')
    group = parser.add_argument_group('optional arguments')
    group.add_argument('-ma', '--max_val', type=int, default=4000,
        help='fluorescence intensity value that triggers cleaning')
    group.add_argument('-mn', '--min_val', type=int, default=1500,
        help='fluorescence value that triggers end of cleaning')

    # print help if no arguments provided
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    # parse command line arguments
    args = parser.parse_args()

    # Pre-defined variables, constants, and settings
    image_extension = 'tif'

    # Check directories and output files
    image_dir = args.image_directory
    if not os.path.isdir(image_dir):
        print("Error: invalid image directory selection. Exiting...")
        sys.exit()

    cleaned_image_dir = args.cleaned_image_directory
    if not os.path.isdir(image_dir):
        print("Error: invalid cleaned image directory selection. Exiting...")
        sys.exit()


    # Run Clean image script

    logstring = cleanImages(image_dir, cleaned_image_dir, args.max_val, args.min_val, args.global_vars_path)


def cleanImages(image_dir, cleaned_image_dir, max_val, min_val, global_vars_path):
    """
    Run the matlab script 'CleanImageScript.m'
    """
    matlabFunctionCallString = "CleanImageScript('{0}','{1}',{2},{3});".format(image_dir, cleaned_image_dir, max_val, min_val)
    logstring = spawnMatlabJob(matlabFunctionCallString, global_vars_path)
    return logstring



def spawnMatlabJob(matlabFunctionCallString,globalVarsPath):
    """
    Adapted from CPlibs.py 

    """
    try:
        #construct the command-line matlab call 
        functionCallString =                      "try,"
        functionCallString = functionCallString +     "addpath('{0}');".format(globalVarsPath) #placeholder TEMP DEBUG CHANGE
        functionCallString = functionCallString +     matlabFunctionCallString + ';'
        functionCallString = functionCallString + "catch e,"
        functionCallString = functionCallString +     "disp(getReport(e,'extended'));"
        functionCallString = functionCallString + "end,"
        functionCallString = functionCallString + "quit;"
    
        logFilename = 'matlabProcess_' + str(uuid.uuid4()) + str(time.time()) + '.tempLog' #timestamped logfile filename
    
        cmdString ='matlab -nodesktop -nosplash -singleCompThread -r "{0}"'.format(functionCallString)
        cmdString = cmdString + ' 1>> {0}'.format(logFilename)
        cmdString = cmdString + ' 2>> {0}'.format(logFilename)
       
        print('issuing subprocess shell command: ' + cmdString)
       
        returnCode = subprocess.call(cmdString,shell=True) #execute the command in the shell
        returnCode2 = subprocess.call('stty sane',shell=True) #matlab messes up the terminal in a weird way--this fixes it 
    
        #read log file into a string
        try:
            with open(logFilename) as logFilehandle:
                logString = logFilehandle.read()
            # delete logfile
            try:
                os.unlink(logFilename)
            except OSError:
                pass
        except IOError:
            logString = 'Log file not generated for command "' + functionCallString + '".'
    
        # return log
        return logString
    except Exception as e:
        return 'Python exception generated in spawnMatlabJob: ' + e



if __name__ == '__main__':
    main()
