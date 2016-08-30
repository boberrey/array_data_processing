
"""
Tools for working with 'CP__' files

Some functions taken from Sarah's IMlibs.py
Others taken from Curtis' CPlibs and/or CPscripts
Some made by me

Ben Ober-Reynolds, boberrey@stanford.edu


Last edited 20160802
"""

import os
import sys
import time
import re
import numpy as np
import pandas as pd
import datetime
import uuid
import subprocess


def get_tile_number_from_filename(inFilename):
    """
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


def find_files_in_directory(dirPath, extensionList=None, 
                            excludedExtensionList=None):
    """
    Locate files in a given directory path. Optionally, desired files are 
    identified as matching one of the extension types provided in 
    'extensionList'
    Input: directory path, list of approved extensions, (list of excluded extensions)
    Output: List of found files 
    """
    def extension_match(filename, extensionList=None):
        # from CPlibs
        if extensionList is not None:
            for currExt in extensionList:
                if filename.lower().endswith(currExt.lower()):
                    return True
        return False

    dirList = os.listdir(dirPath)
    fileList = []
    for currFilename in dirList:
        if (extension_match(currFilename, extensionList) 
            and not extension_match(currFilename, excludedExtensionList)):
            fileList.append(currFilename)
    if len(dirList) == 0:
        print '\tNONE FOUND'
    else:
        #for filename in fileList:
        #    print "found:\t\t{}".format(filename)
        return fileList

def make_tile_dict(fileList, directory):
    """
    Make a dictionary of files keyed by tile number.  
    Input: list of files containing tile numbers
    Output: dictionary of file names keyed by tile number
    """
    fileDict = {}
    for f in fileList:
        tile = get_tile_number_from_filename(f)
        if tile == '':
            print "Error: no tile number in file: "+ f
            sys.exit()
        else:
            if tile in fileDict:
                print "Error: multiple files per tile"
                sys.exit()
            fileDict[tile] = os.path.join(directory, f)
    return fileDict

def make_tile_dict_multiple(fileList, directory):
    """
    Make a dictionary of files keyed by tile number, where each tile has multiple files.  
    Input: list of files containing tile numbers
    Output: dictionary of file names keyed by tile number
    """
    fileDict = {}
    for f in fileList:
        tile = get_tile_number_from_filename(f)
        if tile == '':
            print "Error: no tile number in file: "+ f
            sys.exit()
        else:
            if tile in fileDict:
                fileDict[tile].append(os.path.join(directory, f))
            else:
                fileDict[tile] = [os.path.join(directory, f)]
    for key in fileDict.keys():
        # automatically sort binding series by timestamp:
        fileDict[key].sort(key = lambda x: x.split('_')[-1])
    return fileDict


def generate_CPseries_files(cpSeqFilename, allRNA, bindingSeries, 
                            CPseriesFilename, tile):
    """
    Generate a CPseries file from the appropriate CPseq data, binding series, 
    and optionally all RNA data.
    Inputs: CP seq filename, list of binding series CPfluors, 
    CPseries filename, tile, all RNA CPfluor filename (optional)
    Output: writes the CPseries file
    """
    # Get the number of lines in the CPseq file. Currently assumes all files
    # will have the same number of lines
    numLines = int(subprocess.check_output(
        ("wc -l {} | ".format(cpSeqFilename)
        +" awk \'{print $1}\'"), shell=True).strip())
    # First get the all RNA values, if present
    if os.path.isfile(allRNA):
        allRNAsignal = get_signal_from_CPFluor(allRNA)
    else:
        allRNAsignal = np.ones(numLines)*np.nan

    # Prepare to calculate all signal for binding series:
    bindingSeriesSignal = np.zeros((numLines, len(bindingSeries)))

    for i, fluor in enumerate(bindingSeries):
        # Check to see if the fluor file is a 'real' file by looking
        # at the filesize. Fluor files will never be less that 1000 bytes.
        if os.path.isfile(fluor) and os.path.getsize(fluor) > 1000:
            bindingSeriesSignal[:,i] = get_signal_from_CPFluor(fluor)
        else:
            bindingSeriesSignal[:,i] = np.ones(numLines)*np.nan
    # Combine binding series signals:
    bs_comma_format = np.array(
        [','.join(bindingSeriesSignal[i].astype(str)) for i in range(numLines)])
    
    # read in CPseq file as a pandas dataframe
    cp_seq = pd.read_table(cpSeqFilename, header=None)
    # Here is where we can decide which components of the CPseq file we want
    # to carry over into the CPseries file:
    """
    CPseq columns:
    0: ClusterID
    1: filter
    2: read1_seq
    3: read1_quality
    4: read2_seq
    5: read2_quality
    6: index1_seq
    7: index1_quality
    8: index2_seq
    9: index2_quality
    """
    CPseriesframe = cp_seq.iloc[:,[0,1,2,4,6,8]]

    # create new dataframe with signal data, 
    # then concatenate the two dataframes together
    signal_df = pd.DataFrame({'allRNA': allRNAsignal, 
        'bindingSeriesSignal': bs_comma_format}, dtype=str)
    final_df = pd.concat([CPseriesframe, signal_df], axis=1)
    # Save data
    np.savetxt(CPseriesFilename, final_df.values, fmt='%s', delimiter='\t')
    print "Successfully made file : {}".format(CPseriesFilename)



def get_signal_from_CPFluor(CPfluorfilename):
    """
    From Sarah's IMlibs.py: calculates the signal values from a given 
    CPfluor file.
    Volume under a 2D gaussian function given by 
    2*pi*Amplitude*sigma_x*sigma_y
    Input: CPfluor filename
    Output: np.array of calculated signals (i.e. volume under 2D gaussian) 
    """
    fitResults = pd.read_csv(
        CPfluorfilename, usecols=range(7, 12), 
        sep=':', header=None, 
        names=['success', 'amplitude', 'sigma', 'fit_X', 'fit_Y'] )
    signal = np.array(2*np.pi*fitResults['amplitude']*fitResults['sigma']*fitResults['sigma'])
    signal[np.array(fitResults['success']==0)] = np.nan
    return signal



def parse_timestamp_from_filename(filename):
    """
    Extract the time stamp from a provided filename, assuming the timestamp 
    is at the end of the file and separated from the rest of the filename 
    by a '_'
    Input: filename (string)
    Output: timestamp object
    """
    root, ext = os.path.splitext(filename)
    try: 
        timestamp=filename.strip(ext).split('_')[-1]
        date, time = timestamp.split('-')
        year, month, day = np.array(date.split('.'), dtype=int)
        hour, minute, second, ms = np.array(time.split('.'), dtype=int)
        timestamp_object = datetime.datetime(
            year=year, month=month, day=day, hour=hour, minute=minute, 
            second=second, microsecond=ms*1000)
    except ValueError:
        print "ERROR: no timestamp on file: {}".format(filename)
        sys.exit()
    return timestamp_object

def get_time_delta(timestamp_final, timestamp_initial):
    """
    Get the time delta in seconds
    """
    return (timestamp_final - timestamp_initial).seconds + (timestamp_final
            - timestamp_initial).microseconds/1E6


def spawn_matlab_job(matlabFunctionCallString,tempPaths):
    """
    From CPlibs.py
    Spawn a matlab job by correctly formatting command-line call.
    Inputs: formatted string of matlab function to call, temporary paths to use
    Output: call the matlab function

    *note: this version differs from the original spawnMatlabJob 
    by requiring three new paths

    *Warning: this function probably needs some work*
    """
    try:
        #construct the command-line matlab call 
        functionCallString =                      "try,"
        functionCallString = functionCallString +     "addpath('{0}','{1}','{2}');".format(tempPaths[0], tempPaths[1], tempPaths[2]) #placeholder TEMP DEBUG CHANGE
        functionCallString = functionCallString +     matlabFunctionCallString + ';'
        functionCallString = functionCallString + "catch e,"
        functionCallString = functionCallString +     "disp(getReport(e,'extended'));"
        functionCallString = functionCallString + "end,"
        functionCallString = functionCallString + "quit;"
        #print "function Call string:" + '\n' + functionCallString
    
        logFilename = 'matlabProcess_' + str(uuid.uuid4()) + str(time.time()) + '.tempLog' #timestamped logfile filename
        
        print logFilename

        cmdString ='matlab -nodesktop -nosplash -singleCompThread -r "{0}"'.format(functionCallString)
        cmdString = cmdString + ' 1>> {0}'.format(logFilename)
        cmdString = cmdString + ' 2>> {0}'.format(logFilename)
       
        print 'issuing subprocess shell command: ' + cmdString
       
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
    except Exception, e:
        return 'Python exception generated in spawn_matlab_job: ' + e.message

def get_registration_offset(CPseqFilename, imageFilename, dataScaling, 
                            filterSubsets, tempPaths):
    try:
        matlabFunctionCallString = "GenerateRegistrationOffsetMap('{0}','{1}','{2}',{3},'','','');".format(CPseqFilename, imageFilename, dataScaling, filterSubsets)
        #print matlabFunctionCallString
        logString = spawn_matlab_job(matlabFunctionCallString,tempPaths)
        return (CPseqFilename,logString)
    except Exception,e:
        return(CPseqFilename,'Python excpetion generated in getRegistrationOffset: ' + e.message)

def analyse_image(CPseqFilename, imageFilename, dataScaling, filterSubsets, 
                registrationOffsetMapFilename, tempPaths):
    try:
        matlabFunctionCallString = "AnalyseImage('{0}','{1}','{2}', {3}, '{4}','','');".format(CPseqFilename, imageFilename, dataScaling, filterSubsets, registrationOffsetMapFilename)
        #print matlabFunctionCallString
        logString = spawn_matlab_job(matlabFunctionCallString,tempPaths)
        return (CPseqFilename,logString)
    except Exception,e:
        return(CPseqFilename,'Python excpetion generated in quantifyFluorescence: ' + e.message)


def analyse_series(seqDataFilename, imageListFilename, darkImageIndex, 
                registrationImageIndex, dataScaling, filterSubsets, 
                workingPath, tempPaths):
    try:
        matlabFunctionCallString = "AnalyseSeries('{0}','{1}',{2},{3},'{4}',{5},'{6}');".format(
            seqDataFilename, imageListFilename, darkImageIndex, 
            registrationImageIndex, dataScaling, filterSubsets, workingPath)
        #print matlabFunctionCallString
        logString = spawn_matlab_job(matlabFunctionCallString, tempPaths)
        return (imageListFilename,logString)
    except Exception,e:
        return(imageListFilename,'Python excpetion generated in analyseSeries: ' + e.message)

def update_progress(current, total):
    """
    Curtis wrote this cool little function to display a progress bar for 
    long processes
    """
    updateInterval = 50
    if(current % updateInterval == 0) | (current+1 == total):
        barLength = 50 # Modify this to change the length of the progress bar
        progress = float(current+1)/float(total)
        status = ""
        if current+1 == total:
            status = "Done...\r\n"
        block = int(round(barLength*progress))
        text = "\rPercent: [{0}] {1:3.2f}% ({2} sequences) {3}".format(
         "#"*block + "-"*(barLength-block), progress*100, current+1, status)
        sys.stdout.write(text)
        sys.stdout.flush()


def printList(lst):
    for l in lst:
        print "\t{}".format(l)
