#!/usr/bin/env python
""" 
Bootstrap fluorescence values from an array experiment
Note that there are some limitations with the bootstrapped data you will get from 
this script. This script assumes that you have the same number of data points for
each tile. It also will combine variants/filters accross tiles without respect to
'true' times that those data come from. Given the wide variability in the data to 
begin with, this is likely not a huge problem. If you plan on using these bootstrapped 
values for any downstream data analysis in which time is a factor, be aware that you
are bootstrapping values from a ~2-minute window (assuming you imaged all tiles)

 This script requires you to already have quantified images with another pipeline
 and generated CPsignal files.

 Inputs:
   Quantified array data (.CPsignal files)

 Outputs:
   Bootstrapped fluorescence values (.CPdata)

 Ben Ober-Reynolds, boberrey@stanford.edu
 20160811
 """

import sys
import os
import argparse
import cpfiletools
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import time


########################################################
################# Set Global Variables #################
# Set the variables that define what the columns in the CPseries file correspond to.
fluorSeriesCol = 8 # Column containing comma delimited fluorescence series
clusterIdCol = 0 # Column containing the cluster IDs
filterCol = 1 # Column containing the filter annotations
concCol = 10 # Column containig the concentrations
variantIDCol = 2 # Column containing the variant IDs
numRoundingDigits = 10 # number of digits to round bootstrapped values to
tileNumCol = 9 # Column containing the tile number
### MAIN ###

def main():
	################ Parse input parameters ################

	#set up command line argument parser
	parser = argparse.ArgumentParser(description='Script for bootstrapping fluorescence values from CPsignal files')
	group = parser.add_argument_group('required arguments for processing data')
	group.add_argument('-sd', '--CPsignal_dir', required=True,
	                    help='directory that holds the CPsignal files you want to get data from')

	group = parser.add_argument_group('optional arguments for processing data')
	group.add_argument('-bt','--bootstrap_type', default="v",
	                    help='How to subset data for bootstrapping: f = by filter, v = by variant. Default = "v"')
	group.add_argument('-fs','--filter_set', default="all",
	                    help='which filters you want to bootstrap. Default = "all"')
	group.add_argument('-vs','--variant_set', default="all",
	                    help='which variants you want to bootstrap. Default = "all"')
	group.add_argument('-st','--statistic', default='median',
	                    help='statistic you want to bootstrap. Default = "median", Options: "median","mean"')
	group.add_argument('-nb','--num_bootstraps', default=1000,
	                    help='number of times to bootstrap. Default = 1000')
	group.add_argument('-mr','--min_replicates', default=10,
	                    help='minimum number of replicates a variant must have for bootstrapping. Default = 10')
	group.add_argument('-ci','--confidence_interval', default='95',
	                    help='percent confidence interval to provide on bootstrapped statistic. Default = 95')
	group.add_argument('-od','--output_dir', default="CPsignal_dir",
	                    help='save output files to here. Default is provided CPsignal directory')
	group.add_argument('-op','--output_prefix', default="bootstrap_fluorescence",
	                    help='output file prefix. Default = "bootstrap_fluorescence"')
	group = parser.add_argument_group('other settings')
	group.add_argument('-n','--num_cores', type=int, default=1,
	                    help='maximum number of cores to use. default=1')
	
	if not len(sys.argv) > 1:
	    parser.print_help()
	    sys.exit()


	#parse command line arguments and check for problems
	args = parser.parse_args()
	numCores = int(args.num_cores)

	signal_files = cpfiletools.find_files_in_directory(args.CPsignal_dir, ['.CPseries'])
	bootstrap_type = args.bootstrap_type
	if bootstrap_type != 'f' and bootstrap_type != 'v':
		print >> sys.stderr, "Error: bootstrap type invalid (must be either 'f' or 'v'). Exiting..."
		sys.exit()
	filter_set = str.split(args.filter_set, ',')
	variant_set = str.split(args.variant_set, ',')
	statistic = args.statistic
	if statistic != 'median' and statistic != 'mean':
		print >> sys.stderr, "Error: statistic choice invalid. Exiting..."
		sys.exit()
	num_bootstraps = int(args.num_bootstraps)
	min_replicates = int(args.min_replicates)
	if int(args.confidence_interval) < 100:
		confidence_interval = [(100-float(args.confidence_interval))/2, 
								100-((100-float(args.confidence_interval))/2)]
	else:
		print >> sys.stderr, "Error: confidence interval must be between 0 and 100. Exiting..."
		sys.exit()
	
	if args.output_dir == "CPsignal_dir":
		output_directory = args.CPsignal_dir
	if not os.path.isdir(output_directory):
		print >> sys.stderr, "Error: output directory is not a valid directory. Exiting..."
		sys.exit()


	# Read in the CPseries files:
	print "Reading in data and subsetting if necessary..."
	start = time.time()
	series = loadAndConcatAllTiles(signal_files, args.CPsignal_dir)


	# Subset data of interest:
	# (If you don't reset the index here, pandas gives you a hard time concatenating the two data 
	# frames in the next step)
	series = selectData(filter_set, variant_set, filterCol, variantIDCol, series).reset_index(drop=True)
	print "\nStructuring data for bootstrapping..."

	### Restructure data frame such that fluorescence values are in their own columns
	all_fluor_series = []
	indexes = range(len(series.iloc[0, fluorSeriesCol].split(',')))

	# Pull out the fluorescence series and put into a data frame
	for i in xrange(len(series)):
		fluorescence_series = np.array([float(j) for j in series.iloc[i, fluorSeriesCol].split(',')])
		# Take the time to label unquantified clusters now, since it allows for fast removal later
		if all(np.isnan(fluorescence_series)):
			fluorescence_series = np.append(fluorescence_series, 0)
		else:
			fluorescence_series = np.append(fluorescence_series, 1)
		all_fluor_series.append(fluorescence_series)
	fluor_data_df = pd.DataFrame(all_fluor_series)
	# Quantified clusters get a '1'
	fluor_data_df.columns = indexes + ['Quantified']

	# separate out the ID columns from the working series and give them names
	id_cols = series.iloc[:,[clusterIdCol, filterCol, variantIDCol]]
	id_cols.columns = ["clusterID", "filterID", "variantID"]
	# Create the new working series
	frames = [id_cols, fluor_data_df]
	series = pd.concat(frames, axis=1)
	print "Done: {0:.2f} seconds".format(time.time() - start)


	# Remove all clusters that have no associated values
	print "\nRemoving unquantified clusters..."
	start = time.time()
	count = len(series.index)
	series = series.loc[series["Quantified"] == 1]
	series.drop("Quantified", axis=1, inplace=True)
	count = count - len(series.index)
	print "Removed "+str(count)+" unquantified clusters: {0:.2f} seconds".format(time.time() - start)



	### Perform Bootstrapping ###
	print "\nPerforming bootstrapping..."
	start = time.time()

	if bootstrap_type == 'v':
		allVariants = set(series.iloc[:,variantIDCol])
		namesToBootstrap = list(allVariants)
		label = "variantID"
	if bootstrap_type == 'f':
		allFilters = set(series.iloc[:,filterCol])
		namesToBootstrap = list(allFilters)
		label = "filterID"

	print "bootstrapping {} unique variants...".format(len(namesToBootstrap))
	# bootstrapOneVariant(variantSeries, indexes, variantName, numBootstraps, minReplicates, statistic, confidence_interval):
	if numCores > 1:
		allBootstrappedValues = (Parallel(n_jobs=numCores, verbose = 10)(delayed(bootstrapOneVariant)(series.loc[series[label] == name,:], 
							indexes, name, num_bootstraps, min_replicates, statistic, 
							confidence_interval) for name in namesToBootstrap))
	else:
		allBootstrappedValues = [bootstrapOneVariant(series.loc[series[label] == name,:], 
								indexes, name, num_bootstraps, min_replicates, statistic, 
								confidence_interval) for name in namesToBootstrap]
	allBootstrappedValues = filter(None, allBootstrappedValues)
	print "Done: {0:.2f} seconds".format(time.time() - start)
	print "{} variants passed minimum replicate cutoff of {}".format(len(allBootstrappedValues), min_replicates)
	

	### Write to file ###
	with open(output_directory+args.output_prefix+".CPdata", 'w') as f:
		for variant in allBootstrappedValues:
			for line in variant:
				for i in line:
					f.write(str(i)+'\t')
				f.write('\n')
		




def loadAndConcatAllTiles(tilenames, dirPath):
	# From 'fitSingleClusters.py' (Winston)
	#
    # This function loads the tiles listed in the tilenames variable and concatenates them all into 
    # a single pandas DataFrame.
    # Inputs: tilenames (the locations of all the tiles to be loaded)
    # Outputs: series (a single DataFrame containing information from all of the CPseries files)
    series = pd.DataFrame.from_csv(dirPath+tilenames[0], sep='\t', header=None, index_col=False)

    if len(tilenames) > 1:
        for i in xrange(len(tilenames)-1):
            series = pd.concat([series, pd.DataFrame.from_csv(dirPath+tilenames[i+1], sep='\t', 
            												header=None, index_col=False)])

    return series


def selectData(filtersToFit, variantsToFit, filterCol, variantIDCol, series):
    # Edited from 'fitSingleClusters.py' (Winston)
    # This function prunes data based on what the user wants fit and collects data about the type 
    # of fit to be performed
    # Input:
    # filtersToFit--List of filters applied to the clusters that the user wants to fit
    # filtersToFit--List of variant IDs that the user wants to fit
    # filterCol--the column in CPseries containing the filterID
    # variantIDCol--the column in CPseries containing the variant ID
    # series--the data frame containing the CPseries data
    # Output:
    # series--the pruned data frame containing the CPseries data that the user wants to fit

    if filtersToFit != ['all']:
        # Subset the current data based on what filters the user wants to fit
        series = series.loc[series[filterCol].isin(filtersToFit)]
    
    if variantsToFit != ['all']:
        # Subset the current data based on what variants the user wants to fit
        series = series.loc[series[variantIDCol].isin(variantsToFit)]
    
    # Print a message to the user about how many clusters are being fit
    print 'Found ' + str(len(series)) + ' single clusters...'
    
    return series


def bootstrapOneVariant(variantSeries, indexes, variantName, numBootstraps, minReplicates,
						 statistic, confidence_interval):
    # Inputs: 
    # series- data to be used for bootstrapping: pandas data frame
    # indexes- column indexes that contain data to bootstrap: list
    # variantName- variant/filter name that you are bootstrapping: string
    # Output:
    # bootstrapCI--the bootstrapped confidence intervals: list of lists

    median_statistic = []
    low_ci = []
    high_ci = []
    sampleSize = len(variantSeries.index)
    if sampleSize >= minReplicates:
    	for index in indexes:
    		bootstrappedValues = bootstrap(variantSeries.loc[:, index], 
    										numBootstraps, statistic, 
    										confidence_interval)
    		median_statistic.append(bootstrappedValues[0])
    		low_ci.append(bootstrappedValues[1])
    		high_ci.append(bootstrappedValues[2])
    	low_ci_list = [variantName, sampleSize, "{} percentile".format(confidence_interval[0])] + low_ci
    	median_list = [variantName, sampleSize, "50 percentile"] + median_statistic
    	high_ci_list = [variantName, sampleSize, "{} percentile".format(confidence_interval[1])] + high_ci
    	return [low_ci_list, median_list, high_ci_list]


def bootstrap(parameter, numBootstraps, statistic, confidence_interval):
    # Inputs:
    # parameter--the kinetic parameter computed for all of the single cluster fits--numpy array
    # numBootstraps--number of times to bootstrap the samples
    # Outputs:
    # bootstrappedConfidenceIntervals--the 95% confidence intervals of the bootstrapped medians
    

    # Remove nan's from the parameter array 
    parameter = parameter[~np.isnan(parameter)]
    # Calculate the number of parameter there are to bootstrap
    numFits = len(parameter[parameter])
    if numFits < 1:
    	return [float('nan')]*3
    
    # Initialize data structure to store bootstrapped values
    bootstrappedValues = []
    if statistic == 'median':
    	# Boostrap numBootstraps times
    	for j in xrange(numBootstraps):
    		# During each instance of bootstrapping, sample the data with replacement
    		bootstrappedParams = np.random.choice(parameter, size = (len(parameter)), replace=True)
    		# Compute the medians of the re-sampled parameters and save
    		bootstrappedValues.append(np.median(bootstrappedParams))
    if statistic == 'mean':
    	# Boostrap numBootstraps times
    	for j in xrange(numBootstraps):
    		# During each instance of bootstrapping, sample the data with replacement
    		bootstrappedParams = np.random.choice(parameter, size = (len(parameter)), replace=True)
    		# Compute the medians of the re-sampled parameters and save
    		bootstrappedValues.append(np.mean(bootstrappedParams))
    
    # Define the confidence intervals from the bootstraped medians
    bootstrappedConfidenceIntervals = np.percentile(bootstrappedValues, [50, confidence_interval[0], 
    												confidence_interval[1]])
    
    # Return the bootstraped confidence intervals to the user
    return bootstrappedConfidenceIntervals




if __name__ == '__main__':
    main()