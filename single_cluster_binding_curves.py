#!/usr/bin/env python
""" 
Fit binding curves across a CPseries file

 Ben Ober-Reynolds, boberrey@stanford.edu
 20190620
"""

import pandas as pd
import numpy as np
import argparse
import sys
import os
import random
from lmfit import Parameters, minimize, Model, report_fit
from joblib import Parallel, delayed
import time


### MAIN ###


def main():
    ################ Parse input parameters ################

    #set up command line argument parser
    parser = argparse.ArgumentParser(description='Script for fitting equilibrium binding curves for a CPseries file')
    group = parser.add_argument_group('required arguments:')
    group.add_argument('-cs', '--cpseries', required=True,
                        help='CPseries.pkl file')
    group.add_argument('-c', '--concentrations', required=True,
                        help='flat file of concentrations in matched order to fluorescence data in CPseries')
    group.add_argument('-fa', '--fid_cpannot', required=True,
                        help='CPannot file containing only fiducial clusters')
    group.add_argument('-tf', '--tight_fmaxes', required=True,
                        help='df of previously identified tight binder fits')

    group = parser.add_argument_group('optional arguments for processing data')
    group.add_argument('-n','--num_cores', type=int, default=10,
                        help='number of cores to use')


    if not len(sys.argv) > 1:
        parser.print_help()
        sys.exit()


    #parse command line arguments
    args = parser.parse_args()

    cpseries_file = args.cpseries
    fid_cpannot_file = args.fid_cpannot
    concentration_file = args.concentrations
    tight_fmax_file = args.tight_fmaxes

    num_cores = args.num_cores

    ##############################
    # Read in and preprocess data
    ##############################

    print "Reading in data..."
    binding_df = pd.read_pickle(cpseries_file)
    # drop any rows with all nans (This drops more than half of the clusters? Do they just have low/no signal?)
    nclust = len(binding_df.index)
    print "{} clusters in CPseries file.".format(nclust)
    binding_df.dropna(axis=0, thresh=4, inplace=True) # Clusters must have at least 4 points
    print "{} clusters that have usable data ({}%)".format(
        len(binding_df.index), round((len(binding_df.index)/float(nclust))*100, 3))

    concentrations = np.recfromtxt(concentration_file)

    annot_df = pd.read_pickle(fid_cpannot_file)

    # Get fmax values from library fits
    tight_fmaxes = pd.read_table(tight_fmax_file)['fmax'].values.tolist()

    # Merge annot and cpseries (Fiducial marks only)
    fid_df = annot_df.merge(binding_df, left_index=True, right_index=True, how='inner')

    # Normalize everything by fiducial signal
    print "Normalizing to fiducial... ({} fiducial clusters)".format(len(fid_df.index))
    fiducial_meds = fid_df.groupby('variant_ID').get_group('11111111').iloc[:,1:].median().values
    binding_df = binding_df / fiducial_meds

    # Force fmin to simply be the median signal across all clusters
    fmin_val = float(binding_df.iloc[:,0].median())

    ###########################################
    # Perform single cluster fits
    ###########################################

    single_cluster_filename = 'single_cluster_fits.txt'
    
    # Split data into n groups for fitting
    df_chunks = split(binding_df, chunk_size=int(np.ceil(binding_df.shape[0]/float(num_cores))))

    # Perform single cluster fits in parallel
    print "Fitting single clusters on {} cores...".format(num_cores)
    start = time.time()
    if num_cores > 1:
        fit_df_list = (Parallel(n_jobs=num_cores, verbose=10)\
            (delayed(single_cluster_fits)(
                chunk, concentrations, tight_fmaxes, fmin_val) for chunk in df_chunks))
    else:
        fit_df_list = [single_cluster_fits(
            chunk, concentrations, tight_fmaxes, fmin_val) for chunk in df_chunks]

    print "Single cluster fitting finished, {} minutes.".format(round((time.time() - start)/60.0, 3))
    sc_fit_df = pd.concat(fit_df_list)
    #sc_fit_df.sort_values('variant_ID', inplace=True)
    sc_fit_df.to_csv('single_cluster_fits.txt', sep='\t')


# Splitting df:
# (Adapted from: http://yaoyao.codes/pandas/2018/01/23/pandas-split-a-dataframe-into-chunks)
def index_marks(nrows, chunk_size):
    return range(chunk_size, int(np.ceil(nrows / chunk_size)) * chunk_size, chunk_size)

def split(df, chunk_size):
    indices = index_marks(df.shape[0], chunk_size)
    return np.split(df, indices)


# Fitting functions:

def hill_equation(x, fmin, fmax, Kd, n=1):
    """
    Hill-langmuir equation for equilibrium binding
    x = concentrations (in nM)
    fmin = minimum signal (at 0 concentration)
    fmax = maximum concentration (at infinite concentration)
    return = signal value
    """
    return fmin + fmax*(x**n /(Kd + x**n))


def hill_equation_params(fmax=None, fmin=None, Kd=None, n=None):
    # Define parameters object
    params = Parameters()
    default_params = {
        "fmax":{"value": 1.0, "vary": True, "min": 0.0, "max": np.inf},
        "fmin":{"value": 0.0, "vary": True, "min":0.0, "max":np.inf},
        "Kd":{"value": 1.0, "vary": True, "min": 0.0, "max": np.inf},
        "n":{"value": 1.0, "vary": False, "min": -np.inf, "max": np.inf}

    }
    if fmax:
        for opt, val in fmax.items():
            default_params["fmax"][opt] = val
    if fmin:
        for opt, val in fmin.items():
            default_params["fmin"][opt] = val
    if Kd:
        for opt, val in Kd.items():
            default_params["Kd"][opt] = val
    if n:
        for opt, val in n.items():
            default_params["n"][opt] = val

    for p, dct in default_params.items():
        params.add(p, **dct)


    return params


def single_cluster_fits(df, x, tight_fmaxes, fmin_val=0.0):
    """
    Fit single clusters for cluster in df
    Input:
        df = dataframe of clusters
        x = concentrations
    """
    clusterIDs = df.index.values.tolist()
    data = df.values
    results_dict = {}
    fit_model = Model(hill_equation)
    tight_fmax_2p5, med_fmax, tight_fmax_97p5 = np.nanpercentile(tight_fmaxes, q=[2.5, 50.0, 97.5])
    

    for i, clust in enumerate(clusterIDs):
        
        fluorescence = data[i,:]
        #non_binder = 0.0
        params = hill_equation_params(
                # Force fmax to be at least 2.5th percentile of previous fmax distribution
                fmax={"value":med_fmax, "vary":True, "min": tight_fmax_2p5},
                Kd={"value":max(x)},
                fmin={"value":fmin_val, "vary":False})
        """
        non_binder = 0.0
        if all(fluorescence[~np.isnan(fluorescence)] < tight_fmax_2p5*0.2):
            non_binder = 1.0
            params = hill_equation_params(
                    fmax={"value": med_fmax, "vary": False}, 
                    Kd={"value":max(x)},
                    fmin={"value":fmin_val, "vary":False})
        else:
            params = hill_equation_params(
                fmax={"value":med_fmax, "vary":True},
                Kd={"value":max(x)},
                fmin={"value":fmin_val, "vary":False})
        """

        try:
            fit = fit_model.fit(fluorescence, params, x=x)
        except Exception as e:
            print "Error while fitting {}".format(vID)
            print str(e)
            continue
        ## Things we want to report
        # quality of fit:
        ss_error = np.sum((fit.residual)**2)
        ss_total = np.sum((fluorescence - np.nanmean(fluorescence))**2)
        rsq = 1 - ss_error/ss_total
        rmse = np.sqrt(np.nanmean((fit.residual)**2))

        results_dict[clust] = {
            'Kd': fit.params['Kd'].value, 
            'fmax': fit.params['fmax'].value, 
            'fmin': fit.params['fmin'].value, 
            'rsq': rsq, 
            'rmse': rmse,
            'ier': fit.ier,
            'Kd_stderr': fit.params['Kd'].stderr,
            'fmax_stderr': fit.params['fmax'].stderr,
            'fmin_stderr': fit.params['fmin'].stderr
            #'non_binder': non_binder
        }
            
        
    return pd.DataFrame(results_dict).T



if __name__ == '__main__':
    main()
