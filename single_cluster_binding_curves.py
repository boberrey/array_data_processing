#!/usr/bin/env python
""" 
Fit binding curves across a CPseries file

 Ben Ober-Reynolds, boberrey@stanford.edu
 20190620
"""

import matplotlib as mpl
mpl.use('Agg') # Don't display plots
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy as sp
import argparse
import sys
import os
import random
from lmfit import Parameters, minimize, Model, report_fit
from joblib import Parallel, delayed
from itertools import compress
import time



### MAIN ###


def main():
    ################ Parse input parameters ################

    #set up command line argument parser
    parser = argparse.ArgumentParser(description='Script for fitting equilibrium binding curves for a CPseries file')
    group = parser.add_argument_group('required arguments:')
    group.add_argument('-cs', '--cpseries', required=True,
                        help='CPseries.pkl file')
    group.add_argument('-ca', '--cpannot', required=True,
                        help='CPannot.pkl file')
    group.add_argument('-c', '--concentrations', required=True,
                        help='flat file of concentrations in matched order to fluorescence data in CPseries')
    group.add_argument('-vt', '--variant_table', required=True,
                        help='experiment variant table')

    group = parser.add_argument_group('optional arguments for processing data')
    group.add_argument('-nt', '--norm_table',
                        help='normalization variant table')
    group.add_argument('-n','--num_cores', type=int, default=10,
                        help='number of cores to use')
    group.add_argument('-rf','--rsq_filter', type=float, default=0.7,
                        help='R squared filter for removing clusters based on single cluster fits')
    group.add_argument('-mg','--multiguide_group', default=None,
                        help='Which multiguide group to restrict to.')
    group.add_argument('-pl','--plots', action="store_true",
                        help='Use flag to generate plots of fits.')
    group.add_argument('-pd','--plot_dir', default="curve_plots/",
                        help='Use flag to generate plots of fits.')


    if not len(sys.argv) > 1:
        parser.print_help()
        sys.exit()


    #parse command line arguments
    args = parser.parse_args()

    cpseries_file = args.cpseries
    cpannot_file = args.cpannot
    concentration_file = args.concentrations
    variant_file = args.variant_table
    norm_file = args.norm_table
    num_cores = args.num_cores
    rsq_filter = float(args.rsq_filter)

    # Make the plot directory if plotting
    if args.plots:
        if not os.path.isdir(args.plot_dir):
            print "Plot directory {} does not exist. Creating...".format(args.plot_dir)
            os.mkdir(args.plot_dir)
        plot_dir = args.plot_dir
    else:
        plot_dir = None

    ##############################
    # Read in and preprocess data
    ##############################

    print "Reading in data..."
    binding_df = pd.read_pickle(cpseries_file)
    # drop any rows with all nans
    binding_df.dropna(axis=0, how='all', inplace=True)

    concentrations = np.recfromtxt(concentration_file)

    annot_df = pd.read_pickle(cpannot_file)
    variant_df = pd.read_table(variant_file, header=None)

    # Construct dict of labels for plots
    label_dict = {}
    if len(variant_df.columns) == 6:
        variant_df.columns = ['Sequence', 'variant_ID', 'group_name', 'mut_annotation', 'wt_seq', 'guide_number']
        for idx, row in variant_df.iterrows():
            ID = str(row.variant_ID)
            group = str(row.group_name)
            annot = str(row.mut_annotation)
            guide = "p".join(str(row.guide_number).split('.'))
            label_dict[ID] = "_".join([guide, group, annot])

    else:
        variant_df.columns = ['Sequence', 'variant_ID', 'group_name', 'mut_annotation']
        for idx, row in variant_df.iterrows():
            ID = str(row.variant_ID)
            group = str(row.group_name)
            annot = str(row.mut_annotation)
            label_dict[ID] = "_".join([group, annot])

    if args.multiguide_group:
        print "Only fitting variants in multiguide group {}".format(args.multiguide_group)
        group_ids = list(variant_df[variant_df.guide_number.apply(guide_group) == args.multiguide_group].variant_ID.values)
        group_ids.append('11111111')
        groups = annot_df.groupby('variant_ID').groups.keys()
        # Need to do this elaborate filtering procedure to keep variants that have multiple variant IDs
        valid = [contains_variant_ID(g, group_ids) for g in groups]
        valid_groups = list(compress(groups, valid))
        annot_df = annot_df[annot_df.variant_ID.isin(valid_groups)].copy()


    # Merge annot and cpseries
    merged_df = annot_df.merge(binding_df, left_index=True, right_index=True, how='inner')

    # Normalize everything by fiducial signal
    #print "Normalizing to fiducial..."
    #fiducial_meds = merged_df.groupby('variant_ID').get_group('11111111').iloc[:,1:].median().values
    #merged_df.iloc[:,1:] = merged_df.iloc[:,1:] / fiducial_meds

    # Normalize everything by normalization variants if provided
    if norm_file:
        print "Normalizing to variants in norm table..."
        norm_df = pd.read_table(norm_file, header=None)
        if len(variant_df.columns) == 6:
            norm_df.columns = ['Sequence', 'variant_ID', 'group_name', 'mut_annotation', 'wt_seq', 'guide_number']
        else:
            norm_df.columns = ['Sequence', 'variant_ID', 'group_name', 'mut_annotation']
        norm_var_IDs = [str(x) for x in norm_df.variant_ID.values]
        norm_series = merged_df[merged_df.variant_ID.isin(norm_var_IDs)]
        norm_meds = norm_series.iloc[:,1:].median().values

        merged_df.iloc[:,1:] = merged_df.iloc[:,1:] / norm_meds

    # Force fmin to simply be the median signal across all variants
    fmin_val = float(merged_df.iloc[:,1].median())

    ###########################################
    # Perform single cluster fits
    ###########################################

    # Check to see if single clusters have already been fit:
    single_cluster_filename = 'single_cluster_fits.txt'
    if os.path.exists(single_cluster_filename):
        print "Single cluster file already exists. Loading..."
        sc_fit_df = pd.read_table(single_cluster_filename, index_col=0)
    else:
        print "Splitting data into {} chunks and fitting...".format(num_cores)
        # Now perform the actual fits
        variant_IDs = list(set(merged_df.variant_ID))

        # shuffle list and then make (num_cores) chunks
        random.shuffle(variant_IDs)
        chunk_list = list(make_chunks(variant_IDs, num_cores))

        # Create a new list of dataframes containing chunked variant_IDs
        grouped_list = [merged_df[merged_df.variant_ID.isin(chunk)].copy().groupby('variant_ID') for chunk in chunk_list]

        # Perform single cluster fits in parallel
        start = time.time()
        if num_cores > 1:
            fit_df_list = (Parallel(n_jobs=num_cores, verbose=10)\
                (delayed(single_cluster_fits)(
                    sub_grouped, concentrations, fmin_val) for sub_grouped in grouped_list))
        else:
            fit_df_list = [single_cluster_fits(
                sub_grouped, concentrations, fmin_val) for sub_grouped in grouped_list]

        print "Single cluster fitting finished, {} minutes.".format(round((time.time() - start)/60.0, 3))
        sc_fit_df = pd.concat(fit_df_list)
        sc_fit_df.sort_values('variant_ID', inplace=True)
        sc_fit_df.to_csv('single_cluster_fits.txt', sep='\t')


    #####################################################
    # Filter clusters based on their single cluster fit
    #####################################################
    print "Filtering clusters with fit rsq of less than {}".format(rsq_filter)
    good_sc_df = sc_fit_df[sc_fit_df.rsq > rsq_filter]
    print "{} of {} ({}%) clusters passed rsq cutoff".format(len(good_sc_df.index), len(sc_fit_df.index), (len(good_sc_df.index)/float(len(sc_fit_df.index)))*100)

    #####################################################
    # Get empirical fmax distribution of 'tight binders'
    #####################################################

    # Ensure certain columns are floats:
    good_sc_df['fmax'] = good_sc_df['fmax'].astype(float)

    # Define 'tight binders' as those that are 95% saturated at second to last point
    tight_Kd_cutoff = concentrations[-2]/0.95 - concentrations[-2]
    print "Identifying 'tight binders' with Kd < {}".format(tight_Kd_cutoff)
    tight_binders = good_sc_df[(good_sc_df.Kd < tight_Kd_cutoff)].copy()
    tight_binders.to_csv('tight_binder_fits.txt', sep='\t')
    tight_fmaxes = tight_binders.groupby('variant_ID')['fmax'].median()
    tight_fmax_2p5, tight_fmax_97p5 = np.nanpercentile(tight_fmaxes, q=[2.5,97.5])
    conf_fmaxes = tight_fmaxes[(tight_fmaxes > tight_fmax_2p5) & (tight_fmaxes < tight_fmax_97p5)]


    #####################################################
    # Bootstrap fits, enforcing fmax when necessary
    #####################################################

    # First, resplit data after filtering
    print "Splitting data into {} chunks and fitting...".format(num_cores)
    filtered_df = merged_df.loc[good_sc_df.index,:].copy()

    variant_IDs = list(set(filtered_df.variant_ID))

    # shuffle list and then make (num_cores) chunks
    random.shuffle(variant_IDs)
    chunk_list = list(make_chunks(variant_IDs, num_cores))

    # Create a new list of dataframes containing chunked variant_IDs
    grouped_list = [filtered_df[filtered_df.variant_ID.isin(chunk)].copy().groupby('variant_ID') for chunk in chunk_list]

    nboot = 1000
    print "Bootstrapping medians to get confidence intervals on fit parameters..."
    start = time.time()
    if num_cores > 1:
        fit_df_list = (Parallel(n_jobs=num_cores, verbose=10)\
            (delayed(bootstrap_fits)(
                sub_grouped, concentrations, tight_fmaxes, label_dict, fmin_val=fmin_val, nboot=nboot, plot_dir=plot_dir) for sub_grouped in grouped_list))
    else:
        fit_df_list = [bootstrap_fits(
            sub_grouped, concentrations, tight_fmaxes, label_dict, fmin_val=fmin_val, nboot=nboot, plot_dir=plot_dir) for sub_grouped in grouped_list]

    print "Fitting finished, {} minutes.".format(round((time.time() - start)/60.0, 3))
    full_fit_df = pd.concat(fit_df_list)
    full_fit_df.to_csv('bootstrapped_fits.txt', sep='\t')




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


def single_cluster_fits(grouped, x, fmin_val=0.0):
    """
    Fit single clusters for each group of variants.
    Save results of fit for use in constraining lower affinity fits.
    Input:
        grouped = pandas groupby object
        x = concentrations
    """
    results_dict = {}
    fit_model = Model(hill_equation)
    group_IDs = grouped.groups.keys()

    for vID in group_IDs:
        df = grouped.get_group(vID)
        clusterIDs = df.index.values
        data = df.iloc[:,1:].values
        nclust = data.shape[0]
        
        for c in range(nclust):
            fluorescence = data[c,:]
            # Initialize parameters as such:
            # fmax = max median fluorescence observed (minimum is fmin)
            # Kd = max protein concentration used
            # fmin = minimum median fluorescence observed (don't let it float)
            params = hill_equation_params(
                fmax={"value":max(fluorescence), "min":fmin_val},
                Kd={"value":max(x)},
                fmin={"value":fmin_val, "vary":False})
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

            results_dict[clusterIDs[c]] = {
                'variant_ID': vID,
                'Kd': fit.params['Kd'].value, 
                'fmax': fit.params['fmax'].value, 
                'fmin': fit.params['fmin'].value, 
                'rsq': rsq, 
                'rmse': rmse,
                'ier': fit.ier,
                'Kd_stderr': fit.params['Kd'].stderr,
                'fmax_stderr': fit.params['fmax'].stderr,
                'fmin_stderr': fit.params['fmin'].stderr,
            }
        
    return pd.DataFrame(results_dict).T




def bootstrap_fits(grouped, x, tight_fmaxes, label_dict, fmin_val=0.0, nboot=1000, ci=[2.5,97.5], plot_dir=None):
    """
    Fit every group in grouped using the indicated params.
    The median fit is the actually reported fit. Estimate error by 
    resampling clusters and refitting the medians, reporting 
    the 95CI of the bootstrapped parameters.
    """
    results_dict = {}
    group_IDs = grouped.groups.keys()
    tight_fmax_2p5, tight_fmax_97p5 = np.nanpercentile(tight_fmaxes, q=[2.5,97.5])

    for vID in group_IDs:
        data = grouped.get_group(vID).iloc[:,1:].values
        nclust = data.shape[0]
        median_fluorescence = np.nanmedian(data, axis=0)

        # If last median fluorescence point is below 2.5 percentile of tight binder
        # fmaxes, then enforce fmax distribution
        fmax_force = 0.0
        if median_fluorescence[-1] < tight_fmax_2p5:
            fmax_force = 1.0
            sampled_fmax = np.random.choice(tight_fmaxes, size=1)[0]
            params = hill_equation_params(
                fmax={"value": sampled_fmax, "vary": False}, 
                Kd={"value":max(x)},
                fmin={"value":fmin_val, "vary":False})
        else:
            params = hill_equation_params(
                fmax={"value":max(median_fluorescence), "min":fmin_val},
                Kd={"value":max(x)},
                fmin={"value":fmin_val, "vary":False})

        fit_model = Model(hill_equation)

        try:
            fit = fit_model.fit(median_fluorescence, params, x=x)
        except:
            print "Error while fitting {}".format(vID)
            continue
        
        ## Things we want to report
        # quality of fit:
        ss_error = np.sum((fit.residual)**2)
        ss_total = np.sum((median_fluorescence - np.nanmean(median_fluorescence))**2)
        rsq = 1 - ss_error/ss_total
        rmse = np.sqrt(np.nanmean((fit.residual)**2))

        results_dict[vID] = {
            'Kd': fit.params['Kd'].value, 
            'fmax': fit.params['fmax'].value, 
            'fmin': fit.params['fmin'].value, 
            'rsq': rsq, 
            'rmse': rmse,
            'ier': fit.ier,
            'Kd_stderr': fit.params['Kd'].stderr,
            'fmax_stderr': fit.params['fmax'].stderr,
            'fmin_stderr': fit.params['fmin'].stderr,
            'nclust': nclust,
            'forced_fmax': fmax_force
        }
        
        # Now bootstrap parameters
        med_array = np.empty((nboot, len(x)))
        Kd_array = np.empty(nboot)
        fmax_array = np.empty(nboot)
        fmin_array = np.empty(nboot)
        for b in range(nboot):
            #meds = data.sample(n=nclust, replace=True).median().values
            meds = np.nanmedian(data[np.random.choice(nclust, size=nclust, replace=True)], axis=0)
            if fmax_force:
                params = hill_equation_params(
                    fmax={"value": sampled_fmax, "vary": False}, 
                    Kd={"value":max(x)},
                    fmin={"value":fmin_val, "vary":False})
            else:
                params = hill_equation_params(
                    fmax={"value":max(meds), "min":fmin_val},
                    Kd={"value":max(x)},
                    fmin={"value":fmin_val, "vary":False})
            
            try:
                fit = fit_model.fit(meds, params, x=x)
            except:
                print "Error while bootstrap fitting {}".format(vID)
                continue
            med_array[b] = meds
            Kd_array[b] = fit.params['Kd'].value
            fmax_array[b] = fit.params['fmax'].value
            fmin_array[b] = fit.params['fmin'].value
        
        # Get confidence intervals
        # Take the 95% ci around the Kd, and whatever the corresponding fmax and fmin values are
        #ci_1_idx = np.where(Kd_array == np.nanpercentile(Kd_array, ci[0], interpolation='nearest'))[0][0]
        #ci_2_idx = np.where(Kd_array == np.nanpercentile(Kd_array, ci[1], interpolation='nearest'))[0][0]
        #ci_pos = [ci_1_idx, ci_2_idx]
        #Kd_2p5, Kd_97p5 = Kd_array[ci_pos]
        #fmax_2p5, fmax_97p5 = fmax_array[ci_pos]
        #fmin_2p5, fmin_97p5 = fmin_array[ci_pos]
        #results_dict[vID]['Kd_2p5'] = Kd_2p5
        #results_dict[vID]['Kd_97p5'] = Kd_97p5
        #results_dict[vID]['fmax_2p5'] = fmax_2p5
        #results_dict[vID]['fmax_97p5'] = fmax_97p5
        #results_dict[vID]['fmin_2p5'] = fmin_2p5
        #results_dict[vID]['fmin_97p5'] = fmin_97p5

        # Doing it the (less correct) way just looks better
        results_dict[vID]['Kd_2p5'] = np.nanpercentile(Kd_array, ci[0], interpolation='nearest')
        results_dict[vID]['Kd_97p5'] = np.nanpercentile(Kd_array, ci[1], interpolation='nearest')
        results_dict[vID]['fmax_2p5'] = np.nanpercentile(fmax_array, ci[0], interpolation='nearest')
        results_dict[vID]['fmax_97p5'] = np.nanpercentile(fmax_array, ci[1], interpolation='nearest')
        results_dict[vID]['fmin_2p5'] = np.nanpercentile(fmin_array, ci[0], interpolation='nearest')
        results_dict[vID]['fmin_97p5'] = np.nanpercentile(fmin_array, ci[1], interpolation='nearest')
        
        # Get median confidence intervals for plotting
        med_ci = np.nanpercentile(med_array, q=ci, axis=0)
        yerr = abs(median_fluorescence - med_ci)
        
        # Plot fit
        if plot_dir:
            fig, ax = plt.subplots()
            ax = plot_bootstrapped_Kd_fit(ax, x, median_fluorescence, yerr, 
                                             results_dict[vID]['Kd'], 
                                             [results_dict[vID]['Kd_2p5'], results_dict[vID]['Kd_97p5']], 
                                             results_dict[vID]['fmin'], 
                                             [results_dict[vID]['fmin_2p5'], results_dict[vID]['fmin_97p5']], 
                                             results_dict[vID]['fmax'], 
                                             [results_dict[vID]['fmax_2p5'], results_dict[vID]['fmax_97p5']], nclust, tight_fmax_97p5)
            for v in vID.split(';'):
                file_name = "/{}_{}.pdf".format(v,label_dict[v])
                plt.savefig(plot_dir+file_name, dpi=300)
            plt.close()
        
    return pd.DataFrame(results_dict).T


def plot_bootstrapped_Kd_fit(ax, x, y, y_ci, 
                                     Kd, Kd_ci,
                                     fmin, fmin_ci, 
                                     fmax, fmax_ci, nclust, max_fmax,
                                     showParams=True, showR2=True):
    """
    Plot the bootstrapped double exponential decay
    """
    x = np.array(x)
    y = np.array(y)
    if any([c == 0 for c in x]):
        zero_pos = np.where(x == 0)[0][0]
        plot_x = np.delete(x, zero_pos)
        plot_y = np.delete(y, zero_pos)
    ax.errorbar(x=x,y=y, yerr=y_ci, fmt='o', c='black', ecolor='black', elinewidth=0.8, ms=4)

    extend_min = min(plot_x) / 2.0
    extend_max = max(plot_x) * 2.0
    ax.set_xlim(extend_min,extend_max)
    fit_x = np.logspace(np.log10(extend_min), np.log10(extend_max), num=100)
    fit_y = hill_equation(fit_x, fmin, fmax, Kd)
    ax.plot(fit_x, fit_y, c="black", linestyle="--")
    ax.set_xscale('log')

    ymin, ymax = ax.get_ylim()
    if ymax > max_fmax:
        ax.set_ylim(ymin, max_fmax)
    
    rsq = 1 - np.var(hill_equation(x, fmin, fmax, Kd) - y) / np.var(y)
    
    y_lower = hill_equation(fit_x, fmin_ci[0], fmax_ci[0], Kd_ci[0])
    y_upper = hill_equation(fit_x, fmin_ci[1], fmax_ci[1], Kd_ci[1])
    ax.fill_between(fit_x, y_lower, y_upper, alpha=0.2)
    ax.set_xlabel("Concentration (nM)")
    ax.set_ylabel("Normalized Fluorescence (a.u.)")
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(0.0, ymax)

    # Make plot square
    x0,x1 = ax.get_xlim()
    y0,y1 = ax.get_ylim()
    ax.set_aspect(abs(x1-x0)/abs(y1-y0))
    
    if showParams and showR2:
        if Kd > 0.001:
            label_txt = "$K_{{D}} = {:0.3f}$ $nM$\n$R^2 = {:0.3f}$\nclusters = {}".format(
                Kd, rsq, int(nclust))
        else:
            label_txt = "$K_{{D}} = {:0.3e}$ $nM$\n$R^2 = {:0.3f}$\nclusters = {}".format(
                Kd, rsq, int(nclust))
        if Kd > np.median(plot_x):
            ax.text(0.05, 0.95, label_txt, transform=ax.transAxes, 
                    verticalalignment='top', horizontalalignment='left', fontsize=12, 
                    bbox={'facecolor': ax.get_facecolor(), 'alpha': 1.0, 'pad': 10, 'edgecolor':'none'})
        else:
            ax.text(0.95, 0.05, label_txt, transform=ax.transAxes, 
                    verticalalignment='bottom', horizontalalignment='right', fontsize=12, 
                    bbox={'facecolor': ax.get_facecolor(), 'alpha': 1.0, 'pad': 10, 'edgecolor':'none'})
    
    return ax



def make_chunks(l, n):
    # Make n chunks of list l
    ll = len(l)
    chunk_size = int(np.ceil(ll / float(n)))
    for i in range(0, ll, chunk_size):
        yield l[i:i+chunk_size]

def guide_group(guide_label):
    if isinstance(guide_label, basestring):
        group, num = guide_label.split('.')
        return group
    else:
        return "no group"


def contains_variant_ID(varID, valid_IDs):
    ids = varID.split(';')
    if any([v in valid_IDs for v in ids]):
        return True
    return False



if __name__ == '__main__':
    main()
