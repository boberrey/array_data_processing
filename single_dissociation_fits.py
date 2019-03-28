#!/usr/bin/env python
""" 
Perform single phase dissociation curve fits across a CPseries file

 Ben Ober-Reynolds, boberrey@stanford.edu
 20190325
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
import time



### MAIN ###


def main():
	################ Parse input parameters ################

	#set up command line argument parser
	parser = argparse.ArgumentParser(description='Script for fitting two-phase dissociation for a CPseries file')
	group = parser.add_argument_group('required arguments:')
	group.add_argument('-cs', '--cpseries', required=True,
						help='CPseries.pkl file')
	group.add_argument('-ca', '--cpannot', required=True,
						help='CPannot.pkl file')
	group.add_argument('-td', '--time_dict', required=True,
						help='pickled dictionary of image times')
	group.add_argument('-vt', '--variant_table', required=True,
						help='experiment variant table')

	group = parser.add_argument_group('optional arguments for processing data')
	group.add_argument('-n','--num_cores', type=int, default=10,
						help='number of cores to use')
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
	diss_time_file = args.time_dict
	variant_file = args.variant_table
	num_cores = args.num_cores

	# Make the plot directory if plotting
	if args.plots:
		if not os.path.isdir(args.plot_dir):
			print "Plot directory {} does not exist. Creating...".format(args.plot_dir)
			os.mkdir(args.plot_dir)
		plot_dir = args.plot_dir
	else:
		plot_dir = None

	# Read in data
	print "Reading in data..."
	diss_df = pd.read_pickle(cpseries_file)
	# drop any rows with all nans
	diss_df.dropna(axis=0, how='all', inplace=True)

	diss_times = pd.read_pickle(diss_time_file)
	# Just use the 9th image times
	diss_times = diss_times['009']

	# If series length is one more than diss_times, assume first point is the baseline
	series_points = len(diss_df.columns)
	if series_points > len(diss_times):
		if series_points == len(diss_times) + 1:
			print "Series length and times don't match. Setting first point in series to 0.0 seconds..."
			diss_times = [0.0] + diss_times
		else:
			print "Series length and times differ by more than 1. Are these the right files?"
			sys.exit()


	annot_df = pd.read_pickle(cpannot_file)
	variant_df = pd.read_table(variant_file, header=None)
	variant_df.columns = ['Sequence', 'variant_ID', 'group_name', 'mut_annotation', 'wt_seq', 'guide_number']

	if args.multiguide_group:
		print "Only fitting variants in multiguide group {}".format(args.multiguide_group)
		group_ids = variant_df[variant_df.guide_number.apply(guide_group) == args.multiguide_group].variant_ID.values
		annot_df = annot_df[annot_df.variant_ID.isin(group_ids)].copy()


	# Merge annot and cpseries
	merged_df = annot_df.merge(diss_df, left_index=True, right_index=True, how='inner')

	# Normalize everything by fiducial signal
	print "Normalizing to fiducial..."
	fiducial_meds = merged_df.groupby('variant_ID').get_group('11111111').iloc[:,1:].median().values
	merged_df.iloc[:,1:] = merged_df.iloc[:,1:] / fiducial_meds


	# groupby variant_ID
	grouped = merged_df.groupby('variant_ID')

	
	# Create params for fitting
	params = single_exp_decay_params()

	print "Splitting data into {} chunks and fitting...".format(num_cores)
	# Now perform the actual fits
	variant_IDs = grouped.groups.keys()

	# shuffle list and then make (num_cores) chunks
	random.shuffle(variant_IDs)
	chunk_list = list(make_chunks(variant_IDs, num_cores))

	# Create a new list of dataframes containing chunked variant_IDs
	grouped_list = [merged_df[merged_df.variant_ID.isin(chunk)].copy().groupby('variant_ID') for chunk in chunk_list]

	# Now fit variants in parallel
	nboot = 1000
	start = time.time()
	if num_cores > 1:
		fit_df_list = (Parallel(n_jobs=num_cores, verbose=10)\
			(delayed(bootstrap_fits)(
				sub_grouped, diss_times, params, nboot=nboot, plot_dir=plot_dir) for sub_grouped in grouped_list))
	else:
		fit_df_list = [bootstrap_fits(
			sub_grouped, diss_times, params, nboot=nboot, plot_dir=plot_dir) for sub_grouped in grouped_list]

	print "Fitting finished, {} minutes.".format(round((time.time() - start)/60.0, 3))
	full_fit_df = pd.concat(fit_df_list)
	full_fit_df.to_csv('double_exp_koff_fits_uptoa2.txt', sep='\t')




# Fitting functions:

def single_exp_decay(x, fmin, fmax, koff):
	"""
	Single exponential decay function. 
	Functions provided to lmfit Model must have dependent variable as first argument.
	x = time points 
	fmin = decay floor
	fmax = max signal
	koff = decay rate
	"""
	return fmin + (fmax - fmin)*np.exp(-koff*x)


def single_exp_decay_params(fmax=None, koff=None):
	# Define parameters object
	params = Parameters()
	default_params = {
		"fmax":{"value": 1.0, "vary": True, "min": -np.inf, "max": np.inf},
		"span":{"value": 1.0, "vary": True, "min":0.0, "max":np.inf},
		"koff":{"value": 0.001, "vary": True, "min": -np.inf, "max": np.inf}

	}
	if fmax:
		for opt, val in fmax.items():
			default_params["fmax"][opt] = val
	if koff:
		for opt, val in koff.items():
			default_params["koff"][opt] = val

	for p, dct in default_params.items():
		params.add(p, **dct)

	# Enforce that fmax > fmin
	params.add("fmin", value=0.0, expr="fmax - span")
	return params



def bootstrap_fits(grouped, x, params, nboot=1000, ci=[2.5,97.5], plot_dir=None):
	"""
	Fit every group in grouped using the indicated params.
	The median fit is the actually reported fit. Estimate error by 
	resampling clusters and refitting the medians, reporting 
	the 95CI of the bootstrapped parameters.
	"""
	results_dict = {}
	group_IDs = grouped.groups.keys()
	for vID in group_IDs:
		data = grouped.get_group(vID).iloc[:,1:].values
		nclust = data.shape[0]
		median_fluorescence = np.nanmedian(data, axis=0)

		# Heuristics to decide if fit should happen
		# (add here)
		fit_model = Model(single_exp_decay)

		try:
			fit = fit_model.fit(median_fluorescence, params, x=x)
		except:
			print "Error while fitting {}".format(vID)
			continue
		
		# Things we want to report
		results_dict[vID] = {
			'koff': fit.params['koff'].value, 
			'fmax': fit.params['fmax'].value, 
			'fmin': fit.params['fmin'].value, 
			'diss_rsq': 1 - fit.residual.var() / np.var(median_fluorescence), 
			'ier': fit.ier,
			'koff_stderr': fit.params['koff'].stderr,
			'fmax_stderr': fit.params['fmax'].stderr,
			'fmin_stderr': fit.params['fmin'].stderr,
			'nclust': nclust,
		}
		
		# Now bootstrap parameters
		med_array = np.empty((nboot, len(x)))
		koff_lst = []
		fmax_lst = []
		fmin_lst = []
		for b in range(nboot):
			#meds = data.sample(n=nclust, replace=True).median().values
			meds = np.nanmedian(data[np.random.choice(nclust, size=nclust, replace=True)], axis=0)
			try:
				fit = fit_model.fit(meds, params, x=x)
			except:
				print "Error while bootstrap fitting {}".format(vID)
				continue
			med_array[b] = meds
			koff_lst.append(fit.params['koff'].value)
			fmax_lst.append(fit.params['fmax'].value)
			fmin_lst.append(fit.params['fmin'].value)
		
		# Get confidence intervals
		koff_2p5, koff_97p5 = np.nanpercentile(koff_lst, q=ci)
		fmax_2p5, fmax_97p5 = np.nanpercentile(fmax_lst, q=ci)
		fmin_2p5, fmin_97p5 = np.nanpercentile(fmin_lst, q=ci)
		results_dict[vID]['koff_2p5'] = koff_2p5
		results_dict[vID]['koff_97p5'] = koff_97p5
		results_dict[vID]['fmax_2p5'] = fmax_2p5
		results_dict[vID]['fmax_97p5'] = fmax_97p5
		results_dict[vID]['fmin_2p5'] = fmin_2p5
		results_dict[vID]['fmin_97p5'] = fmin_97p5
		
		# Get median confidence intervals for plotting
		med_ci = np.nanpercentile(med_array, q=ci, axis=0)
		yerr = abs(median_fluorescence - med_ci)
		
		# Plot fit
		if plot_dir:
			fig, ax = plt.subplots()
			ax = plot_bootstrapped_single_exp_fit(ax, x, median_fluorescence, yerr, 
											 results_dict[vID]['koff'], 
											 [koff_2p5, koff_97p5], 
											 results_dict[vID]['fmin'], 
											 [fmin_2p5, fmin_97p5], 
											 results_dict[vID]['fmax'], 
											 [fmax_2p5, fmax_97p5])
			plt.savefig(plot_dir+'/{}_full.pdf'.format(vID), dpi=300)
			if results_dict[vID]['koff'] > 0.001:
			   ax.set_xlim(0,7200)
			   plt.savefig(plot_dir+'/{}_2h.pdf'.format(vID), dpi=300)
			plt.close()
		
	return pd.DataFrame(results_dict).T


def plot_bootstrapped_single_exp_fit(ax, x, y, y_ci, 
									 koff, koff_ci,
									 fmin, fmin_ci, 
									 fmax, fmax_ci, 
									 showParams=True, showR2=True):
	"""
	Plot the bootstrapped double exponential decay
	"""
	x = np.array(x)
	y = np.array(y)
	ax.errorbar(x=x,y=y, yerr=y_ci, fmt='o', c='black', ecolor='black', elinewidth=0.8, ms=4)
	xmin, xmax = ax.get_xlim()
	ax.set_xlim(-10.0,xmax)
	fit_x = np.linspace(0.0, xmax, 1000)
	fit_y = single_exp_decay(fit_x, fmin, fmax, koff)
	ax.plot(fit_x, fit_y, c="black", linestyle="--")
	
	rsq = 1 - np.var(single_exp_decay(x, fmin, fmax, koff) - y) / np.var(y)
	
	y_lower = single_exp_decay(fit_x, fmin_ci[0], fmax_ci[0], koff_ci[0])
	y_upper = single_exp_decay(fit_x, fmin_ci[1], fmax_ci[1], koff_ci[1])
	ax.fill_between(fit_x, y_lower, y_upper, alpha=0.2)
	ax.set_xlabel("Time (s)")
	ax.set_ylabel("Normalized Fluorescence (a.u.)")
	ymin, ymax = ax.get_ylim()
	ax.set_ylim(0.0, ymax)
	half_life = np.log(2.0)/koff
	
	if showParams and showR2:
		label_txt = "$k_{{off}} = {:0.3e}$ $s^{{-1}}$\n$t_{{1/2}} = {:0.3f}$ $min$\n$R^2 = {:0.3f}$".format(
				koff, half_life/60.0, rsq)
		ax.text(0.95, 0.95, label_txt, transform=ax.transAxes, 
				verticalalignment='top', horizontalalignment='right', fontsize=12, 
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


if __name__ == '__main__':
	main()