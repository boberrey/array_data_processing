#!/usr/bin/env python
""" 
Perform two-phase dissociation curve fits across a CPseries file

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
	group.add_argument('-ck2','--constrain_k2', action="store_true",
						help='Constrain the second rate to be a constant value for all variants')
	group.add_argument('-ca2','--constrain_a2', action="store_true",
						help='Constrain the amplitude of the second rate to be a constant value for all variants')
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
		group_ids = list(variant_df[variant_df.guide_number.apply(guide_group) == args.multiguide_group].variant_ID.values)
		group_ids.append('11111111')
		annot_df = annot_df[annot_df.variant_ID.isin(group_ids)].copy()

	# Merge annot and cpseries
	merged_df = annot_df.merge(diss_df, left_index=True, right_index=True, how='inner')

	# Normalize everything by fiducial signal
	print "Normalizing to fiducial..."
	fiducial_meds = merged_df.groupby('variant_ID').get_group('11111111').iloc[:,1:].median().values
	merged_df.iloc[:,1:] = merged_df.iloc[:,1:] / fiducial_meds


	# groupby variant_ID
	grouped = merged_df.groupby('variant_ID')

	k2 = None
	a2 = None

	if args.constrain_k2 or args.constrain_a2:
		start = time.time()
		print "Estimating constrained parameters..."
		to_constrain = []
		if args.constrain_k2:
			to_constrain.append("k2")
		if args.constrain_a2:
			to_constrain.append("a2")
		# If constraining one or more parameters, bootstrap global fits to get estimate for constrained value
		nboot = 500
		fit_size = 250

		# Only use variants that have reasonably high binding to estimate constrained parameters
		median_df = grouped.agg(np.median)
		median_df.drop('11111111', inplace=True)
		valid_data = median_df[median_df.iloc[:,0] > 0.5].copy()

		# sample as many datasets as required
		datasets = []
		for i in range(nboot):
			data_ids = np.random.choice(valid_data.index.values, size=fit_size, replace=False)
			datasets.append(np.array(valid_data.loc[data_ids,:].copy()))
			#print valid_data.loc[data_ids,:].copy().head()

		if num_cores > 1:
			print "Estimating constrained parameters on {} cores...".format(num_cores)
			param_estimates = (Parallel(n_jobs=num_cores, verbose=10)\
				(delayed(estimate_params_global)(
					ds, diss_times, to_constrain) for ds in datasets))
		else:
			print "Estimating constrained parameters on a single core"
			param_estimates = [estimate_params_global(
				ds, diss_times, to_constrain) for ds in datasets]

	
		param_df = pd.DataFrame(param_estimates)
		param_df.columns = to_constrain
		param_df.to_csv('param_estimates.txt', sep='\t', index=False)
		if args.constrain_k2:
			k2_low, k2_high = np.nanpercentile(param_df.k2.values, q=[5, 95])
			k2 = {"value":param_df.k2.median(), "vary":True, "min":k2_low, "max": k2_high}
		if args.constrain_a2:
			#a2 = {"value":param_df.a2.median(), "vary":False}
			a2_low, a2_high = np.nanpercentile(param_df.k2.values, q=[5, 95])
			a2 = {"value":param_df.a2.median(), "vary":False, "min": a2_low, "max": a2_high}
		print "Constrained parameters fit, {} minutes".format(round((time.time() - start)/60.0, 3))
	
	# Create params for fitting
	params = double_exp_decay_params(k2=k2, a2=a2)

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
			(delayed(bootstrap_two_phase_fits)(
				sub_grouped, diss_times, params, nboot=nboot, plot_dir=plot_dir) for sub_grouped in grouped_list))
	else:
		fit_df_list = [bootstrap_two_phase_fits(
			sub_grouped, diss_times, params, nboot=nboot, plot_dir=plot_dir) for sub_grouped in grouped_list]

	print "Fitting finished, {} minutes.".format(round((time.time() - start)/60.0, 3))
	full_fit_df = pd.concat(fit_df_list)
	full_fit_df.to_csv('double_exp_koff_fits_uptoa2.txt', sep='\t')




# Fitting functions:
def estimate_params_global(data, x_vals, param_list=["k2", "a2"]):
	"""
	Estimate params in params list by fitting them globally using a random
	sample of variants
	"""
	nvar = data.shape[0]
	# Create as many sets of parameters as you have datasets
	fit_params = Parameters()
				
	for i in range(nvar):
		fit_params.add('fmin_{}'.format(i+1), value=0.01, min=-np.inf, max=np.inf)
		fit_params.add('fmax_{}'.format(i+1), value=1.0, min=-np.inf, max=np.inf)
		fit_params.add('k1_{}'.format(i+1), value=0.001, min=-np.inf, max=np.inf)
		fit_params.add('k2_{}'.format(i+1), value=0.001, min=-np.inf, max=np.inf)
		fit_params.add('a2_{}'.format(i+1), value=0.1, min=0.0, max=1.0)

	# Now constrain certain parameters
	for i in range(1, nvar):
		for p in param_list:
			fit_params[p + '_{}'.format(i+1)].expr=p+'_1'

	# Run global fit
	x = np.array(x_vals)

	try:
		result = minimize(global_double_exp_objective, fit_params, args=(x, data), 
			method='leastsq', maxfev=data.size*3)
	except:
		# Problem in fitting?
		to_return = []
		for p in param_list:
			to_return.append(np.nan)
			return to_return
	# Extract and return the globally constrained params
	to_return = []
	for p in param_list:
		to_return.append(result.params['{}_1'.format(p)].value)
	return to_return




def double_exp_decay(x, fmin, fmax, k1, k2, a2):
	"""
	Double exponential decay function. 
	x = timepoints (s)
	fmin = minimum signal
	fmax = max signal
	k1 = specific rate
	k2 = non-specific rate
	a2 = amplitude (percent of total signal) of non-specific rate
	"""
	amp = fmax - fmin
	a1 = 1 - a2
	return fmin + a1*amp*np.exp(-k1*x) + a2*amp*np.exp(-k2*x)


def double_exp_decay_params(k2=None, a2=None):
	# Define parameters object
	params = Parameters()
	default_params = {
		"fmax":{"value": 1.1, "vary": True, "min": -np.inf, "max": np.inf},
		"k1":{"value": 0.0001, "vary": True, "min": -np.inf, "max": np.inf},
		"k2":{"value": 0.0001, "vary": True, "min": -np.inf, "max": np.inf},
		"a2":{"value": 0.1, "vary": True, "min": 0.0, "max": 1.0},
		"span":{"value": 1.0, "vary": True, "min":0.0, "max":np.inf}
	}
	if k2:
		for opt, val in k2.items():
			default_params["k2"][opt] = val
	if a2:
		for opt, val in a2.items():
			default_params["a2"][opt] = val
	
	for p, dct in default_params.items():
		params.add(p, **dct)
	# Enforce that fmax > fmin
	params.add("fmin", value=0.0, expr="fmax - span")
	return params


def double_exp_fit_data(params, i, x):
	"""
	Retrieve the parameter set for variant i
	"""
	fmin = params['fmin_{}'.format(1 + i)].value
	fmax = params['fmax_{}'.format(1 + i)].value
	k1 = params['k1_{}'.format(1 + i)].value
	k2 = params['k2_{}'.format(1 + i)].value
	a1 = params['a2_{}'.format(1 + i)].value
	return double_exp_decay(x, fmin, fmax, k1, k2, a1)


def global_double_exp_objective(params, x, y):
	"""
	Objective function for fitting multiple curves with some shared parameters
	y is n x m dimensional, where n is the number of variants, and 
	m is the number of points per variant.
	x is 1 x m dimensional
	"""
	ndata, nx = y.shape
	resid = 0.0*y[:]
	# Calculate residuals per dataset
	for i in range(ndata):
		resid[i,:] = (y[i,:] - double_exp_fit_data(params, i, x))
	# Residuals must be returned as a 1D array
	return resid.flatten()


def bootstrap_two_phase_fits(grouped, x, params, nboot=1000, ci=[2.5,97.5], plot_dir=None):
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
		# (Add here?)

		fit_model = Model(double_exp_decay)
		fit = fit_model.fit(median_fluorescence, params, x=x)
		
		# Things we want to report
		results_dict[vID] = {
			'k1': fit.params['k1'].value,
			'k2': fit.params['k2'].value,
			'a2': fit.params['a2'].value,
			'fmax': fit.params['fmax'].value, 
			'fmin': fit.params['fmin'].value, 
			'diss_rsq': 1 - fit.residual.var() / np.var(median_fluorescence), 
			'ier': fit.ier,
			'k1_stderr': fit.params['k1'].stderr,
			'k2_stderr': fit.params['k2'].stderr,
			'a2_stderr': fit.params['a2'].stderr,
			'fmax_stderr': fit.params['fmax'].stderr,
			'fmin_stderr': fit.params['fmin'].stderr,
			'nclust': nclust,
		}
		
		# Now bootstrap parameters
		med_array = np.empty((nboot, len(x)))
		k1_lst = []
		k2_lst = []
		a2_lst = []
		fmax_lst = []
		fmin_lst = []
		for b in range(nboot):
			meds = np.nanmedian(data[np.random.choice(nclust, size=nclust, replace=True)], axis=0)
			try:
				fit = fit_model.fit(meds, params, x=x)
			except:
				print "Error while bootstrap fitting {}".format(vID)
				continue
			med_array[b] = meds
			k1_lst.append(fit.params['k1'].value)
			k2_lst.append(fit.params['k2'].value)
			a2_lst.append(fit.params['a2'].value)
			fmax_lst.append(fit.params['fmax'].value)
			fmin_lst.append(fit.params['fmin'].value)
		
		# Get confidence intervals
		k1_2p5, k1_97p5 = np.nanpercentile(k1_lst, q=ci)
		k2_2p5, k2_97p5 = np.nanpercentile(k2_lst, q=ci)
		a2_2p5, a2_97p5 = np.nanpercentile(a2_lst, q=ci)
		fmax_2p5, fmax_97p5 = np.nanpercentile(fmax_lst, q=ci)
		fmin_2p5, fmin_97p5 = np.nanpercentile(fmin_lst, q=ci)
		results_dict[vID]['k1_2p5'] = k1_2p5
		results_dict[vID]['k1_97p5'] = k1_97p5
		results_dict[vID]['k2_2p5'] = k2_2p5
		results_dict[vID]['k2_97p5'] = k2_97p5
		results_dict[vID]['a2_2p5'] = a2_2p5
		results_dict[vID]['a2_97p5'] = a2_97p5
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
			ax = plot_bootstrapped_double_exp_fit(ax, x, median_fluorescence, yerr, 
											 results_dict[vID]['k1'], 
											 [k1_2p5, k1_97p5], 
											 fit.params['k2'].value, 
											 [k2_2p5, k2_97p5],
											 fit.params['a2'].value, 
											 [a2_2p5, a2_97p5],
											 results_dict[vID]['fmin'], 
											 [fmin_2p5, fmin_97p5], 
											 results_dict[vID]['fmax'], [fmax_2p5, fmax_97p5])
			plt.savefig(plot_dir+'/{}_full.pdf'.format(vID), dpi=300)
			if results_dict[vID]['k1'] > 0.001:
				ax.set_xlim(0,7200)
				plt.savefig(plot_dir+'/{}_2h.pdf'.format(vID), dpi=300)
			plt.close()
		
	return pd.DataFrame(results_dict).T


def plot_bootstrapped_double_exp_fit(ax, x, y, y_ci, 
									 k1, k1_ci, 
									 k2, k2_ci,
									 a2, a2_ci,
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
	ax.set_xlim(0.0,xmax)
	fit_x = np.linspace(0.0, xmax, 1000)
	fit_y = double_exp_decay(fit_x, fmin, fmax, k1, k2, a2)
	ax.plot(fit_x, fit_y, c="black", linestyle="--")
	
	rsq = 1 - np.var(double_exp_decay(x, fmin, fmax, k1, k2, a2) - y) / np.var(y)
	
	y_lower = double_exp_decay(fit_x, fmin_ci[0], fmax_ci[0], k1_ci[0], k2_ci[0], a2_ci[0])
	y_upper = double_exp_decay(fit_x, fmin_ci[1], fmax_ci[1], k1_ci[1], k2_ci[1], a2_ci[1])
	ax.fill_between(fit_x, y_lower, y_upper, alpha=0.2)
	ax.set_xlabel("Time (s)")
	ax.set_ylabel("Normalized Fluorescence (a.u.)")
	ymin, ymax = ax.get_ylim()
	ax.set_ylim(0.0, ymax)
	
	if showParams and showR2:
		label_txt = "$k_{{off,1}} = {:0.3e}$ $s^{{-1}}$\n$k_{{off,2}} = {:0.3e}$ $s^{{-1}}$ \na = {:0.3f}\n$R^2 = {:0.3f}$".format(
				k1, k2, a2, rsq)
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