"""
Plotting functions and other utilities for having and easier time with
matplotlib and seaborn
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import scipy as sp
import re


# Utility functions

def rev_comp(seq, complement_dict={'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}, missing='N'):
	# Reverse complement a sequence
	return "".join(complement_dict.get(base, missing) for base in reversed(seq))

def texttoint(text):
	# Convert string integers to true integers, if possible
	return int(text) if text.isdigit() else text

def natural_keys(text):
	# Custom key for sorting values that are strings, but contain integers
	# Sorts based on the integer identified in a larger string
	return [texttoint(c) for c in re.split('(\d+)', text)]

def invert_annot_numbering(annot, max_num):
	"""
	Specific to Ago libraries. Original library annotations were target-focused.
	Often times, you want to think about things from the guide persepctive.
	This function alters the mutation position numbering to be guide-focused.
	"""
	context, swaps = annot.split('_')
	all_swaps = swaps.split(',')
	new_swaps = []
	for swap in all_swaps:
		num, text = re.split('(\d+)', swap)[1:]
		new_swaps.append(str(max_num + 1 - int(num)) + text)
	return context + '_' + ','.join(new_swaps)

# R in kcal / mol*K

def dG_to_Kd(dG, temp=37, R=0.0019858775):
	# Returns the Kd value in M
	return (1/np.exp(-dG/(R*(temp + 273))))


def Kd_to_dG(Kd, temp=37, R=0.0019858775):
	# Returns the dG value in kcal/mol, assuming Kd is in M
	return (R*(temp + 273)*np.log(Kd))

def convert_dG_to_kT(dG, temp=37):
	# k is the boltzman constant:
	# k = 0.0019872 kcal / mol * K
	# by default, dG is in kcal/mol
	return dG / (0.0019872*(273 + temp))

# Plotting functions




"""
Ago project plotting

These plotting functions are more specific to the Ago projects Winston and I 
have been working on. Be warned that some of them may have functionality that
is fairly specific to those kinds of data.
"""

### Double Mutant Heatmap Plotting ###

def parse_point_mut_annot(annotation):
	"""
	Pull out the individual mutations from a point mutant annotation
	"""
	swaps = annotation.split('_')[-1]
	individual_swaps = swaps.split(',')
	return individual_swaps


def construct_double_mut_df(dm_df, value):
	"""
	Construct a square data frame for a set of double mutants

	20180224 needs fixing
	"""
	
	# split mutations annotations and get value of interest:
	list_data = dm_df.apply(lambda x: \
		parse_point_mut_annot(x, value) + list(x.value), axis=1).tolist()
	df_data = pd.DataFrame(list_data)
	df_data.columns = ['first_mut', 'second_mut', 'value']
	
	# Get a sorted list of all possible base swaps:
	keys = sorted(list(set([x[0] for x in list_data] + [x[1] for x in list_data])), 
		key=natural_keys)

	# Fill in square dataframe with desired double mutant values
	new_frame = pd.DataFrame(index=keys[::-1], columns=keys)

	def fill_in_dm_frame(row, df_to_fill):
		# Apply function for filling in a square dataframe
		mut1 = row['first_mut']
		mut2 = row['second_mut']
		value = row['value']
		df_to_fill.loc[mut1, mut2] = value
		df_to_fill.loc[mut2, mut1] = value

	df_data.apply(lambda x: fill_in_frame(x, new_frame), axis=1)      
	
	return new_frame


def combine_double_mut_dfs(upper_df, lower_df):
    """
    If you have two square dataframes, this function lets you merge them
    into one dataframe (across a BL to TR diagonal)
    Maintains the indices of the 'upper_df' and the columns of the 'lower_df'
    """
    nrows = len(upper_df.index)
    ncols = len(upper_df.columns)
    combined_df = upper_df.copy()
    
    for i in range(nrows):
        for j in range(ncols):
            if nrows - i < j:
                combined_df.iloc[i,j] = lower_df.iloc[i,j]
    
    combined_df.columns = lower_df.columns
    
    return combined_df
