#!/usr/bin/env python
""" 
Fit association models for TtAgo
"""

import pandas as pd
import numpy as np
import scipy as sp
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import mpl_cookbook as mplcb
import nupack_wrappers as npwp
import re
import time
from collections import OrderedDict
from sklearn.mixture import GaussianMixture
from lmfit import Parameters, minimize, Model, report_fit
from joblib import Parallel, delayed
import pickle



### MAIN ###


def main():
    ################ Parse input parameters ################

    # Read in filtered data
    print("Reading and preparing data...")

    all_data_df = pd.read_csv("/raid/USRdirs/ago/TtAgo/paper_figures/final_filtered_data.tsv", sep='\t')
    outdir = "/raid/USRdirs/ago/TtAgo/paper_figures/figure_2/model_figures/"

    group_names = [
        #'G1.1', 
        #'G1.2', 
        #'G1.3', 
        #'G1.4', 
        #'G1.5', 
        'G2.1', 
        'G2.2', 
        'G2.3', 
        'G2.4', 
        'G2.5'
    ]

    wt_var_ids = {
    #'G1.1': 'FQ8N317Z',
    #'G1.2': 'P6AIWCVY',
    #'G1.3': '29312X1M',
    #'G1.4': '1FV2VMAZ',
    #'G1.5': 'TASZ696K',
    'G2.1': 'QCUCZ0A5',
    'G2.2': 'YUQBTUFX',
    'G2.3': 'CRMJY4Y2',
    'G2.4': '3D2ZBOOW',
    'G2.5': 'GVPYRD08'
    }

    colors = ['#d51f27', "#262f68","#238a43","#89298d","#f37c29"]
    color_dict = dict(zip(group_names, colors))

    # Fix G2.3 to be the correct sequence
    wrong_G2p3 = 'TGATGCCAGGTAATTT'
    fixed_G2p3 = 'TGATGCCAGGTAATTG'
    all_data_df.replace(wrong_G2p3, fixed_G2p3, inplace=True)

    # Run filters
    all_data_df.kon_oligo_37 = all_data_df.apply(lambda x: x.kon_oligo_37 if x.kon_filter_oligo_37 else np.nan, axis=1)
    all_data_df.kon_ttago_37 = all_data_df.apply(lambda x: x.kon_ttago_37 if x.kon_filter_ttago_37 else np.nan, axis=1)

    # Get nupack data

    # Convert kcal/mol to kT
    nupack_variant_file = "/raid/USRdirs/ago/TtAgo/paper_figures/variant_nupack_df.tsv"
    nupack_kmer_file = "/raid/USRdirs/ago/TtAgo/paper_figures/kmer_nupack_df.tsv"

    nupack_var_df = pd.read_table(nupack_variant_file)
    nupack_kmer_df = pd.read_table(nupack_kmer_file)

    all_data_df = all_data_df.merge(nupack_var_df[['variant_ID', 'target_ensemble_37', 
                                                   'rnafold_full_37', 'rnafold_8mer_block_37', 
                                                   'rnafold_6mer_block_37', 'rnafold_4mer_block_37', 
                                                   'seed_4mer_ensemble_37', 'seed_5mer_ensemble_37', 
                                                   'seed_6mer_ensemble_37', 'seed_7mer_ensemble_37']], on='variant_ID', how='inner')

    # Convert energy values to kT
    all_data_df.target_ensemble_37 = all_data_df.target_ensemble_37.apply(lambda x: convert_dG_to_kT(x, temp=37))
    nupack_kmer_df.kmer_ensemble = nupack_kmer_df.kmer_ensemble.apply(lambda x: convert_dG_to_kT(x, temp=37))
    nupack_kmer_df.kmer_mfe = nupack_kmer_df.kmer_mfe.apply(lambda x: convert_dG_to_kT(x, temp=37))

    wt_ttago_kon_37 = get_wt_params(all_data_df, wt_var_ids, 'kon_ttago_37')
    wt_oligo_kon_37 = get_wt_params(all_data_df, wt_var_ids, 'kon_oligo_37')

    wt_oligo_dG_37 = get_wt_params(all_data_df, wt_var_ids, 'dG_kt_oligo_37')
    wt_ttago_dG_37 = get_wt_params(all_data_df, wt_var_ids, 'dG_kt_ttago_37')
    wt_ttago_dG_55 = get_wt_params(all_data_df, wt_var_ids, 'dG_kt_ttago_55')

    wt_kc_green_55 = get_wt_params(all_data_df, wt_var_ids, 'kc_green_55')
    wt_kc_red_55 = get_wt_params(all_data_df, wt_var_ids, 'kc_red_55')
    wt_kc_green_65 = get_wt_params(all_data_df, wt_var_ids, 'kc_green_65')
    wt_kc_red_65 = get_wt_params(all_data_df, wt_var_ids, 'kc_red_65')

    wt_seqs = get_wt_params(all_data_df, wt_var_ids, 'sequence')
    wt_struct = get_wt_params(all_data_df, wt_var_ids, 'target_ensemble_37')

    # Delta kon rate from WT
    all_data_df['delta_kon_ttago_wt'] = all_data_df.apply(lambda x: delta_wt_param(x, 'kon_ttago_37', wt_ttago_kon_37), axis=1)

    # Get oligo/TtAgo delta kon
    all_data_df['delta_kon_oligo'] = all_data_df.kon_ttago_37 - all_data_df.kon_oligo_37

    # Delta predicted structure from WT
    all_data_df['delta_target_ensemble_37'] = all_data_df.apply(lambda x: delta_wt_param(x, 'target_ensemble_37', wt_struct), axis=1)

    point_mut_groups = ['single_mutants', 'double_mutants', 'tandem_triples', 
                     'transition_triples', 'complement_triples']

    non_tandem_pm_groups = ['single_mutants', 'double_mutants', 
                        'transition_triples', 'complement_triples']

    point_mutant_df = all_data_df[all_data_df.mut_group.isin(point_mut_groups)].copy()
    nt_point_mutant_df = all_data_df[all_data_df.mut_group.isin(non_tandem_pm_groups)].copy()

    # Perform X-fold cross validated fitting holding out 1 guide each time

    # Restrict to point mutants that have a measured TtAgo association rate
    mpm_df = point_mutant_df[~np.isnan(point_mutant_df.kon_ttago_37)].copy()

    # Remove parameters at extremes
    qlow, qhigh = np.percentile(mpm_df.delta_kon_ttago_wt.values, q=[1,99])
    mpm_df = mpm_df[(mpm_df.delta_kon_ttago_wt > qlow) & (mpm_df.delta_kon_ttago_wt < qhigh)].copy()

    # Required columns:
    required_cols = ['sequence', 'guide_seq', 'delta_target_ensemble_37', 'delta_kon_ttago_wt']


    # Construct kmer dict
    kd = dict(zip(nupack_kmer_df.kmer_duplex.tolist(), nupack_kmer_df.kmer_ensemble.tolist()))

    # All possible guides:
    all_guides = ['G2.1', 'G2.2', 'G2.3', 'G2.4', 'G2.5']
    #all_guides = ['G2.1', 'G2.2', 'G2.3', 'G2.5']

    kmer_len = 6

    print("Using the following guides for fitting:")
    for g in all_guides:
        print(g)

    lmda = [0.0, 0.5, 1.0, 2.5, 5.0, 7.5, 10, 15, 20, 25]
    guide_combos = [all_guides]
    combo_names = ['all_guides']

    for i, g in enumerate(all_guides):
        guide_combos.append(all_guides[:i] + all_guides[i+1:])
        combo_names.append('{}_held_out'.format(g))

    train_idxes = [mpm_df.guide_num.isin(x).values for x in guide_combos]
    train_sets = [{'lmda':l,'combo_name':cn,'train_idx':ti, 'data_df': mpm_df[required_cols].copy()} for l in lmda for cn,ti in zip(combo_names, train_idxes)]
    set_names = ["{}_lmda{}".format(l['combo_name'], l['lmda']) for l in train_sets]

    print("Will fit the following datasets:")
    for s in set_names:
        print(s)

    num_cores = min(25, len(set_names))

    print("Fitting {} models on {} cores...".format(len(set_names), num_cores))
    fit_list = (Parallel(n_jobs=num_cores, verbose=10)\
            (delayed(pos_by_swap_with_kmers)(
                x['data_df'], train_idx=x['train_idx'], kmer_lst=[(1,kmer_len)], kmer_dict=kd, lmda=x['lmda']) for x in train_sets))

    fit_dict = dict(zip(set_names, fit_list))

    print("Finished fitting, saving fit models as pickled dict...")
    with open('model_fits.pkl', 'wb') as f:
        pickle.dump(fit_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

    


# Misc functions
def convert_dG_to_kT(dG, temp=37, R=0.0019872):
        # k is the boltzman constant:
        # k = 0.0019872 kcal / mol * K
        # by default, dG is in kcal/mol
        return dG / (R*(273 + temp))


def get_wt_params(df, wt_ids, param):
    d = {}
    for nm, var_id in wt_ids.items():
        d[nm] = df[df.variant_ID == var_id][param].values[0]
    return d


def delta_wt_param(row, param, wt_dict):
    """
    Subtract 'param' from the corresponding wt value in wt_dict
    """
    return row[param] - wt_dict[row.guide_num]


# Functions for position by swap model fitting

def vec_translate(a, my_dict):
    """
    Vectorized dictionary lookup.
    See here:
    https://stackoverflow.com/questions/16992713/translate-every-element-in-numpy-array-according-to-key
    """
    return np.vectorize(my_dict.__getitem__)(a)


def parameterize_basepairs(gs, ts):
    """
    return list of base pair identities for a guide sequence and a target sequence
    Expects to receive sequences in guide-paired orientation, i.e.
    guide:    5' CA...CT 3'
    target:   3' GT...GA 3'
    Returns the identity of t1G and then the remaining basepairs for g2:t2 through gn:tn
    """
    l = len(gs)
    if l != len(ts):
        print("Error: lengths must be equal!")
        return []
    params = ["t1{}".format(ts[0])]
    for g, t in zip(gs[1:],ts[1:]):
        #params.append("{}{}".format(*sorted([g,t])))
        params.append("g{}t{}".format(g,t))
    return params


def parameterize_kmers_exact(gs, ts, kmer_dict, kmer_lst):
    """
    Return list of kmer energies and identity of t1
    kmers are provided for a list of tuples specifying the starting position and kmer length
    """
    l = len(gs)
    if l != len(ts):
        print("Error: lengths must be equal!")
        return []
    params = ["t1{}".format(ts[0])]
    # Get kmer indicees
    energies = []
    for st, length in kmer_lst:
        k1 = gs[st:st+length]
        k2 = ts[st:st+length][::-1]
        kpair = "{}:{}".format(k1, k2)
        energies.append(kmer_dict[kpair])

    return params, energies


def logistic(x, C=1, k=1, a=1, A=1):
    return A + (C - A)/(1+a*np.exp(k*x))


def pbs_association_function(param_arr, pos_arr, struct_vec, kmer_vec, params):
    """
    Calculate the association rate given a list of params and the matching param values
    """
    v = params.valuesdict()
    # Scale internal structure prediction
    scaled_sv = struct_vec * v['struct']
    # look up params in order
    bp_arr = vec_translate(param_arr, v)
    # get positional params
    p_arr = vec_translate(pos_arr, v)
    bp_by_pos = bp_arr*p_arr

    # Dampen effect of kmer energy:
    C = v['C']
    k = v['k']
    a = v['a']
    A = v['A']
    scaled_sk = logistic(kmer_vec, C,k,a,A)

    kon_arr = np.sum(bp_by_pos, axis=1)  + scaled_sv + scaled_sk
    return kon_arr


def pbs_association_objective(params, param_arr, pos_arr, struct_vec, kmer_vec, reg_params, lbda, data):
    """
    Objective function for fitting association rates
    Regularize with ridge regression
    """
    v = params.valuesdict()
    bp_v = vec_translate(reg_params, v)
    regularize = lbda * sum(bp_v**2) # Ridge regression
    resid = data - pbs_association_function(param_arr, pos_arr, struct_vec, kmer_vec, params)
    loss = sum(resid**2) + regularize

    # When using Levenberg-Marquardt (leastsq), you MUST return an array with more elements than variables
    #resid = data - pbs_association_function(param_arr, pos_arr, struct_vec, kmer_vec, params)
    #return resid
    
    return loss


def pos_by_swap_with_kmers(data_df, train_idx, kmer_lst, kmer_dict, lmda=0.5):
    """
    Fit the position by swap with kmer energy model
    """
    # Time fitting
    start = time.time()
    
    flatten = lambda l: [item for sublist in l for item in sublist]
    
    # Parameterize position and swap
    mpm_params = data_df.apply(lambda x: parameterize_basepairs(
        mplcb.rev_comp(x.guide_seq), x.sequence[5:-5][::-1]), axis=1).tolist()
    
    # All possible position and swap params
    all_params = np.array(sorted(list(set(flatten(mpm_params)))))
    
    # Parameterize kmer
    kmer_params = data_df.apply(lambda x: parameterize_kmers_exact(
        mplcb.rev_comp(x.guide_seq), x.sequence[5:-5][::-1], 
        kmer_dict, kmer_lst), axis=1).tolist()
    
    # Unpack kmer parameters
    t1_lst, kmer_lst = zip(*kmer_params)
    t1_arr = np.array(t1_lst).flatten()
    kmer_arr = np.array(kmer_lst).flatten()
    
    # Define param object
    fit_params = Parameters()

    for p in all_params:
        # bp params cannot be positive
        p = str(p)
        if p[0] == 'g':
            # Should we constrain these to be negative?
            fit_params.add(p, value=-0.1, max=0.0)
            #fit_params.add(p, value=-0.1)
        else:
            fit_params.add(p, value=-0.1)
    
    # Constrain base pairs to be 0.0
    #fit_params['AT'].set(value=0.0, vary=False)
    #fit_params['CG'].set(value=0.0, vary=False)
    fit_params['gAtT'].set(value=0.0, vary=False)
    fit_params['gTtA'].set(value=0.0, vary=False)
    fit_params['gGtC'].set(value=0.0, vary=False)
    fit_params['gCtG'].set(value=0.0, vary=False)

    # Constrain t1G to be 0.0
    fit_params['t1G'].set(value=0.0, vary=False)
    
    # Constrain p1 and p2 to be 1.0
    fit_params.add("p1", value=1.0, vary=False)
    fit_params.add("p2", value=1.0, vary=False)
    for i in range(3,17):
        # Should we constrain these to be negative?
        fit_params.add("p{}".format(i), value=1.0, min=0.0)
        #fit_params.add("p{}".format(i), value=1.0)
        
    # Convert param list into 2d array:
    param_arr = np.array(mpm_params)
    # Get matching param positions:
    pos_vec = np.array(['p{}'.format(x) for x in range(1,17)])
    pos_arr = np.tile(pos_vec, (len(mpm_params),1))
    

    # Dampening parameters for kmer energy
    fit_params.add("C", value=0.0, vary=False)
    fit_params.add("k", value=1.0, min=0.0, vary=True)
    fit_params.add("a", value=1.0, vary=True)
    fit_params.add("A", value=-5.0, max=0.0, vary=True)
    
    # Parameters we want to regularize:
    reg_params = list(all_params) + ['p{}'.format(x) for x in range(1,17)]
    
    
    # Get internal structure
    struct_vec = data_df.delta_target_ensemble_37.values
    #struct_vec = data_df.target_ensemble_37.values
    
    fit_params.add("struct", value=1.0)
    
    # Get data
    data = data_df.delta_kon_ttago_wt.values
    
    # Isolate training data
    fit_param_arr = param_arr[train_idx]
    fit_pos_arr = pos_arr[train_idx]
    fit_struct_vec = struct_vec[train_idx]
    fit_kmer_arr = kmer_arr[train_idx]
    fit_data = data[train_idx]
    
    # Perform fit
    result = minimize(pbs_association_objective, fit_params, 
                  args=(fit_param_arr, fit_pos_arr, fit_struct_vec, fit_kmer_arr, 
                        reg_params, lmda,
                        fit_data), method='L-BFGS-B')
    
    #print("Fitting complete: {} minutes".format(round((time.time() - start)/60.0), 3))
    result_dict = {
        'model': result,
        'param_arr': param_arr,
        'pos_arr': pos_arr,
        'struct_vec': struct_vec,
        'kmer_arr': kmer_arr,
        'data': data
    }
    return result_dict




if __name__ == '__main__':
    main()
