"""
Plotting functions and other utilities for having and easier time with
matplotlib and seaborn
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as mplcm
import seaborn as sns
import pandas as pd
import numpy as np
import scipy as sp
import re

from scipy.stats import gaussian_kde
from Bio import pairwise2
from collections import OrderedDict


# Custom colormaps

# Solar extra (aka the Jason classic)
solar_extra_colors = ['#3361A5', '#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC', '#EAD397', '#FDB31A','#E42A2A','#A31D1D']
solar_extra = mpl.colors.LinearSegmentedColormap.from_list('CustomMap', solar_extra_colors)
solar_extra_r = mpl.colors.LinearSegmentedColormap.from_list('CustomMap', solar_extra_colors[::-1])

# 'Sunrise' (aka chromatin viridis)
sunrise_colors = ['#352A86', '#343DAE', '#0262E0', '#1389D2', '#2DB7A3', '#A5BE6A', '#F8BA43','#F6DA23','#F8FA0D']
sunrise = mpl.colors.LinearSegmentedColormap.from_list('CustomMap', sunrise_colors)
sunrise_r = mpl.colors.LinearSegmentedColormap.from_list('CustomMap', sunrise_colors[::-1])

# Samba Night
samba_colors = ['#1873CC', '#1798E5', '#00BFFF', '#4AC596', '#00CC00', '#A2E700', '#FFFF00','#FFD200','#FFA500']
samba_night = mpl.colors.LinearSegmentedColormap.from_list('CustomMap', samba_colors)
samba_night_r = mpl.colors.LinearSegmentedColormap.from_list('CustomMap', samba_colors[::-1])

# Fireworks
fireworks_colors = ['#ffffff', '#2488F0', '#7F3F98', '#E22929', '#FCB31A']
fireworks = mpl.colors.LinearSegmentedColormap.from_list('CustomMap', fireworks_colors)
fireworks_r = mpl.colors.LinearSegmentedColormap.from_list('CustomMap', fireworks_colors[::-1])




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
Generic plotting

These plotting functions are more general for exploring new data
"""

def show_colormaps(category=None, cmap_to_show=None):
    """
    ==================
    Colormap reference
    ==================

    Reference for colormaps included with Matplotlib.

    This reference example shows all colormaps included with Matplotlib. Note that
    any colormap listed here can be reversed by appending "_r" (e.g., "pink_r").
    These colormaps are divided into the following categories:

    Sequential:
        These colormaps are approximately monochromatic colormaps varying smoothly
        between two color tones---usually from low saturation (e.g. white) to high
        saturation (e.g. a bright blue). Sequential colormaps are ideal for
        representing most scientific data since they show a clear progression from
        low-to-high values.

    Diverging:
        These colormaps have a median value (usually light in color) and vary
        smoothly to two different color tones at high and low values. Diverging
        colormaps are ideal when your data has a median value that is significant
        (e.g.  0, such that positive and negative values are represented by
        different colors of the colormap).

    Qualitative:
        These colormaps vary rapidly in color. Qualitative colormaps are useful for
        choosing a set of discrete colors. For example::

            color_list = plt.cm.Set3(np.linspace(0, 1, 12))

        gives a list of RGB colors that are good for plotting a series of lines on
        a dark background.

    Miscellaneous:
        Colormaps that don't fit into the categories above.

    """


    # Have colormaps separated into categories:
    # http://matplotlib.org/examples/color/colormaps_reference.html
    cmaps = [('Perceptually Uniform Sequential', [
                'viridis', 'plasma', 'inferno', 'magma']),
             ('Sequential', [
                'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']),
             ('Sequential (2)', [
                'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
                'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
                'hot', 'afmhot', 'gist_heat', 'copper']),
             ('Diverging', [
                'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
                'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']),
             ('Qualitative', [
                'Pastel1', 'Pastel2', 'Paired', 'Accent',
                'Dark2', 'Set1', 'Set2', 'Set3',
                'tab10', 'tab20', 'tab20b', 'tab20c']),
             ('Miscellaneous', [
                'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
                'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'hsv',
                'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar'])]


    nrows = max(len(cmap_list) for cmap_category, cmap_list in cmaps)
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))


    def plot_color_gradients(cmap_category, cmap_list, nrows):
        fig, axes = plt.subplots(nrows=nrows)
        fig.subplots_adjust(top=0.95, bottom=0.01, left=0.2, right=0.99)
        axes[0].set_title(cmap_category + ' colormaps', fontsize=14)

        for ax, name in zip(axes, cmap_list):
            ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(name))
            pos = list(ax.get_position().bounds)
            x_text = pos[0] - 0.01
            y_text = pos[1] + pos[3]/2.
            fig.text(x_text, y_text, name, va='center', ha='right', fontsize=10)

        # Turn off *all* ticks & spines, not just the ones with colormaps.
        for ax in axes:
            ax.set_axis_off()

    if category:
        for cmap_category, cmap_list in cmaps:
            if cmap_category == category:
                plot_color_gradients(cmap_category, cmap_list, nrows)
    else:
        for cmap_category, cmap_list in cmaps:
            plot_color_gradients(cmap_category, cmap_list, nrows)

    plt.show()


def get_colors(cmap_name, num_colors):
    """
    Get a number of evenly spaced colors from a provided colormap
    Inputs:
        cmap (str) -- name of colormap
        num_colors (int) -- number of colors to get
    Outputs:
        colors (list) -- list of RGBA tuples
    """
    cmap = mplcm.get_cmap(cmap_name)
    return [cmap(x) for x in np.linspace(0,1,num_colors)]


def plot_cdf(ax, x, kwargs={}):
    """
    Plot a CDF of array-like data x
    Inputs:
        ax (matplotlib axes class) -- The Axes instance on which to plot
        x (array-like) -- values from which to generate CDF
    Output:
        ax (matplotlib axes class) -- returns the modified axes object
    """
    cleaned_x = [i for i in x if np.isfinite(i)]
    sorted_x = sorted(cleaned_x)
    cumulative_fraction = np.linspace(0,1,len(sorted_x))
    ax.plot(sorted_x, cumulative_fraction, **kwargs)
    ax.set_ylim(0, 1.05)
    return ax


def plot_histogram(ax, x, bins=50, log=False, remove_outliers=False, kwargs={}):
    """
    Plot a histogram of values provided values
    Inputs:
        ax (matplotlib axes class) -- The Axes instance on which to plot
        x (array-like) -- values from which to plot histogram
        bins (int) -- number of bins for histogram
        log (bool) -- should the histogram be plotted in log space? (log10)
        remove_outliers (bool or numeric) -- Should outliers be removed? (uses
            a sigma of 3 by default, but other values may be provided)
    Outputs:
        ax (matplotlib axes class) -- The modified axes object
    """
    cleaned_x = [i for i in x if np.isfinite(i)]
    if remove_outliers:
        # Remove outliers by sigma cutoff
        if isinstance(remove_outliers, int) or isinstance(remove_outliers, float):
            cutoff = remove_outliers
        else:
            cutoff = 3.0
        if log:
            mean = np.mean(np.log10(cleaned_x))
            std = np.std(np.log10(cleaned_x))
        else:
            mean = np.mean(cleaned_x)
            std = np.std(np.log10(cleaned_x))
        
        cut_high = mean + cutoff*std
        cut_low = mean - cutoff*std
        cleaned_x = [i for i in x if (i < cut_high) & (i > cut_low)]
    
    if log:
        ax.hist(x=cleaned_x, bins=np.logspace(np.log10(min(cleaned_x)), np.log10(max(cleaned_x)), bins), **kwargs)
        ax.set_xscale("log")
    else:
        ax.hist(x=cleaned_x, bins=bins, **kwargs)
    
    return ax


def reg_scatter(ax, x, y, 
    showR2=False, showRMSE=False, square=False,
    eqLine=None, linFit=None, splineFit=None, label_loc=None, kwargs={}):
    """
    Create a linear scatter plot with points colored by gaussian density
    Inputs:
        ax (matplotlib axes class) -- The Axes instance on which to plot
        x (array-like) -- values to plot on x axis
        y (array-like) -- values to plot on y axis (must be same length as x)
        square (bool) -- should axis be forced to be equal?
        eqLine (Dict/bool) -- plotting parameters for a one to one line (i.e. x=y line)
        linFit (Dict/bool) -- plotting parameters for a linear regression line
        splineFit (Dict/bool) -- plotting parameters for a smoothing spline
        kwargs (Dict) -- keyword arguments to be passed to the scatter plot function
    Output:
        ax (matplotlib axes class) -- returns the modified axes object
    """
    # Adjust line parameters if necessary
    default_line_params = {"ls": "--", "c": "black", "linewidth": 3}
    eq_params = default_line_params.copy()
    lin_params = default_line_params.copy()
    spl_params = default_line_params.copy()
    for curr_params, new_params in zip([eq_params, lin_params, spl_params],[eqLine, linFit, splineFit]):
        if type(new_params) is dict:
            for key, val in new_params.items():
                curr_params[key] = val
            
    
    # Check data first for correct data
    if len(x) != len(y):
        raise ValueError("x and y are not of same length")
    
    # Scrub missing values
    nonans = [xy for xy in zip(x, y) if all(np.isfinite(xy))]
    x, y = [np.array(i) for i in zip(*nonans)]

    # Get pearson correlation
    xy_pcorr = sp.stats.pearsonr(x, y)
    
    # default R2
    r2 = xy_pcorr[0]**2

    # default RMSE:
    xy_rmse = np.sqrt(np.mean((x - y)**2))
    
    # Make plot:
    ax.scatter(x, y, **kwargs)
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    # Optional plot additions
    if square:
        xmin, ymin = [min([xmin, ymin])]*2
        xmax, ymax = [max([xmax, ymax])]*2
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        # make x ticks and y ticks equal
        xticks = ax.get_xticks()
        ax.set_yticks(xticks[1:-1])
        ax.set_aspect(1)
    
    if eqLine:
        # draw a diagonal line:
        ax.plot([xmin, xmax], [ymin, ymax], ls=eq_params["ls"], c=eq_params["c"], linewidth=eq_params["linewidth"])
        # R2 doesn't change
        # xy_rmse doesn't change
    
    if linFit:
        # get linear regression of y on x
        linreg = sp.stats.linregress(x, y)
        plotx = np.array(ax.get_xlim())
        ploty = plotx*linreg[0] + linreg[1]
        newy = x*linreg[0] + linreg[1]
        ax.plot(plotx, ploty, ls=lin_params["ls"], c=lin_params["c"], linewidth=lin_params["linewidth"])
        r2 = linreg[2]**2
        xy_rmse = np.sqrt(np.mean((y - newy)**2))
    
    if splineFit:
        # get natural cubic spline
        # (spline function requires sorted x values)
        x, y = zip(*sorted(zip(x,y), key=lambda tup: tup[0]))
        # Need to also handle duplicate x-values?
        prev_x = None
        new_x = []
        new_y = []
        for x, y in zip(x,y):
            if x == prev_x:
                continue
            else:
                new_x.append(x)
                new_y.append(y)
                prev_x = x
                
        x = np.array(new_x)
        y = np.array(new_y)

        cs = sp.interpolate.UnivariateSpline(x, y)
        xs = np.arange(min(x), max(x), (max(x)-min(x))/100)
        newy = cs(xs)
        ax.plot(xs, newy, ls=spl_params["ls"], c=spl_params["c"], linewidth=spl_params["linewidth"])
        # R^2 = 1 - RSS/TSS
        RSS = sum((y - cs(x))**2)
        TSS = sum((y - np.mean(y))**2)
        r2 = 1 - RSS/TSS
        xy_rmse = np.sqrt(np.mean((y - cs(x))**2))

    label_txt = None
    if showR2 and showRMSE:
        label_txt = "$R^2 = {:0.3f} $\n$ RMSE = {:0.3f}$".format(r2, xy_rmse)   
    elif showR2:
        label_txt = "$R^2 = {:0.3f}$".format(r2)
    elif showRMSE:
        label_txt = "$RMSE = {:0.3f}$".format(xy_rmse)

    if label_txt:
        if label_loc=="bottom_right":
            ax.text(0.95, 0.05, label_txt, transform=ax.transAxes, 
                verticalalignment='bottom', horizontalalignment='right', 
                bbox={'facecolor': ax.get_facecolor(), 'alpha': 1.0, 'pad': 10, 'edgecolor':'none'})
        else:
            ax.text(0.05, 0.95, label_txt, transform=ax.transAxes, 
                verticalalignment='top', horizontalalignment='left', 
                bbox={'facecolor': ax.get_facecolor(), 'alpha': 1.0, 'pad': 10, 'edgecolor':'none'})

    
    return ax


def density_scatter(ax, x, y, 
    showR2=False, showRMSE=False, square=False,
    eqLine=None, linFit=None, splineFit=None, cmap="viridis", kwargs={}):
    """
    Create a linear scatter plot with points colored by gaussian density
    Inputs:
        ax (matplotlib axes class) -- The Axes instance on which to plot
        x (array-like) -- values to plot on x axis
        y (array-like) -- values to plot on y axis (must be same length as x)
        square (bool) -- should axis be forced to be equal?
        eqLine (Dict/bool) -- plotting parameters for a one to one line (i.e. x=y line)
        linFit (Dict/bool) -- plotting parameters for a linear regression line
        splineFit (Dict/bool) -- plotting parameters for a smoothing spline
        kwargs (Dict) -- keyword arguments to be passed to the scatter plot function
    Output:
        ax (matplotlib axes class) -- returns the modified axes object
    """
    # Adjust line parameters if necessary
    default_line_params = {"ls": "--", "c": "black", "linewidth": 3}
    eq_params = default_line_params.copy()
    lin_params = default_line_params.copy()
    spl_params = default_line_params.copy()
    for curr_params, new_params in zip([eq_params, lin_params, spl_params],[eqLine, linFit, splineFit]):
        if type(new_params) is dict:
            for key, val in new_params.items():
                curr_params[key] = val
            
    
    # Check data first for correct data
    if len(x) != len(y):
        raise ValueError("x and y are not of same length")
    
    # Scrub missing values
    nonans = [xy for xy in zip(x, y) if all(np.isfinite(xy))]
    x, y = [np.array(i) for i in zip(*nonans)]

    # Calculate point density
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    # Note that gaussian_kde returns a kernel function-- 
    # calling it as above evaluates that kernel function (the estimated pdf) for a set of provided points
    
    # sort points by density
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    # Get pearson correlation
    xy_pcorr = sp.stats.pearsonr(x, y)
    
    # default R2
    r2 = xy_pcorr[0]**2

    # default RMSE:
    xy_rmse = np.sqrt(np.mean((x - y)**2))
    
    # Make plot:
    ax.scatter(x, y, c=z, cmap=cmap, edgecolors='', **kwargs)
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    # Optional plot additions
    if square:
        xmin, ymin = [min([xmin, ymin])]*2
        xmax, ymax = [max([xmax, ymax])]*2
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect(1)
    
    if eqLine:
        # draw a diagonal line:
        ax.plot([xmin, xmax], [ymin, ymax], ls=eq_params["ls"], c=eq_params["c"], linewidth=eq_params["linewidth"])
        # R2 doesn't change
        # xy_rmse doesn't change
    
    if linFit:
        # get linear regression of y on x
        linreg = sp.stats.linregress(x, y)
        plotx = np.array(ax.get_xlim())
        ploty = plotx*linreg[0] + linreg[1]
        newy = x*linreg[0] + linreg[1]
        ax.plot(plotx, ploty, ls=lin_params["ls"], c=lin_params["c"], linewidth=lin_params["linewidth"])
        r2 = linreg[2]**2
        xy_rmse = np.sqrt(np.mean((y - newy)**2))
    
    if splineFit:
        # get natural cubic spline
        # (spline function requires sorted x values)
        x, y = zip(*sorted(zip(x,y), key=lambda tup: tup[0]))
        # Need to also handle duplicate x-values?
        prev_x = None
        new_x = []
        new_y = []
        for x, y in zip(x,y):
            if x == prev_x:
                continue
            else:
                new_x.append(x)
                new_y.append(y)
                prev_x = x
                
        x = np.array(new_x)
        y = np.array(new_y)

        cs = sp.interpolate.UnivariateSpline(x, y)
        xs = np.arange(min(x), max(x), (max(x)-min(x))/100)
        newy = cs(xs)
        ax.plot(xs, newy, ls=spl_params["ls"], c=spl_params["c"], linewidth=spl_params["linewidth"])
        # R^2 = 1 - RSS/TSS
        RSS = sum((y - cs(x))**2)
        TSS = sum((y - np.mean(y))**2)
        r2 = 1 - RSS/TSS
        xy_rmse = np.sqrt(np.mean((y - cs(x))**2))

    if showR2 and showRMSE:
        label_txt = "$R^2 = {:0.3f} $\n$ RMSE = {:0.3f}$".format(r2, xy_rmse)
        ax.text(0.7, 0.15, label_txt, transform=ax.transAxes, verticalalignment='top', fontsize=14, 
                bbox={'facecolor': ax.get_facecolor(), 'alpha': 1.0, 'pad': 10, 'edgecolor':'none'})
    elif showR2:
        label_txt = "$R^2 = {:0.3f}$".format(r2)
        ax.text(0.05, 0.95, label_txt, transform=ax.transAxes, verticalalignment='top', fontsize=12, 
                bbox={'facecolor': ax.get_facecolor(), 'alpha': 1.0, 'pad': 10, 'edgecolor':'none'})
    elif showRMSE:
        label_txt = "$RMSE = {:0.3f}$".format(xy_rmse)
        ax.text(0.05, 0.95, label_txt, transform=ax.transAxes, verticalalignment='top', fontsize=12, 
                bbox={'facecolor': ax.get_facecolor(), 'alpha': 1.0, 'pad': 10, 'edgecolor':'none'})
    
    return ax


"""
Ago project plotting

These plotting functions are more specific to the Ago projects Winston and I 
have been working on. Be warned that some of them may have functionality that
is fairly specific to those kinds of data.
"""

### Annotation Corrections ###
def get_inverted_annotations(annot_series, wt_seq, annot_col="annotation"):
    """
    Wrapper to get an inverted annotation series from an uninverted annotation series.
    (Right now only works for point-mutant annotation styles)
    """
    wt_base_labels = list(wt_seq)
    return annot_series.apply(lambda x: invert_annot_numbering(x, len(wt_base_labels)))


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

######################################
### Double Mutant Heatmap Plotting ###
######################################

def parse_point_mut_annot(annotation, out='swaps', get_flank=False):
    """
    Pull out the individual mutations from a point mutant annotation.
    Return these point mutations as a list.
    """
    flank, swaps = annotation.split('_')
    individual_swaps = swaps.split(',')
    if out=='swaps':
        if get_flank:
            return flank, individual_swaps
        return individual_swaps
    if out=='pos':
        positions = []
        for swap in individual_swaps:
            num, text = re.split('(\d+)', swap)[1:]
            positions.append(int(num))
        if get_flank:
            return flank, sorted(positions)
        return sorted(positions)


def construct_double_mut_df(dm_df, value, annot_col="annotation"):
    """
    Construct a square data frame for a set of double mutants
    Inputs:
        dm_df (DataFrame) -- A dataframe of double mutants to be formatted into
            the square dataframe. Must contain a column of 'Ago-style' double
            mutant annotations.
        value (str) -- the column name indicating which value should be put in
            the square double mutant dataframe.
        annot_col (str) -- the annotation column name
    Outputs:
        new_frame (DataFrame) -- a square dataframe containing the indicated
            value at double mutant coordinate positions. This dataframe will
            be symmetric across the BL-TR diagonal.
    """
    
    # split mutations annotations and get value of interest:
    list_data = dm_df.apply(lambda x: \
                            parse_point_mut_annot(x[annot_col])  + [x[value]], 
                            axis=1).tolist()
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

    df_data.apply(lambda x: fill_in_dm_frame(x, new_frame), axis=1)      
    
    return new_frame


def add_singles_to_dm_df(sm_df, dm_df, value, annot_col="annotation"):
    """
    Add single mutant values to the diagonal of a prepared double mutant
    heatmap
    Inputs:
        sm_df (DataFrame) -- A dataframe of single mutants
        dm_df (DataFrame) -- A square dataframe of previously formated double mutants
        value (str) -- column name indicating which value should be put in 
            the dataframe
        annot_col (str) -- the annotation column name
    Outputs:
        new_frame (DataFrame) -- modified dm dataframe with single mutants
            in diagonal
    """
    # split mutations annotations and get value of interest:
    list_data = sm_df.apply(lambda x: parse_point_mut_annot(x[annot_col])  + [x[value]], 
                            axis=1).tolist()
    
    sm_data = pd.DataFrame(list_data)
    sm_data.columns = ['first_mut', 'value']
    #print sm_data.head()
    
    # Put single mutants in dm df
    new_frame = dm_df.copy()
    
    for ix, row in sm_data.iterrows():
        pos = row['first_mut']
        val = row['value']
        new_frame.loc[pos,pos] = val
    
    return new_frame



def construct_position_median_dm_df(dm_df, value, annot_col="mut_annotation"):
    """
    Construct a square data frame of position-median double mutant values 
    (i.e. median of 9 for a single dataset)
    Inputs:
        dm_df (DataFrame) -- A dataframe of double mutants to be formatted into
            the square dataframe. Must contain a column of 'Ago-style' double
            mutant annotations.
        value (str) -- the column name indicating which value should be put in
            the square double mutant dataframe.
        annot_col (str) -- the annotation column name
    Outputs:
        new_frame (DataFrame) -- a square dataframe containing the indicated
            value at double mutant coordinate positions. This dataframe will
            be symmetric across the BL-TR diagonal.
    """
    
    # split mutations annotations and get value of interest:
    list_data = dm_df.apply(lambda x: \
                            parse_point_mut_annot(x[annot_col], out='pos')  + [x[value]], 
                            axis=1).tolist()
    df_data = pd.DataFrame(list_data)
    df_data.columns = ['first_mut', 'second_mut', 'value']
    
    # Get a sorted list of all possible base swaps:
    keys = sorted(list(set([x[0] for x in list_data] + [x[1] for x in list_data])))

    # Fill in square dataframe with desired double mutant values
    new_frame = pd.DataFrame(index=keys[::-1], columns=keys)

    for i in range(1, len(keys)+1,1):
        for j in range(1, i, 1):
            if i == j:
                continue
            med_dat = df_data[(df_data.first_mut == j) & (df_data.second_mut == i)]['value'].median()
            new_frame.loc[i,j] = med_dat
            new_frame.loc[j,i] = med_dat
    
    return new_frame


def combine_double_mut_dfs(upper_df, lower_df):
    """
    If you have two square dataframes, this function lets you merge them
    into one dataframe (across a BL to TR diagonal)
    Maintains the indices of the 'upper_df' and the column names of the 'lower_df'
    """
    nrows = len(upper_df.index)
    ncols = len(upper_df.columns)
    combined_df = upper_df.copy()
    
    for i in range(nrows):
        for j in range(ncols):
            if nrows - i <= j:
                combined_df.iloc[i,j] = lower_df.iloc[i,j]
    
    combined_df.columns = lower_df.columns
    
    return combined_df


def double_mutant_heatmap(ax, dm_df, mask_color="lightgrey", cmap="viridis", cbar_label="", 
    splitCells=None, splitColor="lightgrey", slim_ticks=True, kwargs={}):
    """
    Plot a double mutant heatmap.
    Inputs:
        ax (matplotlib axes class) -- The Axes instance on which to plot
        dm_df (DataFrame) -- The rectangular dataframe containing data to plot
        cbar_label (str) -- The label for the colorbar (color bar limits can be passed to kwargs under 'vmax'/'vmin')
        splitCells (tuple of ints) -- should divisions between groups of cells be highlighted? Send a 2-member tuple
            where the first int is the split frequency for rows and the second int is the split frequency for columns
        splitColor (str) -- The color to be used between cells that are split
    Output:
        hm (matplotlib axes class) -- returns the modified axes object
    """
    # cast all data to float (to prevent sns.heatmap from being a baby)
    plot_df = dm_df.astype(float)
    
    # this is how you change the 'mask' color
    ax.set_facecolor(mask_color)
    
    hm = sns.heatmap(plot_df, ax=ax, cmap=cmap, square=True, cbar_kws={'label':cbar_label}, **kwargs)
    
    if splitCells:
        ax.hlines([splitCells[0]*i for i in range(len(dm_df.index))], *ax.get_xlim(), color=splitColor, linewidth=1)
        ax.vlines([splitCells[1]*i for i in range(len(dm_df.columns))], *ax.get_ylim(), color=splitColor, linewidth=1)

    if slim_ticks:
        xlabels = [item.get_text().split(">")[-1] for item in ax.get_xticklabels()]
        ylabels = [item.get_text().split(">")[-1] for item in ax.get_yticklabels()]
        ax.set_xticklabels(xlabels)
        ax.set_yticklabels(ylabels)
        xtick_pos = ax.get_xticks()
        ytick_pos = ax.get_yticks()
        x_pos_step = xtick_pos[1] - xtick_pos[0]
        y_pos_step = ytick_pos[1] - ytick_pos[0]
        new_xpos = []
        new_ypos = []
        for i, pos in enumerate(xtick_pos):
            if i % 3 == 0:
                new_xpos.append(pos + x_pos_step*0.1)
            elif i % 3 == 2:
                new_xpos.append(pos - x_pos_step*0.1)
            else:
                new_xpos.append(pos)
        for i, pos in enumerate(ytick_pos):
            if i % 3 == 0:
                new_ypos.append(pos + y_pos_step*0.1)
            elif i % 3 == 2:
                new_ypos.append(pos - y_pos_step*0.1)
            else:
                new_ypos.append(pos)
        ax.set_xticks(new_xpos)
        ax.set_yticks(new_ypos)
        plt.xticks(rotation=0)
    
    return hm



def single_mutant_heatmap(ax, sm_df, mask_color="lightgrey", cmap="viridis", cbar_label="", 
    splitCells=None, splitColor="lightgrey", slim_ticks=True, kwargs={}):
    """
    Plot a single mutant heatmap.
    Inputs:
        ax (matplotlib axes class) -- The Axes instance on which to plot
        sm_df (DataFrame) -- The dataframe containing data to plot
        cbar_label (str) -- The label for the colorbar (color bar limits can be passed to kwargs under 'vmax'/'vmin')
        splitCells (tuple of ints) -- should divisions between groups of cells be highlighted? Send a 2-member tuple
            where the first int is the split frequency for rows and the second int is the split frequency for columns
        splitColor (str) -- The color to be used between cells that are split
    Output:
        hm (matplotlib axes class) -- returns the modified axes object
    """
    # cast all data to float (to prevent sns.heatmap from being a baby)
    plot_df = sm_df.astype(float)
    
    # this is how you change the 'mask' color
    ax.set_facecolor(mask_color)
    
    hm = sns.heatmap(plot_df, ax=ax, cmap=cmap, cbar_kws={'label':cbar_label}, **kwargs)
    
    if splitCells:
        ax.hlines([splitCells[0]*i for i in range(len(sm_df.index))], *ax.get_xlim(), color=splitColor, linewidth=1)
        ax.vlines([splitCells[1]*i for i in range(len(sm_df.columns))], *ax.get_ylim(), color=splitColor, linewidth=1)

    if slim_ticks:
        xlabels = [item.get_text().split(">")[-1] for item in ax.get_xticklabels()]
        ax.set_xticklabels(xlabels)
        xtick_pos = ax.get_xticks()
        x_pos_step = xtick_pos[1] - xtick_pos[0]
        new_xpos = []
        for i, pos in enumerate(xtick_pos):
            if i % 3 == 0:
                new_xpos.append(pos + x_pos_step*0.1)
            elif i % 3 == 2:
                new_xpos.append(pos - x_pos_step*0.1)
            else:
                new_xpos.append(pos)

        ax.set_xticks(new_xpos)
        plt.xticks(rotation=0)
    
    return hm



def construct_single_mut_df(sm_df, value, annot_col="annotation"):
    """
    Construct a data frame for a set of single mutants, potentially
    with multiple base contexts
    Inputs:
        sm_df (DataFrame) -- A dataframe of single mutants to be formatted into
            the (context) dataframe. Must contain a column of 'Ago-style' single
            mutant annotations.
        value (str) -- the column name indicating which value should be put in
            the single mutant dataframe.
        annot_col (str) -- the annotation column name
    Outputs:
        new_frame (DataFrame) -- a context x mutant dataframe containing the indicated
            value at each context and mutant position.
    """
    
    # split mutations annotations and get value of interest:
    flanks_and_swaps = sm_df.apply(lambda x: \
                            parse_point_mut_annot(x[annot_col], get_flank=True), 
                            axis=1).tolist()
    flanks, swaps = zip(*flanks_and_swaps)
    # Reformat
    flanks = list(flanks)
    swaps = [x[0] for x in swaps]
    vals = sm_df[value].tolist()

    df_data = pd.DataFrame([flanks, swaps, vals]).T
    df_data.columns = ['flank', 'swap', 'value']
    
    # Get a sorted list of all possible base swaps:
    keys = sorted(list(set(swaps)), 
        key=natural_keys)

    # Fill in square dataframe with desired double mutant values
    new_frame = pd.DataFrame(index=list(set(flanks)), columns=keys)

    def fill_in_sm_frame(row, df_to_fill):
        # Apply function for filling in a square dataframe
        flank = row['flank']
        swap = row['swap']
        value = row['value']
        df_to_fill.loc[flank, swap] = value

    df_data.apply(lambda x: fill_in_sm_frame(x, new_frame), axis=1)

    
    return new_frame



##########################################
### Cumulative Mutant Heatmap Plotting ###
##########################################

def is_cumulative_mut(row, wt_seq, flank = "AAAAA",
                      cm_dict={"A": "T", "T": "A", "G": "C", "C": "G"}):
    """
    Manually determine if sequence is cumulative mutant
    """
    if row.sequence[:5] != flank or row.sequence[-5:] != flank:
        # Make sure flank is correct
        return "NA"
    seq = row.sequence[5:-5]
    if len(seq) != len(wt_seq):
        # If seq lengths are not equal, not a cumulative mutant
        return "NA"
    start = 1
    end = 1
    in_mut = False
    found_end = False
    for i, (a, b) in enumerate(zip(seq, wt_seq)):
        if a == b:
            if in_mut:
                # End of mutant stretch, mark end
                end = i
                found_end = True
                in_mut = False
        elif a == cm_dict[b]:
            if found_end:
                return "NA" # Already had one mutant stretch. Not continuous
            if not in_mut:
                start = i + 1
                in_mut = True
        else:
            return "NA" # Not a complement mismatch
    # If got all the way to end of sequence, set end
    if in_mut and not found_end:
        end = len(seq)
    if not in_mut and not found_end:
        # no mutations found
        return "NA"
    return "_".join(["polyA", "comp", str(start), "to", str(end)])


def parse_cumulative_mut_annot(annotation, get_flank=False):
    """
    Pull out the beginning and end point of cumulative mutant.
    Return these as a list.
    """
    flank, comp, start, to, end = annotation.split('_')
    return [int(x) for x in [start, end]]

def construct_cumulative_mut_df(cm_df, value, annot_col="annotation", reverse=True):
    """
    Construct a square data frame for cumulative mutants
    Inputs:
        cm_df (DataFrame) -- A dataframe of cumulative mutants to be formatted into
            the square dataframe. Must contain a column of 'Ago-style' double
            mutant annotations.
        value (str) -- the column name indicating which value should be put in
            the square double mutant dataframe.
        annot_col (str) -- the annotation column name
    Outputs:
        new_frame (DataFrame) -- a square dataframe containing the indicated
            value at double mutant coordinate positions. This dataframe will
            be symmetric across the BL-TR diagonal.
    """
    
    # split mutations annotations and get value of interest:
    list_data = cm_df.apply(lambda x: \
                            parse_cumulative_mut_annot(x[annot_col])  + [x[value]], 
                            axis=1).tolist()
    df_data = pd.DataFrame(list_data)
    df_data.columns = ['start_pos', 'end_pos', 'value']
    
    # Get range of positions
    first_pos = min(df_data.start_pos.values)
    last_pos = max(df_data.end_pos.values)
    all_pos = range(first_pos, last_pos + 1)
    
    # Reverse orientation to be guide-centric
    if reverse:
        new_start = last_pos - df_data.end_pos + 1
        new_end = last_pos - df_data.start_pos + 1
        df_data.start_pos = new_start
        df_data.end_pos = new_end

    # Fill in square dataframe with desired double mutant values
    new_frame = pd.DataFrame(index=all_pos[::-1], columns=all_pos)

    def fill_in_cm_frame(row, df_to_fill):
        # Apply function for filling in a square dataframe
        p1 = row['start_pos']
        p2 = row['end_pos']
        value = row['value']
        df_to_fill.loc[p1, p2] = value
        df_to_fill.loc[p2, p1] = value

    df_data.apply(lambda x: fill_in_cm_frame(x, new_frame), axis=1)      
    
    return new_frame

def cumulative_mutant_heatmap(ax, cm_df, mask_color="lightgrey", cmap="viridis", cbar_label="",  kwargs={}):
    """
    Plot a cumulative mutant heatmap.
    Inputs:
        ax (matplotlib axes class) -- The Axes instance on which to plot
        cm_df (DataFrame) -- The square dataframe containing data to plot
        cbar_label (str) -- The label for the colorbar (color bar limits can be passed to kwargs under 'vmax'/'vmin')
    Output:
        hm (matplotlib axes class) -- returns the modified axes object
    """
    # cast all data to float (to prevent sns.heatmap from being a baby)
    plot_df = cm_df.astype(float)
    
    # this is how you change the 'mask' color
    ax.set_facecolor(mask_color)
    
    hm = sns.heatmap(plot_df, ax=ax, cmap=cmap, square=True, cbar_kws={'label':cbar_label}, **kwargs)
    
    return hm





############################
### Tandem mismatch Plots ##
############################

def make_tandem_mutant_df(mut_df):
    """
    Given a dataframe already filtered to the mutant class you're interested in,
    generate a tandem mutant df
    """
    orig_col_names = mut_df.columns.tolist()
    tandem_mutants = []
    for index, row in mut_df.iterrows():
        mutation_positions = sorted([int(re.findall('\d+', x)[0]) for x in row.mut_annotation.split('_')[-1].split(',')])
        mut_intervals = [mutation_positions[i+1] - mutation_positions[i] for i in range(len(mutation_positions)-1)]
        if all([x == 1 for x in mut_intervals]):
            tandem_mutants.append(row.tolist() + [mutation_positions[0]])
    tandem_df = pd.DataFrame(tandem_mutants)
    tandem_df.columns = orig_col_names + ['mut_pos']
    return tandem_df.sort_values(by='mut_pos')


def get_missing_vals(tandem_df, value):
    """
    Get the missing values for each position from a tandem mutant df
    """
    # Get tandem data
    groups = tandem_df.groupby('mut_pos')
    plot_data = [group[value].tolist() for name, group in groups]
    missing_vals = []
    for pos_list in plot_data:
        missing_vals.append(sum([1 for x in pos_list if np.isnan(x)]))
    
    if 1 not in groups.groups.keys():
        missing_vals = [[np.nan]] + missing_vals
    return missing_vals


def get_out_of_bounds_vals(tandem_df, value, oob, direction='high'):
    """
    Get the missing values for each position from a tandem mutant df
    (direction = {high, low})
    """
    # Get tandem data
    groups = tandem_df.groupby('mut_pos')
    plot_data = [group[value].tolist() for name, group in groups]
    oob_vals = []
    for pos_list in plot_data:
        if direction == 'high':
            oob_vals.append(sum([1 for x in pos_list if x >= oob]))
        elif direction == 'low':
            oob_vals.append(sum([1 for x in pos_list if x <= oob]))
        else:
            print("'direction' must either be 'high' or 'low'!")
            return None
    
    if 1 not in groups.groups.keys():
        oob_vals = [[np.nan]] + oob_vals
    return oob_vals


def plot_tandem_mutant_boxplot(ax, tandem_df, value, wt_rate, num_bases, wt_line_color="black", kwargs={}):
    """
    Plot a tandem mutant boxplot
    """
    # Get tandem data
    groups = tandem_df.groupby('mut_pos')
    plot_data = [group[value].tolist() for name, group in groups]
    cleaned_plot_data =[]
    for pos_list in plot_data:
        cleaned_plot_data.append([x for x in pos_list if not np.isnan(x)])
    
    if 1 not in groups.groups.keys():
        cleaned_plot_data = [[np.nan]] + cleaned_plot_data
    if 1 not in groups.groups.keys():
        plot_data = [[np.nan]] + plot_data

    
    # Make plot
    sns.set_style('ticks')
    
    bp = ax.boxplot(cleaned_plot_data, patch_artist=True, **kwargs)
    """
    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(bp[element], color='black')
    for patch in bp['boxes']:
        patch.set(facecolor='#B39E9C')
    """
    ax.axhline(y=wt_rate, ls='--', color=wt_line_color,alpha=0.6)
        
    ax.tick_params(labelsize=12)
    xlabels = ax.get_xticklabels()

    tick_labels = [",".join([str(int(x.get_text()) + i) for i in range(num_bases)]) for x in xlabels]

    ax.set_xticklabels(tick_labels)
    ax.tick_params(axis='x', rotation=90)
    
    #sns.despine()
    
    return ax





##########################
### Single Indel Plots ###
##########################

def parse_single_insertion(row, wt_seq, val_col):
    flank, p1, p2etc = row.mut_annotation.split("_")
    p2, base = p2etc.split("ins")
    p1 = len(wt_seq) - float(p1) + 1
    p2 = len(wt_seq) - float(p2) + 1
    return pd.Series([(p1 + p2)/2, base, row[val_col]])

def parse_single_deletion(row, wt_seq, val_col):
    no_flank = row.mut_annotation.split("_")[-1]
    pos = int(no_flank.split("del")[0])
    base = wt_seq[pos-1]
    pos = len(wt_seq) - pos + 1
    return pd.Series([pos, base, row[val_col]])


def plot_single_insertions_and_deletions(ax, full_df, val_col, wt_seq, lines=False, spanning_ins=True, log=True,
                                         base_dict=OrderedDict([("A", "#E24E42"), ("G", "#E9B000"), 
                                                               ("C", "#008F95"), ("T", "#30415D"), 
                                                               ("Deletion", "white")]), 
                                         ins_kwargs={}, del_kwargs={}, ins_line_kwargs={}, del_line_kwargs={}):
    """
    plot insertions and deletions together for a given parameter
    """
    single_deletions = full_df[full_df.mutant_group == "single_deletions"]
    single_insertions = full_df[full_df.mutant_group == "single_insertions"]
    
    del_data = single_deletions.apply(lambda x: parse_single_deletion(x, wt_seq, val_col), axis=1)
    del_data.columns = ["pos", "base", "val"]
    ins_data = single_insertions.apply(lambda x: parse_single_insertion(x, wt_seq, val_col), axis=1)
    ins_data.columns = ["pos", "base", "val"]
    
    # Log 
    if log:
        del_data["val"] = np.log(del_data["val"])
        ins_data["val"] = np.log(ins_data["val"])
    
    # Sort positions
    del_data.sort_values("pos", inplace=True)
    ins_data.sort_values("pos", inplace=True)
    
    # fill in missing deletion data
    positions = del_data["pos"].tolist()
    for i, p in enumerate(positions[1:]):
        if p - positions[i] > 1:
            dup_data = del_data[del_data.pos == positions[i + 1]].copy()
            dup_data["pos"] = positions[i] + 1
            del_data = del_data.append(dup_data, ignore_index=True)
    
    if spanning_ins:
        # fill in missing insertion data (spanning insertions)
        spanners = []
        ins_pos = list(set(ins_data["pos"]))
        bases = set(ins_data["base"])
        for i, pos in enumerate(ins_pos):
            missing_base = bases - set(ins_data[ins_data["pos"] == pos]["base"])
            if missing_base:
                dup_data = ins_data[(ins_data["pos"] == pos + 1) & (ins_data["base"].isin(missing_base))].copy()
                if dup_data.empty:
                    continue
                dup_data["pos"] = pos
                ins_data = ins_data.append(dup_data, ignore_index=True)
                b, v = dup_data.values[0][1:]
                spanners.append([pos, pos+1, b, v])
        
        #for sp in spanners:
        #    ax.plot([sp[0], sp[1]], [sp[3], sp[3]], c=base_dict[sp[2]])
    
    
    # resort
    del_data.sort_values("pos", inplace=True)
    ins_data.sort_values("pos", inplace=True)
    
    # do the plotting
    
    del_pos = del_data["pos"].tolist()
    del_y = del_data["val"].tolist()
    del_bases = del_data["base"].tolist()
    
    ins_pos = ins_data["pos"].tolist()
    ins_y = ins_data["val"].tolist()
    ins_bases = ins_data["base"].tolist()
    ins_colors = [base_dict[x] for x in ins_bases]
    
    if lines:
        ax.plot(del_pos, del_y, alpha=0.5, **del_line_kwargs)
        ins_mean = ins_data.groupby("pos").mean()
        ax.plot(ins_mean.index, ins_mean, alpha=0.5, **ins_line_kwargs)

    ax.scatter(del_pos, del_y, c=base_dict["Deletion"], edgecolors="black", 
               facecolors=base_dict["Deletion"], zorder=10, **del_kwargs)
    
    ax.scatter(ins_pos, ins_y, c=ins_colors, zorder=11, **ins_kwargs)
                
    return ax



##################################
### Homopolymer insertion plot ###
##################################

def parse_ins_mut(annotation, get_flank=False):
    """
    Retrieve the position, size, and identity of an insertion or bulge mutant
    Assumes either single insertion or bulge annotation styles!
    (single ins: polyA_8_9insG)
    (bulge: polyA_PC_10insAAAAA)
    """
    # First determine if single insertion or bulge mutant
    flank, PC, etc = annotation.split('_')
    if PC == 'PC':
        # bulge mutant
        pos, ins_str = etc.split('ins')
        base = ins_str[0]
        size = len(ins_str)
        return [int(pos), size, base]
    elif PC.isdigit():
        pos = int(PC)
        _, base = etc.split('ins')
        return [pos, 1, base]
    # Incorectly formatted annotation
    return [np.nan, np.nan, np.nan]


def construct_ins_df(ins_df, value, annot_col="annotation", reverse=True):
    """
    Construct a homopolymer insertion dataframe for plotting
    """
    list_data = ins_df.apply(lambda x: parse_ins_mut(x[annot_col]) + [x[value]], axis=1).tolist()
    df_data = pd.DataFrame(list_data)
    df_data.columns = ['ins_pos', 'ins_size', 'ins_base', 'value']
    
    if reverse:
        max_len = max(df_data.ins_pos.values)
        new_pos = max_len - df_data.ins_pos + 1
        df_data.ins_pos = new_pos
    
    return df_data


def plot_homopolymer_insertions(ax, ins_df, sScale=5, pScale=0.15, shade_bkg=True,
                                base_dict={'A': "#E24E42", 'G': "#E9B000", 'C': "#008F95", 'T': "#30415D"}):
    """
    Plot homopolymer insertions
    ax = ax used for plotting
    ins_df = DataFrame with insertion information
    sScale = scale factor for point size
    pScale = scale factor for position spread
    shade_bkg = Should background be shaded in alternating colors
    base_dict = colors to be used for each base
    """
    # Get required information in lists
    ins_pos = ins_df.ins_pos.tolist()
    ins_y = ins_df.value.tolist()
    ins_bases = ins_df.ins_base.tolist()
    ins_size = ins_df.ins_size.tolist()
    ins_colors = [base_dict[x] for x in ins_bases]
    # Adjust x positions
    unique_pos = sorted(list(set(ins_pos)))
    unique_bases = list(set(ins_bases))
    nbases = len(unique_bases)
    shift_pos = range(nbases) - np.mean(range(nbases))
    shift_dict = dict(zip(unique_bases, [pScale*x for x in shift_pos]))
    new_pos = []
    for p, b in zip(ins_pos, ins_bases):
        new_pos.append(p + shift_dict[b])
        
    # Now plot:
    ax.scatter(new_pos, ins_y, c=ins_colors, s=[x*sScale for x in ins_size])
    ax.set_xticks(ins_pos)
    ax.set_xticklabels(ins_pos)
    
    # Shade, if necessary
    if shade_bkg:
        for i in unique_pos:
            if i%2:
                plt.axvspan((i) - 0.5, (i) + 0.5, facecolor='grey', alpha=0.1)
    
    ax.set_xlim(0.5, max(unique_pos) + 0.5)
    return ax


