"""
Curve fitting functions and related utilities
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy as sp

from lmfit import Parameters, Model


"""
Most functions here will use the lmfit library:
    "Non-linear Least-squares Minimization and Curve-fitting for Python"

I will be making use of two lmfit classes: 'Parameters' and 'Model'.
The 'Parameters' class to lets you define the parameters of a model and 
optionally provide some information about their allowed behavior during fitting. 
Most models in here will have an associated function to pre-define their parameters.

The 'Model' class uses a 'model function', some provided data, and a 'Parameters'
instance to store information about a curve-fitting problem.
"""


"""
How to add/define parameters for lmfit:
    
add(name, value=None, vary=True, min=-inf, max=inf, expr=None, brute_step=None)

Parameters:
    name (str) -- Name of parameter. Cannot be a Python reserved word.
    value (float, optional) -- Numerical Parameter value, typically the initial value.
    vary (bool, optional) -- Whether the Parameter is varied during a fit (default is True).
    min (float, optional) -- Lower bound for value (default is -numpy.inf, no lower bound).
    max (float, optional) -- Upper bound for value (default is numpy.inf, no upper bound).
    expr (str, optional) -- Mathematical expression used to constrain the value during the fit.
    brute_step (float, optional) -- Step size for grid points in the brute method.
"""


#################################
### General fitting functions ###
#################################

"""
These functions can be used for a variety of objective functions
"""


def fit_curve(x, y, params, objective_function):
    """
    Fit an objective function using the provided params.
    Inputs:
        x (array-like) -- independent variable for fit
        y (array-like) -- observed data to be fit
        params (Parameters) -- lmfit 'Parameters' object matching objective_function
        objective_function (function) -- the model to be fit
    Output:
        (ModelResult) -- a lmfit 'ModelResult' object
    """
    fit_model = Model(objective_function)
    # x is explicitly passed to overwrite the independent parameter x
    return fit_model.fit(y, params, x=x)


################################
### Single Exponential Decay ###
################################


def single_exp_decay_params(fmin=None, fmax=None, kobs=None):

    # Define parameters object
    params = Parameters()
    default_params = {
        "fmin":{"value": 0.0, "vary": True, "min": -np.inf, "max": np.inf}, 
        "fmax":{"value": 1.0, "vary": True, "min": -np.inf, "max": np.inf},
        "kobs":{"value": 0.5, "vary": True, "min": -np.inf, "max": np.inf}
    }
    if fmin:
        for opt, val in fmin.items():
            default_params["fmin"][opt] = val
    if fmax:
        for opt, val in fmax.items():
            default_params["fmax"][opt] = val
    if kobs:
        for opt, val in kobs.items():
            default_params["kobs"][opt] = val

    for p, dct in default_params.items():
        params.add(p, value=dct["value"], vary=dct["vary"], min=dct["min"], max=dct["max"])
    return params


def single_exp_decay(x, fmin, fmax, kobs):
    """
    Single exponential decay function. 
    Functions provided to lmfit Model must have dependent variable as first argument.
    """
    return fmin + (fmax - fmin)*np.exp(-kobs*x)


def plot_single_exp_data_and_curve(ax, x, y, fit_model, showParams=True, showR2=True):
    """
    Plot a fit curve over some provided data
    """
    x = np.array(x)
    y = np.array(y)
    ax.scatter(x=x, y=y, c="red", s=50)
    model_x = np.linspace(min(x), max(x), 50)
    model_y = fit_model.eval(x=model_x)
    ax.plot(model_x, model_y, c="black", linestyle="--")
    
    # Show parameters if requested
    params_dict = fit_model.params.valuesdict()
    ss_total = np.sum((y - y.mean())**2)
    ss_resid = np.sum((fit_model.residual)**2)
    rsq = 1 - ss_resid/ss_total
    half_life = np.log(2.0)/params_dict["kobs"]
    
    if showParams and showR2:
        if params_dict["kobs"] < 0.001:
            label_txt = "$k_c = {:0.3e}$ $h^{{-1}}$\n$t_{{1/2}} = {:0.3f}$ $h$\n$R^2 = {:0.3f}$".format(
                params_dict["kobs"], half_life, rsq)
        else:
            label_txt = "$k_c = {:0.3f}$ $h^{{-1}}$\n$t_{{1/2}} = {:0.3f}$ $h$\n$R^2 = {:0.3f}$".format(
                params_dict["kobs"], half_life, rsq)
        ax.text(0.95, 0.95, label_txt, transform=ax.transAxes, 
                verticalalignment='top', horizontalalignment='right', fontsize=12, 
                bbox={'facecolor': ax.get_facecolor(), 'alpha': 1.0, 'pad': 10, 'edgecolor':'none'})
    return ax



###################################
### Single Exponentail Function ###
###################################


def exponential_params(a=None, b=None):
    """
    Simple exponential equation parameters
    """

    # Define parameters object

    params = Parameters()
    default_params = {
        "a":{"value": 0.5, "vary": True, "min": -np.inf, "max": np.inf}, 
        "b":{"value": 1.0, "vary": True, "min": -np.inf, "max": np.inf}
    }
    if a:
        for opt, val in a.items():
            default_params["a"][opt] = val
    if b:
        for opt, val in b.items():
            default_params["b"][opt] = val
    
    for p, dct in default_params.items():
        params.add(p, value=dct["value"], vary=dct["vary"], min=dct["min"], max=dct["max"])
    return params


def exp_function(x, a, b):
    """
    Basic exponential function 
    Functions provided to lmfit Model must have dependent variable as first argument.
    """
    return a*b**x


def plot_exp_data_and_curve(ax, x, y, fit_model, showParams=True, showR2=True):
    """
    Plot a fit curve over some provided data
    """
    x = np.array(x)
    y = np.array(y)
    ax.scatter(x=x, y=y, c="red", s=50)
    model_x = np.linspace(min(x), max(x), 50)
    model_y = fit_model.eval(x=model_x)
    ax.plot(model_x, model_y, c="black", linestyle="--")
    
    # Show parameters if requested
    params_dict = fit_model.params.valuesdict()
    ss_total = np.sum((y - y.mean())**2)
    ss_resid = np.sum((fit_model.residual)**2)
    rsq = 1 - ss_resid/ss_total
    
    if showParams and showR2:
        label_txt = "$k_c = ab^T$\n$a = {:0.3e}$\n$b = {:0.3f}$\n$R^2 = {:0.3f}$".format(
            params_dict["a"], params_dict["b"], rsq)
        ax.text(0.05, 0.95, label_txt, transform=ax.transAxes, 
                verticalalignment='top', horizontalalignment='left', fontsize=12, 
                bbox={'facecolor': ax.get_facecolor(), 'alpha': 1.0, 'pad': 10, 'edgecolor':'none'})
    return ax


##########################
### Arrhenius equation ###
##########################

def arrhenius_params(A=None, Ea=None):
    """
    Parameters for the Arrhenius equation
    """

    # Define parameters object
    params = Parameters()
    default_params = {
        "A":{"value": 0.1, "vary": True, "min": -np.inf, "max": np.inf}, 
        "Ea":{"value": 1.0, "vary": True, "min": -np.inf, "max": np.inf}
    }
    if A:
        for opt, val in a.items():
            default_params["A"][opt] = val
    if Ea:
        for opt, val in b.items():
            default_params["Ea"][opt] = val
    
    for p, dct in default_params.items():
        params.add(p, value=dct["value"], vary=dct["vary"], min=dct["min"], max=dct["max"])
    return params


def arrhenius_equation(x, A, Ea):
    """
    The Arrhenius equation relating temperature to reaction rate. x is temp. 
    Functions provided to lmfit Model must have dependent variable as first argument.
    """
    #R = 1.987e-3 # kcal/K*mol
    R = 8.314e-3 # kJ / K*mol
    return A*np.exp(-Ea / (R*(273 + x)))


def plot_arrhenius_data_and_curve(ax, x, y, fit_model, showParams=True, showR2=True):
    """
    Plot a fit curve over some provided data
    """
    x = np.array(x)
    y = np.array(y)
    ax.scatter(x=x, y=y, c="red", s=50)
    model_x = np.linspace(min(x), max(x), 50)
    model_y = fit_model.eval(x=model_x)
    ax.plot(model_x, model_y, c="black", linestyle="--")
    
    # Show parameters if requested
    params_dict = fit_model.params.valuesdict()
    ss_total = np.sum((y - y.mean())**2)
    ss_resid = np.sum((fit_model.residual)**2)
    rsq = 1 - ss_resid/ss_total
    
    if showParams and showR2:
        #eq_form = "$k_c = Ae^{\\frac{-E_a}{RT}}$"
        label_txt = "$A = {:0.3e}$ $h^{{-1}}$\n$E_a = {:0.3f}$ $kJ/mol$\n$R^2 = {:0.3f}$".format(
            params_dict["A"], params_dict["Ea"], rsq)
        ax.text(0.05, 0.95, label_txt, transform=ax.transAxes, 
                verticalalignment='top', horizontalalignment='left', fontsize=12, 
                bbox={'facecolor': ax.get_facecolor(), 'alpha': 1.0, 'pad': 10, 'edgecolor':'none'})
    return ax


