import json
import glob
import os
import sys 
import subprocess
from tracemalloc import start
import numpy as np
from astropy.io import fits
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import shutil
import math
import time
from scipy import stats
import csv
import re
import astropy.units as u
from astropy.time import Time, TimeDelta
import multiprocessing as mp
from xspec import Xset
import signal


from xspec import * 
from get_results_data import get_swift_xrt_counts


import sys
import os
sys.path.append(os.path.abspath("../"))
from input_parameters import *


##########################################################################################


def iso2mjd(iso_dates):
    """
    Convert ISO dates to MJD using astropy Time.
    """
    times = Time(iso_dates, format='isot', scale='utc')
    return times.mjd



def swift_met_to_mjd(met):
    """
    Convert Swift MET to MJD -- approximate way.
    """
    swift_epoch = Time("2001-01-01T00:00:00", scale="tt")
    return (swift_epoch + TimeDelta(met, format="sec")).utc.mjd



def extract_file_info(spectrum_file):
    """
    Get the required information from a spEctral file.
    """

    try:
        with fits.open(spectrum_file) as hdul:
            hdr = hdul[1].header

            bin_counts = hdr.get("COUNTGRP") # min number of counts in the bins
            tot_counts = hdr.get("COUNTS") #... should be same as s.rate[0] * s.exposure
            exposure_sec = hdr.get("EXPOSURE", 0) # exposure length... same as s.exposure
            exp_days = exposure_sec / 86400.0 # 24*60*60
            
            # Extract observation ID
            obs_id = hdr.get("OBS_ID", "UNKNOWN")

            ## Get the middle MJD, WAY 1
            # Extract UTC-based info 
            date_start_str = hdr.get("DATE-OBS") # start date
            date_end_str = hdr.get("DATE-END") # start date
            t_start_1 = Time(date_start_str, format='isot', scale='utc').mjd
            t_end_1   = Time(date_end_str, format='isot', scale='utc').mjd
            mid_mjd_1 = 0.5 * (t_start_1 + t_end_1)
            dt_1 = t_end_1 - t_start_1
            
            
            ## Get the middle MJD, WAY 2
            # Extract MET-based info
            tstart_met = hdr.get("TSTART")
            tstop_met = hdr.get("TSTOP")
            # Convert MET to MJD
            try:
                mjd_start = swift_met_to_mjd(tstart_met)
                mjd_end = swift_met_to_mjd(tstop_met)
                mid_mjd_2 = 0.5*(mjd_start + mjd_end)
                dt_2 = (mjd_end - mjd_start)
            except Exception as e:
                print(f"met conversion error for {spectrum_file}: {e}")
                mid_mjd_2 = None
                dt_2 = None

        
    except Exception as e:
        print(f"Error reading {spectrum_file}: {e}")


    return obs_id, bin_counts, tot_counts, date_start_str, date_end_str, mid_mjd_1, dt_1, mid_mjd_2, dt_2, exposure_sec
    




def f_test(chi1, dof1, chi2, dof2):
    """
    Run a statistical test to compare two models. Number 2 is the new model
    """

    chidiff = chi1 - chi2
    dofdiff = dof1 - dof2

    new_model_preferred = False

    # Deal with the case where the two models have the same number of parameters, and thus the same dof, in which case we can just compare the chi-squared values
    if dofdiff == 0:
        if chi2 < chi1: new_model_preferred = True
        else: new_model_preferred = False

    else: 
        F = (chidiff ) / chi2 * (dof2 / dofdiff)
        p_value = 1 - stats.f.cdf(F, dofdiff, dof2)
        if p_value < 0.001: new_model_preferred = True 
        else: new_model_preferred = False 
     
    if new_model_preferred: 
        #print("New model preferred with p-value: ", p_value)
        return True
    else: 
        #print("Previous model is preferred with p-value: ", p_value)
        return False




##########################################################################################


def plot_resid(spectrum_name, mod_name, mod, save=True, log=True, fit_with_binning=True, setplot_rebin_mincounts= None, setplot_rebin_maxbins = None):
    """
    Helper method to plot the residuals for the fits
    """

    Plot.device = '/null'
    Plot.xAxis = "KeV"  # x-axis for plotting set to energy instead of channel

    if setplot_rebin_mincounts is not None and setplot_rebin_maxbins is not None:
        Plot.setRebin(setplot_rebin_mincounts, setplot_rebin_maxbins)

    mpl.rcParams['xtick.labelbottom'] = True

    plots_dir = "./spectral_fit_residuals"
    if fit_with_binning: plots_dir+="_bin/"
    else: plots_dir+="_bin1/"

    Plot('data resid') # data from most recent Fit that was performed

    # Get coordinates from plot:
    x = Plot.x()
    rates = Plot.y()
    yer_counts = Plot.yErr()
    xer_counts = Plot.xErr()
    folded = Plot.model()
    #resids = np.array(rates) - np.array(folded)
    resids = Plot.y(1,2)

    dataLabels = Plot.labels(1)
    residLabels = Plot.labels(2)

    fig, ax = plt.subplots(3, 1, figsize=(12, 9), sharex=True, gridspec_kw={'height_ratios': [3, 1, 1]})

    ax[0].errorbar(x, rates,xerr=xer_counts, yerr=yer_counts, fmt='ro',ms=1, label="data", elinewidth=0.2)
    ax[0].plot(x, folded, 'b', label="model", linewidth=1, zorder=100)
    ax[0].set_ylabel(dataLabels[1]) #ax[0].set_ylabel(r'counts/cm$^2$/sec/keV')
    ax[0].legend(fontsize=11)

    ax[1].plot(x, resids, 'g', label="residuals (data-model)", linewidth=1)
    ax[1].set_ylabel(residLabels[1])
    ax[1].set_xlabel(residLabels[0]) #ax[1].set_xlabel('Energy [keV]')
    ax[1].legend(fontsize=11)

    ax[2].plot(x, np.array(resids)/yer_counts, 'g', label=r"(data-model)/$\sigma_{\text{data}}$", linewidth=1)
    ax[2].set_ylabel(r'')
    ax[2].set_xlabel(residLabels[0])
    ax[2].legend()

    title = "Spectrum: " + spectrum_name + " & model: " + mod_name
    if setplot_rebin_mincounts is not None and setplot_rebin_maxbins is not None: title+= + "\n bin mincounts: " + str(setplot_rebin_mincounts) + "; bin maxbins: " + str(setplot_rebin_maxbins)
    ax[0].set_title(title)

    # Save the results
    if save: plt.savefig(plots_dir+ spectrum_name+"_"+mod_name)
    
    if not save and not log: plt.show()

    # Save results on log scales
    ax[0].set_yscale("log")
    ax[0].set_xscale("log")
    if save: plt.savefig(plots_dir+ spectrum_name+"_"+mod_name+"_log")
    
    if not save and log: 
        plt.draw()  # force redraw
        plt.show()

    # Clear the plot
    #plt.clf()
    plt.close()


    #################

    ## Plot the unfolded spectrum
    Plot.add = True
    Plot("eeufspec resid")

    x = np.array(Plot.x())
    y = np.array(Plot.y())
    yer = np.array(Plot.yErr())
    xer = np.array(Plot.xErr())
    unfolded = np.array(Plot.model())

    # Build additive component label list (excluding multiplicative like TBabs)
    MULTIPLICATIVE = {'tbabs', 'phabs', 'wabs', 'tbnew', 'zphabs', 'cabs',
                      'pcfabs', 'pwab', 'gabs', 'edge', 'notch', 'smedge',
                      'constant', 'cflux', 'zashift', 'zwabs'}
    additive_comp_names = [
        name for name in mod.componentNames
        if name.lower().split('_')[0] not in MULTIPLICATIVE
    ]

    fig, ax = plt.subplots(2, 1, figsize=(12, 9), sharex=True,
                           gridspec_kw={'height_ratios': [1, 1]})

    ax[0].errorbar(x, y, xerr=xer, yerr=yer, fmt='ro', ms=1,
                   label='Data (unfolded)', elinewidth=0.2)
    ax[0].plot(x, unfolded, color='black', lw=2, label='Total model')


    n_add = Plot.nAddComps()

    linestyles = ['--', '-.', ':', (0, (3, 1, 1, 1))]
    for i in range(n_add):
        try:
            comp_arr = np.array(Plot.addComp(i + 1, plotGroup=1))
            label = additive_comp_names[i] if i < len(additive_comp_names) else f"addcomp {i+1}"
            ax[0].plot(x, comp_arr, ls=linestyles[i % len(linestyles)],
                       lw=1.5, label=label)
        except Exception as e:
            print(f"Could not retrieve additive component {i+1}: {e}")

    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[0].set_xlabel("Energy (keV)")
    ax[0].set_ylabel(r"$E^2 F(E)$ (keV$^2$ photons cm$^{-2}$ s$^{-1}$ keV$^{-1}$)")
    ax[0].legend()
    title = f"{spectrum_name} – {mod_name} (Unfolded)" 
    if setplot_rebin_mincounts is not None and setplot_rebin_maxbins is not None: title+= "\n bin mincounts: " + str(setplot_rebin_mincounts) + "; bin maxbins: " + str(setplot_rebin_maxbins)
    ax[0].set_title(title)

    ax[1].plot(x, np.array(resids)/yer_counts, 'g', label=r"(data-model)/$\sigma_{\text{data}}$", linewidth=1)
    ax[1].set_ylabel(residLabels[1])
    ax[1].set_xlabel(residLabels[0])
    ax[1].legend(fontsize=11)

    if save: plt.savefig(plots_dir+ spectrum_name+"_"+mod_name+"_unfolded_log")
    else: plt.show()
    plt.close(fig)

    Plot.add = False 



##########################################################################################


def initialise_model(model, parameters=None, fix_names=None, fix_values=None, flux_guess=None):
    """
    Initialise different spectral models.
    Note: The model-fitting for the more complex compound models can be tricky, and hands-on fitting would be preferred here. 
    Since we are only interested in the flux (not a detailed fit), the simpler models are likely appropriate, at the cost of worse fits.
    In some cases, emission/absorption lines are visible (especially in NICER data), and we may also need a more manual fit. 
    
    Parameters:
    - model: the spectral model to fit. The options are: ['powerlaw', 'pegged_powerlaw', 'diskbb', 'powerlaw+diskbb', 'pegged_powerlaw+diskbb', 'diskbb+bbodyrad','powerlaw+bbodyrad'].
    - parameters: dictionary of parameter initialisations. If None, default values are used. The keys should be the same as the parameter names in the xspec model (e.g., 'nh', 'gamma', 'Tin', 'kT').
    - fix_names: list of parameter names to fix during the fit. The names should be the same as the keys in the 'parameters' dictionary.
    - fix_values: list of values corresponding to the parameters in 'fix_names' to which they should be fixed during the fit.
    - flux_guess: an initial guess for the flux, used for the pegged powerlaw normalisation and cflux flux parameter. This can be obtained from the count rate and a rough count rate to flux conversion factor (e.g., 1 count/s ~ 1e-11 erg/cm^2/s), or from a previous fit with a simpler model.


    Information on parameter initialisation in xspec/PyXspec:
    - We can enter ranges for the parameter by passing in a string containing <param>,<delta>,<hard min>,<soft min>,<soft max>,<hard max> where:
    - param: trial param value used initially in fit
    - delta: step size used in the numerical determination of the derivatives during the fitting process. This value may be overriden for all parameters by the xset delta command option, which will apply a proportional rather than a fixed delta.
    - Or, if we specify a value and then -1, it is fixes to this value.
    - In our fitting, if any of the values get pegged at the bounds, it either means; (i) the spectrum is insensitive to the pegged model, and thus it might not be necessary; (ii) we are stuck in a local minimum and may require more hands-on spectral fitting.

    The models we make sure of are:
    - tbabs (see XSPEC manual pp. 348-349): Absorption due to the ISM including molecules and grains. Allows the user just to vary the hydrogen column
    - pgpwrlw (see XSPEC manual pg. 297): Power law with pegged normalisation. A(E) = K E^{-alpha}
    - cflux (see XSPEC manual pp. 359-360): A convolution model to calculate the flux of other model components.
    - diskbb (see XSPEC manual pg. 249): Multiple blackbody accretion disk model.
    - powerlaw (see XSPEC manual pg. 303)
    
    """

    Emin, Emax = EMIN, EMAX
    nH_init, Tin_init, gamma_init, gamma_init_IMS, kT_init = NH_INIT, DISKBB_TIN_INIT, PLAW_GAMMA_INIT, PLAW_GAMMA_INIT_IMS, BBODY_KT_INIT
    
    # Set the initialisations
    if parameters==None:
        parameters = {
        "nh":nH_init,  # absorption
        "Tin": Tin_init,  # temperature
        "norm1": '1.0',
        "norm2": '1.0', 
        "gamma":gamma_init, # spectral index
        "kT": kT_init
        }
        # The IMS fits can be a bit sensitive to initialisation.
        # So, we may want to, for example, set the diskbb+powerlaw initial gamma to be slightly higher.
        if model==  'powerlaw+diskbb' or model== 'pegged_powerlaw+diskbb': 
            parameters["gamma"] = gamma_init_IMS

    # Set parameters to the values specified in 'fix'
    for name, value in zip(fix_names, fix_values):
        if name in parameters: parameters[name] = f"{value} -1"

    # Get the parameters
    nh = parameters["nh"]
    try: gamma = parameters["gamma"]
    except: gamma=None
    try: Tin = parameters["Tin"]
    except: Tin=None
    norm1 = parameters["norm1"]
    norm2 = parameters["norm2"]
    try: kT = parameters["kT"]
    except: kT=None

    # Distance information
    #D = DIST_KPC # kpc
    #bbody_norm = 4 * np.pi * (D * 3.086e+21) ** 2 * flux_guess * 1e-12 / 1e39 / (D / 10) ** 2
    #diskbb_norm = (20 / (D / 10)) ** 2 # Assume 20km isco


    flux_guess_logged = np.log10(flux_guess * 1e-12)

    
    ############# Initialise the models

   
    ## Absorbed power law with pegged normalisation
    # par1: nH; equivalent hydrogen column, in units of 10^{22} atoms/cm^2 
    # par2: alpha; photon index of pwrlw (dimensionless)
    # par3: lower peg energy range (keV)
    # par4: upper peg energy range (keV)
    # par5: norm; flux (in units of 10^{-12} ergs/cm^2/s over the energy par2-par3). If par3 = par4, it is the flux in micro-Jy at par3.
    if model == "pegged_powerlaw":
        mod = 'tbabs * pegpwrlw' 
        initial_pars = {1:nh, 2: gamma, 3:Emin, 4:Emax, 5:f'{flux_guess}'}

    ## Absorbed power law 
    # par1: nH; equivalent hydrogen column, in units of 10^{22} atoms/cm^2 
    # par2: Emin; Minimum energy over which flux is calculated.
    # par3: Emax; Maximum energy over which flux is calculated.
    # par4: lg10Flux; log (base 10) flux in erg/cm^2/s
    # par5: alpha; photon index of pwrlw (dimensionless)
    # par6: norm; K, photons/keV/cm2/s at 1 keV.
    elif model == "powerlaw":
        pwrlw_flux = flux_guess_logged
        mod = 'tbabs * cflux * powerlaw' 
        initial_pars = {1:nh, 2: Emin, 3: Emax, 4: f'{pwrlw_flux}', 5:gamma, 6:'1.0 -1'}

    ## Absorbed disk blackbody
    # par1: nH; equivalent hydrogen column, in units of 10^{22} atoms/cm^2 
    # par2: Emin; Minimum energy over which flux is calculated.
    # par3: Emax; Maximum energy over which flux is calculated.
    # par4: lg10Flux; log (base 10) flux in erg/cm^2/s
    # par5: temperature of the inner disk radius (keV)
    # par6: norm; (Rin/ D10)^2 cos0, where where Rin is an apparent inner disk radius in km, D10 the distance to the source in units of 10 kpc, and the angle of the disk ( = 0 is face-on). On the correction factor between the apparent inner disk radius and the realistic radius, see e.g. Kubota et al. 1998.
    elif model == "diskbb":
        bb_flux = flux_guess_logged
        mod = 'tbabs * cflux *  diskbb'
        initial_pars = {1:nh, 2:Emin, 3:Emax, 4:f'{bb_flux}', 5: Tin, 6:'1.0 -1'}
 
    ## Absorbed powerlaw + disk blackbody
    # par1: nH; equivalent hydrogen column, in units of 10^{22} atoms/cm^2 
    # par2: alpha; photon index of pwrlw (dimensionless)
    # par3: lower peg energy range (keV) for pwrlw
    # par4: upper peg energy range (keV) for pwrlw
    # par5: norm; pwrlw flux (in units of 10^{-12} ergs/cm^2/s over the energy par2-par3). If par3 = par4, it is the flux in micro-Jy at par3.
    # par6: Emin; Minimum energy over which bb flux is calculated.
    # par7: Emax; Maximum energy over which bb flux is calculated.
    # par8: lg10Flux; log (base 10) bb flux in erg/cm^2/s
    # par9: temperature of the inner disk radius (keV)
    # par10: norm for bb; (Rin/ D10)^2 cos0, where where Rin is an apparent inner disk radius in km, D10 the distance to the source in units of 10 kpc, and the angle of the disk ( = 0 is face-on). On the correction factor between the apparent inner disk radius and the realistic radius, see e.g. Kubota et al. 1998.
    elif model == "pegged_powerlaw+diskbb":
        mod = 'tbabs * (pegpwrlw + cflux *  diskbb)'
        bb_flux = np.log10(flux_guess * 0.1 * 1e-12)
        pwrlw_flux = flux_guess*0.9
        initial_pars = {1:nh, 2: gamma, 3:Emin, 4:Emax, 5:f'{pwrlw_flux}',
                        6:Emin, 7:Emax, 8:f'{bb_flux}', 9: Tin, 10:'1.0 -1'}


    ## Absorbed powerlaw + disk blackbody
    # par1: nH; equivalent hydrogen column, in units of 10^{22} atoms/cm^2 
    # par2: alpha; photon index of pwrlw (dimensionless)
    # par3: norm; K, photons/keV/cm2/s at 1 keV.
    # par4: temperature of the inner disk radius (keV)
    # par5: norm for bb; (Rin/ D10)^2 cos0, where where Rin is an apparent inner disk radius in km, D10 the distance to the source in units of 10 kpc, and the angle of the disk ( = 0 is face-on). On the correction factor between the apparent inner disk radius and the realistic radius, see e.g. Kubota et al. 1998.
    elif model=="powerlaw+diskbb":
        mod = 'tbabs*(powerlaw + diskbb)'
        initial_pars = {1:nh, 2:gamma, 3:'1.0', 4: Tin, 5:'1.0'}


    ## Absorbed powerlaw + disk blackbody
    # par1: nH; equivalent hydrogen column, in units of 10^{22} atoms/cm^2 
    # par2: Emin; Minimum energy over which flux is calculated.
    # par3: Emax; Maximum energy over which flux is calculated.
    # par4: lg10Flux; log (base 10) flux in erg/cm^2/s
    # par5: alpha; photon index of pwrlw (dimensionless)
    # par6: norm; K, photons/keV/cm2/s at 1 keV.
    # par7: temperature of the inner disk radius (keV)
    # par8: norm for bb; (Rin/ D10)^2 cos0, where where Rin is an apparent inner disk radius in km, D10 the distance to the source in units of 10 kpc, and the angle of the disk ( = 0 is face-on). On the correction factor between the apparent inner disk radius and the realistic radius, see e.g. Kubota et al. 1998.
    # IMPORTANT: When we pass parameters to this model, norm1 is always fixed
    elif model=="cflux_(powerlaw+diskbb)":
        mod = 'tbabs*cflux(powerlaw + diskbb)'
        flux = flux_guess_logged
        initial_pars = {1:nh, 2: Emin, 3: Emax, 4: f'{flux}', 5:gamma, 6:norm1, 7: Tin, 8: norm2}
        

    elif model=="diskbb+bbodyrad":
        mod = 'tbabs * (diskbb + bbodyrad)'
        initial_pars = {1:nh, 2:Tin, 3:norm1, 4: kT, 5:norm2}
        
        
    # IMPORTANT: When we pass parameters to this model, norm1 is always fixed
    elif model=="cflux_(diskbb+bbodyrad)":
        mod = 'tbabs * cflux* (diskbb + bbodyrad)'
        flux = flux_guess_logged
        initial_pars = {1:nh, 2: Emin, 3: Emax, 4: f'{flux}', 5:Tin, 6:norm1, 7: kT, 8: norm2}
        
    # IMPORTANT: When we pass parameters to this model, norm1 is always fixed
    elif model=="powerlaw+bbodyrad":
        mod = 'tbabs * (powerlaw + bbodyrad)'
        initial_pars = {1:nh, 2:gamma, 3:norm1, 4: kT, 5:norm2}
        

    # IMPORTANT: When we pass parameters to this model, norm1 is always fixed
    elif model=="cflux_(powerlaw+bbodyrad)":
        mod = 'tbabs * cflux* (powerlaw + bbodyrad)'
        flux = flux_guess_logged
        initial_pars = {1:nh, 2: Emin, 3: Emax, 4: f'{flux}', 5:gamma, 6:norm1, 7: kT, 8: norm2}

    else:
        sys.exit('Invalid model.')

    print(initial_pars)

    return mod, initial_pars # return the model and parameters



########################################################################################## 



def run_spectral_fit( spectral_folder = "./spectra_swift_xrt/" ):
    """
    Run the spectral fitting for the Swift XRT data.

    Parameters:
    - mjd_min, mjd_max: the MJD range over which to perform the spectral fitting. Only spectra with mid-MJD within this range will be fitted.
    - models: a list of spectral models to fit to the data. The options are: ['powerlaw', 'pegged_powerlaw', 'diskbb', 'powerlaw+diskbb', 'pegged_powerlaw+diskbb', 'diskbb+bbodyrad','powerlaw+bbodyrad'].
    - spectral_folder: the folder containing the spectral files to fit. The code expects the files to be in the format '*final_bin.pi' or '*final_bin1.pi', where 'bin' or 'bin1' corresponds to whether the spectra have been grouped with a minimum number of counts per bin or not, respectively. 

    
    Notes on xspec/PyXspec:
    - AllData: container for all loaded data sets (objects of class Spectrum).
    - AllModels: container for all Model objects.
    - XSPEC models are defined by creating a Model object.
    Note: getattr(object, 'x') is completely equivalent to object.x, but used in the case we don't know exactly what 'x' is.


    IMPORTANT NOTE:
    - cflux is found to somtimes behave weirdly when two-component fits are used. 
    - So in these cases, we first fit without it, and then add the cflux convolution after. 
    - Since pyXspec does not have 'editmod', we have to do a workaround to do this. 
    - After adding cflux, the norm of one of the components needs to be fixed to prevent degeneracy with the lg10flux.
    """

    #############################################
    ## GET THE PARAMETERS

    models = MODELS

    incbad, fit_with_binning, renorm, add_systematic, counts_threshold = INCBAD, FIT_WITH_BINNING, RENORM, ADD_SYS_ERR, COUNTS_THRESHOLD
    min_E_keV, Emin, Emax = MIN_E_KEV, EMIN, EMAX
    nH, nH_fix_all_epochs, nH_counts_threshold = NH, NH_FIX_ALL_EPOCHS, NH_COUNTS_THRESHOLD
    low_count_threshold, Gamma_low_count = LOW_COUNT_THRESHOLD, PLAW_GAMMA_LOW_COUNT
    nH_init, Tin_init, gamma_init, gamma_init_IMS, kT_init = NH_INIT, DISKBB_TIN_INIT, PLAW_GAMMA_INIT, PLAW_GAMMA_INIT_IMS, BBODY_KT_INIT
    min_mjd, max_mjd = MJD_MIN, MJD_MAX


    if any(model not in ['powerlaw', 'pegged_powerlaw', 'diskbb', 'powerlaw+diskbb', 'pegged_powerlaw+diskbb', 'diskbb+bbodyrad','powerlaw+bbodyrad'] for model in models):
        raise ValueError(f"Error: Invalid model(s) found in array. Allowed values are: {['powerlaw', 'pegged_powerlaw', 'diskbb', 'powerlaw+diskbb', 'pegged_powerlaw+diskbb', 'diskbb+bbodyrad','powerlaw+bbodyrad']}")




    #############################################
    ## SET UP THE FOLDERS FOR THE RESULTS

    if fit_with_binning: name = "bin"
    else: name = "bin1"

    # Folder for the results
    folder = './spectral_fit_results'
    folder+= "_"+name+"/"
    if os.path.exists(folder):
        shutil.rmtree(folder)  
    os.makedirs(folder)  
    


    # Folder for the residuals
    plots_dir = "./spectral_fit_residuals"
    plots_dir+= "_"+name+"/"
    if os.path.exists(plots_dir):
        shutil.rmtree(plots_dir)  
    os.makedirs(plots_dir) 


    print(f"Spectral results will be saved to: {folder}")
    print(f"Spectral fit residuals will be saved to: {plots_dir}\n")
    print()


    print("Xspec version: ", Xset.version)
    if fit_with_binning: print("Fitting with binning...\n")
    else: print("Fitting with bin1...\n") 
      


    #############################################
    ## INITIALISE XSPEC
    
    print('Configuring XSPEC settings')

    # Xset: storage class for XSPEC settings
    Xset.abund = "wilm"    # set the abundance table used in the plasma emission and photoelectric absorption models; 'wilm' is an XSPEC built-in abundance table
    Xset.xsect = "vern"    # set the photoelectric absorption cross-sections
    Xset.chatter = 2 # set the console chatter level... change this to 10 for more chatter
    Xset.logChatter = 20
    Xset.openLog(folder+"xspec.log") 
    #Xset.parallel.error = 1 # previously 10, i.e. use up to 10 parallel processes during the Fit.error() command runs... however, this caused issues with the parameter saving being done before the values were updated during the error calculation

    # Fit: manager class for setting properties and running a fit
    Fit.query = "yes"      # when nIterations is reached, continue the fit without stopping to query
    Fit.nIterations = 100  # fit iterations for each attempt
    Fit.statTest = "chi" # default 

    # Plot: manager class for performing XSPEC plots
    Plot.xAxis = "KeV"     # x-axis for plotting set to energy instead of channel
    Plot.device = '/null'
    Plot.yLog = True  # Use logarithmic y-axis for plots


    #############################################
    ## GET THE SPECTRAL FILES 

    print('Reading spectral files')

    # Load in the spectral files
    try:
        spectra_pattern = spectral_folder+f'*final_{name}.pi'
        spectral_files = glob.glob(spectra_pattern) 
        n_spectra = len(spectral_files) 
        print("Number of files: ", n_spectra)
    except Exception as e:
        print(f"Error reading spectral files: {e}")
        return
    
    # Initialise variables
    start_obs_isots  = [] # start observation times (ISO format)
    end_obs_isots = [] # end observation times (ISO format)
    mid_mjds = [] # mid observation MJD (ISO format)
    dt_mjds = [] # observation durations in MJD
    bin_counts = [] # to hold number of counts in grouped bin 
    tot_counts = [] # counts in the spectrum over the fitted energy range
    exp = [] # exposure time

    
    # Loop through the spectral files and extract the relevant information for each spectrum.
    for i, spectrum_file in enumerate(spectral_files):

        # Get the required information from the header.
        obs_id, bin_count, tot_count, date_start_str, date_end_str, mid_mjd_1, dt_1, mid_mjd_2, dt_2, exp_days = extract_file_info(spectrum_file)
        start_obs_isots.append(date_start_str) # start date
        end_obs_isots.append(date_end_str) # end date
        mid_mjds.append(mid_mjd_2)
        dt_mjds.append(dt_2)
        bin_counts.append(bin_count) # min number of counts in the bins
        tot_counts.append(tot_count) #... should be same as s.rate[0] * s.exposure
        exp.append(exp_days) # exposure length... same as s.exposure


    # Convert dates to MJDs
    start_obs_mjds = iso2mjd(np.array(start_obs_isots)) 
    end_obs_mjds = iso2mjd(np.array(end_obs_isots)) 


    # Order lists by time
    sort_index = np.argsort(mid_mjds)
    spectral_files = np.array(spectral_files)[sort_index]
    start_obs_isots = np.array(start_obs_isots)[sort_index]
    start_obs_mjds = start_obs_mjds[sort_index]
    end_obs_isots = np.array(end_obs_isots)[sort_index]
    end_obs_mjds = end_obs_mjds[sort_index]
    mid_mjds = np.array(mid_mjds)[sort_index]
    dt_mjds = np.array(dt_mjds)[sort_index]
    bin_counts = np.array(bin_counts)[sort_index]
    tot_counts = np.array(tot_counts)[sort_index]
    exp = np.array(exp)[sort_index]
    
    # Also get the ID and mode (WT/PC)
    ID_LEN = 13  # ID length, including the mode
    IDs = []
    modes = []
    for sf in spectral_files:
        name = sf.split("Obs_")[1]   # everything after "Obs_"
        ID = name[:ID_LEN]
        mode = name[ID_LEN-2:ID_LEN]
        IDs.append(ID)
        modes.append(mode.lower())
    IDs = np.array(IDs)
    modes = np.array(modes)


    # Boolean stating whether the data was binned to >20 counts per bin
    # We opt to do all the fits with C-stat, as the chi-squared fits can be biased. 
    binned = np.array([True if bin_count>=20 else False for bin_count in bin_counts])


    #############################################
    ## FILTERING THE SPECTRAL FILES 
    # We ignore the following: 
    # (1) Observations identified as uplims from the lightcurve generation. 
    # (2) PC observations with fewer than 10 counts. These are too few to reliably fit the spectra.
    # (3) WT observations with fewer than 15 counts. These are likely spurious.
    # (4) Filter based on the MJD range specified for the fit.
 

    # Get the IDs of the uplims from the lightcurve generation, and ignore these spectra.
    uplim_ids = get_swift_xrt_counts(verbose=False, incbad = incbad, just_ids = True)
    mask_detections = np.array([False if ID in uplim_ids else True for ID in IDs])
    # Also ignore spectra with too few counts, based on the mode of the observation (WT vs. PC)
    mask_valid = ((modes == 'pc') & (tot_counts >= 10)) | ((modes == 'wt') & (tot_counts >= 15))
    

    # Filter based on the MJD range specified for the fit.
    mjd_mask = np.ones(len(mid_mjds), dtype=bool)
    if min_mjd is not None or max_mjd is not None:
        if min_mjd is not None:
            mjd_mask &= mid_mjds >= min_mjd
        if max_mjd is not None:
            mjd_mask &= mid_mjds <= max_mjd

    # Final mask
    mask = mask_detections & mask_valid & mjd_mask

    # Print to file the obs_isot, obs_mjds, IDs, counts, exp [s] for the ignored spectra
    xrt_dict_no_fit = {"start_obs_isot": start_obs_isots[~mask],"mid_mjds": mid_mjds[~mask],"IDs": IDs[~mask],"counts": tot_counts[~mask],"exp": exp[~mask]}
    df = pd.DataFrame(xrt_dict_no_fit)
    df.to_string(folder+"not_fit.txt", index=False, justify="left")

    # Get the parameters for the observations that we will fit
    spectral_files = spectral_files[mask]
    start_obs_isots = start_obs_isots[mask]
    mid_mjds = mid_mjds[mask]
    dt_mjds = dt_mjds[mask]
    bin_counts = bin_counts[mask]
    tot_counts = tot_counts[mask]
    exp = exp[mask]
    IDs = IDs[mask]
    modes = modes[mask]
    binned = binned[mask]
    print("Number of files to be fit (after filtering): ", len(spectral_files))

    # Initialise a dictionary to hold the results
    xrt_dict = {model: {"IDs": IDs.tolist(), "isot_i": start_obs_isots.tolist(), "mjd_mid": mid_mjds.tolist(), "dt_mjd": dt_mjds.tolist(), "exp [s]": exp.tolist(), "counts": tot_counts.tolist(), "binned?": binned.tolist(), "chi2": [], "dof": [], "redchi2": []} for model in models}

    print("------------------------------------------------------------")


    #############################################

    # Outputs from here on will be saved to file
    sys.stdout = open(folder + "output.txt", "w", buffering=1) 


    #############################################
    ## Define functions to be used during the fit

    class FitTimeout(Exception):
        pass

    def _timeout_handler(signum, frame):
        raise FitTimeout("Fit.error timed out")
    

    #############################################
    
    ## Iterate through each spectrum (i.e. each observation), fit it, get the fit parameters, and then save these.
    for k, spectrum_file in enumerate(spectral_files):

        print('\n\n', spectrum_file , " MJD: ", mid_mjds[k], " counts: ", tot_counts[k], " exp [s]: ", exp[k], " mode: ", modes[k], " binned? ", binned[k])


        #################
        ## INITIAL SET-UP FOR THE SPECTRUM

        # Load the spectrum into the AllData object, removing any previously loaded data sets... 
        # We also run AllData.clear() before, just to be sure.
        # If the rmf and arf files are located in the same folder (which they are), these will automatically be loaded. In our case, they are also listed in the header to be safe.
        AllModels.clear() # clear any previously defined models
        AllData.clear()
        AllData(spectrum_file)
        # AllData.show() # check current state of the AllData container


        # Define range of energy to use -- i.e. ignore certain channels. This is done for all the loaded spectra.
        # The floating point values are assumed to be in keV, as this is what was set for the Plot.xAxis value above.
        AllData.ignore('*:10.0-**')
        AllData.ignore(f'*:**-{min_E_keV}')
        AllData.notice(f'*:{min_E_keV}-10.0') # so that edge points are used
        AllData.ignore('bad') # ignore data that is bad -- using the quality column

        # If too much data is ignored, the spectrum won't fit 
        # Also print out a warning 
        #if (len(AllData(1).ignored) == len(AllData(1).noticed)): print("ERROR: too much bad data")

        # We choose to always use C-stat to fit, as chi^2 fits can result in biases. 
        # Note that when a background is loaded, this actually uses the W-statistic, which is a modified version of the C-statistic that accounts for the background.
        Fit.statMethod = 'cstat' # Poisson data
        if bin_counts[k] == 1: print("Bin1. Using C-statistics.")
        else: print("Binned. Using C-statistics.")

        # If there are more than 20 counts per bin, it is still reasonable to use chi-squared statistics to evaluate the fit.
        if bin_counts[k] >= 20: chi_stat_test = True # equivalent to if binned[k]==True
        else: chi_stat_test = False
        


        #################
        ## DETERMINE WHAT PARAMETERS TO FIX DURING THE FIT

        fix_names = [] # names of parameters to fix for this spectrum
        fix_values = [] # values to fix these parameters to for this spectrum
        models_to_fit = models # models to fit for this spectrum

        # If it was specified that nH should be fixed for all epochs 
        if nH_fix_all_epochs: 
            fix_names = ["nh"]
            fix_values = [nH] 

        # If it was specified that nH should be fixed for spectra with fewer than nH_counts_threshold counts
        if  nH_counts_threshold!=None and tot_counts[k] < nH_counts_threshold:  
            fix_names = ["nh"]
            fix_values = [nH] 
        
        # If this spectrum has fewer than low_count_threshold counts, we only fit with a powerlaw
        # and we fix both nH and gamma to reasonable values for faint spectra
        if tot_counts[k] < low_count_threshold:
            fix_names = ["nh", "gamma"] # names of parameters to fix for this spectrum
            fix_values = [nH, Gamma_low_count] # values to fix these parameters to for this spectrum
            models_to_fit = ["pegged_powerlaw"]


        print("Fixed parameters: ", fix_names)
        print("Fixed parameter values: ", fix_values)


        #################
        ## GET FLUX ESTIMATE
        # Perform an initial power-law fit to estimate the flux (flux_guess)
        # flux_guess will be used in the spectral fitting routine. 
        

        # Define the spectral model -- absorbed power law
        # 1: nH (10^{22} atoms/cm^2); 2: photon index; 3: lower peg energy range (keV); 4: upper peg energy range (keV); 5: norm i.e. flux (10^{-12} ergs/cm^2/s)
        AllModels.clear()
        mod1 = Model('tbabs * pegpwrlw')
        mod1.setPars({1:f'{nH} -1', 2:'1.7 -1', 3:Emin, 4:Emax, 5:'1000.0'}) 
        Fit.delta = 1e-2 # controls the convergence criterion for the fitting algorithm
        try: 
            Fit.perform()
            flux_guess = mod1.pegpwrlw.norm.values[0] # flux (in units of 10^{-12} ergs/cm^2/s)
        except: flux_guess = 1 # faint data... i.e. assume flux = 1 * 10^{-12} ergs/cm^2/s


        #################
        ## ITERATE THROUGH THE MODELS
        # Iterate through the models, fitting, and appending the results to an output dictionary
        # If chi_stat_test, then we can compare the chi-squared values of the fits

        # Parameters when comparing models for when the spectra are sufficiently binned.
        best_model = None
        best_model_chi2 = np.inf
        best_model_dof = np.inf   
        
        
        mod_2 = False # whether the second model succeeded in the 2-step model approach
        for mod_name in models:
            print('\n', mod_name)

            # Get the model
            AllModels.clear() # AllModels represents the collection of all currently defined models in the XSPEC session
            parameters = None
            mod, initial_pars = initialise_model(mod_name, parameters, fix_names, fix_values, flux_guess)
            mod_obj = Model(mod)
            

            if mod_name in models_to_fit: fit = True
            else: 
                print("Fit will not be performed.")
                fit = False


            if fit:

                # Add systematic error, if required
                # Note that this needs to be before the model is initialised
                # And it does not have an effect when C-statistics is used, but it does when chi-squared statistics is used.
                if add_systematic: AllModels.systematic = 0.03 
                
                # Initialise the model
                mod_obj.setPars(initial_pars) # set the parameters to the initial parameters specified in the initialisation
                AllModels.show() # for checking


                # Renormalisation can sometimes help
                # However, this only changes the values when the stat method is chi
                # So change back to C-stat fitting afterwards
                if renorm:
                    # Renormalisation helps sometimes -- this only works when the stats method is set to chi-squared, so we need to change to this first 
                    Fit.statMethod = 'chi' 
                    Fit.renorm()
                    print("Renormalised")
                    Fit.statMethod = 'cstat' # Poisson data
        
                             
                ## Try to perform the fitting
                try:

                    # Adjust the fit precision, and perform the fit of the spectrum. 
                    # The first pass allows the fit to reach a close approximation of the best-fit parameters quickly. 
                    # Each subsequent pass with a smaller delta refines the fit further.
                    for delt in [1e-2, 1e-3, 1e-4, 1e-5]:
                        Fit.delta = delt
                        Fit.perform()
                    time.sleep(1)

          
                    #######

                    ## For the two-component models, we add cflux afterwards
                    # The first normalisation needs to be fixed: https://heasarc.gsfc.nasa.gov/docs/software/xspec/manual//node303.html

                    # Two-step powerlaw+diskbb
                    if mod_name == "powerlaw+diskbb":
                        print("Adding cflux for powerlaw+diskbb...")
                        # Get best-fit parameters and pass to initialise_model
                        parameters = {
                        "nh": f'{mod_obj.TBabs.nH.values[0]},'+ nH_init.split(',', 1)[1],  
                        "gamma": f'{mod_obj.powerlaw.PhoIndex.values[0]},'+ gamma_init_IMS.split(',',1)[1],  
                        "Tin": f'{mod_obj.diskbb.Tin.values[0]},'+ Tin_init.split(',',1)[1], 
                        "norm1": f'{mod_obj.powerlaw.norm.values[0]} -1', # fix the first normalisation
                        "norm2": f'{mod_obj.diskbb.norm.values[0]}'
                        }
                        AllModels.clear()
                        mod, initial_pars = initialise_model("cflux_(powerlaw+diskbb)", parameters= parameters, fix_names = fix_names, fix_values = fix_values, flux_guess = flux_guess)
                        if add_systematic: AllModels.systematic = 0.03 # add systematic error
                        mod_obj = Model(mod)
                        mod_obj.setPars(initial_pars) # set the parameters to the initial parameters specified in the initialisation
                        mod_2 = True 
                        AllModels.show() # for checking
                        Fit.statMethod = 'cstat'
                        Fit.nIterations = 100 
                        Fit.delta = 1e-5
                        Fit.perform()


                    elif mod_name =='diskbb+bbodyrad':
                        print("Adding cflux for diskbb+bbodyrad...")
                        # Get best-fit parameters and pass to initialise_model
                        parameters = {
                        "nh": f'{mod_obj.TBabs.nH.values[0]},'+ nH_init.split(',',1)[1],  
                        "Tin": f'{mod_obj.diskbb.Tin.values[0]},'+ Tin_init.split(',',1)[1], 
                        "kT": f'{mod_obj.bbodyrad.kT.values[0]},'+ kT_init.split(',',1)[1], 
                        "norm1": f'{mod_obj.diskbb.norm.values[0]} -1', # fix the first normalisation
                        "norm2": f'{mod_obj.bbodyrad.norm.values[0]}'
                        }
                        #print(parameters)
                        AllModels.clear()
                        mod, initial_pars = initialise_model("cflux_(diskbb+bbodyrad)",  parameters= parameters, fix_names = fix_names, fix_values = fix_values, flux_guess = flux_guess)
                        if add_systematic: AllModels.systematic = 0.03 # add systematic error
                        mod_obj = Model(mod)
                        print( "Initial parameters for cflux_(diskbb+bbodyrad): ", initial_pars)
                        mod_obj.setPars(initial_pars) # set the parameters to the initial parameters specified in the initialisation
                        mod_2 = True 
                        AllModels.show() # for checking
                        Fit.statMethod = 'cstat'
                        Fit.nIterations = 100 
                        Fit.delta = 1e-5
                        Fit.perform()

                    elif mod_name =='powerlaw+bbodyrad':
                        print("Adding cflux for powerlaw+bbodyrad...")
                        # Get best-fit parameters and pass to initialise_model
                        parameters = {
                        "nh": f'{mod_obj.TBabs.nH.values[0]},'+ nH_init.split(',',1)[1],   
                        "gamma": f'{mod_obj.powerlaw.PhoIndex.values[0]},'+ gamma_init.split(',',1)[1], 
                        "kT": f'{mod_obj.bbodyrad.kT.values[0]},'+ kT_init.split(',',1)[1],  
                        "norm1": f'{mod_obj.powerlaw.norm.values[0]} -1', # fix the first normalisation
                        "norm2": f'{mod_obj.bbodyrad.norm.values[0]}'
                        }
                        AllModels.clear()
                        mod, initial_pars = initialise_model("cflux_(powerlaw+bbodyrad)",  parameters= parameters, fix_names = fix_names, fix_values = fix_values, flux_guess = flux_guess)
                        if add_systematic: AllModels.systematic = 0.03 # add systematic error
                        mod_obj = Model(mod)
                        mod_obj.setPars(initial_pars) # set the parameters to the initial parameters specified in the initialisation
                        mod_2 = True 
                        AllModels.show() # for checking
                        Fit.statMethod = 'cstat'
                        Fit.nIterations = 100 
                        Fit.delta = 1e-5
                        Fit.perform()


                    #######
                    # Store the test statistics
                    # Note that this is only reliable to intepret when the data were binned (i.e. when binned[k] == True). 

                    xrt_dict[mod_name]['chi2'].append(Fit.testStatistic) # test statistic value from the most recent fit
                    xrt_dict[mod_name]['dof'].append(Fit.dof) # the degrees of freedom from the fit
                    try: xrt_dict[mod_name]['redchi2'].append(Fit.testStatistic/Fit.dof)
                    except: xrt_dict[mod_name]['redchi2'].append(-1) # dof =0 
                    print("Fit succeeded.")


                    #######
                    # Print the parameter results

                    string_to_print = ""
                    for comp_name in mod_obj.componentNames: 
                        comp = getattr(mod_obj, comp_name)
                        for par_name in comp.parameterNames:
                            param = getattr(comp, par_name)
                            string_to_print += f'{par_name} : {param.values[0]}; '
                    print("RESULT: ", string_to_print)

                
                except: # Fit performing failed
                    print("Fit failed.")
                    fit = False

            if not fit: # Fit failed or not performed
                xrt_dict[mod_name]['chi2'].append(-1) 
                xrt_dict[mod_name]['dof'].append(-1) 
                xrt_dict[mod_name]['redchi2'].append(-1)

            
            else: # The fit succeeded

                # If the spectrum is sufficiently binned, chi^2 test statistics are applicable, so compare models
                if chi_stat_test: 
                    print("chi2: ", Fit.testStatistic, " dof: ", Fit.dof, " redchi2: ", Fit.testStatistic/Fit.dof)
                    if best_model is None: # if this is the first model, set it as the best model by default
                        best_model = mod_name
                        best_model_chi2 = Fit.testStatistic
                        best_model_dof = Fit.dof
                    bool_better_model = f_test(chi1=best_model_chi2, dof1=best_model_dof, chi2=Fit.testStatistic, dof2=Fit.dof)
                    if bool_better_model:
                        best_model = mod_name
                        best_model_chi2 = Fit.testStatistic
                        best_model_dof = Fit.dof


                # Get the fit errors
                try: 
                
                    # Determine the confidence intervals of a fit (FitManager), i.e. get the parameter errors
                    # Note, the default is to estimate the 90% confidence ranges for all parameters. 
                    # If we put "1.0", this indicates that we want 68% uncertainty.
                    # The maximum keyword ensures that error will not be run if the reduced chi-squared of the best fit exceeds <redchi>. The default value for <redchi> is 2.0. We set it a bit higher.
                    n_params = mod_obj.nParameters

                    # Fit.error(f'maximum 3.0 nonew 1.0 1-{n_params}') 

                    signal.signal(signal.SIGALRM, _timeout_handler)
                    signal.alarm(60)
                    try:
                        Fit.error(f"maximum 3.0 1.0 1-{n_params}")
                    finally:
                        signal.alarm(0)
                   
                    print("Fit error calculation succeeded.")


                    # Get the results
                    for comp_name in mod_obj.componentNames: # Iterate through components, e.g. diskbb
                        comp = getattr(mod_obj, comp_name)
                        for par_name in comp.parameterNames: # Iterate through parameters, e.g. Tin
                            if par_name not in ["Emin", "Emax", "eMin", "eMax"]: # Skip fixed parameters (like norm for diskbb, and Emin/Emax)
                                if par_name == 'norm' and comp_name != 'pegpwrlw': pass 
                                else:
                                    param = getattr(comp, par_name)
                                    xrt_dict[mod_name].setdefault(par_name, []).append(param.values[0]) # best value
                                    # error[0] is lower bound and error[1] is upper bound
                                    if param.error[0]==0 and param.error[1]==0: neg_er, pos_er = 0, 0 # parameter was fixed during fitting
                                    else: neg_er, pos_er = abs(param.values[0] - param.error[0]), abs(param.error[1] - param.values[0]) 
                                    xrt_dict[mod_name].setdefault(par_name + '_neg', []).append(neg_er) 
                                    xrt_dict[mod_name].setdefault(par_name + '_pos', []).append(pos_er) 
                
                except: 
                    print("Fit error calculation failed.")
                    fit = False


            # If the fit was not performed or failed, or the error calculation failed
            # Fill in -1s for the results
            if not fit: 
                for comp_name in mod_obj.componentNames:
                    comp = getattr(mod_obj, comp_name)
                    for par_name in comp.parameterNames:
                        if par_name not in ["Emin", "Emax", "eMin", "eMax"]:
                            if par_name == 'norm' and comp_name != 'pegpwrlw': pass 
                            else: # append -1s to indicate failure
                                xrt_dict[mod_name].setdefault(par_name, []).append(-1)
                                xrt_dict[mod_name].setdefault(par_name + '_neg', []).append(-1)
                                xrt_dict[mod_name].setdefault(par_name + '_pos', []).append(-1)
                # For the 2-step diskbb+powerlaw, if the fit failed before the cflux model was initialised:
                if mod_name in ["powerlaw+diskbb","diskbb+bbodyrad","powerlaw+bbodyrad"] and mod_2==False: 
                    xrt_dict[mod_name].setdefault('lg10Flux', []).append(-1)
                    xrt_dict[mod_name].setdefault('lg10Flux' + '_neg', []).append(-1)
                    xrt_dict[mod_name].setdefault('lg10Flux' + '_pos', []).append(-1)


            # Plot the residuals
            if fit: 
                date = f"{math.floor(mid_mjds[k])}"
                base = os.path.splitext(os.path.basename(spectrum_file))[0]
                start = base.find("Obs")
                spectrum_file_name = base[start:]
                spectrum_name = "MJD"+date+"_"+spectrum_file_name # spectrum name includes date


                # When no binning was conducted for the fitting, when plotting the residuals, for visualisation (only), it is useful to rebin the spectra. 
                rebin_mincounts, rebin_maxbins = None, None
                if bin_counts[k] == 1:
                    if tot_counts[k] > counts_threshold:
                        rebin_mincounts = 20
                        rebin_maxbins = 5
                    elif tot_counts[k] > low_count_threshold: 
                        rebin_mincounts = 10
                        rebin_maxbins = 3
                    else: 
                        rebin_mincounts = 2
                        rebin_maxbins = 3
                
                # The function below plots the residuals from the last fit 
                plot_resid(spectrum_name, mod_name, mod_obj, fit_with_binning=fit_with_binning, setplot_rebin_mincounts= rebin_mincounts, setplot_rebin_maxbins = rebin_maxbins)
            
            print()



        # If we have done the model comparison
        if chi_stat_test: 
            print(f"Best model for {spectrum_file}: {best_model}")
        
        print("------------------------------------------------------------")


    #################
    ## SAVE THE RESULTS DICTIONARY

    with open(folder+'xrt_spectral_dict.json', 'w') as j:
        json.dump(xrt_dict, j, indent = 4)   
    print(f"Saved the results dictionary to {folder+'xrt_spectral_dict.json'}")
       

    # Close the log file 
    Xset.closeLog()
    
    


