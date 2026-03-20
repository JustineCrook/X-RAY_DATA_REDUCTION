import subprocess
import glob
import json
import os
import numpy as np
np.set_printoptions(threshold=np.inf, precision=16, linewidth=100)
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pickle as p
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition,mark_inset)
from astropy.time import Time
import scipy.special 
#from scipy.special import iv
from astropy.constants import c
from astropy.table import Table
from astropy.io import fits

#from PyDynamic import interp1d_unc
#rom scipy.interpolate import interp1d
import pandas as pd
pd.set_option('display.width', 1000)  
pd.set_option('display.max_columns', None)  
import warnings
mpl.pyplot.close()
import shutil
from matplotlib.lines import Line2D
import tempfile
import re


from plotting_helpers import *
from get_results_data import get_swift_xrt_counts


import sys
import os
sys.path.append(os.path.abspath("../"))
from input_parameters import *




####################################################################################################################
## SPECTRAL RESULTS FUNCTIONS
####################################################################################################################



####################################################################################################################
## GET RESULTS


## Print the spectral fit results (json file) to a txt file in table format.
# Once we have an idea which model we would like to use for each observation, we can run this with the models specified
# The model_indexes are an inclusive range and correspond the time-ordered spectra that are fit
def get_results():

    models = MODELS
    fit_with_binning = FIT_WITH_BINNING
    low_count_threshold = LOW_COUNT_THRESHOLD
    nH_counts_threshold = NH_COUNTS_THRESHOLD
    models_indexes = MODELS_INDEXES


    if any(model not in ['powerlaw', 'pegged_powerlaw', 'diskbb', 'powerlaw+diskbb', 'pegged_powerlaw+diskbb','diskbb+bbodyrad','powerlaw+bbodyrad']  for model in models):
        raise ValueError(f"Error: Invalid model(s) found in array. Allowed values are: {['powerlaw', 'pegged_powerlaw', 'diskbb', 'powerlaw+diskbb', 'pegged_powerlaw+diskbb','diskbb+bbodyrad','powerlaw+bbodyrad'] }")


    # Load JSON data from file
    filename = './spectral_fit_results'
    if fit_with_binning: filename += "_bin/"
    else: filename += "_bin1/"
    json_path = filename+"xrt_spectral_dict.json"
    with open(json_path, 'r') as file:
        data = json.load(file)

    # Output text file
    f = open(filename+"fit_outputs.txt", "w")

    # Initialise an empty list to store DataFrames
    dataframes = []

    nH_all_valid = []
    nH_all_valid_unc_l = []
    nH_all_valid_unc_u = []

    # Loop through each entry in the dictionary (i.e. for each model)
    for i, model_name in enumerate(models): 

        df = pd.DataFrame(data[model_name])
        df = df[['isot_i'] + ['mjd_mid']+ [col for col in df.columns if col != 'isot_i' and col!='mjd_mid']]
        dataframes.append(df)

        if i==0: # i.e. first model, print the low-count spectra
            counts = df['counts'].to_numpy()
            mask = counts < low_count_threshold
            obs = df['IDs'].to_numpy()[mask]
            indexes = np.where(mask)[0] 
            f.write(f"The following observations have counts < {low_count_threshold}:\n")
            f.write("IDs: " + np.array2string(obs, separator=",") +"\n")
            f.write("Indexes: " + np.array2string(indexes, separator=",") + "\n")
            f.write("\n\n\n")


            if nH_counts_threshold!=None:
                mask_new = counts < nH_counts_threshold 
                mask = (~mask) & mask_new
                obs = df['IDs'].to_numpy()[mask]
                indexes = np.where(mask)[0] 
                f.write(f"The following observations have counts > {low_count_threshold} but < {nH_counts_threshold}:\n")
                f.write("IDs: " + np.array2string(obs, separator=",") +"\n")
                f.write("Indexes: " + np.array2string(indexes, separator=",") + "\n")
                f.write("\n\n\n")

                mask = counts < nH_counts_threshold
                obs = df['IDs'].to_numpy()[mask]
                indexes = np.where(mask)[0] 
                f.write(f"The following observations have counts < {nH_counts_threshold}: \n")
                f.write("IDs: " + np.array2string(obs, separator=",") +"\n")
                f.write("Indexes: " + np.array2string(indexes, separator=",") + "\n")
                f.write("\n\n\n")


        f.write("All results: \n")
        f.write(model_name +"\n") # model name
        f.write(df.to_string(index=True))
        f.write("\n")

        # Mask when index ranges are specified for this model
        length = len(df['isot_i'])
        try: 
            model_range = models_indexes[i]
            mask_filter = np.full(length, False)
            for period in model_range:
                if period==[]: continue
                start_index, end_index = period[0], period[1]   
                end_index = min(end_index, length - 1)  # Ensure end_index doesn't exceed the bounds of the array 
                mask_filter[start_index: end_index+1] = True # Replace the False with True in those positions in mask
        except: # there are no model_indexes specified
            mask_filter = np.full(length, True)
        
        ## Calculate the average fit parameters -- nH / gamma / Tin
        # To do this, filter the rows where nH is not -1 (i.e. fit was successful), and redchi2 is between 0.8 and 1.2
        # In other words, we only use high-quality fits when determining the average of these parameters
        mask_valid = (df['nH'] != -1) & (df['redchi2'] >= 0.85) & (df['redchi2'] <= 1.15) & (df['counts'] > low_count_threshold)# & (df['cstat?'] == False) 
        
        ## Filter the dataframe
        # Only average when the parameter was not fixed
        mask = mask_filter & mask_valid
        filtered_df = df[mask]
        t = filtered_df['mjd_mid'].to_numpy()
        parameters = ['nH'] + [col for col in ['PhoIndex', 'Tin'] if col in filtered_df.columns] # list of relevant parameter names

        true_indexes = np.where(mask)[0]
        f.write("Indexes used for calculating parameters: \n" + np.array2string(true_indexes, separator=",") + "\n")

        
        for parameter in parameters:
            
            # Calculate the average of the column for the filtered rows
            values = filtered_df[parameter].to_numpy()
            unc_neg = filtered_df[f'{parameter}_neg']
            unc_pos = filtered_df[f'{parameter}_pos']
            
            # Remove the ones where a parameter was fixed
            mask = (unc_neg != 0)
            values = values[mask]
            unc_neg = unc_neg[mask]
            unc_pos = unc_pos[mask]

            if len(values>0):
                avg = values.mean()
                f.write(f"Average {parameter}: {avg}\n") 
            else: f.write(f"Average {parameter} not calculated.\n") 

            # For the overall average, only use the results from powerlaw / diskbb, not the 2-component models as these are tricky to constrain
            if parameter=="nH" and model_name in ['powerlaw', 'pegged_powerlaw', 'diskbb' ]:
                nH_all_valid.append(values)
                nH_all_valid_unc_l.append(unc_neg)
                nH_all_valid_unc_u.append(unc_pos)

            # Do weighted average
            weighted_avg = None
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try: # will only work if some of the parameter values haven't been fixed
                    # Compute the weighted average of nh, using the inverse square of the maximum uncertainty values as the weights
                    weighted_avg = np.average(values, weights = np.amax((unc_neg, unc_pos), axis=0).astype(float) **(-2)) 
                    neg_er = np.sum(unc_neg.astype(float) ** (-2)) ** (-0.5)
                    pos_er = np.sum(unc_pos.astype(float) ** (-2)) ** (-0.5)
                    # (Note: I could also propagate into this the uncertainty due to the variance of the results)
                    f.write(f"Weighted average {parameter}: {weighted_avg} + {pos_er} - {neg_er}\n") 
                except:
                    f.write(f"Weighted average of {parameter} not calculated.\n")

        
        f.write("\n\n") 


    # Overall averages of nH
    f.write("\n\n")
    f.write("\n\nUSING ALL NH RESULTS\n")
    nH_all_valid = np.array([item for sublist in nH_all_valid for item in sublist])
    nH_all_valid_unc_l = np.array([item for sublist in nH_all_valid_unc_l for item in sublist])
    nH_all_valid_unc_u = np.array([item for sublist in nH_all_valid_unc_u for item in sublist])
    if len(nH_all_valid) > 0:
        avg = nH_all_valid.mean()
        f.write(f"Average all nH: {avg}\n") 
    else: f.write(f"Average all nH not calculated.\n") 
    try: # will only work if some of the parameter values haven't been fixed
        # Compute the weighted average of nh, using the inverse square of the maximum uncertainty values as the weights
        weighted_avg = np.average(nH_all_valid, weights = np.amax((nH_all_valid_unc_l, nH_all_valid_unc_u), axis=0).astype(float) **(-2)) 
        print()
        print(weighted_avg)
        neg_er = np.sum(nH_all_valid_unc_l.astype(float) ** (-2)) ** (-0.5)
        pos_er = np.sum(nH_all_valid_unc_u.astype(float) ** (-2)) ** (-0.5)
        # (Note: I could also propagate into this the uncertainty due to the variance of the results)
        f.write(f"Weighted average all nH: {weighted_avg} + {pos_er} - {neg_er}\n") 
    except:
        f.write(f"Weighted average of all nH not calculated.")

    f.close()
    # return tuple(dataframes)  



####################################################################################################################
## PLOT RESULTS


def plot_all_spectral_fit_results():

    models = MODELS

    if any(model not in ['powerlaw', 'pegged_powerlaw', 'diskbb', 'powerlaw+diskbb', 'pegged_powerlaw+diskbb','diskbb+bbodyrad','powerlaw+bbodyrad']  for model in models):
        raise ValueError(f"Error: Invalid model(s) found in array. Allowed values are: {['powerlaw', 'pegged_powerlaw', 'diskbb', 'powerlaw+diskbb', 'pegged_powerlaw+diskbb','diskbb+bbodyrad','powerlaw+bbodyrad'] }")


    for model in models: plot_spectral_results( [model], models_indexes= [] ) # for each model individually
    plot_spectral_results( models, models_indexes= []) # for all the models together




def plot_spectral_results(models = None, models_indexes = None):

    if models is None: models =  MODELS
    if models_indexes is None: models_indexes = MODELS_INDEXES

    fit_with_binning = FIT_WITH_BINNING
    transitions = TRANSITIONS
    xrt_soft_lo, xrt_soft_hi, xrt_hard_lo, xrt_hard_hi = SOFT_LO, SOFT_HI, HARD_LO, HARD_HI

    if models_indexes!=[]: 
        print("Getting final results...")

        # Store final fits in a separate folder, for evaluation
        final_dir = "./spectral_fit_final_results"
        if fit_with_binning: final_dir += "_bin/"
        else: final_dir += "_bin1/"

        if os.path.exists(final_dir):
            shutil.rmtree(final_dir)  
        os.makedirs(final_dir)  


    plots_dir = "./spectral_fit_residuals"
    if fit_with_binning: plots_dir += "_bin/"
    else: plots_dir += "_bin1/"

    mpl.rcParams['xtick.labelbottom'] = False
    colours = ['blue', 'red', 'green', 'purple']

    # Check that none of the defined state ranges overlap
    ranges_list = [np.array(value) for value_list in models_indexes for value in value_list]
    ranges_list = [x for x in ranges_list if x.size > 0] # Remove empty arrays
    if not ranges_okay(ranges_list):
        print(f"The ranges have overlaps.")
        return

    ## Load in the results file
    filename = './spectral_fit_results'
    if fit_with_binning: filename+= "_bin/"
    else: filename+= "_bin1/"
    with open(filename+'xrt_spectral_dict.json', 'r') as j:
        xrt_fit_dict = json.load(j)
    
    # Convert each entry in the nester dictionary to a numpy array
    for key in xrt_fit_dict.keys():
        for keyi in xrt_fit_dict[key].keys():
            xrt_fit_dict[key][keyi] = np.array(xrt_fit_dict[key][keyi])

    
    all_dt_MJD = xrt_fit_dict[models[0]]['dt_mjd'] # all the models have the same MJDs
    all_dates_MJD =  xrt_fit_dict[models[0]]['mjd_mid'] # middle MJD



    # Split into contiguous segments when gap > 60 days
    _gap_days = 60.0
    _mjd_all = all_dates_MJD.copy()
    _diffs = np.diff(_mjd_all)
    _break_idxs = np.where(_diffs > _gap_days)[0]
    _segments = []
    _start_idx = 0
    for _bi in _break_idxs:
        _segments.append((_start_idx, _bi))
        _start_idx = _bi + 1
    _segments.append((_start_idx, len(_mjd_all)-1))
    # _segments now contains (start_index, end_index) tuples in indices of _mjd_all


    if models_indexes!=[]: # Initialise dataframe to store final results
        df = pd.DataFrame(columns=["obs_id", "middle_mjds", "mjd_range", "flux", "flux_er_neg", "flux_er_pos", "uplims", "model", "cstat_bool", "redchi2"])

    
    # Make the figure
    fig, ax = plt.subplots(6, figsize=(30,18), sharex='col', gridspec_kw={'hspace': 0.05})
    fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)


    # HR
    # Filter the array to only include values where the error is low
    lc_mjd, lc_cps, lc_cps_nerr, lc_cps_perr, uplims_bool, obs_type, hr_mjds, hr, hr_err, cps_low, cps_low_err, cps_high, cps_high_err = get_swift_xrt_counts(verbose=False)
    mask= hr_err <= 0.5
    ax[1].errorbar(Time(hr_mjds[mask], format='mjd').datetime, hr[mask], [hr_err[mask], hr_err[mask]], fmt='o', color='k',mfc='black')


    length = len(all_dates_MJD)    
    for i, model in enumerate(models): # make a mask for all the models, just for easier tracking of indexes

        print("Model: ", model)

        mask_valid = xrt_fit_dict[model]['nH']!=-1 # -1 indicates the fit was unsuccessful
        
        # Mask when index ranges are specified
        try: 
            model_range = models_indexes[i]
            mask_filter = np.full(length, False)
            for period in model_range:
                if period==[]: continue
                start_index, end_index = period[0], period[1]   
                end_index = min(end_index, length - 1)  # Ensure end_index doesn't exceed the bounds of the array 
                mask_filter[start_index: end_index+1] = True # Replace the False with True in those positions in mask
        except: # there are no model_indexes specified
            mask_filter = np.full(length, True)
        
        mask = mask_valid & mask_filter
        print("Mask: ", mask)
        print()

        dates_MJD = all_dates_MJD[mask] # middle MJD
        exposures = all_dt_MJD[mask]
        
        data = xrt_fit_dict[model]
        IDs = data['IDs'][mask]
        
        # Flux
        # For the log values, we need to propagate uncertainties: if y = log_10(x), then x_unc = x * ln(10) * y_unc = 10**y * ln(10) * y_unc
        # HOWEVER, if the SNR is low, we cannot use the simple propagation of uncertainties for the log. 
        # In that case, rather use x_unc_l = 10**y - 10**(y-y_unc_l)
        # x_unc_u = 10**(y+y_unc_u) - 10**y 
        if model=='pegged_powerlaw': flux, flux_neg, flux_pos = 1e-12*data['norm'][mask], 1e-12*data['norm_neg'][mask], 1e-12*data['norm_pos'][mask]
        
        
        elif model=='powerlaw' or model=='diskbb' or model=='powerlaw+diskbb' or model=='diskbb+bbodyrad' or model=='powerlaw+bbodyrad':
            #flux, flux_neg, flux_pos = 10**data['lg10Flux'][mask], 10**data['lg10Flux'][mask]*np.log(10)*data['lg10Flux_neg'][mask], 10**data['lg10Flux'][mask]*np.log(10)*data['lg10Flux_neg'][mask]
            flux = 10**data['lg10Flux'][mask]
            flux_neg = 10**data['lg10Flux'][mask] - 10**(data['lg10Flux'][mask] - data['lg10Flux_neg'][mask])
            flux_pos =  10**(data['lg10Flux'][mask] + data['lg10Flux_pos'][mask]) - 10**data['lg10Flux'][mask] 
        
        
        elif model == "pegged_powerlaw+diskbb": 
            #flux, flux_neg, flux_pos = 1e-12*data['norm'][mask] + 10**data['lg10Flux'][mask], np.sqrt ( (1e-12*data['norm_neg'][mask])**2 + (10**data['lg10Flux'][mask]*np.log(10)*data['lg10Flux_neg'][mask])**2 ) , np.sqrt( (1e-12*data['norm_pos'][mask])**2 + (10**data['lg10Flux'][mask]*np.log(10)*data['lg10Flux_neg'][mask])**2 )
            flux = 1e-12*data['norm'][mask] + 10**data['lg10Flux'][mask]
            flux_neg_norm_part = 1e-12*data['norm_neg'][mask]
            flux_neg_log_part = 10**data['lg10Flux'][mask] - 10**(data['lg10Flux'][mask] - data['lg10Flux_neg'][mask])
            flux_neg = np.sqrt(flux_neg_norm_part**2 + flux_neg_log_part**2)
            flux_pos_norm_part =  1e-12*data['norm_pos'][mask]
            flux_pos_log_part = 10**(data['lg10Flux'][mask] + data['lg10Flux_pos'][mask]) - 10**data['lg10Flux'][mask]
            flux_pos = np.sqrt(flux_pos_norm_part**2 + flux_pos_log_part**2)
        
        
        ax[0].errorbar(Time(dates_MJD, format='mjd').datetime, flux, [flux_neg, flux_pos], fmt='o',color='k', mfc=colours[i])



        if models_indexes!=[]:

            df = pd.concat([df, pd.DataFrame({
            "obs_id": IDs,
            "middle_mjds": dates_MJD,
            "mjd_range": exposures,
            "flux": flux,
            "flux_er_neg": flux_neg,
            "flux_er_pos": flux_pos,
            "uplims": np.nan,
            "model": model,
            "binning": data['binned?'][mask],
            "redchi2": data['redchi2'][mask]
            })], ignore_index=True)


            for k, ID in enumerate(IDs):

                search_pattern = os.path.join(plots_dir, f"*{ID}*_{model}_log.png")
                if k==0: print(search_pattern )

                # Find matching files
                matching_files = glob.glob(search_pattern)

                if matching_files:
                    shutil.copy(matching_files[0], final_dir)
                    print(f"Copied {matching_files[0]} to {final_dir}")
                else:
                    print("No matching file found.")

        

        # nH
        nH, nH_neg, nH_pos = data["nH"][mask], data["nH_neg"][mask], data["nH_pos"][mask]
        fixed_mask = (nH_neg == 0) & (nH_pos == 0)  # plot square if parameter was fixed during fitting
        ax[2].errorbar(Time(dates_MJD[fixed_mask], format='mjd').datetime, nH[fixed_mask], [nH_neg[fixed_mask], nH_pos[fixed_mask]], fmt='s',color='k', mfc=colours[i])
        ax[2].errorbar(Time(dates_MJD[~fixed_mask], format='mjd').datetime, nH[~fixed_mask], [nH_neg[~fixed_mask], nH_pos[~fixed_mask]], fmt='o',color='k', mfc=colours[i])



        # Tin; only for models 'diskbb' and 'pegged_powerlaw+diskbb' and 'powerlaw+diskbb'
        if model=="diskbb" or model=="pegged_powerlaw+diskbb" or model=="powerlaw+diskbb" or model=='diskbb+bbodyrad': 
            Tin, Tin_neg, Tin_pos = data['Tin'][mask], data['Tin_neg'][mask], data['Tin_pos'][mask]
            fixed_mask = (Tin_neg == 0) & (Tin_pos == 0)  # plot square if parameter was fixed during fitting
            ax[3].errorbar(Time(dates_MJD[fixed_mask], format='mjd').datetime, Tin[fixed_mask], [Tin_neg[fixed_mask], Tin_pos[fixed_mask]], fmt='s', color='k', mfc=colours[i])
            ax[3].errorbar(Time(dates_MJD[~fixed_mask], format='mjd').datetime, Tin[~fixed_mask], [Tin_neg[~fixed_mask], Tin_pos[~fixed_mask]], fmt='o', color='k', mfc=colours[i])

        # Gamma; only for models 'pegged_powerlaw', 'powerlaw', 'powerlaw+diskbb', 'pegged_powerlaw+diskbb'
        if model=="powerlaw" or model=="pegged_powerlaw" or model=="powerlaw+diskbb" or model=="pegged_powerlaw+diskbb" or model=='powerlaw+bbodyrad': 
            gamma, gamma_neg, gamma_pos = data['PhoIndex'][mask], data['PhoIndex_neg'][mask], data['PhoIndex_pos'][mask]
            fixed_mask = (gamma_neg == 0) & (gamma_pos == 0) # plot square if parameter was fixed during fitting
            ax[4].errorbar(Time(dates_MJD[fixed_mask], format='mjd').datetime, gamma[fixed_mask], [gamma_neg[fixed_mask], gamma_pos[fixed_mask]], fmt='s', color='k', mfc=colours[i])
            ax[4].errorbar(Time(dates_MJD[~fixed_mask], format='mjd').datetime, gamma[~fixed_mask], [gamma_neg[~fixed_mask], gamma_pos[~fixed_mask]], fmt='o', color='k', mfc=colours[i])

        # Chi^2
        if models_indexes==[]: mask = np.full(length, True) # show all chi^2 results, even when it is above the threshold for the fit and error calculation to be considered successful
        #mask_no_cstat = xrt_fit_dict[model]['cstat?'] == False # exclude cstat points
        #mask = mask & mask_no_cstat
        chi, dates_MJD = xrt_fit_dict[model]['redchi2'][mask], all_dates_MJD[mask]
        ax[5].errorbar(Time(dates_MJD, format='mjd').datetime, chi, 0.0, fmt='o', color='k', mfc=colours[i], label=model)



        # Set plot constraints
        ax[0].set_yscale('log')
        ax[0].set_ylabel('Flux [1$-$10 keV]\n(erg s$^{-1}$ cm$^{-1}$)')
        ax[1].set_ylabel(rf'HR $\left(\frac{{[{xrt_hard_lo}-{xrt_hard_hi}\,\text{{keV}}]}}{{[{xrt_soft_lo}-{xrt_soft_hi}\,\text{{keV}}]}}\right)$')
        ax[2].set_ylabel(r'$n_\text{H}$ ($\times 10^{22}$)')
        ax[3].set_ylabel(r'$k_B T_\text{in}$ (keV)')
        ax[4].set_ylabel('$\Gamma$')
        ax[5].set_ylabel(r'$\chi^2_\text{red}$')
        ax[5].set_yscale('log')
        ax[5].legend(fontsize=11)

        for i in range(6):
            if transitions is not None:
                for t in transitions: # transition points
                    ax[i].axvline(Time(t, format='mjd').datetime, color='yellow', linestyle='--', linewidth=1.5)



    top_ax = ax[0]    
    ax_main = ax[-1]      
    orig_xlim = ax_main.get_xlim()
    orig_ylim = ax_main.get_ylim()
    orig_top_locator = top_ax.xaxis.get_major_locator()
    orig_top_formatter = top_ax.xaxis.get_major_formatter()
    orig_top_xlabel = top_ax.get_xlabel()
    orig_axes = list(fig.axes)   # snapshot of axes that exist now
    # Save separate images for each contiguous segment 
    for seg_i, (sidx, eidx) in enumerate(_segments, start=1):
        
        seg_mjd_min = _mjd_all[sidx]
        seg_mjd_max = _mjd_all[eidx]
        seg_xmin = seg_mjd_min - 5.0
        seg_xmax = seg_mjd_max + 5.0

        # convert to datetime for set_xlim (your plots use Time(...).datetime)
        seg_xlim_dt = Time([seg_xmin, seg_xmax], format='mjd').datetime

        ax[-1].set_xlim(seg_xlim_dt[0], seg_xlim_dt[1])

        # Optionally adjust formatting ticks for this narrower range:
        _all_dates_seg = _mjd_all[(_mjd_all >= seg_mjd_min) & (_mjd_all <= seg_mjd_max)]
        if _all_dates_seg.size > 1:
            T_seg = _all_dates_seg[-1] - _all_dates_seg[0]
            dt_seg = int(max(1, T_seg // 4))
        else:
            dt_seg = 1
        FormatAxis(ax, _all_dates_seg, interval=dt_seg)


        # Save segment image
        name = "_".join(models)
        if models_indexes!=[]: seg_name = f"{final_dir}final_fit_selection_{name}_segment{seg_i}.png"
        else: seg_name = f"{filename}all_fits_{name}_segment{seg_i}.png"
        print(f"Saving segment {seg_i}: MJD {seg_mjd_min:.2f}--{seg_mjd_max:.2f} -> {seg_name}")
        plt.savefig(seg_name)

        sec = getattr(ax_main, "_mjd_secondary_axis", None)
        if sec is not None:
            try:
                sec.remove()
            except Exception:
                pass
            ax_main._mjd_secondary_axis = None

        # restore top axis locator/formatter etc if you saved them earlier
        top_ax.xaxis.set_major_locator(orig_top_locator)
        top_ax.xaxis.set_major_formatter(orig_top_formatter)
        top_ax.set_xlabel(orig_top_xlabel)
        top_ax.tick_params(axis='x', labeltop=True, labelbottom=False)

        # restore main axis limits or autoscale back
        ax_main.set_xlim(orig_xlim)
        ax_main.relim()
        ax_main.autoscale_view()


    
    if models_indexes!=[]: # i.e. getting final results

        # Tabulate the results
        df = df.sort_values(by="middle_mjds").reset_index(drop=True)
        with open(final_dir + "final_fit_selection.txt", "w") as f:
            f.write(df.to_string(index=False))  
            f.write("\n")
        # Save as a CSV file
        df.to_csv(final_dir + "final_fit_selection.csv", index=False)
        



####################################################################################################################

