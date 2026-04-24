from fileinput import filename
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
import scipy.stats 
import scipy.special 
#from scipy.special import iv
from astropy.constants import c
from astropy.table import Table
from astropy.io import fits
from scipy import stats
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
from astropy.time import Time
from pathlib import Path

from plotting_helpers import *


import sys
import os
sys.path.append(os.path.abspath("../"))
from input_parameters import QS, HS, SS, IMS, SOFT_LO, SOFT_HI, HARD_LO, HARD_HI, MIN_E_LC


####################################################################################################################
## GET DATA FUNCTIONS
####################################################################################################################


colours_states = {
    "QS":  "tab:blue",
    "HS":  "tab:orange",
    "SS":  "tab:green",
    "IMS": "tab:purple",
}




####################################################################################################################
## GET MEERKAT RESULTS


def read_radio_file():
    """
    Reads the radio text file using the header line (starting with '#').
    Returns cleaned numpy arrays with NaN rows removed.
    """

    filename = Path(__file__).resolve().parent.parent / "other_lightcurve_data" / "meerkat.txt"


    # Read header line manually
    with open(filename, "r") as f:
        header_line = f.readline()
        for line in f:
            if line.startswith("#"):
                header_line = line
                break

    if header_line is None:
        raise ValueError("Header line starting with '#' not found.")

    columns = [col.strip() for col in header_line[1:].split(",")]

    df = pd.read_csv(
        filename,
        comment="#",
        names=columns,
        usecols=range(len(columns)),
        skip_blank_lines=True
    )

    # Replace empty strings with NaN
    df = df.replace(r'^\s*$', pd.NA, regex=True)

    # Drop rows containing NaNs in relevant columns
    df = df.dropna(subset=['mid_mjd'])

    # Create uplim flag (True if 3sigma_uplim exists and is not NaN)
    df['uplim_bool'] = df['3sigma_uplim'].notna()
    # Where uplim exists: move value to Fr and set Fr_unc = Fr/3
    uplim_mask = df['uplim_bool']
    df.loc[uplim_mask, 'Fr'] = df.loc[uplim_mask, '3sigma_uplim']
    df.loc[uplim_mask, 'Fr_unc'] = df.loc[uplim_mask, 'Fr'] / 3.0

    # Sort by MJD
    df = df.sort_values('mid_mjd').reset_index(drop=True)

    return (
        df['mid_mjd'].to_numpy(),
        df['Fr'].to_numpy(),
        df['Fr_unc'].to_numpy(),
        df['uplim_bool'].to_numpy()
    )




####################################################################################################################
## GET SWIFT/XRT RESULTS


def get_incbad_ids():
    """ 
    Prints which observations were labelled as "bad" from the light curve generation, menaing no centroid could be found. 
    Note that the light curve generation excludes observations shorter than 7 seconds, while the spectal generator includes 
    them, so these observations (if any) are missing from the lists.
    """

    all_ids = []
    all_obs_type = []
    incbad_ids = []
    incbad_obs_type = []

    file = Path(__file__).resolve().parent.parent / "lightcurves_swift_xrt" 


    # The routine creates separate sub-folders for each ID
    for lc_file in glob.glob(f'{file}/00*'):

        for wt_pc in ["wt", "pc"]:
            pc_wt_cap = wt_pc.upper()
            
            try: lc_name = glob.glob(f'{lc_file}/{pc_wt_cap}.qdp')[0]
            except: continue # no WC/PC data for this ID, so skip

            lc_name_incbad = None
            try: lc_name_incbad = glob.glob(f'{lc_file}/{pc_wt_cap}_incbad.qdp')[0]
            except: pass

            clean_lines = []
            with open(lc_name) as f:
                for line in f:

                    if line.startswith('READ'):
                        continue
                    if line.startswith('!'):
                        continue
                    if line.startswith('NO'):
                        continue
                    
                    # Remove ObsID part
                    clean_line = line.split('::')[0]
                    clean_lines.append(clean_line)
                    obsid = line.partition('ObsID=')[2].strip()
                    id = obsid + wt_pc
                    all_ids.append(id)
                    all_obs_type.append(wt_pc)

            if lc_name_incbad is not None:
                with open(lc_name_incbad) as f:
                    for line in f:

                        if line.startswith('READ'):
                            continue
                        if line.startswith('!'):
                            continue
                        if line.startswith('NO'):
                            continue
                        
                        # Remove ObsID part
                        clean_line = line.split('::')[0]
                        clean_lines.append(clean_line)
                        obsid = line.partition('ObsID=')[2].strip()
                        id = obsid + wt_pc
                        incbad_ids.append(id)
                        incbad_obs_type.append(wt_pc)

    # Determine which IDs are in incbad but not in the main list (i.e. which IDs are only in incbad)
    incbad_only_ids = set(incbad_ids) - set(all_ids)
    incbad_only_obs_type = [incbad_obs_type[incbad_ids.index(id)] for id in incbad_only_ids]

    print("IDs in incbad but not in main list: \n" + np.array2string(np.array(list(incbad_only_ids)), separator=","))
    print("Observation types for IDs in incbad but not in main list: \n" + np.array2string(np.array(incbad_only_obs_type), separator=","))

    #return np.array(all_ids), np.array(all_obs_type), np.array(incbad_ids), np.array(incbad_obs_type)




def get_swift_xrt_counts(verbose=True, incbad = True, just_ids = False):
    """
    Retrieve Swift XRT light curve data.

    Note that the IDs are only in the WT.qdp and PC.qdp files, not in the curve.qdp files.
    """
 

    file = Path(__file__).resolve().parent.parent / "lightcurves_swift_xrt" 


    # Initialise arrays
    lc_mjd      = [] # date
    lc_cps      = [] # count rate
    lc_cps_nerr = [] # count rate lower error
    lc_cps_perr = [] # count rate upper error
    obs_type = []
    ids = []


    # If downloaded via the API, the different IDs will be in different folders
    # If copied from the website, all the results will be in one folder
    sub_files = glob.glob(f'{file}/00*')
    if not sub_files:
        sub_files = [file]

    # The routine creates separate sub-folders for each ID
    for lc_file in sub_files:

        """"
        if incbad: 
            try: lc_name = glob.glob(f'{lc_file}/curve_incbad.qdp')[0]
            except: lc_name = glob.glob(f'{lc_file}/curve.qdp')[0]
        else: lc_name = glob.glob(f'{lc_file}/curve.qdp')[0]


        headers = []
        with open(lc_name, 'r') as f:
            for line in f: 
                if line.startswith('! '): headers.append(line[1:].strip()) 

        # Data is stored in different table_id values in the .qdp table, so iterate through this
        i = 0
        while True:
            try: 
                data = Table.read(lc_name, format='ascii.qdp', table_id=i)
                lc_mjd = np.append(lc_mjd, data['col1'].data) # creates a copy of the array, so we need to re-assign it to the original variable
                lc_mjd_exp = np.append(lc_mjd_exp, abs(data['col1_perr'].data) + abs(data['col1_nerr'].data))
                lc_cps = np.append(lc_cps, data['col2'].data)
                lc_cps_nerr = np.append(lc_cps_nerr, abs(data['col2_nerr'].data))
                lc_cps_perr = np.append(lc_cps_perr, abs(data['col2_perr'].data))
                obs_type.extend([headers[i]] * len(data))
                i+=1
            except IndexError:
                break
        """

        for wt_pc in ["wt", "pc"]:
            pc_wt_cap = wt_pc.upper()
            
    
            try: lc_name = glob.glob(f'{lc_file}/{pc_wt_cap}.qdp')[0]
            except: continue

            if incbad: 
                try: lc_name = glob.glob(f'{lc_file}/{pc_wt_cap}_incbad.qdp')[0]
                except: pass

            clean_lines = []
            with open(lc_name) as f:
                for line in f:

                    if line.startswith('READ'):
                        continue
                    if line.startswith('!'):
                        continue
                    if line.startswith('NO'):
                        continue
                    
                    # Remove ObsID part
                    clean_line = line.split('::')[0]
                    clean_lines.append(clean_line)
                    obsid = line.partition('ObsID=')[2].strip()
                    id = obsid + wt_pc
                    ids.append(id)
                    obs_type.append(wt_pc)

                    split_line = clean_line.split()
                    mjd = float(split_line[0])
                    lc_mjd.append(mjd)
                    cps = float(split_line[3])
                    lc_cps.append(cps)
                    cps_perr = float(split_line[4])
                    lc_cps_perr.append(cps_perr)
                    cps_nerr = abs(float(split_line[5]))
                    lc_cps_nerr.append(cps_nerr)



    # Order light curve arrays in time, converting first to numpy arrays
    index = np.argsort(lc_mjd)
    lc_mjd = np.array(lc_mjd)[index]
    lc_cps = np.array(lc_cps)[index]
    lc_cps_nerr = np.array(lc_cps_nerr)[index]
    lc_cps_perr = np.array(lc_cps_perr)[index]
    obs_type = np.array(obs_type)[index]
    ids = np.array(ids)[index]

    ## Find upper limits
    # These are located at the indexes where the count error is zero
    index_uplims = np.where(lc_cps_nerr == 0.0)
    time_uplims = lc_mjd[index_uplims]
    time_uplims_utc = mjd2utc(time_uplims)
    cps_uplims = lc_cps[index_uplims]
    obs_type_uplims = obs_type[index_uplims]
    ids_uplims = ids[index_uplims]

    uplims_bool = np.zeros_like(lc_cps_nerr, dtype=bool)
    uplims_bool[index_uplims] = True

  
    if len(time_uplims)!=0 and verbose: 
        print()
        print("Only upper limits were obtained for some observations:")
        print("IDs for upper limits: \n" + np.array2string(ids_uplims, separator=","))
        print("MJDs for upper limits: \n" +  np.array2string(time_uplims, separator=","))
        print("UTC for upper limits: \n" + np.array2string(time_uplims_utc, separator=","))
        print("Observation type for the upper limits: \n" + np.array2string(obs_type_uplims, separator=","))
        print("Count rates [counts/s] for the upper limits: \n" + np.array2string(cps_uplims, separator=","))
        

    if just_ids: 
        #index_det = np.where(lc_cps_nerr != 0.0)
        #detected_ids = ids[index_det]
        #return detected_ids
        return ids_uplims
    

    ###########################
    ## Get:
    # - hardness ratios (HR)
    # - 0.2-2 keV and 2-10 keV count rates 


    # Initialise arrays
    hr_mjds  = np.array([])
    hr      = np.array([])
    hr_err = np.array([])

    mjd_high  = np.array([])
    cps_high      = np.array([])
    cps_high_err = np.array([])

    mjd_low  = np.array([])
    cps_low      = np.array([])
    cps_low_err = np.array([])

    
    for lc_file in sub_files:

        hr_name = glob.glob(f'{lc_file}/hardrat.qdp')[0] 

        pattern = re.compile(r"\s*!::.*$")
        with open(hr_name) as fin, tempfile.NamedTemporaryFile(
            mode="w+", delete=False
        ) as fout:
            for line in fin:
                fout.write(pattern.sub("", line))
            clean_name = fout.name
    
        # The hardness ratio is only in every third table -- the other two are the hard-band data and soft-band data 
        for i in [2,5]:
            try: data = Table.read(clean_name, format='ascii.qdp', table_id=i)
            except: continue
            hr_mjds = np.append(hr_mjds, data['col1'].data)
            hr = np.append(hr, data['col2'].data)
            hr_err= np.append(hr_err, abs(data['col2_err'].data))


        for i in [0,3]:
            try: data = Table.read(clean_name, format='ascii.qdp', table_id=i)
            except: continue
            mjd_high = np.append(mjd_high, data['col1'].data)
            cps_high = np.append(cps_high, data['col2'].data)
            cps_high_err = np.append(cps_high_err, abs(data['col2_err'].data))

        for i in [1,4]:
            try: data = Table.read(clean_name, format='ascii.qdp', table_id=i)
            except: continue
            mjd_low = np.append(mjd_low, data['col1'].data)
            cps_low = np.append(cps_low, data['col2'].data)
            cps_low_err = np.append(cps_low_err, abs(data['col2_err'].data))



    # Order HR arrays in time
    index = np.argsort(hr_mjds)
    hr_mjds = hr_mjds[index]
    hr = hr[index]
    hr_err = hr_err[index]

    # Order low and high arrays in time
    index = np.argsort(mjd_high)
    mjd_high = mjd_high[index]
    cps_high = cps_high[index]
    cps_high_err = cps_high_err[index]
    index = np.argsort(mjd_low)
    mjd_low = mjd_low[index]
    cps_low = cps_low[index]
    cps_low_err = cps_low_err[index]


    return lc_mjd, lc_cps, lc_cps_nerr, lc_cps_perr, uplims_bool, obs_type, hr_mjds, hr, hr_err, cps_low, cps_low_err, cps_high, cps_high_err



####################################################################################################################
## GET MAXI RESULTS


def get_maxi_counts():

    path = Path(__file__).resolve().parent.parent / "other_lightcurve_data" / "maxi_daily_ave.txt"
    

    # Format: MJDcenter, 2-20keV[ph/s/cm2], err, 2-4keV, err, 4-10keV, err, 10-20keV, err 
    df = pd.read_csv(path, delim_whitespace=True, header=None)

    # Extract columns
    mjd = df[0].to_numpy()
    count_rate_2_20 = df[1].to_numpy() # 2-20 keV
    count_rate_2_20_err = df[2].to_numpy() 
    count_rate_2_4 = df[3].to_numpy() # 2-4 keV ... low
    count_rate_2_4_err = df[4].to_numpy()
    count_rate_4_10 = df[5].to_numpy() # 4-10 keV ... high
    count_rate_4_10_err = df[6].to_numpy()
    count_rate_10_20 = df[7].to_numpy() # 10-20 keV ... high
    count_rate_10_20_err = df[8].to_numpy()

    # Mask where rate is zero or negative
    valid_mask = (count_rate_2_4 > 0) & (count_rate_4_10 > 0) & (count_rate_2_20 > 0) & (count_rate_10_20 > 0)

    # Apply the mask to all arrays
    mjd = mjd[valid_mask]
    count_rate_2_4 = count_rate_2_4[valid_mask]
    count_rate_2_4_err = count_rate_2_4_err[valid_mask]
    count_rate_4_10= count_rate_4_10[valid_mask]
    count_rate_4_10_err = count_rate_4_10_err[valid_mask]
    count_rate_2_20  = count_rate_2_20 [valid_mask]
    count_rate_2_20_err = count_rate_2_20_err[valid_mask]
    count_rate_10_20 = count_rate_10_20[valid_mask]
    count_rate_10_20_err = count_rate_10_20_err[valid_mask]


    # Compute the ratio and propagate uncertainties
    hardness = count_rate_4_10 / count_rate_2_4
    hardness_err = hardness * np.sqrt((count_rate_4_10_err / count_rate_4_10)**2 + (count_rate_2_4_err / count_rate_2_4)**2)

    # Compute log hardness
    #log_hardness = np.log10(hardness)
    # Error propagation for log10(hardness)
    # d(log(x/y)) = 1/ln(10) * sqrt[(err_x / x)^2 + (err_y / y)^2] = 1/ln(10) * hardness_err/hardness
    log_hardness_err = (1 / np.log(10)) * np.sqrt((count_rate_4_10_err / count_rate_4_10) ** 2 + (count_rate_2_4_err / count_rate_2_4) ** 2)


    # Filter where error is too large
    #error_mask = np.ones_like(hardness_err, dtype=bool)  # No filtering for now
    error_mask = log_hardness_err < 0.2

    mjd = mjd[error_mask]
    hardness = hardness[error_mask]
    hardness_err = hardness_err[error_mask]
    count_rate_2_4 = count_rate_2_4[error_mask]
    count_rate_2_4_err = count_rate_2_4_err[error_mask]
    count_rate_4_10 = count_rate_4_10[error_mask]
    count_rate_4_10_err = count_rate_4_10_err[error_mask]
    count_rate_2_20 = count_rate_2_20[error_mask]
    count_rate_2_20_err = count_rate_2_20_err[error_mask]
    count_rate_10_20 = count_rate_10_20[error_mask]
    count_rate_10_20_err = count_rate_10_20_err[error_mask]


    return mjd, hardness,hardness_err , count_rate_2_4 ,count_rate_2_4_err ,count_rate_4_10 , count_rate_4_10_err ,count_rate_2_20 , count_rate_2_20_err , count_rate_10_20 , count_rate_10_20_err







####################################################################################################################
## GET SWIFT/BAT RESULTS



def get_bat_counts():


    filename = Path(__file__).resolve().parent.parent / "other_lightcurve_data" / "bat_daily_ave.txt"
    #print("BAT filename: " + str(filename))

    data = np.loadtxt(filename, comments="#")

    # Columns (based on your file structure)
    mjd = data[:, 0]
    rate = data[:, 1]
    error = data[:, 2]

    # Remove negative count rates and when error is too large
    mask = (rate >= 0) & (error < 0.8*rate)  # Example: keep points where error is less than 80% of the rate
    mjd = mjd[mask]
    rate = rate[mask]
    error = error[mask]

    return mjd, rate, error





####################################################################################################################
## COUNT RATE LIGHT CURVE PLOTTING FUNCTIONS
###################################################################################################################



def plot_count_rates_and_hr(gap_days=60):
    """
    Uses MeerKAT time coverage to define the plotting window(s). If MeerKAT has a gap > gap_days,
    creates one separate figure per contiguous MeerKAT segment. Uses only raw data (no interpolation).
    """
    mpl.rcParams['xtick.labelbottom'] = False


    qs, hs, ss, ims = QS, HS, SS, IMS
    xrt_soft_lo, xrt_soft_hi, xrt_hard_lo, xrt_hard_hi, xrt_min_E_lc = SOFT_LO, SOFT_HI, HARD_LO, HARD_HI, MIN_E_LC

    # -------------------------
    # Load Swift/XRT
    # -------------------------
    
    # Swift/XRT
    xrt_mjds, xrt_cps, xrt_cps_nunc, xrt_cps_punc, xrt_uplims_bool, xrt_obs_type, xrt_hr_mjds, xrt_hr, xrt_hr_unc, xrt_low_cps, xrt_low_cps_unc, xrt_high_cps, xrt_high_cps_unc = get_swift_xrt_counts()



    # -------------------------
    # Load MeerKAT data -- use this to define the plotting windows, if it exists
    # -------------------------
    
    segments = [] # contains tuples of (start_idx, end_idx) based on either MeerKAT or XRT time axis.
    incl_mkt = True
    try: 
        mkt_mjds, mkt_Fr, mkt_Fr_unc, mkt_uplims_bool = read_radio_file() # This data is already sorted in time
        mjds_for_segments = mkt_mjds
    except:
        incl_mkt = False
        mjds_for_segments = xrt_mjds
        print("No MeerKAT data")
    
    # Define segments based on gaps in the data
    # find breaks > gap_days
    diffs = np.diff(mjds_for_segments)
    break_idxs = np.where(diffs > gap_days)[0]
    start = 0
    for bi in break_idxs:
        segments.append((start, bi))
        start = bi + 1
    segments.append((start, len(mjds_for_segments) - 1))


    # -------------------------
    # Load other X-ray data
    # -------------------------

    incl_maxi = True
    try: 
        # MAXI 
        (maxi_mjds, maxi_hr, maxi_hr_unc,
        maxi_2_4_cps, maxi_2_4_cps_unc,
        maxi_4_10_cps, maxi_4_10_cps_unc,
        maxi_2_20_cps, maxi_2_20_cps_unc,
        maxi_10_20_cps, maxi_10_20_cps_unc) = get_maxi_counts()
    except: 
        incl_maxi = False
        print("No MAXI data")
    
    incl_bat = True
    try: 
        # BAT
        bat_mjds, bat_cps, bat_unc = get_bat_counts()
    except: 
        incl_bat = False
        print("No BAT data")



    # -------------------------
    # Prepare output dir
    # -------------------------
    plots_dir = "./final_lightcurves/"
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)
        print(f"Created {plots_dir}")


    # -------------------------
    # Loop segments and plot
    # -------------------------
    for part, (s, e) in enumerate(segments, start=1):

        seg_mjds = mjds_for_segments[s:e + 1]
        seg_tmin = seg_mjds.min() - 20 
        seg_tmax = seg_mjds.max() + 20
        seg_t = np.linspace(seg_tmin, seg_tmax, 10)

        print(f"Processing segment {part}: {seg_tmin} to {seg_tmax}")
    

        # Masks selecting data within this window
        xrt_mask = (xrt_mjds >= seg_tmin) & (xrt_mjds <= seg_tmax)
        xrt_hr_mask = (xrt_hr_mjds >= seg_tmin) & (xrt_hr_mjds <= seg_tmax)
        if incl_bat: bat_mask = (bat_mjds >= seg_tmin) & (bat_mjds <= seg_tmax)
        if incl_maxi: maxi_mask = (maxi_mjds >= seg_tmin) & (maxi_mjds <= seg_tmax)
        if incl_mkt: mkt_mask = (mkt_mjds >= seg_tmin) & (mkt_mjds <= seg_tmax)


        # Figure with stacked panels
        n_panels = 10
        fig_height = 2.8 * n_panels
        fig, axes = plt.subplots(n_panels, 1, figsize=(14, fig_height), sharex=True, gridspec_kw={'hspace': 0.05, 'height_ratios': [1] * n_panels})
        fig.set_facecolor('white')


        # ---------- Panel 0: BAT 15-50 ----------
        if incl_bat:
            ax = axes[0]
            if np.any(bat_mask):
                ax.errorbar(Time(bat_mjds[bat_mask], format='mjd').datetime, bat_cps[bat_mask],
                            yerr=bat_unc[bat_mask], fmt='o', ms=5, color='orange')
                ax.set_yscale('log')
            ax.set_ylabel('BAT\n15-50 keV')
            #ax.grid(True, which='both', ls=':')


        # ---------- Panel 1: MAXI 10-20 ----------
        if incl_maxi:
            ax = axes[1]
            if np.any(maxi_mask):
                valid = ~np.isnan(maxi_10_20_cps[maxi_mask])
                if np.any(valid):
                    ax.errorbar(Time(maxi_mjds[maxi_mask][valid], format='mjd').datetime,
                                maxi_10_20_cps[maxi_mask][valid],
                                yerr=maxi_10_20_cps_unc[maxi_mask][valid],
                                fmt='o', ms=5, color='blue')
                    ax.set_yscale('log')
            ax.set_ylabel('MAXI\n10-20 keV')
            #ax.grid(True, which='both', ls=':')

        # ---------- Panel 2: MAXI 4-10 ----------
        if incl_maxi:
            ax = axes[2]
            if np.any(maxi_mask):
                valid = ~np.isnan(maxi_4_10_cps[maxi_mask])
                if np.any(valid):
                    ax.errorbar(Time(maxi_mjds[maxi_mask][valid], format='mjd').datetime,
                                maxi_4_10_cps[maxi_mask][valid],
                                yerr=maxi_4_10_cps_unc[maxi_mask][valid],
                                fmt='o', ms=5, color='blue')
                    ax.set_yscale('log')
            ax.set_ylabel('MAXI\n4-10 keV')
            #ax.grid(True, which='both', ls=':')

        # ---------- Panel 3: MAXI 2-4 ----------
        if incl_maxi:
            ax = axes[3]
            if np.any(maxi_mask):
                ax.errorbar(Time(maxi_mjds[maxi_mask], format='mjd').datetime, maxi_2_4_cps[maxi_mask],
                            yerr=maxi_2_4_cps_unc[maxi_mask], fmt='o', ms=5, color='blue')
                ax.set_yscale('log')
            ax.set_ylabel('MAXI\n2-4 keV')
        #ax.grid(True, which='both', ls=':')

        # ---------- Panel 4: MAXI HR ----------
        if incl_maxi:
            ax = axes[4]
            if np.any(maxi_mask):
                ax.errorbar(Time(maxi_mjds[maxi_mask], format='mjd').datetime, maxi_hr[maxi_mask],
                            yerr=maxi_hr_unc[maxi_mask], fmt='o', ms=5, color="blue")
            ax.set_ylabel('MAXI\nHR\n(4-10 / 2-4)')
            ymin = max(0.0, np.min(maxi_hr[maxi_mask] - maxi_hr_unc[maxi_mask]))
            ymax = min(2.0, np.max(maxi_hr[maxi_mask] + maxi_hr_unc[maxi_mask]))
        ax.set_ylim(ymin, ymax)
        #ax.grid(True, which='both', ls=':') 


        # ---------- Panel 5: XRT hard ----------
        ax = axes[5]
        if np.any(xrt_hr_mask):
            ax.errorbar(Time(xrt_hr_mjds[xrt_hr_mask], format='mjd').datetime,
                        xrt_high_cps[xrt_hr_mask],
                        yerr=[xrt_high_cps_unc[xrt_hr_mask], xrt_high_cps_unc[xrt_hr_mask]],
                        fmt='o', ms=5, color='red')
            ax.set_yscale('log')
        ax.set_ylabel(f'XRT\n({xrt_hard_lo}-{xrt_hard_hi}) keV c/s')
        #ax.grid(True, which='both', ls=':')


        # ---------- Panel 6: XRT soft ----------
        ax = axes[6]
        if np.any(xrt_hr_mask):
            ax.errorbar(Time(xrt_hr_mjds[xrt_hr_mask], format='mjd').datetime,
                        xrt_low_cps[xrt_hr_mask],
                        yerr=[xrt_low_cps_unc[xrt_hr_mask], xrt_low_cps_unc[xrt_hr_mask]], 
                        fmt='o', ms=5, color='red')
            ax.set_yscale('log')
        ax.set_ylabel(f'XRT\n({xrt_soft_lo}-{xrt_soft_hi}) keV c/s')
        #ax.grid(True, which='both', ls=':')


        # ---------- Panel 7: XRT ----------
        # Colour by obs type (WT/PC) 
        ax = axes[7]
        if np.any(xrt_mask):
            # Colour by obs type (WT/PC) if possible
            xrt_obs_type_mask = xrt_obs_type[xrt_mask]
            xrt_cps_mask = xrt_cps[xrt_mask]
            xrt_cps_nunc_mask = xrt_cps_nunc[xrt_mask]
            xrt_cps_punc_mask = xrt_cps_punc[xrt_mask]
            xrt_uplims_bool_mask = xrt_uplims_bool[xrt_mask]

            # Create a mask for each observation type
            wt_mask = (xrt_obs_type_mask == 'wt')
            pc_mask = (xrt_obs_type_mask == 'pc')

            # Plot WT points in blue
            if np.any(wt_mask):
                ax.errorbar(Time(xrt_mjds[xrt_mask][wt_mask], format='mjd').datetime,
                            xrt_cps_mask[wt_mask],
                            yerr=[xrt_cps_nunc_mask[wt_mask], xrt_cps_punc_mask[wt_mask]], uplims = xrt_uplims_bool_mask[wt_mask], 
                            fmt='o', ms=5, color='red', label='WT')

            # Plot PC points in red
            if np.any(pc_mask):
                ax.errorbar(Time(xrt_mjds[xrt_mask][pc_mask], format='mjd').datetime,
                            xrt_cps_mask[pc_mask],
                            yerr=[xrt_cps_nunc_mask[pc_mask], xrt_cps_punc_mask[pc_mask]], uplims = xrt_uplims_bool_mask[pc_mask], 
                            fmt='o', ms=5, color='red', alpha = 0.5, label='PC')

            ax.set_yscale('log')
        ax.legend()
        ax.set_ylabel(f'XRT\n({xrt_min_E_lc}-10) keV c/s')
        #ax.grid(True, which='both', ls=':')


        # ---------- Panel 8: XRT HR ----------
        ax = axes[8]
        if np.any(xrt_hr_mask):
            ax.errorbar(Time(xrt_hr_mjds[xrt_hr_mask], format='mjd').datetime, xrt_hr[xrt_hr_mask],
                        yerr=xrt_hr_unc[xrt_hr_mask], fmt='o', ms=5, color='red')
        ax.set_ylabel(f'XRT\nHR\n({xrt_hard_lo}-{xrt_hard_hi} / {xrt_soft_lo}-{xrt_soft_hi})')
        ymin = max(0.0, np.min(xrt_hr[xrt_hr_mask] - xrt_hr_unc[xrt_hr_mask]))
        ymax = min(2.0, np.max(xrt_hr[xrt_hr_mask] + xrt_hr_unc[xrt_hr_mask]))
        ax.set_ylim(ymin, ymax)
        #ax.grid(True, which='both', ls=':')


        # ---------- Panel 9: MeerKAT (radio) ----------
        if incl_mkt:
            ax = axes[9]
            ax.errorbar(Time(mkt_mjds[mkt_mask], format='mjd').datetime, mkt_Fr[mkt_mask], yerr=mkt_Fr_unc[mkt_mask], fmt='o', ms=5, color='green', uplims= mkt_uplims_bool[mkt_mask])
            # choose linear for radio unless all >0 and user prefers log
            if np.all(mkt_Fr[mkt_mask] > 0):
                ax.set_yscale('log')
            ax.set_ylabel('MeerKAT\n(mJy)')
            #ax.grid(True, which='both', ls=':')



        # ---------- vertical transition lines if present ----------
        for i, ax in enumerate(axes):
            for (xk, label) in [(qs, "QS"), (hs, "HS"), (ss, "SS"), (ims, "IMS")]:
                    for k, (x0, x1) in enumerate(xk):
                        ax.axvspan(
                            Time(x0, format='mjd').datetime, Time(x1, format='mjd').datetime,
                            color=colours_states[label],
                            alpha=0.2,
                            label=label if (k == 0 and i == len(axes) - 1)  else None  # only label once for legend
                        )
     

        axes[9].legend(fontsize=11)
     
        # ---------- X-axis formatting ----------
        plt.xlim([seg_tmin, seg_tmax ])  
        # use the FormatAxis helper you already have; pick interval from duration
        T_seg = seg_tmax - seg_tmin
        dt = int(max(10, T_seg // 4))
        # FormatAxis was used previously with (ax, mjd_array, interval=...), so call similarly:
        FormatAxis(axes, seg_t, interval=dt)
        

        outname = os.path.join(plots_dir, f'xray_cps_and_mkt_part_{part}.png')
        plt.savefig(outname, bbox_inches='tight', dpi=150)
        print(f"Saved {outname} (MeerKAT MJD range: {seg_tmin:.2f} - {seg_tmax:.2f})")
        plt.close(fig)


