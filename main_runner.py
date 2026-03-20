# Before running, activate your venv environment:
# e.g., source /mnt/users/crookmansourj/swiftenv/bin/activate
# This environment has swiftools 3.0.22

import sys
import os
sys.path.append(os.path.abspath("./CODE/"))

from get_xrt_from_pipeline import get_xrt_prods_runner, group_spectra
from get_results_data import plot_count_rates_and_hr
from fit_xrt_spectra import run_spectral_fit
from get_results_spec import plot_all_spectral_fit_results, get_results, plot_spectral_results
from input_parameters import ANALYSIS_TYPE


if __name__ in "__main__": 

    print("Analysis type: ",  ANALYSIS_TYPE)


    ## STEP 1: GET LIGHTCURVES & SPECTRA
    if ANALYSIS_TYPE=="get_swift_xrt_data": 
        get_xrt_prods_runner()
        

    ## STEP 2: RE-BIN THE SPECTRA
    elif ANALYSIS_TYPE=="group_spectra":
        # Group the spectra -- required before running spectral fits
        group_spectra()
        
  
    ## STEP 3: PLOT THE LIGHT CURVE AND HARDNESS RATIO
    elif ANALYSIS_TYPE=="plot_count_rates_and_hr": 
        plot_count_rates_and_hr() 


    ## STEP 4: FIT THE SPECTRA WITH NO FIXED PARAMETERS, AND GET INITIAL SPECTRAL FIT RESULTS
    elif ANALYSIS_TYPE=="fit": 
        run_spectral_fit() # fit the spectra
        get_results() # output the fit results
        plot_all_spectral_fit_results() # plot all the fit results


    ## GET THE FINAL RESULTS, BASED ON CHOSEN MODELS
    elif ANALYSIS_TYPE=="get_final_results": 
        plot_spectral_results()


    else:
        print("Analysis_type not valid.")
