* Swift/XRT data reduction for the ThunderKAT & X-KAT programme:

In this folder, I share code to automatically extract and process Swift/XRT data. This code is based on an early version of code by Andrew Hughes.


Given the fact that fits are very sensitive to initial parameters, the ideal scenario would be to conduct the fits by hand. 
However, in some cases, we have too much data for this to be possible, so automating this is necessary. 
Even in cases where we do not have many epochs, running these automated fits can be useful as a first 'quick-look' to try out different spectral models.  


----------------------------------------------

* GENERAL NOTES ABOUT SPECTRAL FITTING:

Further details regarding the strategy for spectral fitting can be found here: https://www.overleaf.com/read/zqgkvmdbrffj#5921f9


MOST NOTEWORTHY POINTS:

(1) For ThunderKAT and X-KAT, we are only interested in the flux, so I make use of simple spectral models for the fits.


(2) Research has shown (e.g., Humphrey+2009) that fitting with $\chi^2$ statistics can lead to biased fits. 
Indeed, in our tests, we found that there are sometimes large differences in the extracted flux. 
As a result, we conduct all our fitting using C-statistics (which is actually W-statistics, since we specify a background file). 
C-statistic fitting can be conducted (as is best performed) on data that has been been binned to one count per bin (bin1). 
However, model comparison and interpretation are more involved in this case. 
As such, it is also useful to bin the data (that have sufficient counts) to >=20 counts per bin. Even if we fit these with C-statistics, we can nevertheless evaluate the models using $\chi$^2 statistics to get an idea of how different models compare. 
In the code in this package, I have structure it such that you can run either with bin1 data or binned data. The former results are saved to files called *_bin1 while the latter are saved to *bin, so that the results can be compared. 
Generally, the flux differences between both cases are negligible.     


----------------------------------------------

* REQUIREMENTS:
- The swiftools package is required. 
- Additionally, for some of the notebooks in ./NOTEBOOKS/, the package plotly is required.  



* OVERALL STRUCTURE OF THIS CODE BASE:
- All the parameters that the user needs to enter are in input_parameters.py. 
- To run the code, run main_runner.py.
- In input_parameters.py, change the ANALYSIS_TYPE to be one of the following, which corresponds to the sub-sections below:
* get_swift_xrt_data
* plot_count_rates_and_hr
* group_spectra
* fit
* get_final_results





----------------------------------------------

* GETTING THE SWIFT/XRT DATA -- get_swift_xrt_data: 
- First, sign up for the Swift product generator here, if you have not already: https://www.swift.ac.uk/user_objects/register.php

Then, in input_parameters.py, change the following:
- EMAIL: The email address you used to register for the Swift service. 
- TARGET_NAMES, TARGET_IDS, SEGMENTS_LC, SEGMENTS_SPEC: These can be obtained from https://www.swift.psu.edu/operations/obsSchedule.php?. 
- TARGET_COORDS: You can use the coordinates listed after selecting "Resolve coordinates" on the link above, or coodinates from simbad. These are RA and dec (J2000).
- GRADES: Which grades to use. Due to systematics that differ between grade selection, it is generally preferable to use the same grade selection for all observations of a particular source. 'all' uses all grades and is generally recommended. 
If there is heavy pile-up, using the restricted grade selection ('4') can help; however, the difference in results between this and using all grades is generally negligible for the precision required.
- MIN_E_LC: Mininum energy for the light curves. Usually 0.6 or 0.5 is recommended, as there can be some artefacts in high-count WT spectra at lower energies. The maximum energy is set to 10 keV in the code. 
- SOFT_LO, SOFT_HI: Respectively the low and high energies to use for the soft band when calculating hardness ratios. 
- HARD_LO, HARD_HI: Respectively the low and high energies to use for the hard band when calculating hardness ratios. 


NOTEBOOKS/HID_diagrams.ipynb:
- In this notebook, you can plot the hardness ratios as a function of time, as well as hardness intensity diagrams. 


If you instead choose to use the Swift tools GUI, make sure to place the resultant spectra in ./spectra_swift_xrt/. 
The spectral files (Obs*) should be in that folder, not in any subfolders. 
Also, place the light curve results in ./lightcurves_swift_xrt/. Again, there is no need to make any sub-folders. 
But if you have the results from different IDs separately, you can place each in subfolders named by the ID (00*).   


----------------------------------------------

* PLOTTING THE COUNT RATE LIGHT CURVES -- plot_count_rates_and_hr:

In input_parameters.py, change the following, if required:
- TRANSITIONS: If you want to plot any vertical lines at specific times, specify these here. 


If you have MeerKAT data, place this in ./other_lightcurve_data/meerkat.txt using the template. 
If included, these data will be plotted.


Also, fetch other X-ray data for the source, as these are also useful to plot, and place them in ./other_lightcurve_data/:
- BAT daily data: Get it at https://swift.gsfc.nasa.gov/results/transients/ and name the file bat_daily_ave.txt
- MAXI daily data: Get it at https://maxi.riken.jp/star_data/ and name the file maxi_daily_ave.txt



----------------------------------------------

* GROUPING THE DATA -- group_spectra:

The code is set up that for each spectra it:
- Groups to one count per bin and calls the resultant spectrum *final_bin1.pi
- If the number of counts are high enough (>=COUNTS_THRESHOLD), groups to MIN_COUNTS_CHI per bin, and calls the resultant spectrum *final_bin.pi

In input_parameters.py, there are the following parameters:
- MIN_E_KEV: By default, this is the same as MIN_E_LC
- MIN_COUNTS_CHI: When binning, how many counts should there be per bin. For $\chi$^2 statistics, this should be a minimum of 20.
- COUNTS_THRESHOLD: Set the minimum counts threshold above which we bin to MIN_COUNTS_CHI per bin.


----------------------------------------------

* FITTING THE DATA -- fit:

In input_parameters.py, there are the following parameters:
- INCBAD: Whether to include observations labelled as "bad" from the light curve generator (see note below).
- MJD_MIN, MJD_MAX: The MJD range (inclusive) for the spectra to fit. If None, it will fit all spectra in the folder ./spectra_swift_xrt/. 
- FIT_WITH_BINNING: If True, it will use the spectra *final_bin.pi and output results to ./spectral_fit_results_bin/ and ./spectral_fit_residuals_bin/. 
In this case, it is appropriate to use the $\chi^2$ test statistic to compare models. 
If False, it will it will use the spectra *final_bin1.pi and output results to ./spectral_fit_results_bin1/ and ./spectral_fit_residuals_bin1/.
In this case, it is not appropriate to use the $\chi^2$ test statistic to compare models. 
- ADD_SYS_ERR: Whether to add systematic error to the response. Note that since we use cstat for all the fitting, this actually has no effect. 
- MODELS: Which models to fit for each spectrum. The options are: ['powerlaw', 'pegged_powerlaw', 'diskbb', 'powerlaw+diskbb', 'pegged_powerlaw+diskbb', 'diskbb+bbodyrad','powerlaw+bbodyrad']
- EMIN, EMAX: The energy range over which to extract the unabsorbed flux.


IMPORTANT: 
The fitting results are found within ./spectral_fit_results_bin/ if a binned spectra were used (FIT_WITH_BINNING=true), or ./spectral_fit_results_bin1/ if the bin1 spectra were used (FIT_WITH_BINNING=False). 



NOTEBOOKS/incbad_obs.ipynb:
Products generated by the Swift pipeline for some IDs are labelled "bad". 
This does not actually mean they are bad per se -- just that no centroid was found, so the input position was used. 
If the astronometry is good, the result shouldn't be bad, but is less good statistically. 
We may want to exclude such spectra, or check them manually. 
In this notebook, we show how to print which ones were labelled as "bad" by the light curve generator. 
Note that the light curve generation excludes observations shorter than 7 seconds (for historical reasons, since this was first built for GRBs), while the spectal generator includes them, so these observations (if any) are missing from the lists. 
It is, however, fine to use these spectra (assuming the data are of good quality).


NOTEBOOKS/example_manual_fitting.ipynb
- Since the code is automated, there are cases



----------------------------------------------

* GETTING THE FINAL RESULTS

NOTEBOOKS/binning_comparison


----------------------------------------------

* TO DO:
- Parallise the fitting so that it is done faster. 
- Add the upper limits functionality again?
- Add functionality for models with spectral lines, which may be needed for NICER spectra. 
- Add code for extracting NICER spectra. 




----------------------------------------------

Please send any comments and suggestions to justine.crook-mansour@physics.ox.ac.uk.




