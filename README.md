# Swift/XRT Data Reduction for the ThunderKAT & X-KAT Programme

This repository contains code to automatically extract and process Swift/XRT data.  
The code is based on an early version developed by Andrew Hughes.


---

## Overview

Spectral fitting is highly sensitive to initial parameters, so manual fitting is often ideal.  
However, in many cases the volume of data makes this impractical. 

This pipeline enables:
- Automated fitting across many epochs
- Quick-look analysis for testing different spectral models and/or parameters such as N_h
- Consistent and reproducible processing

---

## General Notes on Spectral Fitting

Further details on the fitting strategy can be found here:  
https://www.overleaf.com/read/zqgkvmdbrffj#5921f9


### Key Points

- For ThunderKAT and X-KAT, the primary quantity of interest is **flux**, so simple spectral models are used.
- Fitting is performed using **C-statistics (W-statistics)** rather than χ², as χ² can produce biased results (e.g. Humphrey et al. 2009, and as we found in our tests).
- Each spectrum is binned in two ways:
  - **bin1**: one count per bin (preferred for fitting with C-statistics)
  - **binned**: ≥20 counts per bin (can still be fit with C-statistics, and is useful for model comparison using the χ² **test** statistic)

Both fitting approaches are supported:
- `*_bin1` → bin1 results  
- `*_bin` → binned results  

In practice, flux differences between these are usually negligible.


---

## Requirements

Notably:
- `swiftools`
- If running the notebooks in `./NOTEBOOKS/`, `plotly` is required


---

## Code Structure

- All user inputs are defined in `input_parameters.py`
- In `input_parameters.py`, change the ANALYSIS_TYPE to be one of the following, which corresponds to the sub-sections below:
    - `get_swift_xrt_data`
    - `plot_count_rates_and_hr`
    - `group_spectra`
    - `fit`
    - `get_final_results`


- Run the pipeline using:
```bash
python main_runner.py
```

---

## 1. Getting Swift/XRT Data (`get_swift_xrt_data`)

### Setup

- Register for the Swift product generator, if you haven't:  
  https://www.swift.ac.uk/user_objects/register.php

- Retrieve observation details from:  
  https://www.swift.psu.edu/operations/obsSchedule.php?


### Parameters 

- `EMAIL` – The email address you used to register for the Swift service. 
- `TARGET_NAMES`, `TARGET_IDS`, `SEGMENTS_LC`, `SEGMENTS_SPEC` – Observation details
- `TARGET_COORDS` – RA/Dec (J2000). You can use the coordinates listed after selecting "Resolve coordinates" on the link above, or coodinates from simbad. 
- `GRADES` – 'all'/'4'/'0'. Which grades to use. Due to systematics that differ between grade selection, it is generally preferable to use the same grade selection for all observations of a particular source. 'all' uses all grades and is generally recommended. 
If there is heavy pile-up, using the restricted grade selection ('4') can help; however, the difference in results between this and using all grades is generally negligible for the precision required.
- `MIN_E_LC` – Mininum energy for the light curves. Usually 0.6 or 0.5 is recommended, as there can be some artefacts in high-count WT spectra at lower energies. The maximum energy is set to 10 keV in the code. 
- `SOFT_LO`, `SOFT_HI` – Respectively the low and high energies to use for the soft band when calculating hardness ratios. 
- `HARD_LO`, `HARD_HI` – Respectively the low and high energies to use for the hard band when calculating hardness ratios. 



### Results
- The extracted light curve data and spectra will respestively be in `./lightcurves_swift_xrt/` and `./spectra_swift_xrt/`.


### Notes

- If using the Swift GUI to retrive the data instead:
  - Place spectra (Obs*) in `./spectra_swift_xrt/`, and do not use sub-folders.
  - Place light curves in `./lightcurves_swift_xrt/`. Do not use subfolders unless grouping by ID (00*).




---

## 2. Plotting Light Curves (`plot_count_rates_and_hr`)


### Optional Data

Place additional datasets in `./other_lightcurve_data/`:
- MeerKAT: `meerkat.txt`, using the template
- BAT: `bat_daily_ave.txt` - https://swift.gsfc.nasa.gov/results/transients/
- MAXI: `maxi_daily_ave.txt` - https://maxi.riken.jp/top/lc.html



### Parameters

- `TRANSITIONS` – Optional vertical lines at specific times (e.g., state transitions) in the plots



### NOTEBOOKS/HID_diagrams.ipynb:
- In this notebook, you can plot the hardness ratios as a function of time, as well as hardness intensity diagrams. 




---

## 3. Grouping Spectra (`group_spectra`)

Each spectrum is processed as follows:

- `*final_bin1.pi` → 1 count per bin
- `*final_bin.pi` → ≥ `MIN_COUNTS_CHI` counts per bin (if counts ≥ `COUNTS_THRESHOLD`)


### Parameters

- `MIN_E_KEV` - Minimum energy below which data is ignored. By default, this is the same as MIN_E_LC.
- `MIN_COUNTS_CHI` - When binning, how many counts should there be per bin. For χ² as a test statistic, ≥20 is recommendeded.
- `COUNTS_THRESHOLD` - Set the minimum count threshold above which we bin to MIN_COUNTS_CHI per bin.

---


## 4. Fitting Spectra (`fit`)

### Key Parameters

- `INCBAD` – Whether to include observations labelled as "bad" from the light curve generator (see note below).
- `MJD_MIN`, `MJD_MAX` – The MJD range (inclusive) for the spectra to fit. If None, it will fit all spectra in the folder ./spectra_swift_xrt/. 
- `FIT_WITH_BINNING`:
  - `True` → uses `*final_bin.pi` and outputs results to `./spectral_fit_results_bin/` and `./spectral_fit_residuals_bin/`. In this case, it is appropriate to use the χ² test statistic to compare models. 
  - `False` → uses `*final_bin1.pi` and outputs results to `./spectral_fit_results_bin1/` and `./spectral_fit_residuals_bin1/`. In this case, it is NOT appropriate to use the χ² test statistic to compare models. 
- `RENORM` – Renormalise models before fitting, which sometimes helps fitting convergence
- `ADD_SYS_ERR` – Whether to add systematic error to the response. Note that since we use cstat for all the fitting, this actually has no effect. 
- `MODELS` – Models to fit for each spectrum. Options are:
  ```python
  ['powerlaw', 'pegged_powerlaw', 'diskbb',
   'powerlaw+diskbb', 'pegged_powerlaw+diskbb',
   'diskbb+bbodyrad', 'powerlaw+bbodyrad']
- `EMIN`, `EMAX` – Energy range over which to extract unabsorbed flux  
- `NH` – The expected column density (N_h) in units of 10^{22} cm^{-2}. Or if it is expected to change, just list a good value for initialisation. 
- `NH_FIX_ALL_EPOCHS` - Whether to fix N_h in all the fits.   
- `LOW_COUNT_THRESHOLD` - For the spectra with counts < `LOW_COUNT_THRESHOLD`, we just fit with a powerlaw, with fixed parameters (`PLAW_GAMMA_LOW_COUNT` and `NH`). 
- `PLAW_GAMMA_LOW_COUNT` - The powerlaw photon index to use for spectra with counts < `LOW_COUNT_THRESHOLD`. 
- `NH_COUNTS_THRESHOLD` - This parameter is applicable/only used in cases when `NH_FIX_ALL_EPOCHS=False` and when `NH_COUNTS_THRESHOLD` > `LOW_COUNT_THRESHOLD` (since N_h will always be fixed when counts < `LOW_COUNT_THRESHOLD`).   


### Model Initialisation Format

The model parameter intialisation is a string with the format:
```text
<param>,<delta>,<hard min>,<soft min>,<soft max>,<hard max>
```
where: 
- param: trial param value used initially in fit
- delta: step size used in the numerical determination of the derivatives during the fitting process. This value may be overriden for all parameters by the xset delta command option, which will apply a proportional rather than a fixed delta.
- hard/soft min: minimum for the param. hard/soft max: maximum for the param.
- Note that since we are simply interested in the flux, the important thing is that the fit is good rather than the parameters being physical. 



### Initial Parameters

- `NH_INIT` - For epochs (if any) where NH is not fixed, set the initialisation for the fits.   
- `DISKBB_TIN_INIT` - Diskbb temperature in keV.   
- `PLAW_GAMMA_INIT` - Powerlaw photon index. 
- `BBODY_KT_INIT` - Blackbody temperature in keV.  
- `PLAW_GAMMA_INIT_IMS` -  If we want a different gamma for the powerlaw+diskbb models.  


### Adding Models

- If you wish to add a model, you should do so in `initialise_model`. 
- Also every time the code checks whether all the models specified in `MODELS` are valid models, you should add the name for the new model.



### Outputs

If binned spectra were used (`FIT_WITH_BINNING=True`), the fitting results and residuals are respectively found in:
- `./spectral_fit_results_bin/`  
- `./spectral_fit_residuals_bin/`  

If bin1 spectra were used (`FIT_WITH_BINNING=False`), the fitting results and residuals are respectively found in:
- `./spectral_fit_results_bin1/`  
- `./spectral_fit_residuals_bin1/`  


Note regarding the residuals: 
- Even if the spectra were fit with bin1, it is useful to rebin the spectra just for visualisation -- so that the fit can be better examined by eye. 
- In XSPEC, this can be achieved through `setplot rebin mincounts maxbins` where `mincounts` is the minimum number of counts per bin, and `maxbins` is the maximum number of bins to combine. 
- As such, for each spectrum that was not binned during fitting, if the counts is greater than `COUNTS_THRESHOLD`, I use `mincounts = 20` and `maxbins = 5`. 
If the counts is less than this but greater than `LOW_COUNT_THRESHOLD`, I use `mincounts = 10` and `maxbins = 3`. 
Else, I used `mincounts = 3` and `maxbins = 3`. 
The setting used is shown in the title of the residual plots. 


Each folder `./spectral_fit_results_bin*/ contains:
- Plots showing all results.  
- `fit_outputs.txt` tabulating all the results.  
- `not_fit.txt` are the spectra that were not fit (upper limits or too few counts).  
- `run.txt` shows the input and output parameters for each fit for each spectrum (which is useful for debugging). It also shows the χ² comparisons in the case of FIT_WITH_BINNING=True, which is useful for model selection. 
- `xspec.log` shows a log of the fits. This is useful to see the details of the fits and whether any parameters were pegged, for example. 



### NOTEBOOKS/incbad_obs.ipynb:
- For some IDs, products generated by the Swift pipeline are labelled "bad". 
- This does not actually mean they are bad per se -- just that no centroid was found, so the input position was used. 
- If the astronometry is good, the result shouldn't be bad, but is less good statistically. 
- We may want to exclude such spectra, or check them manually. 
- In this notebook, we show how to print which ones were labelled as "bad" by the light curve generator. 
- Note that the light curve generation excludes observations shorter than 7 seconds (for historical reasons, since this was first built for GRBs), while the spectal generator includes them, so these observations (if any) are missing from the lists. It is, however, fine to use these spectra (assuming the data are of good quality).


### NOTEBOOKS/example_manual_fitting.ipynb:
- Since the code is automated, there are cases where it may be beneficial to fit the data manually. 
- This is mainly because some spectral fits can be extremely sensitive to initial parameters -- especially in cases where more complex models are used. 
- This notebook shows how to do manual fits in PyXspec for different models.



---

## 5. Getting Final Results (`get_final_results`)

After fitting, based on the results (and any other knowledge about spectral states), choose which model should be used for each spectrum. 


### Parameter: `MODELS_INDEXES`
- A nested list defining which model applies to which index range. 
- Each outer element corresponds to the `MODELS` array defined above. The inner elements are the index ranges for each model. 
- Indexes correspond to rows in `fit_outputs.txt` (leftmost number). 



### Example

```python
MODELS = ["powerlaw", "diskbb"]

MODELS_INDEXES = [
    [[0, 5], [10, 15], [20, 20]],
    [[6, 9], [16, 19]]
]
```
This means that we should use the powerlaw model for points 0-5 and 10-15 and 20, and the diskbb model for points 6-9 and 16-20.




### Outputs:
- The results will either be in ./final_spectral_fit_results_bin/ (for FIT_WITH_BINNING=True) or ./final_spectral_fit_results_bin1/ (FIT_WITH_BINNING=False).



### NOTEBOOKS/binning_comparison.ipynb:
- If you ran both with both FIT_WITH_BINNING=True and FIT_WITH_BINNING=False, you can use this notebook to compare these flux results. 




---

## Get Dates 

In the code folder, the script `fit_xrt_spectra.py` is the most accurate way to get the Swift/XRT dates. This uses `swifttime`. 


---

## TO DO:
- Parallelise the fitting so that it is done faster. 
- Add the upper limits functionality again?
- Add functionality for models with spectral lines, which may be needed for NICER spectra. 
- Add code for extracting NICER spectra. 


---

Please send any comments and suggestions to justine.crook-mansour@physics.ox.ac.uk.




