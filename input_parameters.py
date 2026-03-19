######################################################################################################################

## FUNCTION TO RUN
# Change the line below to run different functions
# Options: get_swift_xrt_data, group_spectra, plot_count_rates_and_hr, fit, get_final_results
ANALYSIS_TYPE = "plot_count_rates_and_hr"
 



######################################################################################################################

## DATA RETRIEVAL -- get_swift_xrt_data
# Function: Get the data from the swifttools.ukssdc pipeline.


# Email address to use for the request
# An account needs to be set up first
EMAIL = 'justine.crook-mansour@physics.ox.ac.uk'


# Get the following information from the observation schedule website: https://www.swift.psu.edu/operations/obsSchedule.php?
TARGET_NAMES = ["Aquila X-1", "Aquila X-1", "Aquila X-1"]
TARGET_IDS =["00033665", "00035323", "00089981"]
SEGMENTS_LC =[[136, 137, 138, 139,140,142,143,144,145,146,147,148], [53,55,57,58,59,60,61,62,64,65,66,67,68,69,70,71,72], [1]]
SEGMENTS_SPEC = SEGMENTS_LC

# Target coordinates
# Use, for example, simbad: https://simbad.u-strasbg.fr/simbad/
# Or the swift observations website: https://www.swift.psu.edu/operations/obsSchedule.php?
TARGET_COORDS = [287.8169, 0.5850]


# Light curve and spectrum settings
GRADES = 'all' # options are 'all', '0', '4'

# Light curve settings
MIN_E_LC = 0.6
SOFT_LO = 0.5
SOFT_HI = 2.0
HARD_LO = 2.0
HARD_HI = 10.0
 

######################################################################################################################

## GROUPING SPECTRA -- group_spectra
# Function: Group the Swift/XRT spectral data.
# It will create both binned spectra (with MIN_COUNTS_CHI counts per bin) and unbinned spectra (with 1 count per bin) for each observation, and will use the former for spectra with counts above the COUNTS_THRESHOLD and the latter for spectra with counts below this threshold.

# Minimum energy of the first bin when conducting spectral fitting
MIN_E_KEV = MIN_E_LC


# Set the minimum counts threshold above which we bin to MIN_COUNTS_CHI per bin and below which we only do single-count binning 
MIN_COUNTS_CHI = 20 # should be greater than 20; recommended is 20 or 25
COUNTS_THRESHOLD = 300




######################################################################################################################

## PLOT COUNT RATES AND HARNESS RATIO -- plot_count_rates_and_hr
# Function: Plot the count rate light curve and hardness ratio, as calculated by the swifttools.ukssdc pipeline.


TRANSITIONS = None # MJD of state transitions to mark with vertical lines (if any); e.g. [59715.0, 59730.0]




######################################################################################################################

## Swift/XRT SPECTRAL FITTING -- fit

# Note that the routine does not fit spectra that:
# Are classified as upper limits based on the product generator light curve routine. 
# OR are WT observations with fewer than 15 counts (as these are probably spurious)
# OR PC observations with fewer than 10 counts

# Whether to include observations labelled as incbad from the lightcurve product generator -- i.e. these are ones for which a centroid could not be determined, so the position was fixed.
INCBAD = True

# MJD range for filtering spectra (None = no filtering)
MJD_MIN = None  # Minimum MJD to include (e.g., 59700.0)
MJD_MAX = None # Maximum MJD to include (e.g., 59750.0)


# Whether to run with binning or not
# If False, we use single-count binning
FIT_WITH_BINNING = True

# Sometimes it helps to renormalise the spectra before running the fits
RENORM = False

# When running with chi^2 stats, we can add systematic errors
# However, with C-stat (which is what we adopt), we cannot add systematic errors, so we set this to False
ADD_SYS_ERR = False


# Define models to fit. 
# The options are: ['powerlaw', 'pegged_powerlaw', 'diskbb', 'powerlaw+diskbb', 'pegged_powerlaw+diskbb', 'diskbb+bbodyrad','powerlaw+bbodyrad'].
MODELS= ['pegged_powerlaw', 'diskbb', 'diskbb+bbodyrad', 'powerlaw+bbodyrad'] 


# Energy range to extract the unabsorbed flux for
EMIN = 1.0
EMAX = 10.0


# Galactic absorption, NH (source-specific)
# Set the NH value, in units of 10^{22} cm^{-2}
NH = 0.5
NH_FIX_ALL_EPOCHS = True

# For the spectra with counts < LOW_COUNT_THRESHOLD, we just fit with a powerlaw, with fixed parameters (PLAW_GAMMA_LOW_COUNT and NH).
LOW_COUNT_THRESHOLD = 50 
PLAW_GAMMA_LOW_COUNT = 1.6


# If NH_FIX_ALL_EPOCHS = False and we wish to fix NH just for epochs with counts below a certain threshold, specify this here. 
# Only epochs with counts < nH_counts_threshold will have NH fixed.
NH_COUNTS_THRESHOLD = None # e.g., 100

# For epochs (if any) where NH is not fixed, set the initialisation for the fits
NH_INIT = f'{NH},,0.1,0.1,2.0,2.0' # absorption


# Other parameter initialisations and constraints for the fits
# Note that since we are simply interested in the flux, the important thing is that the fit is good rather than the parameters being physical. 
DISKBB_TIN_INIT = '0.5,,0.05,0.05,4.0,4.0'  # temperature
PLAW_GAMMA_INIT = '1.7,,0.5,0.5,4.0,4.0' # spectral index, range 0-4 
BBODY_KT_INIT = '0.6,,0.05,0.05,4.0,4.0'
# If we want a different gamma for the powerlaw+diskbb models
PLAW_GAMMA_INIT_IMS = '2.0,,0.5,0.5,4.0,4.0'


# Distance to the source in kpc
#DIST_KPC = 6




######################################################################################################################

## FINAL RESULTS -- get_final_results


# Based on the spectral fit results and any prior knowledge regarding spectral state ranges, define which spectral model to use for each point in time.
# MODELS_INDEXES is a multi-dimensional array. Each outer element corresponds to the MODELS defined above. 
# The inner elements are the index ranges for each model. 
# The index ranges are inclusive and correspond to the indexes in the fit_outputs.txt file (leftmost number).
# For example, if MODELS = ["powerlaw", "diskbb"], then MODELS_INDEXES = [ [[0,5], [10,15], [20,20]], [[6,9], [16,19]] ] would mean that we use the powerlaw model for points 0-5 and 10-15 and 20, and the diskbb model for points 6-9 and 16-20.
MODELS_INDEXES =[ [[0,5], [18,18], [20,26]], [], [[6,7], [10,10], [12,16]] , [ [8,9], [11,11], [17,17], [19,19] ]  ]


######################################################################################################################

