import time
import json
import glob
import subprocess
import os
import sys
import numpy as np
from astropy.io import fits
import shutil
from xspec import *

from input_parameters import *

from swifttools.ukssdc.xrt_prods import XRTProductRequest


import sys
import os
sys.path.append(os.path.abspath("../"))
from input_parameters import *


#######################################################################

def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp+txt)


#######################################################################



def get_xrt_prods_runner():

    target_names, target_ids, segment_lc = TARGET_NAMES, TARGET_IDS, SEGMENTS_LC
    
    # Iterate through Target name(s), id(s), and segments and get xrt_prods
    print(target_names, target_ids, segment_lc)

    for target_name, target_id, segment_lc in zip(target_names, target_ids, segment_lc):
        print("TARGET: ", target_name, "ID: ", target_id)
        print("Getting light curve:")
        get_xrt_prods(target_id, target_name, segment_lc, centroid=True, prod_type = 'lightcurve')
        print("-------------------------------------------------------------------------------------------------------")
    print("-------------------------------------------------------------------------------------------------------")
    
    
    for target_name, target_id, segment_spec in zip(TARGET_NAMES, TARGET_IDS, SEGMENTS_SPEC):
        print("TARGET: ", target_name, "ID: ", target_id)
        try:
            print("Getting spectrum:")
            get_xrt_prods(target_id, target_name, segment_spec, centroid=True, prod_type = 'spectrum')
        except:
            print("Error getting spectrum.")
        print("-------------------------------------------------------------------------------------------------------")
        


def get_xrt_prods(target_id, target_name, segments, centroid = True, prod_type = 'spectrum'):

    target_coords = TARGET_COORDS
    email = EMAIL
    min_E_lc, softLo , softHi , hardLo, hardHi, grades = MIN_E_LC, SOFT_LO, SOFT_HI, HARD_LO, HARD_HI, GRADES

    # Make sure its a valid data type
    if prod_type != 'spectrum' and prod_type != 'lightcurve':
        sys.exit('Error: Please specify variable [prod_type] as either "spectrum" or "lightcurve"')

    # Ensure the segments is atleast a 1d array
    segments = np.atleast_1d(segments)

    # Initialise the observation IDS
    obs_ids = ','.join([f'{target_id}{seg:03d}' for seg in segments])
    print("Requested observation IDs: ", obs_ids)

    # Create an XRTProductRequest object
    # silent: whether messages should be printed out to stdout
    myRequest = XRTProductRequest(email, silent=False) 

    # Set the global parameters
    # name: object name
    # targ: comma-separated list of targetIDs to use
    # SinceTO: whether all time variables are relative to T0
    # RA and dec: object coordinates in J2000
    # centroid: boolean indicating whether to try centroid in the XRT coordinate frame
    # useSXPS: boolean indicating whether to use source lists from SXPS where possible. The 2SXPS catalogue is used to identify sources in the field which should be excluded from the background. This is not really necessary for single source pointed observations, especially when our source is the brightest.
    # poserr: how far from the input position the centroid position can be (arcmin)
    myRequest.setGlobalPars(name=target_name,
                        targ=target_id,
                        SinceT0=False,
                        RA=target_coords[0],
                        Dec=target_coords[1],
                        centroid=centroid,
                        useSXPS=True, 
                        posErr=1)

    # ~~~~~~~~ #
    # Spectrum #
    # ~~~~~~~~ # 
    # This code collects the spectra. 
    # For each obsID, it returns the following files:
    #- .areas (extraction areas used)
    #- _pow_fit.fit OR if I turned off fitting _nofit_fit.fit
    #- .arf
    #- .pi
    #- .rmf
    #- back.pi
    #- source.pi
    if prod_type == 'spectrum':

        print("\n Spectrum")

        # Add a spectrum to the request
        # useObs: specific observations to use to create the spectrum
        # hasRedshift: redshift of the source
        # timeslice: what spectra to create, must be one of {'single', 'user', 'snapshot', 'obsid'}
        myRequest.addSpectrum(whichData = 'user',
                        useObs = obs_ids,
                        hasRedshift = False,
                        timeslice = 'obsid', 
                        grades = grades,
                        doNotFit=True)
        print(myRequest.getAllPars())


        if(myRequest.isValid()[0]==True): # Check the request is valid
            # Submit the request and check whether the returned value is False (i.e. the job was not successfully submitted)
            if not myRequest.submit(): 
                 msg(f'I could not submit error:{myRequest.submitError}')
    
            # Wait until the job is complete
            while not myRequest.complete:
                #print("Products not yet complete. Waiting 60 seconds.")
                time.sleep(60)

            # Obtain the spectral data in the standard spectra dict.
            outdict = myRequest.retrieveSpectralFits(returnData=True)
            print(outdict.keys()) # e.g. ['T0', 'GalNH', 'rnames', 'Obs_00016584019']
        
        else: # If the request was invalid, print the problem
            msg('BAD REQUEST: {}'.format(myRequest.isValid()[1]))
            exit()   

        # To download the data, we could alternatively use: myRequest.downloadProduct('myPath', format='zip')
        # Or, using outdict above, we could use SaveSpectralData()

        observations = [key for key in outdict.keys() if 'Obs' in key]
        print("Observations: ", observations) # e.g. ['Obs_00016584019']

        # Ensure the spectra directory exists
        spectra_dir = "../spectra_swift_xrt/"
        if not os.path.exists(spectra_dir):
            os.makedirs(spectra_dir)
            print(f"{spectra_dir} has been created.")

        # Downloading and de-tarring files
        for obs in observations:
            try:
                pipeline_output = outdict[obs]['DataFile'] # get URL of data file
                subprocess.run([f'wget -O {spectra_dir}{obs}.tar.gz {pipeline_output}'], shell=True)
                subprocess.run([f'tar -xf {spectra_dir}{obs}.tar.gz -C {spectra_dir}'], shell=True)
                subprocess.run([f'rm -rf {spectra_dir}{obs}.tar.gz'], shell=True)
            except:
                print("Something went wrong with ", obs)

    

    # ~~~~~~~~~~ #
    # Lightcurve #
    # ~~~~~~~~~~ # 
    
    # We are just interested in the returned:
    # (1) curve_incbad.qdp. Note that curve2_incbad.qdp, is the same but with additional information.
    # (2) hardrat_incbad.qdp
    # (3) WT_incbad.qdp and WT.qdp
    # (4) PC_incbad.qdp and PC.qdp


    else:
        print("\n Light curve")
        
        # Add a light curve to the request
        # useObs: specific observations to use to create the light curve
        # binMeth: which binning method to use. Must be one of {'counts', 'time', 'snapshot', 'obsid'}
        # timeFormat: units to use on the time axis. Must be one of {'s'(=seconds), 'm'(=MJD)}
        # minEnergy & maxEnergy: min and max energy for the main light curve, in keV
        # softLo & softHi: min and max energy for the soft-band, in keV. Must respectively be ≥minEnergy and ≤maxEnergy.
        # hardLo & hardHi: min and max energy for the hard-band, in keV. Must respectively be ≥minEnergy and ≤maxEnergy.
        # allowUL: whether upper limit are allowed. Must be one of {'no' 'pc', 'wt', 'both'}.
        myRequest.addLightCurve(whichData = 'user',
                        useObs = obs_ids,
                        binMeth ='obsid',
                        timeFormat ='m',
                        minEnergy = min_E_lc,
                        maxEnergy = 10.0,
                        softLo = softLo,
                        softHi = softHi,
                        hardLo = hardLo,
                        hardHi = hardHi,
                        grades = grades,
                        allowUL='both')

        # Submit the XRT_prod request
        if(myRequest.isValid()[0]==True):
            if not myRequest.submit():
                 msg(f'I could not submit error:{myRequest.submitError}')
    
            while not myRequest.complete:
                time.sleep(60)

            # Obtain the light curve data in the standard light curve dict.
            # returnData: a bool indicating whether the function should return the data, as well as storing it internally (default: False).
            # nosys: whether to return the data from which WT systematics have been excluded. Can be 'yes', 'no' or 'both' (default: 'no').
            # incbad: whether to return the light curve(s) which include data from times when no centroid could be obtained. Can be 'yes', 'no' or 'both' (default: 'no').
            outdict = myRequest.retrieveLightCurve(returnData = True, incbad = True)
            print(outdict.keys())
        
        else:
            msg('BAD REQUEST: {}'.format(myRequest.isValid()[1]))
            exit()   

        # Check a directory exists to store the light curves
        lightcurve_dir = "../lightcurves_swift_xrt/"
        if not os.path.exists(lightcurve_dir):
            os.makedirs(lightcurve_dir)
            print(f"{lightcurve_dir} has been created.")

        # Check if there is a directory for the target_id
        lightcurve_target_dir = f'{lightcurve_dir}{target_id}'
        #if os.path.exists(lightcurve_target_dir):
        #    shutil.rmtree(lightcurve_target_dir)
        #    print(f"{lightcurve_target_dir} has been removed.")
        #else:
        #    os.makedirs(lightcurve_target_dir)
        #    print(f"{lightcurve_target_dir} has been created.")
        if not os.path.exists(lightcurve_target_dir):
            os.makedirs(lightcurve_target_dir)
            print(f"{lightcurve_target_dir} has been created.")
        
        # Download the light curve products
        # By default, it downloads the products into tar.gz files
        ##TODO Might want to add 'LightCurve' argument so that it only downloads the lightcurves??
        myRequest.downloadProducts(f'{lightcurve_dir}{target_id}')

        # De-tar the files 
        subprocess.run([f'tar -xf {lightcurve_dir}{target_id}/lc.tar.gz -C {lightcurve_dir}{target_id}/'], shell = True)
        # Move files around
        subprocess.run(f'mv {lightcurve_dir}{target_id}/USERPROD*/lc/* {lightcurve_dir}{target_id}/.', shell=True)
        # Clean up
        subprocess.run([f'rm -rf {lightcurve_dir}{target_id}/USERPROD*; rm -rf {lightcurve_dir}{target_id}/lc.tar.gz'], shell=True) 

        # Remove nonsense string
        # Finds all .qdp files in the specified directory, and in each file, it removes any text following !:: on every line
        #for qdp_name in glob.glob(f'{lightcurve_dir}*/*.qdp'):
        #    subprocess.run([f"sed -i 's/!::.*//g' {qdp_name}"], shell=True)




#######################################################################

def group_spectra():

    min_E_keV = MIN_E_KEV
    min_counts_chi = MIN_COUNTS_CHI
    counts_threshold = COUNTS_THRESHOLD

    print("min_E_keV: ", min_E_keV)
    print("min_counts_chi: ", min_counts_chi)
    print("counts_threshold: ", counts_threshold)


    # Get the current working directory
    spectrum_directory = os.getcwd()+ '/spectra_swift_xrt/'
    print("Spectrum directory: ", spectrum_directory)


    # Sorts all the source spectrum files alphabetically
    source_spectra = sorted(glob.glob(f'{spectrum_directory}*source.pi'))
    print("Source spectra: ", source_spectra)
    print()
    print()


    # Iterate through spectral files, i.e. for every obsID
    for source_spectrum in source_spectra[:]:
        print()
        print("Spectrum: ", source_spectrum)
        print('\n', fits.getheader(source_spectrum)['DATE-OBS'])


        source_spectrum_with_header = source_spectrum.replace('source.pi', 'source_with_header.pi')
        background_spectrum = source_spectrum.replace('source.pi', 'back.pi')
        rmf = source_spectrum.replace('source.pi', '.rmf')
        arf = source_spectrum.replace('source.pi', '.arf')


        syscall = f'grppha {source_spectrum} !{source_spectrum_with_header} comm="chkey backfile {background_spectrum} & chkey respfile {rmf} & chkey ancrfile {arf} & exit"'
        subprocess.run([syscall], shell = True)  


        # Get the number of counts
        AllData.clear()
        # Load or select the spectrum
        s = Spectrum(source_spectrum_with_header)
        AllData.ignore(f"**-{min_E_keV}")  # Ignore below lower bound
        AllData.ignore("10.0-**")                          # Ignore above 10 keV
        AllData.notice(f"{min_E_keV}-10.0") # Notice desired energy range
        AllData.ignore("bad") 
        time.sleep(2)
        print("With background subtraction:", s.rate[0]*s.exposure)
        print("Without background subtraction:", s.rate[2]*s.exposure)
        total_counts = s.rate[2]*s.exposure


        if total_counts >= counts_threshold:
            binning = min_counts_chi
            print(f"Number of counts: {total_counts}. Binning one file to {min_counts_chi} per bin, and another to one count per bin.")
        else:
            binning = 1
            print(f"Number of counts: {total_counts}. Binning just to one count per bin.")

        final_spectrum = source_spectrum.replace('source.pi', 'final_bin.pi')
        final_spectrum_bin1 = source_spectrum.replace('source.pi', 'final_bin1.pi')
        
        ## With binning
        ignore_max = int(min_E_keV/0.01) - 1 # the number of channels to ignore above the lower bound, given that the channel energy binning is 10 eV
        syscall = f'grppha {source_spectrum} !{final_spectrum} comm="bad 0-{ignore_max} & group min {binning} & chkey backfile {background_spectrum} & chkey respfile {rmf} & chkey ancrfile {arf} & exit"'
        #syscall = f'grppha {source_spectrum} !{final_spectrum} comm="group min {binning} & chkey backfile {background_spectrum} & chkey respfile {rmf} & chkey ancrfile {arf} & exit"'
        subprocess.run([syscall], shell = True)  


        ## Also create a bin1 file
        syscall = f'grppha {source_spectrum} !{final_spectrum_bin1} comm="bad 0-{ignore_max} & group min 1 & chkey backfile {background_spectrum} & chkey respfile {rmf} & chkey ancrfile {arf} & exit"'
        subprocess.run([syscall], shell = True)  


        # Add results to header
        with fits.open(final_spectrum, mode="update") as hdul:
            print("Updating the header")
            hdr = hdul[1].header  # Extension 1 is typically the PHA data
            hdr["COUNTS"] = int(total_counts)
            hdr["COUNTGRP"] = binning
            hdul.flush()  # Save changes


        # Add results to header
        with fits.open(final_spectrum_bin1, mode="update") as hdul:
            print("Updating the header")
            hdr = hdul[1].header  # Extension 1 is typically the PHA data
            hdr["COUNTS"] = int(total_counts)
            hdr["COUNTGRP"] = 1
            hdul.flush()  # Save changes


        print("-----------------------------------------------------")



#######################################################################

def group_spectra_custom():
    # Prepare spectral data for spectral fitting using Sivakoff's custom grouping code.

    min_E_keV = MIN_E_KEV
    min_counts_chi = MIN_COUNTS_CHI
    counts_threshold = COUNTS_THRESHOLD

    print("-------------------------------------------------------------------------------------------------------")
    print("Prepare spectral data")  
    print()

    min_E = min_E_keV*1e3 # eV

    # Get the current working directory
    spectrum_directory = os.getcwd()+ '/spectra_swift_xrt/'

    # Sorts all the source spectrum files alphabetically
    source_spectra = sorted(glob.glob(f'{spectrum_directory}/*source.pi'))
    print("Source spectra: ", source_spectra)
    print()
    print()

    # Iterate through spectral files, i.e. for every obsID
    for source_spectrum in source_spectra[:]:
        print()
        print("Spectrum: ", source_spectrum)
        print('\n', fits.getheader(source_spectrum)['DATE-OBS'])

        # Note that the arf and rmf files are returned by the data retrieval
        # Define file names, by replacing certain parts of the source file name for the observation under consideration
        group_spectrum = source_spectrum.replace('source.pi', 'group.pi')
        final_spectrum = source_spectrum.replace('source.pi', 'final.pi')
        background_spectrum = source_spectrum.replace('source.pi', 'back.pi')
        rmf = source_spectrum.replace('source.pi', '.rmf')
        arf = source_spectrum.replace('source.pi', '.arf')

        # The background, rmf, and arf files need to be set in the header, for later use by XSPEC.
        # Call grppha to read the source spectrum file ({source_spectrum}) and write to the grouped spectrum file ({group_spectrum}). 
        # The ! prefix tells grppha to overwrite the file if it already exists.
        # comm="...": Specifies a series of grppha commands to be executed sequentially within the grppha session:
        # chkey backfile {background_spectrum}: Sets the background file key to {background_spectrum} in the spectrum header.
        # chkey respfile {rmf}: Sets the response matrix file (RMF) key to {rmf}.
        # chkey ancrfile {arf}: Sets the ancillary response file (ARF) key to {arf}.
        # I.e. This creates a file Obs{obsID&obsType}group.pi
        syscall = f'grppha {source_spectrum} !{group_spectrum} comm="chkey backfile {background_spectrum} & chkey respfile {rmf} & chkey ancrfile {arf} & exit"'
        subprocess.run([syscall], shell = True)         

        # Run Sivakoff GRPPHA to group together the counts (i.e. group the spectra)
        # i: input pi file
        # o: output (grouped) pi file
        # l: minimum energy of first bin in eV -- default is 500eV
        # u: maximum energy of first bin in eV
        # c: minimum counts per bin, when ncounts> counts_threshold -- default is 25
        # e: minimum energy change per bin in eV or fraction
        # t: counts threshold for switching between chi^2 and cstat binning -- default is 300
        # This creates a file Obs{obsID&obsType}final.pi
        syscall = f'python3 new_grppha.py -i {group_spectrum} -o {final_spectrum} -l {min_E} -u 10000 -c {min_counts_chi} -e 0.0 -t {counts_threshold}'
        subprocess.run([syscall], shell = True)   

        print("----------------------------------------------------------------------------------------------")
   

    