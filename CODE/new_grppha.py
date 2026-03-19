#!/usr/bin/python

import os, optparse, sys
import numpy as np
from astropy.io import fits
np.set_printoptions(threshold=np.inf) 


## Parse the arguments to the python file
def parse_args( argv):

    desc="""%prog creates a grouped pha spectral file suitable for Xspec over a specified energy range. \
    User must specify input and ouput file. User can specifiy minimum binning based on spectral resolution and minimum number of counts per bin."""
    
    parser = optparse.OptionParser(description=desc)

    parser.add_option('-i', '--input', \
                    help='input pi file', \
                    dest='input_pi')
    parser.add_option('-r', '--rmf', \
                        help='rmf file', \
                        dest='rmf')
    parser.add_option('-a', '--arf', \
                        help='arf file', \
                        dest='arf')
    parser.add_option('-b', '--background', \
                        help='background pi file', \
                        dest='background')
    parser.add_option('-o', '--output', \
                        help='output (grouped) pi file', \
                        dest='output_pi')
    parser.add_option('-e', '--energy_bin', \
                        help='minimum energy change per bin in eV or fraction [%default]', \
                        dest='bin_min_energy', \
                        type='float', \
                        default=100.0)
    parser.add_option('-f', '--fractional_energy_binning', \
                        help='switch to set fractional energy binning [%default]', \
                        dest='fraction_switch', \
                        default=False, \
                        action='store_true')
    parser.add_option('-l', '--lower_energy', \
                        help='minimum energy of first bin in eV [%default]', \
                        dest='lower_energy', \
                        type='float', \
                        default=300.0)
    parser.add_option('-u', '--upper_energy', \
                        help='maximum energy of last bin in eV [%default]', \
                        dest='upper_energy', \
                        type='float', \
                        default=10000.0)
    parser.add_option('-c', '--count_bin', \
                        help='minimum counts per bin [%default]', \
                        dest='bin_min_counts', \
                        type='long', \
                        default=25)
    parser.add_option('-t', '--counts_threshold', \
                        help = 'minimum counts threshold for grouping [%default]',\
                        dest='counts_threshold', \
                        type='long', \
                        default=300)

    if (argv == []):
        parser.print_help()
        exit(-1)

    return parser



## Check whether the arguments are valid
def valid_args_check(opts, parser):

    # Check the input file exists and is well specified
    if opts.input_pi is None:
        print("The input pi file is missing\n")
        parser.print_help()
        exit(-1)
    if os.path.isfile(opts.input_pi) == 0:
        print("The specified input pi file ("+opts.input_pi+") does not exist\n")
        exit(-1)

    input_fits   = fits.open(opts.input_pi)
    input_pi     = input_fits[1]
    #  input_pi_hdr = input_pi.header.ascardlist()
    input_pi_hdr = input_pi.header


    # Check the output file exists and is well specified
    if opts.output_pi is None:
        print("The output pi file is missing\n")
        parser.print_help()
        exit(-1)

    # Check the RMF file exists and is well specified
    if opts.rmf is None:
        if 'RESPFILE' in input_pi_hdr:
            opts.rmf = input_pi_hdr['RESPFILE']
            print("Found RESPFILE: {}".format(opts.rmf))
        if opts.rmf is None:
            print("No rmf file is specified in either input pi file or command line\n")
            parser.print_help()
            exit(-1)

    # Check the background file exists and is well specified
    if opts.background is None:
        if 'BACKFILE' in input_pi_hdr:
            opts.background = input_pi_hdr['BACKFILE']
            print("Found BACKFILE: {}".format(opts.background))
        if opts.background is None:
            print("*WARNING*: The grouped spectra has no background file\n")

    # Check the ARF file exists and is well specified
    if opts.arf is None:
        if 'ANCRFILE' in input_pi_hdr:
            opts.arf = input_pi_hdr['ANCRFILE']
            print("Found ANCRFILE: {}".format(opts.arf))
        if opts.arf is None:
            print("*WARNING*: The grouped spectra has no arf file\n")

    input_fits.close()





def count_source_counts(opts, pi_data, ebounds):
    """Count total source counts within the specified energy range without any grouping."""
    count_index    = np.where(np.array(pi_data.columns.names).astype(str) == 'COUNTS')[0][0]
    channels       = ebounds.data.field('CHANNEL')
    e_min          = ebounds.data.field('E_MIN') * 1000
    e_max          = ebounds.data.field('E_MAX') * 1000
    tot_counts     = 0

    for row in pi_data:
        pi_row         = row[0]
        pi_counts      = row[count_index]
        matching_index = (np.where(channels == pi_row))[0][0]
        if (e_min[matching_index] > opts.lower_energy) and (opts.upper_energy > e_max[matching_index]):
            tot_counts += pi_counts

    return tot_counts




def group_pha(opts, pi_data, ebounds):
    """
    Build QUALITY and GROUPING arrays for the spectrum.

    QUALITY: 5 = ignored (outside energy range), 0 = good
    GROUPING: 1 = start of a new bin, -1 = continuation of previous bin

    Bins are merged until both the count threshold (bin_min_counts) and the
    energy-width threshold (bin_min_energy, or fractional) are satisfied.
    The final incomplete bin at the upper boundary is merged back into the
    previous bin if it falls short of the count threshold.
    """

    # print("Pi data shape: ", pi_data.shape) # (1024,)
    pi_rows = pi_data.shape[0] 
    print("Pi rows: ", pi_rows) # there are 1024 detector channels

    # Seed quality from the input file so pipeline flags (bad channels, pile-up,
    # etc.) are preserved.  Fall back to all-5 (ignored) if no column exists.
    if 'QUALITY' in pi_data.columns.names:
        quality = pi_data['QUALITY'].astype(float).copy()
        unique_q = np.unique(quality).astype(int).tolist()
        print(f'Found QUALITY column with unique values: {unique_q}')
    else:
        quality = np.zeros(pi_rows) + 5 # quality of data; set to 5 at the start, which means bad
        print('No QUALITY column found in input; defaulting all channels to 5 (ignored)')
    grouping = np.zeros(pi_rows) + 1  # default: every channel starts its own bin. to track the grouping; set all to 1 at the start
    

    channels            = ebounds.data.field('CHANNEL')
    e_min               = ebounds.data.field('E_MIN')*1000 # shape (1024,)
    e_max               = ebounds.data.field('E_MAX')*1000 # shape (1024,)
    energy              = np.sqrt(e_min*e_max)
    energy_width        = e_max-e_min # width in energy of each each channel
    

    # Minimum energy width required per bin (absolute eV, or fractional if flag set).
    # Always sized to pi_rows so it can be safely indexed by row_counter.
    min_width_at_energy = np.zeros(pi_rows) + opts.bin_min_energy # set minimum energy change per channel
    if (opts.fraction_switch): # using fractional energy binning
        #Peviously: min_width_at_energy = opts.bin_min_energy * energy # scale bin width by energy
        # energy is indexed by matching_index (ebounds channel), not row_counter,
        # so the per-row threshold is looked up via matching_index in the loop.
        pass  # handled below in the loop via energy[matching_index]
        
    # Initiate values
    tot_counts          = 0
    row_counter         = 0
    last_new_row        = 0
    counts_in_bin       = opts.bin_min_counts # initialise as satisfied so first in-range channel starts a new bin
    energy_width_of_bin = np.sum(energy_width) # i.e. total energy of the whole spectrum
    
    all_counts=[] # for testing

    # AKH added this to find the index that corresponded to counts
    count_index = np.where(np.array(pi_data.columns.names).astype(str) == 'COUNTS')[0][0] # = 1
    
    # For each channel in pi_data:
    for row in (pi_data): # row is type FITS_record
        
        pi_counts = row[count_index] # This was originally an index of 2, potentially a difference between chandra v. swift
        all_counts.append(pi_counts) # for testing

        pi_row = row[0] # channel number
        matching_index = (np.where(channels==pi_row))[0][0]

        # If the energy of the current channel falls in the specified energy range
        if ((e_min[matching_index] > opts.lower_energy) and (opts.upper_energy > e_max[matching_index])):
            # In-range: mark good unless the input pipeline already flagged this
            # channel as bad (quality > 0), in which case preserve that flag.
            if quality[row_counter] == 0 or quality[row_counter] == 5:
                quality[row_counter] = 0 # quality is marked as good (0)
            tot_counts += pi_counts # total counts for that channel are added to tot_counts

            # Resolve the per-channel energy threshold here so it is always
            # looked up from ebounds via matching_index, avoiding any size mismatch.
            if opts.fraction_switch:
                min_width_at_energy[row_counter] = opts.bin_min_energy * energy[matching_index]

            # Check whether the counts in the bin and bin width meet the minimum requirements (i.e. previous bin is complete).
            # If so, a new group/bin is started.
            if((counts_in_bin >= opts.bin_min_counts) and (energy_width_of_bin >= min_width_at_energy[row_counter])):
                counts_in_bin         = pi_counts
                energy_width_of_bin   = energy_width[row_counter]
                last_new_row          = row_counter
                grouping[row_counter] =1 # indicates end of one group
            # Otherwise, the counts and energy for the current channel are added to the existing group/bin.
            # i.e. keep accumulating into the current bin
            else:
                counts_in_bin         += pi_counts
                energy_width_of_bin   += energy_width[row_counter]
                grouping[row_counter] =-1 # i.e. this channel is part of an ongoing group/bin
        
        else:
            # Out-of-range: force-ignore regardless of what the input file had
            quality[row_counter] = 5

        row_counter += 1
    
    # If the final bin never reached the count threshold, merge it back into
    # the previous complete bin so Xspec doesn't see an under-filled group.
    if counts_in_bin < opts.bin_min_counts:
        grouping[last_new_row] =-1 # for last channel

    # Testing:
    #print()
    #print("ALL COUNTS:")
    #ar = np.array(all_counts)
    #print(ar)
    #print(len(all_counts))
    #print(np.sum(all_counts))
    #indices = np.where(ar == 1)
    #print(indices[0])
    
    return (quality, grouping, tot_counts)


def write_passthrough(opts, input_fits, tot_counts):
    """
    Write the input spectrum to the output path without altering QUALITY or
    GROUPING.  Only the per-observation header keywords are updated.
    """
    print(f'Total Counts: {tot_counts}  — writing file unchanged (no grouping applied)')

    output_fits = input_fits.copy()
    output_fits[1].header['BACKFILE'] = opts.background
    output_fits[1].header['ANCRFILE'] = opts.arf
    output_fits[1].header['COUNTGRP'] = opts.bin_min_counts
    output_fits[1].header['COUNTTOT'] = tot_counts

    output_fits.writeto(opts.output_pi, overwrite=True)
    output_fits.close()



## Wrapper to group counts and output an updated FITS file
def fitsio_grppha(opts):

    # Open PHA FITS file containing the photon counts
    input_fits   = fits.open(opts.input_pi)
    rmf_fits   = fits.open(opts.rmf)
    # Open RMF which contains energy information
    ebounds = rmf_fits['EBOUNDS']
    counts_threshold = opts.counts_threshold

    # Count source counts first so we can decide whether grouping is worthwhile
    tot_counts = count_source_counts(opts, input_fits[1].data, ebounds)
    print("Total counts: ", tot_counts)


    # If single-count binning was requested, or the spectrum is too faint to
    # group meaningfully, just pass the file through untouched.
    if opts.bin_min_counts == 1 or tot_counts < counts_threshold:
        opts.bin_min_counts = 1  # record what we actually used in the header
        write_passthrough(opts, input_fits, tot_counts)
        rmf_fits.close()
        input_fits.close()
        return


    # Spectrum has enough counts: compute and attach QUALITY / GROUPING columns    
    # Use the group_pha function to group counts
    (quality, grouping, tot_counts) = group_pha(opts, input_fits[1].data, ebounds)


    # Testing:
    # print("Output after binning:")    
    # print("QUALITY: ", quality)
    # print("GROUPING: ", grouping)


    # Remove any pre-existing QUALITY / GROUPING columns before adding new ones
    if 'QUALITY' in input_fits[1].columns.names:
        print('Deleting existing QUALITY column')
        input_fits[1].columns.del_col('QUALITY')
    if 'GROUPING' in input_fits[1].columns.names:
        print('Deleting existing GROUPING column')
        input_fits[1].columns.del_col('GROUPING')

    
    # Create new columns for quality and grouping
    quality = fits.Column(name='QUALITY', format='I', array=quality)
    grouping = fits.Column(name='GROUPING', format='I', array=grouping)


    # Make a copy of the input FITS file for the output
    # this is a reference, not a unique object, despite it claiming to be a unique object
    # input_fits.copy() returns a reference in astropy, not a true deep copy,
    # so we build the output HDU explicitly from the (now-modified) column set.
    output_fits = input_fits.copy() 
    # Add the new quality and grouping columns to the first HDU
    output_fits[1] = fits.BinTableHDU.from_columns(input_fits[1].columns + quality + grouping, header=input_fits[1].header)
    # print("Testing:",  output_fits[1].header['RESPFILE']) # the respfile is already saved in the header
    # Update the FITS headers with additional information
    output_fits[1].header['BACKFILE'] = opts.background
    output_fits[1].header['ANCRFILE'] = opts.arf
    output_fits[1].header['COUNTGRP'] = opts.bin_min_counts
    output_fits[1].header['COUNTS'] = tot_counts
    
    # Write the modified FITS file
    output_fits.writeto(opts.output_pi, overwrite=True)
    output_fits.close()
    rmf_fits.close()
    input_fits.close()


def main(argv):

    # Parse the arguments
    parser = parse_args(argv) # initialise the parser
    (opts, args) = parser.parse_args() # opts: an object that contains the parsed options as attributes; args: a list of positional arguments not associated with options.
    
    # Check that the inputted arguments are valid
    valid_args_check(opts, parser) 
    
    # Run wrapper method to group counts and output file
    fitsio_grppha(opts) 

if __name__ == "__main__":
    main(sys.argv[1:])
