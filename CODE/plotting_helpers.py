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


####################################################################################################################
## Matplotlib formatting 

mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['font.size'] = 14
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['ytick.minor.visible'] = True
mpl.rcParams['ytick.labelsize'] = 'medium'
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['ytick.minor.width'] = 1.5
mpl.rcParams['ytick.major.size'] = 6.0
mpl.rcParams['ytick.minor.size'] = 3.0
mpl.rcParams['xtick.top'] = True
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['xtick.labelsize'] = 'medium'
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['xtick.minor.width'] = 1.5
mpl.rcParams['xtick.major.size'] = 6.0
mpl.rcParams['xtick.minor.size'] = 3.0
mpl.rcParams['axes.linewidth'] = 1.5



####################################################################################################################
## Helper functions

def FormatAxis(ax, mjd, dt=10, interval=60, use_secondary=True):
    """
    Format axes array `ax`:
      - ax[0] is the top datetime axis (UTC)
      - ax[-1] is the main plotting axis (used to create bottom secondary MJD axis)

    This function is idempotent: it removes any previously-created secondary axis
    it stored on ax[-1]._mjd_secondary_axis. It also avoids using plt.gcf().
    """
    if len(mjd) < 2:
        return

    top_ax = ax[0]
    main_ax = ax[-1]

    # Top axis formatting (use local axis only)
    top_ax.set_xlabel('Observing Date (UTC)', fontfamily='serif')
    top_ax.xaxis.set_major_locator(mdates.DayLocator(interval=interval))
    top_ax.set_xlim(Time(mjd[0] - dt, format='mjd').datetime,
                    Time(mjd[-1] + dt, format='mjd').datetime)
    top_ax.xaxis.set_label_position('top')
    xformatter = mdates.DateFormatter('%Y-%m-%d')
    top_ax.xaxis.set_major_formatter(xformatter)
    top_ax.tick_params(axis='x', which='major', rotation=10,
                       labeltop=True, labelbottom=False)

    # Hide main axis bottom labels — we want the secondary axis to display MJD only
    main_ax.tick_params(axis='x', which='both', labelbottom=False)

    # Remove previously-created secondary axis if present
    prev = getattr(main_ax, "_mjd_secondary_axis", None)
    if prev is not None:
        try:
            prev.remove()
        except Exception:
            pass
        main_ax._mjd_secondary_axis = None

    if use_secondary:
        # create the secondary axis and store reference on main axis
        mjd_ax = main_ax.secondary_xaxis('bottom', functions=(plot2mjd, mjd2plot))
        mjd_ax.set_xlabel('Observing Date (MJD)', fontfamily='serif')
        mjd_ax.tick_params(which='major', direction='in', length=0.0, width=0.0)

        # Build tick list from the top axis labels (guard empty labels)
        labels = top_ax.get_xticklabels(which='major')
        mjd_ticks = []
        for lab in labels:
            txt = lab.get_text()
            if txt:
                mjd_ticks.append(txt + 'T00:00:00')

        if mjd_ticks:
            try:
                mjd_ticks_int = (Time(mjd_ticks, format='isot').mjd).astype(int)
                mjd_ax.set_xticks(mjd_ticks_int, labels=mjd_ticks_int)
            except Exception:
                # If parsing fails, just leave auto ticks
                pass

        # store for later removal
        main_ax._mjd_secondary_axis = mjd_ax
    

def plot2mjd(t):
    '''Convert from matplotlib plot date to mjd'''
    return Time(t, format="plot_date", scale='utc').mjd


def mjd2plot(mjd):
    '''Convert from mjd to matplotlib plot'''
    return Time(mjd, format="mjd", scale='utc').plot_date

def mjd2utc(mjd):
    # Convert MJD to UTC using astropy
    t = Time(mjd, format='mjd', scale='utc')
    return t.iso  # Returns in ISO format (YYYY-MM-DD HH:MM:SS.sss)

def iso2mjd(iso_dates):
    # Convert ISO dates to MJD using astropy Time
    times = Time(iso_dates, format='isot', scale='utc')
    return times.mjd


# Input: 2D list -- a list of ranges in the form [start_index, end_index], where the indexes are inclusive
def ranges_okay(ranges):

    # Check validity of each range
    for range_pair in ranges:
        if len(range_pair) != 2:
            return False  # Each range must have exactly 2 values
        if range_pair[0] > range_pair[1]:
            return False  # The first value must be less or equal to than the second

    # Sort ranges based on the start value
    sorted_ranges = sorted(ranges, key=lambda x: x[0])
    #print(sorted_ranges)
    
    for i in range(len(sorted_ranges) - 1):
        # Check if the current range overlaps with the next range
        if sorted_ranges[i][1] > sorted_ranges[i + 1][0]:
            return False  # Overlap found

    return True  # No overlaps and all ranges are valid
