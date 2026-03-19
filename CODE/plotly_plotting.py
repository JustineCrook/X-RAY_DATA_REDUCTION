import plotly.graph_objects as go
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

import sys
import os
sys.path.append(os.path.abspath("./"))
from get_results_data import get_maxi_counts
from astropy.time import Time


####################################################################################################################
## HID PLOTTING FUNCTIONS
###################################################################################################################


def plot_maxi_hid(mjd_low, mjd_high, plotly= True, show_errobars = True, mjds_to_point_to = None, xlim = (1e-1,3e0), ylim = (2e-2,2e-1), hr_unc_cutoff = 0.2):

    # MAXI 
    (mjds, hr, hr_unc,
    cps_2_4, cps_2_4_unc,
    cps_4_10, cps_4_10_unc,
    cps_2_20, cps_2_20_unc,
    cps_10_20, cps_10_20_unc) = get_maxi_counts()

    # Mask for the time range of interest
    time_mask = (mjds >= mjd_low) & (mjds <= mjd_high)

    # Filter high hr
    hr_mask = hr_unc <= hr_unc_cutoff

    mask = time_mask & hr_mask

    # Filter the data based on the time mask
    mjds = mjds[mask]
    hr = hr[mask]
    hr_unc = hr_unc[mask]
    cps_2_4 = cps_2_4[mask]
    cps_2_4_unc = cps_2_4_unc[mask]
    cps_4_10 = cps_4_10[mask]
    cps_4_10_unc = cps_4_10_unc[mask]
    cps_2_20 = cps_2_20[mask]
    cps_2_20_unc = cps_2_20_unc[mask]
    cps_10_20 = cps_10_20[mask]
    cps_10_20_unc = cps_10_20_unc[mask]  


    #########################
    ## Plot HR time series using matplotlib

    plt.figure(figsize=(14, 5))
    plt.errorbar(mjds, hr, yerr=hr_unc, fmt='o', markersize=4,capsize=3, label='Hardness ratio')
    plt.xlabel('MJD')
    plt.ylabel('Hardness ratio [(4-10 keV) / (2-4 keV) count rates]')
    plt.grid(True)
    plt.legend()
    plt.ylim(0.05, 2) 
    plt.yscale('log')
    plt.tight_layout()
    plt.show()



    #########################
    ## Plot HID using matplotlib
    
    plt.figure(figsize=(10, 6))

    # Plot faint dotted line between points
    plt.plot(hr, cps_2_20, linestyle='--', linewidth =0.9, color='gray', alpha=0.5, zorder=1)

    # Scatter plot on top
    sc = plt.scatter(hr, cps_2_20, c=mjds, cmap='viridis', s=60, edgecolor='k', zorder=2)

    # Add error bars
    if show_errobars:
        plt.errorbar(hr,cps_2_20 , xerr=hr_unc, yerr=cps_2_20_unc, fmt='none', ecolor='gray', alpha=0.3, capsize=1)


    # Add arrows pointing to selected MJDs
    if mjds_to_point_to is not None:
        for mjd_target in mjds_to_point_to:
            # Find the index of the closest MJD in the data
            idx = np.argmin(np.abs(mjds - mjd_target))
            x, y = hr[idx], cps_2_20[idx]
            closest_mjd = mjds[idx]

            # Coordinates to start the arrow from (you can tweak offset values)
            x_start = x * 0.85
            y_start = y * 1.5

            # Draw arrow and label
            plt.annotate(f'MJD {closest_mjd:.1f}',
                         xy=(x, y),
                         xytext=(x_start, y_start),
                         arrowprops=dict(arrowstyle='->', color='black', lw=2),
                         fontsize=12, color='black',
                         ha='center')
            

    cbar = plt.colorbar(sc)
    cbar.set_label('MJD', fontsize=15)
    cbar.ax.tick_params(labelsize=15)

    plt.xlabel('Hardness ratio [(4-10 keV) / (2-4 keV) \ncount rates]', fontsize=15)
    plt.ylabel('2–20 keV count rate', fontsize=15)
    #plt.grid(True)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.tight_layout()
    plt.show()



    #########################
    ## Plot HR time series using plotly

    if plotly:

        fig = go.Figure()

        fig.add_trace(go.Scatter(
            x=mjds,
            y=hr,
            mode='markers',
            name='Hardness ratio [(4-10 keV) / (2-4 keV) count rates]',
            error_y=dict(
                type='data',
                array=hr_unc,
                visible=True,
                thickness=1.5,
                width=3,
            ),
            marker=dict(size=6)
        ))

        fig.update_layout(
            xaxis_title='MJD',
            yaxis_title='Hardness ratio [(4-10 keV) / (2-4 keV) count rates]',
            yaxis_type='log',  
            template='simple_white',
            legend=dict(title=''),
            width=1200,
            height=600
        )

        fig.show()




    #########################
    ## Plot HID using plotly 

    if plotly: 

        # Convert MJD to UTC date strings
        dates = Time(mjds, format='mjd').iso  

        trace = go.Scatter(
            x=hr,
            y=cps_2_20,
            mode='markers',
            marker=dict(
                color=mjds,
                colorscale='Viridis',
                colorbar=dict(title='MJD'),
                size=6,
                line=dict(width=0.5, color='black')
            ),
            customdata=np.stack((mjds, dates), axis=-1),
            hovertemplate=(
                "Hardness: %{x:.3f}<br>"
                "Count Rate: %{y:.3f}<br>"
                "MJD: %{customdata[0]}<br>"
                "Date: %{customdata[1]}<extra></extra>"
            ),
            name='Data'
        )

        fig = go.Figure(data=[trace])

        # Add error bars if requested
        if show_errobars:
            fig.data[0].update(
                error_x=dict(
                    type='data',
                    array=hr_unc,
                    visible=True,
                    thickness=1.5,
                    width=2
                ),
                error_y=dict(
                    type='data',
                    array=cps_2_20_unc,
                    visible=True,
                    thickness=1.5,
                    width=2
                )
            )

        fig.update_layout(
            title='Hardness-Intensity Diagram (Coloured by MJD)',
            xaxis=dict(
                title='Hardness ratio [(4–10 keV) / (2–4 keV) count rates]',
                type='log',
                range=[np.log10(1e-1), np.log10(3)]
            ),
            yaxis=dict(
                title='2–20 keV count rate',
                type='log', 
                range=[np.log10(2e-2), np.log10(2e-1)]
            ),
            width=800,
            height=600,
            template='simple_white'
        )

        fig.show()




####################################################################################################################
## BINNING COMPARISON PLOTTING FUNCTIONS
###################################################################################################################


def plot_binning_comparison(plotly = True):
    """
    Function to compare the fluxes obtained from spectral fitting with binning (min_counts_chi = 20) and without binning (min_counts_chi = 1) for the same set of spectra. It reads in the results from the two different fitting runs, extracts the relevant data, and creates both a matplotlib plot and a plotly plot to visualize the comparison.
    """

    file1 = Path(__file__).resolve().parent.parent / "spectral_fit_final_results_bin" / "final_fit_selection.csv"
    
    data1 = pd.read_csv(file1)
    ids1 = data1['obs_id'].to_numpy()
    t1 = data1['middle_mjds'].to_numpy()
    flux1 = data1['flux'].to_numpy()
    flux_er_neg1 = data1['flux_er_neg'].to_numpy()
    flux_er_pos1 = data1['flux_er_pos'].to_numpy()


    file2 = Path(__file__).resolve().parent.parent / "spectral_fit_final_results_bin1" / "final_fit_selection.csv"

    data2 = pd.read_csv(file2)
    ids2 = data2['obs_id'].to_numpy()
    t2 = data2['middle_mjds'].to_numpy()
    flux2 = data2['flux'].to_numpy()
    flux_er_neg2 = data2['flux_er_neg'].to_numpy()
    flux_er_pos2 = data2['flux_er_pos'].to_numpy()


    plt.figure(figsize=(12,5))
    plt.errorbar(t1[:], flux1[:], yerr=[flux_er_neg1[:],flux_er_pos1[:]], fmt='.', label="binning", ms=2)
    plt.errorbar(t2[:], flux2[:], yerr=[flux_er_neg2[:],flux_er_pos2[:]], fmt='.', label="bin1", ms=2)
    plt.yscale('log')
    plt.xlabel("Time")
    plt.ylabel("Flux")
    plt.ylim([1e-13,5e-8])
    plt.legend()
    plt.show()


    ## Plot plotly plot

    if plotly:

        fig = go.Figure()

        # --- all grades ---
        fig.add_trace(
            go.Scatter(
                x=t1,
                y=flux1,
                mode="markers",
                name="binning",
                marker=dict(size=4),
                error_y=dict(
                    type="data",
                    symmetric=False,
                    array=flux_er_pos1,   # positive errors
                    arrayminus=flux_er_neg1,  # negative errors
                    visible=True,
                ),
            )
        )

        # --- restricted grades ---
        fig.add_trace(
            go.Scatter(
                x=t2,
                y=flux2,
                mode="markers",
                name="bin1",
                marker=dict(size=4),
                error_y=dict(
                    type="data",
                    symmetric=False,
                    array=flux_er_pos2,
                    arrayminus=flux_er_neg2,
                    visible=True,
                ),
            )
        )

        # --- layout ---
        fig.update_layout(
            width=1500,
            height=600,
            xaxis_title="Time",
            yaxis_title="Flux",
            yaxis_type="log",
            yaxis_range=[-13, -7.3],  # log10(1e-13) to log10(5e-8)
            legend=dict(x=0.01, y=0.99),
            template="simple_white",
        )

        fig.show()




