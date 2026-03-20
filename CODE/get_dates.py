# module load heasoft/6.33


import os
import csv
import re
import subprocess
from astropy.io import fits
from astropy.time import Time
import numpy as np



    
path = "../spectra_swift_xrt/"
output_csv = "xray_dates.csv"
search_path = path


spectral_file_ending = "source.pi" 


def get_mjd_from_swifttime(met_time):
    """
    Run swifttime interactively and extract converted MJD.
    """
    inputs = f"{met_time}\nMET\ns\nUTC\nm\n"

    process = subprocess.Popen(
        ['swifttime'],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )

    output, _ = process.communicate(inputs)

    for line in output.splitlines():
        if "Converted time:" in line:
            try:
                return float(line.split(":")[1].strip())
            except ValueError:
                raise RuntimeError(f"Could not parse MJD from: {line}")

    raise RuntimeError("Failed to extract MJD from swifttime output.")





output_rows = []
for root, dirs, files in os.walk(search_path):
    for file in files:
        if file.endswith(spectral_file_ending):
            filepath = os.path.join(root, file)

            try:
                with fits.open(filepath) as hdul:
                    hdr = hdul[0].header

                    # --- Extract UTC-based info ---
                    date_start_str = hdr.get("DATE-OBS")
                    date_end_str = hdr.get("DATE-END")
                    exposure_sec = hdr.get("EXPOSURE", 0)


                    if not (date_start_str and date_end_str):
                        continue  # skip files without date strings

                    t_start_1 = Time(date_start_str, format='isot', scale='utc').mjd
                    t_end_1   = Time(date_end_str, format='isot', scale='utc').mjd
                    mid_mjd_1 = 0.5 * (t_start_1 + t_end_1)
                    dt_1 = t_end_1 - t_start_1
                    
                    exp_days = exposure_sec / 86400.0 # 24*60*60

                    # --- Extract MET-based info ---
                    tstart_met = hdr.get("TSTART")
                    tstop_met = hdr.get("TSTOP")

                    if tstart_met is None or tstop_met is None:
                        continue  # skip if TSTART/TSTOP missing

                    #tmid_met = 0.5 * (tstart_met + tstop_met)
                    #dt_met = (tstop_met - tstart_met)/ 86400.0

                    # Use swifttime to convert to MJD
                    try:
                        #mid_mjd_start = get_mjd_from_swifttime(tmid_met)
                        mjd_start = get_mjd_from_swifttime(tstart_met)
                        mjd_end = get_mjd_from_swifttime(tstop_met)
                        mid_mjd_2 = 0.5*(mjd_start + mjd_end)
                        dt_2 = (mjd_end - mjd_start)
                    except Exception as e:
                        print(f"swifttime error for {filepath}: {e}")
                        #mid_mjd_tt = None
                        mid_mjd_2 = None
                        dt_2 = None

                    # Extract observation ID
                    #match = re.search(r"Obs_(\w+)source\.pi$", file)
                    #obs_id = match.group(1) if match else "UNKNOWN"
                    obs_id = hdr.get("OBS_ID", "UNKNOWN")

                    output_rows.append({
                        "filepath": filepath,
                        "obs_id": obs_id,
                        "mid_mjd_1": mid_mjd_1,
                        "dt_days_1": dt_1,
                        "mid_mjd_2": mid_mjd_2,
                        "dt_days_2": dt_2 ,
                        "exp_days": exp_days, 
                        "mid_time_diff_min": np.abs(mid_mjd_1 - mid_mjd_2)*24*60,  # in minutes
                        "start_mjd": mjd_start
                    })

            except Exception as e:
                print(f"Error reading {filepath}: {e}")

# Sort by UTC midpoint for consistency
output_rows = sorted(output_rows, key=lambda row: row["mid_mjd_1"])

# Write to CSV
with open(output_csv, mode='w', newline='') as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "filepath", "obs_id",
            "mid_mjd_1", "dt_days_1",
            "mid_mjd_2", "dt_days_2",
            "exp_days" , "mid_time_diff_min", "start_mjd"
        ]
    )
    writer.writeheader()
    writer.writerows(output_rows)

print(f"Summary written to {output_csv}")


