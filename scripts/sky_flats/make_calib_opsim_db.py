import os
import sqlite3
import click
import numpy as np
import pandas as pd

# Read in the OR3 opsin db file and use the smallest skyBrightness magnitude
# entry as the template for bright sky flats.
opsim_db_file = "/sdf/data/rubin/shared/ops-rehearsal-3/scheduler_sims/ops_rehearsal_apr_2024_v3.db"
assert os.path.isfile(opsim_db_file)
with sqlite3.connect(opsim_db_file) as con:
    df0 = pd.read_sql("select * from observations", con)

# Make a temploate dict with default entries from that bright sky entry.
row = dict(df0.query(f"skyBrightness == {min(df0['skyBrightness'])}").iloc[0])

# Create a new DataFrame with the calibration frame entries.
df = pd.DataFrame(columns=df0.columns)

mjd = row['observationStartMJD']
visit = 4000000

# Flat frames. The exposure time has been adjusted to produce the
# ~desired number of counts per pixel.
exptime = 60.
dt = (exptime + 10.)/86400.
bands = ['g', 'r', 'i']
nflat = 20
for band in bands:
    for i in range(nflat):
        row.update(dict(observationId=visit,
                        observationStartMJD=mjd,
                        numExposures=1,
                        visitExposureTime=exptime,
                        filter=band,
                        target='FLAT'))
        df.loc[len(df)] = row
        visit += 1
        mjd += dt

# Write out the new DataFrame as an sqlite file.
outfile = "bright_sky_calib_frames.db"
if os.path.isfile(outfile):
    if click.confirm(f"Overwrite {outfile}?", default=True):
        os.remove(outfile)

with sqlite3.connect(outfile) as con:
    df.to_sql("observations", con, index=False)
