import os
import sqlite3
import click
import numpy as np
import pandas as pd

# Read in an existing opsim_db file to provide a template for the
# columns in each entry in the calibration exposure db.
opsim_db_file = "/sdf/data/rubin/user/jchiang/imSim/rubin_sim_data/opsim_cadences/baseline_v3.2_10yrs.db"
assert os.path.isfile(opsim_db_file)
with sqlite3.connect(opsim_db_file) as con:
    df0 = pd.read_sql("select * from observations limit 10", con)

# Make a template dict with default entries from the first row.
row = dict(df0.iloc[0])

# Create a new DataFrame with the calibration frame entries.
df = pd.DataFrame(columns=df0.columns)

mjd = 60326.    # 2024-01-17 00:00:00
visit = 3000000

# Bias frames
exptime = 0
dt = (exptime + 10.)/86400.
nbias = 20
for i in range(nbias):
    row.update(dict(observationId=visit,
                    observationStartMJD=mjd,
                    visitExposureTime=exptime,
                    numExposures=1,
                    target='BIAS'))
    df.loc[len(df)] = row
    visit += 1
    mjd += dt

# Dark frames
exptime = 100.
dt = (exptime + 10.)/86400.
ndark = 20
for i in range(ndark):
    row.update(dict(observationId=visit,
                    observationStartMJD=mjd,
                    numExposures=1,
                    visitExposureTime=exptime,
                    target='DARK'))
    df.loc[len(df)] = row
    visit += 1
    mjd += dt

# Flat frames, assume 100 counts/pixel illumination
exptime = 500.
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

# PTC frames, assume 100 counts/pixel illumination
num_pairs = 100
exptimes = np.logspace(np.log10(0.3), np.log10(1e3), num_pairs)
for i, exptime in enumerate(exptimes):
    dt = (exptime + 10.)/86400.
    for j in (0, 1):
        row.update(dict(observationId=visit,
                        observationStartMJD=mjd,
                        numExposures=1,
                        visitExposureTime=exptime,
                        target='PTC'))
        df.loc[len(df)] = row
        visit += 1
        mjd += dt


# Write out the new DataFrame as an sqlite file.
outfile = "calib_frames.db"
if os.path.isfile(outfile):
    if click.confirm(f"Overwrite {outfile}?", default=True):
        os.remove(outfile)

with sqlite3.connect(outfile) as con:
    df.to_sql("observations", con, index=False)
