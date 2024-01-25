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
    df = pd.read_sql("select * from observations limit 10", con)

# Create a template dict with default entries from the first row.
row = dict(df.iloc[0])

# Fill a new DataFrame with the calibration frame entries.
df_new = pd.DataFrame(columns=df.columns)

# Bias frames
bias_id_offset = 3000000
mjd0 = 60326.    # 2024-01-17 00:00:00
dt = 10./86400.
for i in range(20):
    visit = bias_id_offset + i
    mjd = mjd0 + dt*i
    row.update(dict(observationId=visit,
                    observationStartMJD=mjd,
                    visitExposureTime=0.,
                    numExposures=1,
                    target='BIAS'))
    df_new.loc[len(df_new)] = row

# Dark frames
dark_id_offset = 4000000
mjd0 = 60327.
dt = 10./86400.
for i in range(20):
    visit = dark_id_offset + i
    mjd = mjd0 + dt*i
    row.update(dict(observationId=visit,
                    observationStartMJD=mjd,
                    numExposures=1,
                    visitExposureTime=100.,
                    target='DARK'))
    df_new.loc[len(df_new)] = row

# Flat frames
flat_id_offset = 5000000
mjd_last = 60328.
dt = 10./86400.
num_pairs = 20
exptimes = np.logspace(np.log10(0.3), np.log10(1e3), num_pairs)
for i, exptime in enumerate(exptimes):
    for j in (0, 1):
        visit = flat_id_offset + i*2 + j
        mjd = mjd_last + exptime + dt
        row.update(dict(observationId=visit,
                        observationStartMJD=mjd,
                        numExposures=1,
                        visitExposureTime=exptime,
                        target='FLAT'))
        df_new.loc[len(df_new)] = row
        mjd_last = mjd

# Write out the new DataFrame as an sqlite file.

outfile = "calib_frames.db"
if os.path.isfile(outfile):
    if click.confirm(f"Overwrite {outfile}?", default=True):
        os.remove(outfile)

with sqlite3.connect(outfile) as con:
    df_new.to_sql("observations", con, index=False)
