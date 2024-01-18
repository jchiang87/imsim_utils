import os
import sqlite3
import pandas as pd

# Read in an existing opsim_db file to provide a template for the
# columns in each entry in the calibration exposure db.
opsim_db_file = "/home/jchiang/work/DESC/rubin_sim_data/baseline_v3.2_10yrs.db"
assert os.path.isfile(opsim_db_file)
with sqlite3.connect(opsim_db_file) as con:
    df = pd.read_sql("select * from observations limit 10", con)

# Create a template dict with default entries from the first row.
row = dict(df.iloc[0])

# Fill a new DataFrame with the calibration frame entries.
df_new = pd.DataFrame(columns=df.columns)

# Bias frames
bias_id_offset = 3000000
mjd0 = 60326.00000000
dt = 10./86400.
for i in range(20):
    visit = bias_id_offset + i
    mjd = mjd0 + dt*i
    row.update(dict(observationId=visit,
                    observationStartMJD=mjd,
                    visitExposureTime=0.,
                    target='BIAS'))
    df_new.loc[len(df_new)] = row

# Dark frames
dark_id_offset = 4000000
mjd0 = 60327.00000000
dt = 10./86400.
for i in range(20):
    visit = bias_id_offset + i
    mjd = mjd0 + dt*i
    row.update(dict(observationId=visit,
                    observationStartMJD=mjd,
                    visitExposureTime=100.,
                    target='DARK'))
    df_new.loc[len(df_new)] = row

# Write out the new DataFrame as an sqlite file.
outfile = "calib_frames.db"
with sqlite3.connect(outfile) as con:
    df_new.to_sql("observations", con)
