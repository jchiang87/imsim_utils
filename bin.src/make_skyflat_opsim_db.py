#!/usr/bin/env python
import os
import yaml
import sqlite3
import argparse
import click
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("config_file", type=str,
                    help="Config file for containing simulation parameters")
config_file = parser.parse_args().config_file

with open(config_file) as fobj:
    config = yaml.safe_load(fobj)

# Read in the opsim db file and use the smallest skyBrightness magnitude
# entry as the template for bright sky flats.
opsim_db_file = config["template_db_file"]
assert os.path.isfile(opsim_db_file)
with sqlite3.connect(opsim_db_file) as con:
    df0 = pd.read_sql("select * from observations", con)

# Make a template dict with default entries from the brightest sky entry.
row = dict(df0.query(f"skyBrightness == {min(df0['skyBrightness'])}").iloc[0])

# Create a new DataFrame with the calibration frame entries.
df = pd.DataFrame(columns=df0.columns)

mjd = row['observationStartMJD']
visit = 4000000

# Flat frames. The exposure time has been adjusted to produce the
# ~desired number of counts per pixel.
exptime = 60.
dt = (exptime + 10.)/86400.
bands = config['bands']
nflat = config['nflat']
for band in bands:
    for i in range(nflat):
        row.update(dict(observationId=visit,
                        observationStartMJD=mjd,
                        numExposures=1,
                        visitExposureTime=exptime,
                        filter=band,
                        target_name='FLAT'))
        df.loc[len(df)] = row
        visit += 1
        mjd += dt

# Write out the new DataFrame as an sqlite file.
outfile = "sky_flat_calib_frames.db"
if os.path.isfile(outfile):
    if click.confirm(f"Overwrite {outfile}?", default=True):
        os.remove(outfile)

with sqlite3.connect(outfile) as con:
    df.to_sql("observations", con, index=False)
