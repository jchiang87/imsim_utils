#!/usr/bin/env python
import os
import sqlite3
import argparse
import yaml
import click
import numpy as np
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("config_file", type=str,
                    help="Config file containing simulation parameters")
config_file = parser.parse_args().config_file

with open(config_file) as fobj:
    config = yaml.safe_load(fobj)

# Read in an existing opsim_db file to provide a template for the
# columns in each entry in the calibration exposure db.
template_db_file = config['template_db_file']
assert os.path.isfile(template_db_file)
with sqlite3.connect(template_db_file) as con:
    df0 = pd.read_sql("select * from observations limit 10", con)

# Make a template dict with default entries from the first row.
row = dict(df0.iloc[0])

# Create a new DataFrame with the calibration frame entries.
df = pd.DataFrame(columns=df0.columns)

mjd = config['mjd_start']
visit = config['visit_start']

# Bias frames
exptime = 0
dt = (exptime + 10.)/86400.
nbias = config['nbias']
for i in range(nbias):
    row.update(dict(observationId=visit,
                    observationStartMJD=mjd,
                    visitExposureTime=exptime,
                    numExposures=1,
                    target_name='BIAS'))
    df.loc[len(df)] = row
    visit += 1
    mjd += dt

# Dark frames
exptime = config['dark_time']
dt = (exptime + 10.)/86400.
ndark = config['ndark']
for i in range(ndark):
    row.update(dict(observationId=visit,
                    observationStartMJD=mjd,
                    numExposures=1,
                    visitExposureTime=exptime,
                    target_name='DARK'))
    df.loc[len(df)] = row
    visit += 1
    mjd += dt

# Sky flat frames.  The sky level is set in sky_flat.yaml.
exptime = config['flat_exptime']
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

# PTC frames, assuming 100 counts/pixel illumination is set in ptc.yaml.
num_pairs = config['num_pairs']
ptc_exptime_min = float(config['ptc_exptime_min'])
ptc_exptime_max = float(config['ptc_exptime_max'])
exptimes = np.logspace(np.log10(ptc_exptime_min), np.log10(ptc_exptime_max),
                       num_pairs)
band = config['ptc_band']
for i, exptime in enumerate(exptimes):
    dt = (exptime + 10.)/86400.
    for j in (0, 1):
        row.update(dict(observationId=visit,
                        observationStartMJD=mjd,
                        numExposures=1,
                        visitExposureTime=exptime,
                        filter=band,
                        target_name='PTC'))
        df.loc[len(df)] = row
        visit += 1
        mjd += dt


# Write out the new DataFrame as an sqlite file.
outfile = config['outfile']
if os.path.isfile(outfile):
    if click.confirm(f"Overwrite {outfile}?", default=True):
        os.remove(outfile)

with sqlite3.connect(outfile) as con:
    df.to_sql("observations", con, index=False)
