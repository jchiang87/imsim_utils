import os
import sys
import socket
import sqlite3
import pandas as pd
from desc_roman_sims import GalSimJobGenerator
from desc_roman_sims.parsl.parsl_config import load_wq_config

cores_per_node = 54
nproc = 9
max_threads = cores_per_node // nproc  # 1 galsim instance per thread

hostname = socket.gethostname()
username = "jchiang"
run_dir = f'runinfo/{hostname}_{username}'

def full_path(basename):
    return os.path.join("/sdf/data/rubin/user/jchiang/ops_rehearsal_3",
                        "vignetted_flats", basename)

calib_yaml = dict(BIAS="bias.yaml", DARK="dark.yaml", FLAT="sky_flat.yaml",
                  PTC="ptc.yaml")
calib_type = sys.argv[1]

imsim_yaml = full_path(calib_yaml[calib_type])

# Read in the visit list
with sqlite3.connect(full_path("bright_sky_calib_frames.db")) as con:
    df = pd.read_sql(f"select * from observations where target='{calib_type}'", con)
visits = df['observationId'].to_numpy()
print(visits)

target_dets = {_: set(range(9)) for _ in visits}

load_wq_config(max_threads=max_threads, use_work_queue=False, monitor=False,
               run_dir=run_dir)

generator = GalSimJobGenerator(imsim_yaml, visits,
                               nfiles=9, nproc=nproc,
                               target_dets=target_dets,
                               GB_per_CCD=3, GB_per_PSF=8,
                               clean_up_atm_psfs=False,
                               bash_app_executor='thread_pool',
                               log_dir='logging',
                               raw_prefix=calib_type.lower(),
                               verbosity=1)
generator.run()
