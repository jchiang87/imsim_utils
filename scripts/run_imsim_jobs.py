import os
import sys
import socket
import sqlite3
import pandas as pd
from imsim_utils import GalSimJobGenerator, load_wq_config
from desc_roman_sims.parsl.parsl_config import load_wq_config

cores_per_node = 32
nproc = 4
max_threads = cores_per_node // nproc  # 1 galsim instance per thread

hostname = socket.gethostname()
username = "jchiang"
run_dir = f'runinfo/{hostname}_{username}'

load_wq_config(max_threads=max_threads, use_work_queue=False, monitor=False,
               run_dir=run_dir)

def full_path(basename):
    return os.path.join("/sdf/data/rubin/user/jchiang/ops_rehearsal_3/calibs",
                        basename)

calib_yaml = dict(BIAS="bias.yaml", DARK="dark.yaml", FLAT="flat.yaml",
                  PTC="ptc.yaml")
calib_type = sys.argv[1]

# Remaining visits
incomplete_visits = {'FLAT': [3000047, 3000049, 3000051, 3000052, 3000053,
                              3000055, 3000056, 3000057, 3000058, 3000059,
                              3000060, 3000061, 3000062, 3000063, 3000064,
                              3000065, 3000066, 3000067, 3000068, 3000069,
                              3000070, 3000071, 3000072, 3000073],
                     'PTC': []}

imsim_yaml = full_path(calib_yaml[calib_type])

## Read in the visit list
#with sqlite3.connect(full_path("calib_frames.db")) as con:
#    df = pd.read_sql(f"select * from observations where target='{calib_type}'", con)
#visits = df['observationId'].to_numpy()[-4:]
#print(visits)

visits = incomplete_visits[calib_type]

target_dets = {_: set(range(9)) for _ in visits}
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
