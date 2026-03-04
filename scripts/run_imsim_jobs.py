import os
import sys
import socket
import sqlite3
import pandas as pd
from imsim_utils import GalSimJobGenerator, load_wq_config

calib_yaml = dict(BIAS="bias.yaml", DARK="dark.yaml", FLAT="sky_flat.yaml",
                  PTC="ptc.yaml")
calib_type = sys.argv[1]

# For a full perlmutter node, use 120 cores, reserving the remaining 8
# for workflow management and system-related work:
cores_per_node = 120
#cores_per_node = 4

# Number of processes used by each galsim instance:
nproc = 4

# There will be one galsim instance per thread:
max_threads = cores_per_node // nproc

hostname = socket.gethostname()
username = os.environ['USER']
run_dir = f'runinfo/{hostname}_{username}'

memory = 4000*cores_per_node  # Assume 4 GB per core.
load_wq_config(
    memory=memory,
    max_threads=max_threads,
    # Expect all jobs to be the same size so don't need work_queue:
    use_work_queue=False,
    monitor=False,
    run_dir=run_dir
)

def full_path(basename):
    return os.path.join("/pscratch/sd/j/jchiang8/Euclid_Joint_DDP/calibration/",
                        basename)

imsim_yaml = full_path(calib_yaml[calib_type])

# Read in the visit list
if calib_type == "FLAT":
    opsim_db_file = "sky_flat_calib_frames.db"
else:
    opsim_db_file = "calib_frames.db"

assert os.path.isfile(opsim_db_file)
with sqlite3.connect(full_path(opsim_db_file)) as con:
    df = pd.read_sql("select * from observations where "
                     f"target_name='{calib_type}'", con)

visits = df['observationId'].to_numpy()
print("number of visits:", len(visits))

# target_dets specifies the detectors to simulate for each visit.
target_dets = {_: set(range(90, 99)) for _ in visits}
#target_dets = {_: set([94]) for _ in visits}
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
