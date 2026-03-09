#!/usr/bin/env python
import os
import socket
import argparse
import sqlite3
import pandas as pd
from imsim_utils import GalSimJobGenerator, load_wq_config


def config_path(basename):
    config_abspath = os.path.abspath(args.config_path)
    return os.path.join(config_abspath, basename)


parser = argparse.ArgumentParser()
parser.add_argument("calib_type", type=str, help="Calibration frame type")
parser.add_argument("--opsim_db_file", type=str,
                    help="Calibration opsim db file")
parser.add_argument("--cores_per_node", type=int, default=120,
                    help="Number of cores per node to use")
parser.add_argument("--nproc", type=int, default=4,
                    help="number of processes per galsim instance")
parser.add_argument("--mem_per_core", type=int, default=4000,
                    help="memory per core in MB")
parser.add_argument("--target_dets", type=str, default=None,
                    help="Range of detectors to simulate.")
parser.add_argument("--config_path", type=str, default=".",
                    help="Config path to yaml and opsim_db files.")

args = parser.parse_args()

calib_yaml = dict(BIAS="bias.yaml", DARK="dark.yaml", FLAT="sky_flat.yaml",
                  PTC="ptc.yaml")
calib_type = args.calib_type

# For a full perlmutter node, use 120 cores, reserving the remaining 8
# for workflow management and system-related work:
cores_per_node = args.cores_per_node

# Number of processes used by each galsim instance:
nproc = args.nproc

# There will be one galsim instance per thread:
max_threads = cores_per_node // nproc

hostname = socket.gethostname()
username = os.environ['USER']
run_dir = f'runinfo/{hostname}_{username}'

memory = args.mem_per_core * cores_per_node
load_wq_config(
    memory=memory,
    max_threads=max_threads,
    # Expect all jobs to be the same size so don't need work_queue:
    use_work_queue=False,
    monitor=False,
    run_dir=run_dir
)

# Read in the visit list
opsim_db_file = args.opsim_db_file

assert os.path.isfile(opsim_db_file)
with sqlite3.connect(config_path(opsim_db_file)) as con:
    df = pd.read_sql("select * from observations where "
                     f"target_name='{calib_type}'", con)

imsim_yaml = config_path(calib_yaml[calib_type])

visits = df['observationId'].to_numpy()
print("number of visits:", len(visits))

# target_dets specifies the detectors to simulate for each visit.
if args.target_dets is None:
    target_det_range = range(189)
else:
    target_det_range = eval(args.target_dets)
target_dets = {_: set(target_det_range) for _ in visits}

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
