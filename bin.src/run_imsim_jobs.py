#!/usr/bin/env python
import os
import socket
import argparse
from imsim_utils import GalSimJobGenerator, load_wq_config


def config_path(basename):
    config_abspath = os.path.abspath(args.config_path)
    return os.path.join(config_abspath, basename)


parser = argparse.ArgumentParser()
parser.add_argument("imsim_yaml", type=str, help="imSim config yaml file")
parser.add_argument("visit_list_file", type=str,
                    help="file with list of visits to simulate")
parser.add_argument("--cores_per_node", type=int, default=120,
                    help="Number of cores per node to use")
parser.add_argument("--nproc", type=int, default=4,
                    help="number of processes per galsim instance")
parser.add_argument("--mem_per_core", type=int, default=4000,
                    help="memory per core in MB")
parser.add_argument("--target_dets", default=None,
                    help="Detectors to simulate.")
parser.add_argument("--GB_per_CCD", type=int, default=5,
                    help="Memory in GB to reserve for each CCD process")
args = parser.parse_args()

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
# imSim yaml config file
imsim_yaml = args.imsim_yaml

# Read in the visit list
with open(args.visit_list_file) as fobj:
    visits = [int(_.strip()) for _ in fobj.readlines() if not _.startswith("#")]
print("number of visits:", len(visits))


# target_dets specifies the detectors to simulate for each visit.
if args.target_dets is None:
    target_det_range = range(189)
else:
    target_det_range = eval(args.target_dets)
target_dets = {_: set(target_det_range) for _ in visits}
print(target_dets)

GB_per_CCD = args.GB_per_CCD

load_wq_config(
    memory=memory,
    max_threads=max_threads,
    # Expect all jobs to be the same size so don't need work_queue:
    use_work_queue=False,
    monitor=False,
    run_dir=run_dir
)
generator = GalSimJobGenerator(imsim_yaml, visits,
                               nfiles=9, nproc=nproc,
                               target_dets=target_dets,
                               GB_per_CCD=GB_per_CCD, GB_per_PSF=8,
                               clean_up_atm_psfs=False,
                               bash_app_executor='thread_pool',
                               log_dir='logging',
                               raw_prefix="amp",
                               verbosity=1)
generator.run()
