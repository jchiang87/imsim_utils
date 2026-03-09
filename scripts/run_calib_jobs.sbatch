#!/bin/bash

#SBATCH --account=m1727
#SBATCH --time=10:00:00
#SBATCH --qos=regular
##SBATCH --time=30:00
##SBATCH --qos=debug
##SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=120
#SBATCH --module=cvmfs

cd /pscratch/sd/j/jchiang8/Euclid_Joint_DDP/calibration
source ./setup_cvmfs.sh

echo "*** BIAS JOBS ***"
run_calib_jobs.py BIAS --opsim_db_file calib_frames.db

echo "*** DARK JOBS ***"
run_calib_jobs.py DARK --opsim_db_file calib_frames.db

echo "*** FLAT JOBS ***"
run_calib_jobs.py FLAT --opsim_db_file calib_frames.db

echo "*** PTC JOBS ***"
run_calib_jobs.py PTC --opsim_db_file calib_frames.db
