#!/bin/bash

#SBATCH --partition=milano
#SBATCH --account=rubin:commissioning
#SBATCH --time=10:00:00

##SBATCH --nodes=1
##SBATCH --exclusive

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=54
##SBATCH --mem=150GB

cd /sdf/data/rubin/user/jchiang/ops_rehearsal_3/vignetted_flats
source ./setup.sh

echo "*** FLAT JOBS ***"
python run_imsim_jobs.py FLAT
