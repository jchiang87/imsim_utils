#!/bin/bash
set -xe

raw_dir=/sdf/data/rubin/shared/ops-rehearsals/ops-rehearsal-4/image_sims/full_focal_plan_calibs/03*
repo=./repo

butler --log-level DEBUG ingest-raws --transfer direct ${repo} ${raw_dir}/*R22_S11*.fits.fz
#butler define-visits ${repo} LSSTCamSim
