#!/bin/bash
set -xe

REPO=./repo

butler create ${REPO}
butler register-instrument ${REPO} lsst.obs.lsst.LsstCamSim
butler write-curated-calibrations --collection LSSTCamSim/calib ${REPO} 'LSSTCamSim'

# Make a copy of the registry before ingesting any raw data.  We can
# then reset the repo to its pristine state by doing
#
# $ cp ${REPO}/gen3.sqlite3_orig ${REPO}/gen3.sqlite3
# $ rm -rf ${REPO}/u
#
cp ${REPO}/gen3.sqlite3 ${REPO}/gen3.sqlite3_orig
