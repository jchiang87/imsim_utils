set -xe

butler register-skymap --config-file $OBS_LSST_DIR/config/makeSkyMap.py --config name="lsst_cells_v2" ./repo
