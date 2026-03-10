set -xe

REPO=./repo

butler certify-calibrations ${REPO} \
       u/jchiang8/bias_w_cti_90..98_w_2026_09/20260306T162952Z \
       u/jchiang/calib/bias.20260308 \
       --begin-date 2026-03-08 --end-date 2050-01-01 bias

butler certify-calibrations ${REPO} \
       u/jchiang8/dark_w_cti_90..98_w_2026_09/20260306T163543Z \
       u/jchiang/calib/dark.20260308 \
       --begin-date 2026-03-08 --end-date 2050-01-01 dark

butler certify-calibrations ${REPO} \
       u/jchiang8/flat_w_cti_90..98_w_2026_09/20260309T190841Z \
       u/jchiang/calib/flat.20260308 \
       --begin-date 2026-03-08 --end-date 2050-01-01 flat

butler certify-calibrations ${REPO} \
       u/jchiang8/ptc_w_cti_90..98_w_2026_09/20260306T160838Z \
       u/jchiang/calib/ptc.20260308 \
       --begin-date 2026-03-08 --end-date 2050-01-01 ptc

butler certify-calibrations ${REPO} \
       u/jchiang8/cti_90..98_w_2026_09/20260306T155056Z \
       u/jchiang/calib/cti.20260308 \
       --begin-date 2026-03-08 --end-date 2050-01-01 cti

butler certify-calibrations ${REPO} \
       u/jchiang8/linearizer_90..98_w_2026_09/20260304T235407Z \
       u/jchiang/calib/linearizer.20260308 \
       --begin-date 2026-03-08 --end-date 2050-01-01 linearizer

butler certify-calibrations ${REPO} \
       u/jchiang8/defects_90..98_w_2026_09/20260304T211003Z \
       u/jchiang/calib/defects.20260308 \
       --begin-date 2026-03-08 --end-date 2050-01-01 defects

#butler certify-calibrations ${REPO} \
#       u/jchiang8/bf_distortion_matrix_w_cti_90..98_w_2026_09/20260306T184348Z \
#       u/jchiang/calib/electroBFDistortionMatrix.20260308 \
#       --begin-date 2026-03-08 --end-date 2050-01-01 electroBfDistortionMatrix

butler collection-chain ${REPO} u/jchiang8/calib.2026-03-08 \
       u/jchiang/calib/bias.20260308 \
       u/jchiang/calib/dark.20260308 \
       u/jchiang/calib/flat.20260308 \
       u/jchiang/calib/ptc.20260308 \
       u/jchiang/calib/cti.20260308 \
       u/jchiang/calib/linearizer.20260308 \
       u/jchiang/calib/defects.20260308
