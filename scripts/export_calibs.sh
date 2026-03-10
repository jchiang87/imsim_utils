set -xe

REPO=./repo
output_dir=./calib.2026-03-08

butler export-calibs -t direct ${REPO} ${output_dir} u/jchiang8/calib.2026-03-08
