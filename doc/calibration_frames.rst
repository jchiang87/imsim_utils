Simulating Calibration Frames
-----------------------------

1. Create an opsim db file for the calibration frame "pointings":

   .. code-block:: bash

      $ python scripts/make_calib_opsim_db.py

2. Run galsim jobs in slurm using
   * ``config/[bias,dark,flat,ptc].yaml`` as imsim config templates
   * ``scripts/run_imsim_jobs.py`` to manage the galsim instances using
     parsl
   * ``scripts/run_calib_jobs.sbatch`` to submit the jobs to slurm

Generating Calibration Products
-------------------------------

1. Create a data repository for the raw calibration frames:

   .. code-block:: bash

      $ bash scripts/make_repo.sh

2. Ingest the raw files:

   .. code-block:: bash

      $ bash scripts/ingest_raws.sh

3. Run the calibration pipelines using ctrl_bps_parsl, e.g.,

   .. code-block:: bash

      $ bps submit bps/bps_cpBiasBootstrap.yaml
      $ bps submit bps/bps_cpDarkBootstrap.yaml
      $ bps submit bps/bps_cpFlatBootstrap.yaml
      $ bps submit bps/bps_cpDefects.yaml
      $ bps submit bps/bps_cpLinearizer.yaml
      $ bps submit bps/bps_cpPtc.yaml
      $ bps submit bps/bps_cpBias.yaml
      $ bps submit bps/bps_cpDark.yaml
      $ bps submit bps/bps_cpFlat.yaml

4. Apply certification and validity ranges

   .. code-block:: bash

      $ bulter certify-calibrations ...

5. Export the calibs for transfer to production repos.

   .. code-block:: bash

      $ bulter export-calibs ...
