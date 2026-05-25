set -xe

butler register-dataset-type ./repo euclid_stars_20251027 SimpleCatalog htm7
(cd /pscratch/sd/j/jchiang8/Euclid_Joint_DDP/refcat_creation && butler ingest-files -t direct /pscratch/sd/j/jchiang8/Euclid_Joint_DDP/euclid_catalogs/work/repo euclid_stars_20251027 refcats/euclid_stars_20251017 filename_to_htm.ecsv)
butler collection-chain ./repo --mode extend refcats refcats/euclid_stars_20251017
