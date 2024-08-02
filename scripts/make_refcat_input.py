import os
import glob
import sys
import pickle
from collections import defaultdict
import multiprocessing
import warnings
import numpy as np
import pandas as pd
from refcat_generation import MagErrors, RADecErrors
from skycatalogs import skyCatalogs, load_lsst_bandpasses
from skycatalogs.objects.base_object import ObjectList


lsst_bps = load_lsst_bandpasses()


def compute_mags(imin, imax, lsst_bps=lsst_bps, outfile=None):
    global objects
    data = defaultdict(list)
    for i, obj in enumerate(objects[imin:imax]):
        data['object_id'].append(obj.id)
        data['ra'].append(obj.ra)
        data['dec'].append(obj.dec)
        #data['is_variable'].append(obj.get_native_attribute('is_variable'))
        data['is_variable'].append(False)
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')
            sed = obj.get_total_observer_sed()
            for band in "gri":
                bp = lsst_bps[band]
                data[band].append(sed.calculateMagnitude(bp))
    df = pd.DataFrame(data)
    if outfile is None:
        outfile = f"refcat_mags_{imin:08d}_{imax:08d}.parquet"
    df.to_parquet(outfile)


skycatalog_yaml = ('/sdf/data/rubin/shared/ops-rehearsals/'
                   'ops-rehearsal-4/imSim_catalogs/'
                   'skyCatalogs/aos_objects/skyCatalog.yaml')

def read_skyCatalog_hps(hps):
    mjd = 60390.3293054128
    obj_type = 'star'

    skycatalog_root = os.path.dirname(skycatalog_yaml)

    sky_cat = skyCatalogs.open_catalog(skycatalog_yaml,
                                       skycatalog_root=skycatalog_root)
    object_list = ObjectList()
    for hp in hps:
        object_list.append_object_list(
            sky_cat.get_object_type_by_hp(hp, obj_type, mjd=mjd))
    return object_list


def read_skyCatalogs(ra, dec, radius=0.5):
    region = skyCatalogs.Disk(ra, dec, radius*3600.0)
    mjd = 60390.3293054128
    obj_types = {'star'}

    skycatalog_root = os.path.dirname(skycatalog_yaml)

    sky_cat = skyCatalogs.open_catalog(skycatalog_yaml,
                                       skycatalog_root=skycatalog_root)
    obj_list = sky_cat.get_objects_by_region(region, mjd=mjd,
                                             obj_type_set=obj_types)
    return obj_list


fields = {
#    "DEEP_A0": (216.00, -12.50),
#    "DEEP_B0": (310.00, -19.00),
#    "ELAIS_S1": (9.45, -44.00),
#    "Rubin_SV_225_-40": (225.00, -40.00),
#    "Rubin_SV_250_2": (250.00, 2.00),
#    "Rubin_SV_280_-48": (280.00, -48.00),
#    "Rubin_SV_300_-41": (300.00, -41.00),
    "aos_blocks": {7682, 7683, 7684, 7685, 7686,
                   7687, 7688, 7689, 7690, 7691,
                   7811, 7812, 7813, 7814, 7815,
                   7816, 7817, 7818, 7819, 7939}
}

radius = 10.

for field, location_info in fields.items():
    initial_refcat_data = f"initial_refcat_mags_{field}_10deg.parquet"
    if not os.path.isfile(initial_refcat_data):
        if isinstance(location_info, set):
            objects = read_skyCatalog_hps(location_info)
        else:
            objects = read_skyCatalogs(*location_info, radius=radius)
        nobj = len(objects)
        print(field, "# skyCatalog objects:", nobj, flush=True)

        processes = 64
        indexes = np.linspace(0, nobj+1, processes+1, dtype=int)

        with multiprocessing.Pool(processes=processes) as pool:
            workers = []
            for imin, imax in zip(indexes[:-1], indexes[1:]):
                args = (imin, imax)
                workers.append(pool.apply_async(compute_mags, args))
            pool.close()
            pool.join()
            _ = [worker.get() for worker in workers]

        files = sorted(glob.glob('refcat_mags*'))
        dfs = [pd.read_parquet(_) for _ in files]
        df = pd.concat(dfs)
        df.to_parquet(initial_refcat_data)
        for item in files:
            os.remove(item)
    else:
        df = pd.read_parquet(initial_refcat_data)

    # Read object ids for variable stars
    with open("../sky_catalog_generation/variable_stars_set.pickle", "rb") as fobj:
        variable_stars = pickle.load(fobj)

    # Simulate uncertainties and add error info.
    # Baseline r-mag cut to control range over which spline fits are applied.
    r_max = 23.
    df0 = pd.DataFrame(df.query(f"r < {r_max}"))

    mag_errors = MagErrors()
    radec_errors = RADecErrors()
    r_mags = df0['r'].to_numpy()
    ra_err, dec_err = radec_errors(r_mags)
    df0['ra'] += np.random.normal(loc=0, scale=ra_err)
    df0['ra_err'] = ra_err*3600*1000.  # Convert to milliarcseconds
    df0['dec'] += np.random.normal(loc=0, scale=dec_err)
    df0['dec_err'] = dec_err*3600*1000.  # Convert to milliarcseconds

    df0['r_err'] = mag_errors(r_mags, 'r')
    df0['r'] += np.random.normal(loc=0, scale=df0['r_err'])

    df0['g_err'] = mag_errors(df0['g'].to_numpy(), 'g')
    df0['g'] += np.random.normal(loc=0, scale=df0['g_err'])

    df0['i_err'] = mag_errors(df0['i'].to_numpy(), 'i')
    df0['i'] += np.random.normal(loc=0, scale=df0['i_err'])

    # Write csv file with refcat entries after final magnitude cut as been
    # applied.
    r_max_final = 21.
    df = df0.query(f"r < {r_max_final}")

    # Set place holder proper motion and parallax values.
    epoch = "1970-01-01T00:00:00"
    pm_ra = 0.0
    pm_dec = 0.0
    pm_ra_err = 0.0
    pm_dec_err = 0.0
    parallax = 0.0
    parallax_err = 0.0

    output_file = f"uw_stars_{field}_refmags_10deg.csv"
    with open(output_file, 'w') as fobj:
        catalog_line = "object_id,ra,dec,ra_err,dec_err"
        for band in "gri":
            catalog_line += f",lsst_{band},lsst_{band}_mag_err"
        catalog_line += (",epoch,pm_ra,pm_dec,pm_ra_err,pm_dec_err,"
                         "parallax,parallax_err")
        catalog_line += ",is_variable"
        fobj.write(catalog_line + "\n")

        for _, row in df.iterrows():
#            is_variable = row.object_id in variable_stars
            is_variable = False
            catalog_line = (f"{row.object_id},{row.ra},{row.dec},"
                            f"{row.ra_err},{row.dec_err}")
            for band in "gri":
                band_err = band + '_err'
                catalog_line += f",{row[band]},{row[band_err]}"
            catalog_line += (f",{epoch},{pm_ra},{pm_dec},{pm_ra_err},"
                             f"{pm_dec_err},{parallax},{parallax_err}")
            catalog_line += f",{is_variable}"
            fobj.write(catalog_line + "\n")
