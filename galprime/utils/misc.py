from scipy.stats import gaussian_kde

from astropy.io import fits

import pickle
import os

import numpy as np

import datetime


def object_kde(columns):
    """ Generate a gaussian kernel density estimate for the columns of a bin. """
    return gaussian_kde(columns)


def gen_filestructure(outdir):
    os.makedirs(outdir, exist_ok=True)

    file_dict = {"MODEL_PROFS": f"{outdir}model_profiles/",
                 "COADD_PROFS": f"{outdir}coadd_profiles/",
                 "BGSUB_PROFS": f"{outdir}bgsub_profiles/",
                 "MEDIANS": f"{outdir}medians/",
                 "ADDL_DATA": f"{outdir}additional_data/",
                 "TEMP": f"{outdir}tempfiles/"}
    
    for key, value in file_dict.items():
        os.makedirs(value, exist_ok=True)

    return file_dict
    

def flatten_dict(d):
    keys, vals = [], []

    for key, val in d.items():
        if isinstance(d[key], dict):
            for k, v in d[key].items():
                keys.append(f'{key}_{k}')
                vals.append(v)
        else:
            keys.append(key)
            vals.append(val)
    return keys, vals


def header_from_config(config):
    keys, vals = flatten_dict(config)
    header = fits.Header()
    for key, val in zip(keys, vals):
        if isinstance(val, np.ndarray):
            val = str(val)
        header[key] = val
    return header


def save_object(obj, filename):
    with open(filename, "wb") as f:
        pickle.dump(obj, f)


def load_object(filename):
    with open(filename, "rb") as f:
        return pickle.load(f)

def get_dt_intlabel():
    dt = datetime.datetime.now()
    return int(1e6 * dt.month + 1e4 * dt.day + 100 * dt.hour + dt.minute)