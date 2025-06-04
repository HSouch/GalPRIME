from scipy.stats import gaussian_kde

from astropy.io import fits

import pickle
import os

import numpy as np

import datetime

from scipy.interpolate import interp1d


def object_kde(columns):
    """ Generate a gaussian kernel density estimate for the columns of a bin. """
    return gaussian_kde(columns)


def gen_filestructure(outdir, generate=True):
    os.makedirs(outdir, exist_ok=True)

    file_dict = {"MODEL_PROFS": f"{outdir}model_profiles/",
                 "COADD_PROFS": f"{outdir}coadd_profiles/",
                 "BGSUB_PROFS": f"{outdir}bgsub_profiles/",
                 "MEDIANS": f"{outdir}medians/",
                 "ADDL_DATA": f"{outdir}additional_data/",
                 "MODEL_PARAMS": f"{outdir}additional_data/model_params/",
                 "PLOTS": f"{outdir}plots/",
                 "TEMP": f"{outdir}tempfiles/"}
    if generate:
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


    """
    Create a FITS header from a configuration dictionary.
    This function takes a configuration dictionary, flattens it, and creates
    a FITS header using the keys and values from the flattened dictionary.
        :param config: The configuration dictionary to be converted into a FITS header.
        :type config: dict
        :return: A FITS header object created from the configuration dictionary.
        :rtype: astropy.io.fits.Header
    """

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
    """
    Generate an integer label based on the current date and time.
    The label is constructed using the current month, day, hour, and minute,
    formatted as MMDDHHMM and scaled to fit into an integer.
    Returns:
        int: An integer representing the current date and time in the format MMDDHHMM.
    """
    
    dt = datetime.datetime.now()
    return int(1e6 * dt.month + 1e4 * dt.day + 100 * dt.hour + dt.minute)


def get_arcconv(wcs):
    """
    Calculate the arcsecond per pixel conversion factor from a WCS object.
    Parameters:
    wcs (astropy.wcs.WCS): A WCS (World Coordinate System) object containing the pixel scale matrix.
    Returns:
    float: The arcsecond per pixel conversion factor, rounded to 6 decimal places.
    """

    dx = wcs.pixel_scale_matrix[0,0]    
    return np.round(abs(dx * 3600), 6)


def get_output_id(outdir):
    """ Get the unique run IDs of the output files in the specified directory
        It does so by extracting the IDs from any saved config files in the additional_data directory.
    Args:
        outdir (str): The output directory.

    Returns:
        array: An array of the unique run IDs. 
    """
    ids = []
    for fn in os.listdir(f'{outdir}/additional_data/'):
        ids.append(int(fn.split('_')[1].split('.')[0]))
    
    return np.unique(ids)


def fill_nans(x, y):
    """
    Fill NaN values in the y array using linear interpolation based on the x array.
    
    Parameters:
    x (array-like): The x values corresponding to the y values.
    y (array-like): The y values, which may contain NaN values.
    
    Returns:
    array: The y array with NaN values filled using linear interpolation.
    """
    
    mask = np.isnan(y)
    if np.all(mask):
        return y  # If all values are NaN, return the original array
    
    interp_func = interp1d(x[~mask], y[~mask], bounds_error=False, fill_value="extrapolate")
    return interp_func(x)
