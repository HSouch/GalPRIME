""" Binning methods

This module contains all necessary methods for binning catalogues according to various parameters.

"""

from astropy.table import Table
import numpy as np
from matplotlib import pyplot as plt

from ..utils import object_kde


class BinList:
    """ Class to hold a list of bins and rebin efficiently.
    """
    def __init__(self, bins, keys={"mag": "i", "r50": "R_GIM2D", "n": "SERSIC_N_GIM2D", "ellip": "ELL_GIM2D"}, 
                 binning_info = {}):

        if isinstance(bins, Bin):
            bins = [bins]
        self.bins = bins
        self.binning_info = binning_info
        
    def rebin(self, key, lims):
        """ Rebins all bins in the current binlist by the key and specified limits.
        """
        new_bins = []
        for b in self.bins:
            new_bins.extend(b.rebin(key, lims))
        
        self.binning_info[key] = lims
        self.bins = new_bins
    
    def prune_bins(self, min_objects=10, logger=None):
        """ Remove bins with fewer than min_objects objects."""
        current_bin_count = len(self.bins)

        self.bins = [b for b in self.bins if len(b.objects) >= min_objects]

        if logger is not None:
            logger.info(f"Pruned {current_bin_count - len(self.bins)} bins with fewer than {min_objects} objects.")

    def __repr__(self):
        return f'BinList with {len(self.bins)} bins.'


class Bin:
    def __init__(self, objects, params={}, bin_info={}, bin_indices=[]):
        self.objects = objects
        self.params = params
        self.bin_info = bin_info
        self.bin_indices = bin_indices

    def rebin(self, key, lims):
        """ Generate a list of bins by splitting the current bin along the key column by the lims array."""
        new_bins = []
        
        for i in range(len(lims)-1):
            mask = (self.objects[key] > lims[i]) & (self.objects[key] < lims[i+1])
            new_bin = Bin(self.objects[mask], self.params, self.bin_info.copy(), self.bin_indices.copy())
            new_bin.bin_info[key] = (lims[i], lims[i+1])
            new_bin.bin_indices.append(i)
            new_bins.append(new_bin)
        return new_bins
    
    def return_columns(self, structural_params=False):
        values = []
        params = ["MAG", "R50", "N", "ELLIP"] if structural_params else self.params
        for key in params:
            values.append(self.objects[self.params[key]].data)
        return np.array(values)
    
    def to_kde(self):
        return object_kde(self.return_columns(structural_params=True))

    def bin_id(self):
        return "_".join([f"_{i}" for i in self.bin_indices])

    def gen_all_kdes(self, keys=None):
        if keys is None:
            keys = self.params
        return {key: object_kde(self.objects[self.params[key]].data) for key in keys}

    def __repr__(self) -> str:
        return f'Bin with {len(self.objects)} objects.'



def bin_catalogue(table, bin_params = {}, 
                  params={"mag": "i", "r50": "R_GIM2D", "n": "SERSIC_N_GIM2D", "ellip": "ELL_GIM2D"},
                  min_objects=10, logger=None):
    """ Bin a table along a set of parameters.

    Args:
        table (astropy.table.Table): The table to bin into a Binlist.
        bin_params (dict): The binning parameters. Needs to be a dict of the form {key: [lims]}.
        params (dict): The structural parameters

    Returns:
        _type_: _description_

    Usage:
        >>> bin_params = {"Z_BEST": [0.1, 0.3, 0.5, 0.7, 0.9], "MASS_MED": [10, 10.5, 11, 11.5], 
            "sfProb": [0, 0.5, 1.]}
        >>> binned = bin_catalogue(table, bin_params=bin_params)
    """
    binlist = BinList(Bin(table, params=params, bin_info=bin_params))
    
    for key in bin_params:
        binlist.rebin(key, bin_params[key])

    binlist.prune_bins(min_objects=min_objects, logger=logger)
    return binlist


def trim_table(table, config=None, mag_key="i", r50_key="R50", n_key="n", ellip_key="ellip", bin_keys = []):
    if config is None:
        keys = [mag_key, r50_key, n_key, ellip_key] + bin_keys
    else:
        keys = [config["KEYS"][key] for key in config["KEYS"].keys()] + config["BINS"].keys()

    t_trimmed = table[keys]
    return t_trimmed