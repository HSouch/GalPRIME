
from ..core import bootstrap_median, get_run_id
from ..utils import to_sb, load_object, gen_filestructure

import os

import numpy as np

from astropy.io import fits
from astropy.table import Table

from matplotlib import pyplot as plt

class ProfileSet:
    """
    A class to contain and handle sets of surface brightness profiles.
    Parameters
    ----------
    profiles : list, optional
        List of profile tables (e.g., astropy Table objects), by default [].
    headers : list, optional
        List of FITS headers corresponding to each profile, by default [].
    name : str, optional
        Name of the profile set, by default "".
    metadata : dict, optional
        Additional metadata associated with the profile set, by default {}.
    convert_to_sb : bool, optional
        If True, convert intensity profiles to surface brightness upon initialization, by default False.
    Attributes
    ----------
    name : str
        Name of the profile set.
    profiles : list
        List of profile tables.
    headers : list
        List of FITS headers for each profile.
    metadata : dict
        Metadata associated with the profile set.
    Methods
    -------
    gen_medians(samples=1000, dtype='table')
        Generate bootstrapped medians for the profiles.
    to_median_set(samples=1000, dtype='table')
        Create a MedianSet object from the bootstrapped medians.
    to_sb(arcconv=0.168, zpm=27)
        Convert intensity profiles to surface brightness.
    plot_profiles(ax=None, xcol="sma", ycol="intens", **kwargs)
        Plot all profiles in the set.
    rebin(key, values=[])
        Rebin profiles based on a specified key and bin edges.
    gen_hist(key, ax=None, **kwargs)
        Generate a histogram of a specified header key across all profiles.
    from_fits(fitsfile, name="", start_index=1, max_index=None, convert_to_sb=False)
        Create a ProfileSet from a FITS file.
    """
    
    def __init__(self, profiles=[], headers=[], name="", metadata={}, convert_to_sb=False):
        
        self.name = name
        self.profiles = profiles
        self.headers = headers
        self.metadata = metadata

        if convert_to_sb:
            self.to_sb()
    
    def gen_medians(self, samples=1000, dtype='table'):
        return bootstrap_median(self.profiles, samples=samples, dtype=dtype)
    

    def to_median_set(self, samples=1000, dtype='table'):
        med_sma, median, low, high = self.gen_medians(samples=samples, dtype=dtype)
        return MedianSet(med_sma=med_sma, median=median, low=low, high=high, name=self.name)
    

    def to_sb(self, arcconv=0.168, zpm=27):
        for profile in self.profiles:
            profile["intens"] = to_sb(profile["intens"], m_0=zpm, arcconv=arcconv)


    def plot_profiles(self, ax=None, xcol="sma", ycol="intens", **kwargs):
        if ycol == "sb" and "sb" not in self.profiles[0].colnames:
            self.to_sb()
            
        if ax is None:
            fig, ax = plt.subplots()
        for profile in self.profiles:
            ax.plot(profile[xcol], profile[ycol], **kwargs)
        return ax
    

    def rebin(self, key, values=[]):
        values = np.array(values)
        psets = []

        for i in range(len(values) - 1):
            vmin, vmax = values[i], values[i + 1]
            print(vmin, vmax)
    

    def gen_hist(self, key, ax=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots()
        
        vals = [head[key] for head in self.headers]
        ax.hist(vals, **kwargs)
    

    @staticmethod
    def from_fits(fitsfile, name="", start_index=1, max_index=None, convert_to_sb=False):

        with fits.open(fitsfile) as hdul:
            max_index = max_index if max_index is not None else len(hdul)
            profiles = []
            headers = []
            for i in range(start_index, max_index):
                profiles.append(Table(hdul[i].data))
                headers.append(hdul[i].header)
        
        return ProfileSet(profiles=profiles, headers=headers, name=name, convert_to_sb=convert_to_sb)
    

class MedianSet:
    """
    Represents a set of median profile data, including the median, lower, and upper bounds.
    These are based on the galprime routine boostrap_median
    Parameters
    ----------
    med_sma : array-like, optional
        The semi-major axis values corresponding to the profile measurements.
    median : array-like, optional
        The median values of the profile.
    low : array-like, optional
        The lower bound values (e.g., 16th percentile) of the profile.
    high : array-like, optional
        The upper bound values (e.g., 84th percentile) of the profile.
    name : str, optional
        An optional name for the median set.
    Methods
    -------
    from_fits(fitsfile, name="")
        Create a MedianSet instance from a FITS file using a ProfileSet.
    """

    def __init__(self, med_sma=None, median=None, low=None, high=None, name=""):
        self.name = name
        self.med_sma = med_sma
        self.median = median
        self.low = low
        self.high = high
    
    @staticmethod
    def from_fits(fitsfile, name=""):
        pset = ProfileSet.from_fits(fitsfile, )
        med_sma, median, low, high = pset.gen_medians(samples=1000, dtype='table')
        return MedianSet(med_sma=med_sma, median=median, low=low, high=high, name=name)


def retrieve_profile_sets(filedir, run_id=5271237, progress_bar=True, **kwargs):
    output_files = gen_filestructure(filedir, generate=False)
    if run_id is None:
        run_id = get_run_id(filedir)
    config = load_object(f'{output_files["ADDL_DATA"]}/config_{run_id}.pkl')

    mod_dir = output_files["MODEL_PROFS"]
    coadd_dir = output_files["COADD_PROFS"]
    bgsub_dir = output_files["BGSUB_PROFS"]

    qual_dir = f'{filedir}/quality_tables/'
    if not os.path.exists(qual_dir):
        os.makedirs(qual_dir)

    iter = os.listdir(coadd_dir)[:]
    if progress_bar:
        from tqdm import tqdm
        iter = tqdm(iter, desc="Processing profiles", unit="file")
    
    model_sets, coadd_sets, bgsub_sets = {}, {}, {}

    in_sb = kwargs.get('in_sb', False)

    for f in iter:
        model_sets[f] = ProfileSet.from_fits(os.path.join(mod_dir, f), name=f, convert_to_sb=in_sb,
                                              max_index=kwargs.get('max_index', None))
        coadd_sets[f] = ProfileSet.from_fits(os.path.join(coadd_dir, f), name=f, convert_to_sb=in_sb,
                                             max_index=kwargs.get('max_index', None))
        bgsub_sets[f] = ProfileSet.from_fits(os.path.join(bgsub_dir, f), name=f, convert_to_sb=in_sb,
                                             max_index=kwargs.get('max_index', None))

    return config, model_sets, coadd_sets, bgsub_sets
