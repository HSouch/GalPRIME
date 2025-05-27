
import os

from astropy.io import fits
import numpy as np

from scipy.signal import convolve2d

from .. import plotting 

import copy

__all__ = ['Cutouts']

class Cutouts:
    def __init__(self, cutouts=None, cutout_data=[], metadata={}, min_index=0):
        self.cutouts = [] if cutouts is None else cutouts
        self.cutout_data = cutout_data
        self.metadata = metadata
        self.min_index = min_index

        self.ras, self.decs = [], []

    def get_ra_dec(self, ra_key="RA", dec_key="DEC"):

        for i, dataset in enumerate(self.cutout_data):
            try:
                self.ras.append(dataset[ra_key])
                self.decs.append(dataset[dec_key])
            except Exception as e:
                print(f"Failed to get header info at index {i}: {e}")
                continue

    def sample(self):
        index = np.random.randint(0, len(self.cutouts))
        return self.cutouts[index], self.cutout_data[index]
    

    def combine(self, to_add, method="random"):
        """
        Combines the cutouts of two instances of the `Cutouts` class.
        Parameters:
            to_add (Cutouts): The `Cutouts` instance to be combined with.
            method (str, optional): The method used for combining the cutouts. 
                Defaults to "random". Possible values are "random" and "direct".
        Returns:
            Cutouts: A new `Cutouts` instance with the combined cutouts.
        Raises:
            ValueError: If an invalid method is provided.

        """
        out_cutouts = self.copy()
        out_cutouts.cutouts = []

        if method == "random":
            for cutout in self.cutouts:
                out_cutouts.cutouts.append(cutout + to_add.cutouts[np.random.randint(0, len(to_add.cutouts))])
        elif method == "direct":
            for i in range(len(self.cutouts)):
                out_cutouts.cutouts.append(self.cutouts[i] + to_add.cutouts[i])
        else:
            raise ValueError("Invalid method provided. Possible values are 'random' and 'direct'.")

        return out_cutouts
    
    def convolve(self, psf):

        psfs = psf if isinstance(psf, list) else [psf for _ in range(len(self.cutouts))]
        out_cutouts = self.copy()
        out_cutouts.cutouts = [convolve2d(cutout, psf) for cutout, psf in zip(self.cutouts, psfs)]

        return out_cutouts
    

    @staticmethod
    def from_file(filename, logger=None, min_index=0):
        cutouts, cutout_data, metadata = [], [], {}
        with fits.open(filename) as hdul:
            for i in range(min_index, len(hdul)):
                try:
                    data, header = hdul[i].data, hdul[i].header
                    cutouts.append(data)
                    cutout_data.append(header)
                except Exception as e:
                    if logger is not None:
                        logger.warn(f'Failed to recover image at index {i}: {e}')
                    continue
        
        metadata["FILENAME"] = filename
        metadata["N_CUTOUTS"] = len(cutouts)
        metadata["SHAPE"] = cutouts[0].shape

        if logger is not None:
            logger.info(f'Loaded {metadata["N_CUTOUTS"]} cutouts from {metadata["FILENAME"]} with shape {metadata["SHAPE"]}')
        
        return Cutouts(cutouts=cutouts, cutout_data=cutout_data, metadata=metadata, min_index=min_index)
    

    @staticmethod
    def from_directory(directory, logger=None, template=None, cutout_index=0):
        def verify_fits(filename):
            return any(filename.endswith(ext) for ext in (".fits", ".fit", ".fts"))

        suffixes = [".fits", ".fits.gz"]
        
        filenames = [os.path.join(directory, f) for f in os.listdir(directory) if verify_fits(f)]
        
        cutouts = Cutouts()

        for fn in filenames:
            try:
                hdul = fits.open(fn)
                cutout = hdul[cutout_index].data
                header = hdul[cutout_index].header

                # Add the primary header to the cutout header if values not already present
                if cutout_index != 0:
                    primary_header = hdul[0].header
                    for key in primary_header.keys():
                        if key not in header:
                            header[key] = primary_header[key]

                cutouts.cutouts.append(cutout)
                cutouts.cutout_data.append(header)

            except Exception as e:
                if logger:
                    logger.error(f"Error reading {fn}: {e}")
                continue
        
        return cutouts

    
    def save(self, filename="cutouts.fits", overwrite=False):
        if not filename.endswith(".fits"):
            filename += ".fits"
        
        hdul = fits.HDUList()
        for cutout, header in zip(self.cutouts, self.cutout_data):
            hdu = fits.ImageHDU(data=cutout, header=fits.Header(header))
            
            hdul.append(hdu)

        hdul = fits.HDUList(hdul)
        hdul.writeto(filename, overwrite=overwrite)


    def show_cutouts(self, ncols=10, **kwargs):
        nrows = len(self.cutouts) // ncols
        plotting.show_cutouts(self, ncols=ncols, nrows=nrows, **kwargs)
    
    def copy(self):
        return copy.deepcopy(self)
    