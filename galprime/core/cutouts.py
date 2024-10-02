from astropy.io import fits
import numpy as np

from .. import plotting 

import copy

__all__ = ['Cutouts']

class Cutouts:
    def __init__(self, cutouts=[], cutout_data=[], metadata={}, min_index=0):
        self.cutouts = cutouts
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
    

    def show_cutouts(self, ncols=10, **kwargs):
        nrows = len(self.cutouts) // ncols
        plotting.show_cutouts(self, ncols=ncols, nrows=nrows, **kwargs)
    
    def copy(self):
        return copy.deepcopy(self)
