from astropy.io import fits
import numpy as np

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
    
    @staticmethod
    def from_file(filename, verbose=False, min_index=0):
        cutouts, cutout_data, metadata = [], [], {}
        with fits.open(filename) as hdul:
            for i in range(min_index, len(hdul)):
                try:
                    data, header = hdul[i].data, hdul[i].header
                    cutouts.append(data)
                    cutout_data.append(header)
                except Exception as e:
                    if verbose:
                        print(f'Failed to recover image at index {i}: {e}')
                    continue
        
        metadata["FILENAME"] = filename
        metadata["N_CUTOUTS"] = len(cutouts)
        metadata["SHAPE"] = cutouts[0].shape

        if verbose:
            print(f'Loaded {metadata["N_CUTOUTS"]} cutouts from {metadata["FILENAME"]} with shape {metadata["SHAPE"]}')
        
        return Cutouts(cutouts=cutouts, cutout_data=cutout_data, metadata=metadata, min_index=min_index)
    
    
    
