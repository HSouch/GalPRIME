from astropy.io import fits

__all__ = ['Cutouts']

class Cutouts:
    def __init__(self, cutouts=[], cutout_data=[], metadata={}):
        self.cutouts = cutouts
        self.cutout_data = cutout_data
        self.metadata = metadata

    @staticmethod
    def from_file(filename, verbose=False):
        cutouts, cutout_data, metadata = [], [], {}
        with fits.open(filename) as hdul:
            for i in range(0, len(hdul)):
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
            print(f'Loaded {metadata["N_CUTOUTS"]} cutouts from {metadata["FILENAME"]}')
        
        return Cutouts(cutouts=cutouts, cutout_data=cutout_data, metadata=metadata)

