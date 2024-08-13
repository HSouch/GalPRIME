import galprime as gp

from . import config, cutouts, binning, masking, modeling, utils

from scipy.signal import convolve2d
from astropy.table import Table
from astropy.visualization import ZScaleInterval

import numpy as np

import logging

from matplotlib import pyplot as plt

import time

import multiprocessing as mp

from concurrent.futures import ProcessPoolExecutor, as_completed


class GPrime:

    def __init__(self, config_filename, verbose=True, **kwargs):
        self.config = c = config.read_config_file(config_filename)

        self.run_id = kwargs.get("run_id", np.random.randint(1e3, 1e4))
        
        self.outfiles = utils.gen_filestructure(c["DIRS"]["OUTDIR"])
        

        self.log_level = kwargs.get("log_level", 20)
        self.logger = utils.setup_logging(self.run_id, self.log_level, 
                                          log_filename=f"{self.outfiles['ADDL_DATA']}output.log")
        
        self.logger.info(f"Starting run ID:{self.run_id}, GalPRIME Version: {gp.__version__}")

        # Load in all necessary files (backgrounds, psfs, catalogues, etc)
        self.bgs = cutouts.Cutouts.from_file(f'{c["FILE_DIR"]}{c["FILES"]["BACKGROUNDS"]}', 
                                             logger=self.logger)
        
        self.psfs = cutouts.Cutouts.from_file(f'{c["FILE_DIR"]}{c["FILES"]["PSFS"]}', logger=self.logger)
        self.psfs.get_ra_dec(ra_key=c["PSFS"]["PSF_RA"], dec_key=c["PSFS"]["PSF_DEC"])
        
        self.table = Table.read(f'{c["FILE_DIR"]}{c["FILES"]["CATALOGUE"]}')
        self.table = binning.trim_table(self.table, c)
        
        if c["FILES"]["MAG_CATALOGUE"] is not None:
            self.mags = Table.read(f'{c["FILE_DIR"]}{c["FILES"]["MAG_CATALOGUE"]}')
            self.mag_kde = utils.object_kde(self.mags[c["KEYS"]["MAG"]])
        else:
            self.mags = self.mag_kde = None


    def run(self, max_bins=None, verbose=True):
        c = self.config

        binlist = binning.bin_catalogue(self.table, bin_params=c["BINS"], params=c["KEYS"], verbose=True)
        max_bins = len(binlist.bins) if max_bins is None else min(max_bins, len(binlist.bins))

        for i in range(max_bins):

            b = binlist.bins[i]
            containers = self.process_bin(b, bin_id=i)

        return containers


    def process_bin(self, bin, bin_id=0):
        t_start = time.time()
        c = self.config
        kde = bin.to_kde()

        models = modeling.gen_models(c, kde, mag_kde=self.mag_kde, n_models=c["MODEL"]["N_MODELS"])
        containers = gen_containers(c, models, self.psfs, self.bgs)

        self.logger.info(f'Processing bin {bin_id}: {bin} on {c["NTHREADS"]} cores.')

        TIMEOUT = c["TIME_LIMIT"]
        results = []
        with ProcessPoolExecutor(max_workers=c["NTHREADS"]) as executor:
            futures = [executor.submit(self.run_container, container) for container in containers]
            try:
                for future in as_completed(futures):
                    try:
                        data = future.result(timeout=TIMEOUT)
                        results.append(data)
                        logging.debug(f"  Finished processing container {data.id}")
                    except TimeoutError:
                        self.logger.warn(f"Timeout reached")
                    except Exception as e:
                        self.logger.warn(f"Error in container: {e}")
            except TimeoutError:
                self.logger.warn(f"Timeout reached for bin {bin_id}")

        self.logger.info(f"Finished processing bin {bin_id} in {(time.time() - t_start) / 60:.3f} minutes")

        return containers
    

    def run_container(self, container):
        container.run()
        return container



class GPrimeContainer:
    def __init__(self, config, model, psf=None, bg=None, metadata={}, id=0, stop_code=0):
        self.id = id
        
        self.config = config
        self.model = model
        self.convolved_model = self.bgadded = self.bgsub = None
        self.bgadded_masked, self.bgsub_masked = None, None

        self.logger = logging.getLogger(f"CONT_{self.id}")

        self.model_profile = self.bgadded_profile = self.bgsub_profile = None
        self.psf = psf
        self.bg = bg
        self.metadata = metadata
        self.stop_code = stop_code
    

    def run(self):

        self.convolved_model = convolve2d(self.model, self.psf, mode='same')

        self.bgadded = self.convolved_model + self.bg

        self.source_mask, mask_metadata = masking.gen_mask(self.bgadded, self.config)

        self.stop_code=1


    def plot(self, outdir="", **kwargs):
        
        fig, ax = plt.subplots(2, 4, figsize=(10, 5), facecolor="white")

        lims = ZScaleInterval().get_limits(self.bg)

        ax[0][0].imshow(self.model, cmap=kwargs.get("cmap", "Greys"), vmin=lims[0], vmax=lims[1])
        ax[0][1].imshow(self.bg, cmap=kwargs.get("cmap", "Greys"), vmin=lims[0], vmax=lims[1])
        ax[0][2].imshow(self.bgadded, cmap=kwargs.get("cmap", "Greys"), vmin=lims[0], vmax=lims[1])
        # ax[0][1].imshow(self.bg, cmap=kwargs.get("cmap", "Greys"), vmin=lims[0], vmax=lims[1])

        ax[1][0].imshow(self.source_mask, cmap=kwargs.get("cmap", "Greys"))
        

        for axis in ax.flatten():
            axis.set(xticks=[], yticks=[])

        plt.tight_layout()
        plt.savefig(f"{outdir}_{self.id}.pdf", dpi=kwargs.get("dpi", 150))


def gen_containers(config, models, psfs, bgs):

    containers = []
    for i, model in enumerate(models.cutouts):
        psf_index = np.random.randint(len(psfs.cutouts))
        bg_index = np.random.randint(len(bgs.cutouts))
        containers.append(GPrimeContainer(config, 
                                          np.copy(model), 
                                          psf = np.copy(psfs.cutouts[psf_index]), 
                                          bg = np.copy(bgs.cutouts[bg_index]),
                                          id=i))
    return containers
