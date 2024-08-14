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

import warnings



class GPrime:

    def __init__(self, config_filename, verbose=True, **kwargs):
        self.config = c = config.read_config_file(config_filename)

        self.run_id = kwargs.get("run_id", np.random.randint(1e3, 1e4))
        
        self.outfiles = utils.gen_filestructure(c["DIRS"]["OUTDIR"])
        

        self.log_level = kwargs.get("log_level", 20)
        self.logger = utils.setup_logging(self.run_id, self.log_level, 
                                          log_filename=f'{c["DIRS"]["OUTDIR"]}output_{self.run_id}.log')
        
        self.logger.info(f"Starting run ID:{self.run_id}, GalPRIME Version: {gp.__version__}", )
        print(f"Starting run ID:{self.run_id}")

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

        binlist = binning.bin_catalogue(self.table, bin_params=c["BINS"], params=c["KEYS"], logger=self.logger)
        max_bins = len(binlist.bins) if max_bins is None else min(max_bins, len(binlist.bins))

        for i in range(max_bins):

            b = binlist.bins[i]
            containers = self.process_bin(b, bin_id=i)

        return containers


    def process_bin(self, bin, bin_id=0, method=None):
        t_start = time.time()
        c = self.config
        kde = bin.to_kde()

        self.logger.info(f'Processing bin {bin_id}: {bin} on {c["NTHREADS"]} cores.')

        TIMEOUT = c["TIME_LIMIT"]
        results = []

        with ProcessPoolExecutor(max_workers=c["NTHREADS"]) as executor:
            futures = [executor.submit(run_single_sersic, self, kde) for i in range(c["MODEL"]["N_MODELS"])]
            for future in as_completed(futures):
                try:
                    data = future.result(timeout=TIMEOUT)
                    results.append(data)
                except TimeoutError:
                    self.logger.warn(f"Timeout reached")
                except Exception as e:
                    self.logger.warn(f"Error in container: {e}")

        t_finish = time.time()
        self.logger.info(f"Finished processing bin {bin_id} in {(t_finish - t_start) / 60:.3f} minutes")
        self.logger.info(f'Time per object: {(t_finish - t_start) / c["MODEL"]["N_MODELS"]:.3f} seconds')
        self.logger.info(f"{len(results)} successful extractions out of {c['MODEL']['N_MODELS']}")

        print(results)

        return None
    

    def run_container(self, container):
        container.run()
        return container




class GPrimeContainer:
    def __init__(self, config, model, psf=None, bg=None, metadata={}, id=0, stop_code=0):
        self.id = id
        
        self.config = config
        self.model = model
        self.convolved_model = self.bgadded = self.bg_model = None
        self.bgadded_masked, self.bgsub_masked = None, None

        self.logger = logging.getLogger(f"CONT_{self.id}")

        self.model_profile = self.bgadded_profile = self.bgsub_profile = None
        self.psf = psf
        self.bg = bg
        self.metadata = metadata
        self.stop_code = stop_code
    

    def run(self):
        
        # Convolve model with PSF, add background
        self.convolved_model = convolve2d(self.model, self.psf, mode='same')
        self.bgadded = self.convolved_model + self.bg

        # Generate source mask, and background mask/model
        self.source_mask, mask_metadata = masking.gen_mask(self.bgadded, self.config)
        self.metadata.update(mask_metadata)
        self.bg_source_mask, self.bg_model = gp.estimate_background_2D(self.bgadded, self.config)

        # Extract all necessary profiles
        self.model_profile = gp.isophote_fitting(self.model, self.config)

        self.bgmasked = np.copy(self.bgadded)
        self.bgmasked[self.bg_source_mask] = self.bg_model.background[self.bg_source_mask]

        # self.bgmasked = np.ma.masked_array(self.bgadded, mask=self.bg_source_mask)
        self.bgadded_profile = gp.isophote_fitting(self.bgmasked, self.config)

        self.logger.debug(len(self.model_profile), len(self.bgadded_profile))

        self.stop_code=1


    def plot(self, outdir="", **kwargs):
        
        fig, ax = plt.subplots(2, 4, figsize=(10, 5), facecolor="white")

        lims = ZScaleInterval().get_limits(self.bg)

        ax[0][0].imshow(self.model, cmap=kwargs.get("cmap", "Greys"), vmin=lims[0], vmax=lims[1])
        ax[0][1].imshow(self.convolved_model, cmap=kwargs.get("cmap", "Greys"), vmin=lims[0], vmax=lims[1])
        ax[0][2].imshow(self.bg, cmap=kwargs.get("cmap", "Greys"), vmin=lims[0], vmax=lims[1])
        ax[0][3].imshow(self.bgadded, cmap=kwargs.get("cmap", "Greys"), vmin=lims[0], vmax=lims[1])
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


def run_single_sersic(gp_obj: GPrime, kde, plot=False):

    warnings.filterwarnings("ignore")

    psf_index = np.random.randint(len(gp_obj.psfs.cutouts))
    bg_index = np.random.randint(len(gp_obj.bgs.cutouts))

    psf, bg = np.copy(gp_obj.psfs.cutouts[psf_index]), np.copy(gp_obj.bgs.cutouts[bg_index])

    model, model_data = gp.gen_single_sersic(gp_obj.config, kde, gp_obj.mag_kde)
    convolved_model = convolve2d(model, psf, mode='same')
    bgadded = convolved_model + bg

    source_mask, mask_metadata = gp.masking.gen_mask(model, gp_obj.config)
    bg_source_mask, bg_model = gp.estimate_background_2D(bgadded, gp_obj.config)

    model_profile = gp.extraction.isophote_fitting(model, gp_obj.config)

    bgmasked = np.copy(bgadded)
    bgmasked[source_mask] = bg_model.background[source_mask]

    bgadded_profile = gp.extraction.isophote_fitting(bgmasked, gp_obj.config)

    bgsub = bgadded - bg_model.background
    bgsub[source_mask] = 0
    bgsub_profile = gp.extraction.isophote_fitting(bgsub, gp_obj.config)

    return ProfileContainer(model_profile=model_profile, 
                            bgadded_profile=bgadded_profile, 
                            bgsub_profile=bgsub_profile, 
                            metadata=mask_metadata)
    
class ProfileContainer:
    def __init__(self, model_profile=None, bgadded_profile=None, bgsub_profile=None, metadata={}):
        self.model_profile = model_profile
        self.bgadded_profile = bgadded_profile
        self.bgsub_profile = bgsub_profile
        self.metadata = metadata

