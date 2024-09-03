import galprime as gp

from . import config, cutouts, binning, masking, modeling, utils

from scipy.signal import convolve2d
from scipy.interpolate import interp1d

from astropy.io import fits
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

        self.binlist = None

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

        self.binlist = binning.bin_catalogue(self.table, bin_params=c["BINS"], params=c["KEYS"], logger=self.logger)
        max_bins = len(self.binlist.bins) if max_bins is None else min(max_bins, len(self.binlist.bins))

        for i in range(max_bins):

            b = self.binlist.bins[i]
            containers = self.process_bin(b, bin_id=i)

        return containers


    def process_bin(self, bin, bin_id=0, method=None):
        t_start = time.time()
        c = self.config
        kde = bin.to_kde()

        self.logger.info(f'Processing {c["MODEL"]["N_MODELS"]} models for bin {bin_id} on {c["NTHREADS"]} cores.')

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

        self.logger.info(f"Saving results to {self.outfiles['MODEL_PROFS']}bin_{bin_id}.fits")
        container_list = ContainerList(c, results)

        return results


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


class ContainerList:

    def __init__(self, config, containers):
        self.config = config
        self.containers = containers

        self.keys = ["intens", "intens_err", "ellipticity", "pa", "x0", "y0"]


    def shared_rs(self, max_sma, step, linear=True):
        if linear:
            rs = np.linspace(1, max_sma, step)
        else:
            rs = []
            r = 1
            while r < max_sma:
                rs.append(r)
                r *= (1 + step)
        return np.array(rs)[:-1]
    

    def process_profile(self, isolist, rs):
        t = isolist.to_table()
        smas = t["sma"]
        arr_set = np.ndarray((len(self.keys), len(rs)))
        for i, key in enumerate(self.keys):
            interp = interp1d(smas, t[key], kind="cubic")
            arr_set[i] = interp(rs)
        return arr_set

    def combine_containers(self, prof="model"):
        """
        Combines the model profiles from each container and returns the combined result as a numpy ndarray.

        Returns:
            combined (ndarray): The combined model profiles.
        """
        rs = self.shared_rs(self.config["MODEL"]["SIZE"] / 2, self.config["EXTRACTION"]["STEP"], 
                       linear=self.config["EXTRACTION"]["LINEAR"])

        for i, container in enumerate(self.containers):
            try:
                if prof == "model":
                    combined[i] = self.process_profile(container.model_profile, rs)
                elif prof == "bgadded":
                    combined[i] = self.process_profile(container.bgadded_profile, rs)
                elif prof == "bgsub":
                    combined[i] = self.process_profile(container.bgsub_profile, rs)
                combined[i] = self.process_profile(container.model_profile, rs)
            except Exception as e:
                combined[i] = np.nan((len(self.keys), len(rs)))
        combined = np.transpose(combined, (1, 0, 2))
        return combined
    

    def gen_astropy_file(self, outname):
        rs = self.shared_rs(self.config["MODEL"]["SIZE"] / 2, self.config["EXTRACTION"]["STEP"], 
                       linear=self.config["EXTRACTION"]["LINEAR"])
        
        combined = self.combine_containers()

        hdulist = fits.HDUList()
        header = fits.Header()
        
        # Add all config parameters to the header
        header = utils.header_from_config(self.config)

        hdulist.append(fits.PrimaryHDU(header=header))
        hdulist.append(fits.ImageHDU(rs))
        combined_header = fits.Header()
        for i, key in enumerate(self.keys):
            combined_header[f"INDEX_{i}"] = key

        hdulist.append(fits.ImageHDU(combined, header=combined_header))

        hdulist.writeto(outname, overwrite=True)
