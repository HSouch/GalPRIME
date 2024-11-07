import galprime as gp

from . import config, cutouts, binning, masking
from .. import utils

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

import warnings

from joblib import Parallel, delayed



class GalPrimeSingle:
    def __init__(self, config, model, params, bg=None, psf=None):
        self.config = config
        self.model = model
        self.params = params
    
    def process(self):
        pass



class GPrime:
    def __init__(self, config_filename, **kwargs):
        self.config = c = config.read_config_file(config_filename)
        self.binlist = None
        self.outfiles = gp.gen_filestructure(c["DIRS"]["OUTDIR"])

        self.run_id = kwargs.get("run_id", np.random.randint(1e3, 1e4))

        self.model_type = c["MODEL"]["MODEL_TYPE"]

        self.log_level = kwargs.get("log_level", 20)
        self.logger = gp.setup_logging(self.run_id, self.log_level, 
                                          log_filename=f'{c["DIRS"]["OUTDIR"]}output_{self.run_id}.log')
        self.logger.info(f"Starting run ID:{self.run_id}, GalPRIME Version: {gp.__version__}", )


        
        print(f"Starting run ID:{self.run_id}")
        print(f'Logfile saved to: {c["DIRS"]["OUTDIR"]}output_{self.run_id}.log')

        # Load in all necessary files (backgrounds, psfs, catalogues, etc)
        self.bgs = gp.Cutouts.from_file(f'{c["FILE_DIR"]}{c["FILES"]["BACKGROUNDS"]}', 
                                             logger=self.logger)
        
        self.psfs = gp.Cutouts.from_file(f'{c["FILE_DIR"]}{c["FILES"]["PSFS"]}', logger=self.logger)
        self.psfs.get_ra_dec(ra_key=c["PSFS"]["PSF_RA"], dec_key=c["PSFS"]["PSF_DEC"])
        
        self.table = Table.read(f'{c["FILE_DIR"]}{c["FILES"]["CATALOGUE"]}')
        self.table = gp.trim_table(self.table, c)
        self.logger.info(f'Loaded catalogue with {len(self.table)} entries')
        
        if c["FILES"]["MAG_CATALOGUE"] is not None:
            self.mags = Table.read(f'{c["FILE_DIR"]}{c["FILES"]["MAG_CATALOGUE"]}')
            self.mag_kde = gp.object_kde(self.mags[c["KEYS"]["MAG"]])
        else:
            self.mags = self.mag_kde = None
    

    def run(self, max_bins=None):
        c = self.config

        self.binlist = gp.bin_catalogue(self.table, bin_params=c["BINS"], 
                                        params=c["KEYS"], logger=self.logger)
        max_bins = len(self.binlist.bins) if max_bins is None else min(max_bins, len(self.binlist.bins))

        for i in range(max_bins):
            self.process_bin(self.binlist.bins[i])

        


    def process_bin(self, b):

        self.logger.info(f'Processing bin {b.bin_id()}')
        n_objects = self.config["MODEL"]["N_MODELS"]
        cores = self.config["NCORES"]
        time_limit = self.config["TIME_LIMIT"]

        model = gp.galaxy_models[self.model_type]
        keys, kde = gp.setup_kde(model(), self.config, b.objects)

        results = Parallel(n_jobs=cores, prefer="processes",
                           timeout=time_limit)(delayed(self.process_single)(self.config, 
                                                                            model, 
                                                                            kde, 
                                                                            keys) 
                                                                            for _ in range(n_objects))


    def process_single(self, config, model, kde, keys, gpobj=GalPrimeSingle):
        params = gp.sample_kde(config, keys, kde)
        params = gp.update_required(params, config)
        
        np.random.seed()
        bg = self.bgs.cutouts[np.random.randint(0, len(self.bgs.cutouts))]
        psf = self.psfs.cutouts[np.random.randint(0, len(self.psfs.cutouts))]   # TODO replace with ra/dec matching

        gpobj = gpobj(config, model, params, bg, psf)
        gpobj.process()

        return gpobj


