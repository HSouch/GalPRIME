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

from joblib import Parallel


class GPrime:
    def __init__(self, config_filename, **kwargs):
        self.config = c = config.read_config_file(config_filename)
        self.binlist = None
        self.run_id = kwargs.get("run_id", np.random.randint(1e3, 1e4))
        self.outfiles = gp.gen_filestructure(c["DIRS"]["OUTDIR"])

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
        
        if c["FILES"]["MAG_CATALOGUE"] is not None:
            self.mags = Table.read(f'{c["FILE_DIR"]}{c["FILES"]["MAG_CATALOGUE"]}')
            self.mag_kde = gp.object_kde(self.mags[c["KEYS"]["MAG"]])
        else:
            self.mags = self.mag_kde = None
    
    def run(self, max_bins=None):
        c = self.config

        self.binlist = gp.bin_catalogue(self.table, bin_params=c["BINS"], params=c["KEYS"], logger=self.logger)
        max_bins = len(self.binlist.bins) if max_bins is None else min(max_bins, len(self.binlist.bins))

        model = gp.galaxy_models[self.model_type]

        for i in range(max_bins):
            self.process_bin(c, self.binlist.bins[i])


    def process_bin(self, config, b):
        for index in range(config["MODEL"]["N_MODELS"]):
            model = gp.galaxy_models[self.model_type]()
            
            keys, kde = gp.setup_kde(model, config, self.table)
            params = gp.sample_kde(config, keys, kde)


class GalPrimeSingle:
    def __init__(self, config, model, params):
        self.config = config
        self.model = model
        self.params = params
    
    def process():
        pass