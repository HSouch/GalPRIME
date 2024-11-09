import galprime as gp

from scipy.signal import convolve2d
from scipy.interpolate import interp1d

from astropy.io import fits
from astropy.table import Table
from astropy.visualization import ZScaleInterval

import numpy as np

import logging

from matplotlib import pyplot as plt

import time


import warnings


class GPrimeSingle:
    def __init__(self, config, model, params, bg=None, psf=None, logger=None):
        self.config = config
        self.model = model
        self.params = params

        self.id = np.random.randint(1e9, 1e10)

        self.bg = bg
        self.psf = psf

        self.logger = logger

        self.stop_code = 0
        self.isophote_lists = []

    def process(self):
        # Generate model and convolve with PSF
        try:
            self.model_image, self.model_params = self.model.generate(self.params)
            self.convolved_model = gp.convolve_model(self.model_image, self.psf)
        except Exception as e:
            self.stop_code = 1
            return

        # Add model to background, estimate background, and subtract
        try:    
            self.bg_added_model = self.convolved_model + self.bg
            self.source_mask, self.background = gp.estimate_background_2D(self.bg_added_model, self.config)
            
            self.bg_model = self.background.background
            self.bgsub = self.bg_added_model - self.bg_model
        except Exception as e:
            self.stop_code = 2
            return

        # Mask image(s)
        try:    
            self.mask_bgadded, self.mask_data_bgadded = gp.gen_mask(self.bg_added_model, config=self.config)
            self.mask_bgsub, self.mask_data_bgsub = gp.gen_mask(self.bgsub, config=self.config)
        except Exception as e:
            self.stop_code = 3
            return
        
        # Extract profiles
        try:
            # Extract profiles
            for dataset in [self.convolved_model, 
                            np.ma.array(self.bg_added_model, mask=self.mask_bgadded), 
                            np.ma.array(self.bgsub, mask=self.mask_bgsub)]:
                isolist = gp.isophote_fitting(dataset, self.config)
                self.isophote_lists.append(isolist)

        except Exception as e:
            self.stop_code = 4
            return

        self.stop_code = 10
