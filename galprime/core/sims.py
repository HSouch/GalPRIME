import galprime as gp

import numpy as np

import warnings


class GPrimeSingle:
    """ A single instance of a GalPRIME iteration. """
    
    def __init__(self, config, model, params, bg=None, psf=None, 
                 logger=None, id=None, save_output=False, metadata={}):
        self.config = config
        self.model = model
        self.params = params

        self.save_output = save_output

        self.id = id if id is None else np.random.randint(1e9, 1e10)

        self.bg = bg
        self.psf = psf

        self.logger = logger

        self.stop_code = 0
        self.isophote_lists = []

        self.metadata = metadata


    def process(self):
        # Generate model and convolve with PSF
        try:
            self.stop_code = 1
            self.model_image, self.model_params = self.model.generate(self.params)
            self.convolved_model = gp.convolve_model(self.model_image, self.psf)
        except Exception as e:
            raise RuntimeError(f'{self.id} failed convolution: {e}')

        # Add model to background, estimate background, and subtract
        try:    
            self.stop_code = 2
            self.bg_added_model = self.convolved_model + self.bg
            self.source_mask, self.background = gp.estimate_background_2D(self.bg_added_model, self.config)
            
            self.bg_model = self.background.background
            self.bgsub = self.bg_added_model - self.bg_model
        except Exception as e:
            raise RuntimeError(f'{self.id} failed bg-modeling: {e}')

        # Mask image(s)
        try:    
            self.stop_code = 3
            self.mask_bgadded, self.mask_data_bgadded = gp.gen_mask(self.bg_added_model, config=self.config)
            self.mask_bgsub, self.mask_data_bgsub = gp.gen_mask(self.bgsub, config=self.config)
        except Exception as e:
            raise RuntimeError(f'{self.id} failed masking: {e}')
        
        # Extract profiles
        try:
            self.stop_code = 4
            # Extract profiles
            for dataset in [self.convolved_model, 
                            np.ma.array(self.bg_added_model, mask=self.mask_bgadded), 
                            np.ma.array(self.bgsub, mask=self.mask_bgsub)]:
                with warnings.catch_warnings():     # Suppress warnings from astropy fitting
                    warnings.simplefilter("ignore")
                    isolist = gp.isophote_fitting(dataset, self.config)
                    self.isophote_lists.append(isolist)

        except Exception as e:
            raise RuntimeError(f'{self.id} failed extraction: {e}')
        
        self.stop_code = 10

    def condensed_output(self):
        try:
            output = {"ISOLISTS": [n["ISOLIST"] for n in self.isophote_lists],
                "PARAMS": self.model_params, 
                "ITERATION": self.metadata.get("ITERATION", 999),
                "BG_INDEX": self.metadata.get("BG_INDEX", 999),
                "PSF_INDEX": self.metadata.get("PSF_INDEX", 999)}
        except Exception as e:
            raise ValueError(f"Error when creating condensed output: {e}")
        return output

gprime_stop_codes = {
    0: "NOT RUN",
    1: "Started Model Gen",
    2: "Started BG handling",
    3: "Started masking",
    4: "Started profile extraction",

    10: "Finished process"
}