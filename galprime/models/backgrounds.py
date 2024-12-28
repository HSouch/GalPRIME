from ..core import cutouts

from astropy.modeling import models

import numpy as np

# TODO - Implement the BackgroundModels class

class BackgroundModels(cutouts.Cutouts):

    def __init__(self, **kwargs):
        super().__init__(cutouts=[], cutout_data=[], metadata={})

    def populate(self, n_models, **kwargs):
        size = kwargs.get("size", 100)
        for i in range(n_models):
            mod, mod_info = self.gen_single(**kwargs)

            self.cutouts.append(mod)
            self.cutout_data.append(mod_info)
    
    def gen_single(self):
        raise NotImplementedError("Abstract class")
    

class GaussianField(BackgroundModels):

    def __init__(self):
        super().__init__()
    
    def gen_single(self, **kwargs):
        bg_info = {}
        n_objects = kwargs.get('n_objects', 20)
        std_range = kwargs.get('std_range', (1, 5))
        amp_range = kwargs.get('amp_range', (0.1, 10))
        size = kwargs.get('size', 151)
        x0, y0 = kwargs.get('x0', size/2), kwargs.get('y0', size/2)

        min_dist = kwargs.get('min_dist', 5)
        max_attempts = kwargs.get('max_attempts', 100)

        bg = np.zeros((size, size))

        ys, xs =  np.mgrid[0:size, 0:size]

        for i in range(n_objects):
            attempts, good = 0, False
            while not good:
                if attempts > max_attempts:
                    break
                x, y = np.random.uniform(0, size, 2)
                if np.sqrt((x - x0)**2 + (y - y0)**2) < min_dist:
                    attempts += 1
                    continue
                std = np.random.uniform(std_range[0], std_range[1])
                amp = np.random.uniform(amp_range[0], amp_range[1])
                g = models.Gaussian2D(amp, x, y, std, std)
                bg += g(xs, ys)
                good = True
        return bg, bg_info 
            
    

class XGradientBackgrounds(BackgroundModels):
    pass


class YGradientBackgrounds(BackgroundModels):
    pass


class DualGradientBackgrounds(BackgroundModels):
    pass


class RadialGradientBackgrounds(BackgroundModels):
    pass
