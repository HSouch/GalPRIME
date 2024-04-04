import numpy as np

from astropy.modeling.models import Sersic2D

from ..cutouts import Cutouts
from ..utils import I_e


def gen_sersic_models(n_models = 50, width=(81, 81),
                      mag_range = [24, 29], 
                      n_range = [1, 4],
                      re_range = [2, 10],
                      m0=27, **kwargs):
    
    if not isinstance(width, tuple):
        width = (width, width)

    y,x = np.mgrid[:width[0], :width[1]]
    cutouts, cutout_data = [], []
    for i in range(n_models):
        mag = np.random.uniform(*mag_range)
        n = np.random.uniform(*n_range)
        re = np.random.uniform(*re_range)
        ellip = np.random.uniform(0.2, 0.8)
        pa = np.random.uniform(0, 2 * np.pi)
        
        sersic = Sersic2D(amplitude=I_e(mag, re, n), r_eff=re, n=n, x_0=width[1]/2, y_0=width[0]/2, ellip=ellip, 
                          theta=pa)
        z = sersic(x, y)
        cutouts.append(z)
        cutout_data.append({"mag": mag, "n": n, "re": re})
    return Cutouts(cutouts=cutouts, cutout_data=cutout_data)