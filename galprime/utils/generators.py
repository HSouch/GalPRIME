import numpy as np

from astropy.modeling.models import Sersic2D

from ..cutouts import Cutouts
from .fluxes import I_e



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


def add_gaussian_noise(cutouts, mean=0, std=0.5):
    """
    Adds Gaussian noise to the given cutouts.

    Parameters:
    cutouts (array-like): The cutouts to which noise will be added.
    mean (float, optional): The mean of the Gaussian distribution. Default is 0.
    std (float, optional): The standard deviation of the Gaussian distribution. Default is 0.5.

    Returns:
    array-like: The cutouts with added noise.
    """
    cutouts = cutouts.copy()
    for i in range(len(cutouts.cutouts)):
        cutouts.cutouts[i] += np.random.normal(mean, std, cutouts.cutouts[i].shape)
    return cutouts


def add_sersic_objects(cutouts, n_objects=10, reff_range=[1, 5], n_range=[0.5, 4], 
                       mag_range=[20, 25], ellip_range = [0.05, 0.95], min_r=10):
    """
    Add sersic objects to the cutouts
    """
    cutouts = cutouts.copy()

    for i in range(len(cutouts.cutouts)):
        cutout = cutouts.cutouts[i]
        ys, xs = np.mgrid[:cutout.shape[0], :cutout.shape[1]]
        for j in range(n_objects):
            reff = np.random.uniform(reff_range[0], reff_range[1])
            n = np.random.uniform(n_range[0], n_range[1])
            mag = np.random.uniform(mag_range[0], mag_range[1])
            ellip = np.random.uniform(ellip_range[0], ellip_range[1])

            x0, y0 = np.random.uniform(0, cutout.shape[0]), np.random.uniform(0, cutout.shape[1])
            if np.sqrt((x0 - cutout.shape[0] / 2) ** 2 + (y0 - cutout.shape[1] / 2) ** 2) < min_r:
                continue

            i_r50 = I_e(mag, reff, n=n)

            theta = np.random.uniform(0, 2 * np.pi)

            model = Sersic2D(amplitude=i_r50, r_eff=reff, n=n, x_0=x0, y_0=y0, ellip=ellip, theta=theta)

            cutout += model(xs, ys)

        cutouts.cutouts[i] = cutout

    return cutouts