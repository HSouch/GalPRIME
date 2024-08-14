import numpy as np

from .cutouts import Cutouts
from . import utils

from astropy.modeling.models import Sersic2D


def check_sersic(*args):
    mag, r_eff, n, ellip = args
    good = (0.1 < n < 6)
    good = good and (0.1 < ellip < 0.9)
    return good


def gen_single_sersic(config, kde, mag_kde=None, check_method=check_sersic, **kwargs):
    max_tries = kwargs.get("ntries", 100)
    arcconv = float(config["MODEL"]["ARCCONV"])
    size = float(config["MODEL"]["SIZE"])

    zpm = float(config["MODEL"]["ZPM"])
    
    xs, ys = np.mgrid[:size, :size]

    try_index = 0
    while try_index < max_tries:
        mag, r_eff, n, ellip = kde.resample(size=1)
        mag, r_eff, n, ellip = mag[0], r_eff[0], n[0], ellip[0]
        
        if check_method(mag, r_eff, n, ellip):
            break
        else:
            try_index += 1
    if mag_kde is not None:
        mag = mag_kde.resample(size=1)[0]
    
    r_eff_pix = r_eff / arcconv
    
    Ltot = utils.Ltot(mag, zpm)

    theta = np.random.uniform(0, 2 * np.pi)

    model = Sersic2D(amplitude=1, r_eff=r_eff_pix, n=n, x_0=size/2, y_0=size/2, ellip=ellip, theta=theta)
    z = model(xs, ys)

    z *= Ltot / np.nansum(z)

    model_info = {"MAG": mag, "R50": r_eff, "N": n, "ELLIP": ellip, "R50_PIX": r_eff_pix}
    
    return z, model_info



def gen_models(config, kde, mag_kde=None, n_models=100, check_method=check_sersic, **kwargs):
    max_tries = kwargs.get("ntries", 100)
    arcconv = float(config["MODEL"]["ARCCONV"])
    size = float(config["MODEL"]["SIZE"])

    zpm = float(config["MODEL"]["ZPM"])
    
    xs, ys = np.mgrid[:size, :size]

    models, models_info = [], []

    for i in range(n_models):
        z, model_info = gen_single_sersic(config, kde, mag_kde, check_method, ntries=max_tries)
        models.append(z)
        models_info.append(model_info)
    
    
    return Cutouts(cutouts=models, cutout_data=models_info)
