
from astropy.modeling.models import Sersic2D

import numpy as np

from .. import utils
from .. import cutouts


class GalaxyModel:
    def __init__(self, config, params={}, size=151, **kwargs):
        self.config = config
        self.params = params

        self.size = size

    def verify_params(self):
        raise NotImplementedError("Abstract class")

    def generate(self):
        raise NotImplementedError("Abstract class")
    
    def gen_multiple(self, config, kde, n_models=100, **kwargs):
        raise NotImplementedError("Abstract class")


class SingleSersicModel(GalaxyModel):
    def __init__(self, size=151, params={}, **kwargs):
        super().__init__(params, size=size, **kwargs)

    def generate(self):
        mag, r50, n, ellip = self.params["mag"], self.params["r50"], self.params["n"], self.params["ellip"]
        theta = np.random.uniform(0, 2 * np.pi)

        ltot = utils.Ltot(mag, self.config["MODEL"]["ZPM"])
        ys, xs = np.mgrid[:self.size, :self.size]

        z = Sersic2D(amplitude=1, r_eff=r50, n=n, x_0=self.size/2, y_0=self.size/2, ellip=ellip, theta=theta)(xs, ys)
        z *= ltot / np.nansum(z)

        return z
    
    def verify_params(self, mag, r_eff, n, ellip):
        good = (0.1 < n < 6)
        good = good and (0.1 < ellip < 0.9)
        return good
    
    def gen_multiple(self, config, kde, mag_kde=None, n_models=100, **kwargs):
        max_tries = kwargs.get("ntries", 100)
        arcconv = float(config["MODEL"]["ARCCONV"])
        
        zpm = float(config["MODEL"]["ZPM"])

        xs, ys = np.mgrid[:self.size, :self.size]

        models, models_info = [], []
        for _ in range(n_models):
            try_index = 0

            while try_index < max_tries:
                mag, r_eff, n, ellip = kde.resample(size=1)
                mag, r_eff, n, ellip = mag[0], r_eff[0], n[0], ellip[0]

                if self.verify_params(mag, r_eff, n, ellip):
                    break
                else:
                    try_index += 1

            if mag_kde is not None:
                mag = mag_kde.resample(size=1)[0]
            
            r_eff_pix = r_eff / arcconv
            Ltot = utils.Ltot(mag, zpm)

            theta = np.random.uniform(0, 2 * np.pi)

            model = Sersic2D(amplitude=1, r_eff=r_eff_pix, n=n, 
                             x_0=self.size/2, y_0=self.size/2, 
                             ellip=ellip, theta=theta)
            z = model(xs, ys)

            z *= Ltot / np.nansum(z)

            model_info = {"MAG": mag, "R50": r_eff, "N": n, "ELLIP": ellip, "R50_PIX": r_eff_pix}

            models.append(z)
            models_info.append(model_info)
    
        return cutouts.Cutouts(cutouts=models, cutout_data=models_info)
    

class BDSersicModel(GalaxyModel):
    def __init__(self, params={}, **kwargs):
        super().__init__(params, **kwargs)

    def generate(self):
        pass