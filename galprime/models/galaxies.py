
from astropy.modeling.models import Sersic2D

import numpy as np

from .. import utils


class GalaxyModel:
    def __init__(self, config, params={}, size=151, **kwargs):
        self.config = config
        self.params = params

        self.size = size

    def verify_params(self):
        raise NotImplementedError("Abstract class")

    def generate(self):
        raise NotImplementedError("Abstract class")


class SingleSersicModel(GalaxyModel):
    def __init__(self, params={}, **kwargs):
        super().__init__(params, **kwargs)

    def generate(self):
        mag, r50, n, ellip = self.params["mag"], self.params["r50"], self.params["n"], self.params["ellip"]
        theta = np.random.uniform(0, 2 * np.pi)

        ltot = utils.Ltot(mag, self.config["MODEL"]["ZPM"])

        ys, xs = np.mgrid[:self.size, :self.size]

        z = Sersic2D(amplitude=1, r_eff=r50, n=n, x_0=self.size/2, y_0=self.size/2, ellip=ellip, theta=theta)(xs, ys)
        z *= ltot / np.nansum(z)

        return z
    

class BDSersicModel(GalaxyModel):
    def __init__(self, params={}, **kwargs):
        super().__init__(params, **kwargs)

    def generate(self):
        pass