
from astropy.modeling.models import Sersic2D

import numpy as np

from .. import utils
from .. import cutouts


class GalaxyModel:
    def __init__(self, defaults={}):
        self.defaults = defaults

    def generate(self, params):
        raise NotImplementedError("Abstract method")


class SingleSersicModel(GalaxyModel):
    def __init__(self):
        
        self.defaults = {
            "REFF": 1,
            "N": 1,
            "ELLIP": 0.3,
            "PA": np.deg2rad(45),
        }


class BulgeDiskModel:
    def __init__(self):
        pass


models = [SingleSersicModel, BulgeDiskModel]

