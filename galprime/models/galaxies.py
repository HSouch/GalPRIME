
from astropy.modeling.models import Sersic2D

import numpy as np

from .. import utils
from .. import cutouts


class GalaxyModel:
    def __init__(self, defaults={}):
        self.defaults = defaults

    def generate(self, params):
        raise NotImplementedError("Abstract method")


class SingleSersicModel:
    def __init__(self):
        pass


class BulgeDiskModel:
    def __init__(self):
        pass


models = [SingleSersicModel, BulgeDiskModel]

