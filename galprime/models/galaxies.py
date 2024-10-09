
from astropy.modeling.models import Sersic2D

import numpy as np

from .. import utils
from .. import cutouts


class GalaxyModel:
    def __init__(self, defaults={}):
        self.defaults = defaults

    def generate(self, params):
        """ Generate a model, and handle user-inputted and default parameters.

        Args:
            params (_type_): _description_
        """
        for key in self.defaults:
            if key not in params:
                params[key] = self.defaults[key]

        return(self._generate(**params))

    def _generate(self, **params):
        # Subclass-specific implementation of the model generation
        raise NotImplementedError("Abstract class")


class SingleSersicModel(GalaxyModel):
    """
    SingleSersicModel is a class that represents a single Sersic galaxy model.
    Attributes:
        params (dict): A dictionary to store parameters of the model.
        defaults (dict): A dictionary containing default values for the model parameters:
            - "MAG" (int): Magnitude of the galaxy, default is 22.
            - "REFF" (int): Effective radius of the galaxy, default is 1.
            - "N" (int): Sersic index, default is 1.
            - "ELLIP" (float): Ellipticity of the galaxy, default is 0.3.
    Methods:
        __init__():
            Initializes the SingleSersicModel with default parameters.
        _generate(**params):
            Generates a single Sersic model with the given parameters.
            Updates the params attribute with the generated model parameters.
            Args:
                **params: Arbitrary keyword arguments representing model parameters.
            Returns:
                tuple: A tuple containing the generated model and the input parameters.
    """

    def __init__(self):       
        self.params = {}
        self.defaults = {
            "MAG": 22,
            "REFF": 1,
            "N": 1,
            "ELLIP": 0.3,
        }
    
    def _generate(self, **params):
        mod, mod_params = gen_single_sersic(**params)
        self.params.update(mod_params)
        return mod, params


class ExponentialDiskModel(GalaxyModel):
    def __init__(self):
        self.params = {}
        self.defaults = {
            "MAG": 22,
            "REFF": 1,
            "ELLIP": 0.3,
        }

    def _generate(self, **params):
        mod, mod_params = gen_single_sersic(**params)
        self.params.update(mod_params)
        return mod, params



def gen_single_sersic(**kwargs):
    shape = kwargs.get("SHAPE", (101, 101))
    if not isinstance(shape, tuple):
        shape = (shape, shape)
    x_0 = kwargs.get("x_0", shape[0] / 2)
    y_0 = kwargs.get("y_0", shape[1] / 2)


    mod = Sersic2D(amplitude=1, r_eff=kwargs.get("REFF", 1), 
                   n=kwargs.get("N", 1), 
                   x_0=x_0, 
                   y_0=y_0, 
                   ellip=kwargs.get("ELLIP", 0.3), theta=kwargs.get("PA", np.random.uniform(0, np.pi)))
    ys, xs = np.mgrid[:shape[0], :shape[1]]
    z = mod(xs, ys) 

    mag, m0 = kwargs.get("MAG", 22), kwargs.get("M0", 27)

    z *= utils.Ltot(mag, m0=m0) / np.sum(z)

    params = {
        "MAG": mag, "M0": m0,
        "REFF": mod.r_eff.value, "N": mod.n.value,
        "ELLIP": mod.ellip.value, "PA":  mod.theta.value,
        "X0": x_0,  "Y0": y_0,
        "SHAPE": shape,
    }
    
    return z, params
