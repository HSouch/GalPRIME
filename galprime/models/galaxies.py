
from astropy.modeling.models import Sersic2D, Gaussian2D

import numpy as np

from .. import utils

from . import verifiers


class GalaxyModel:
    def __init__(self, defaults={}):
        self.defaults = defaults
        self.verifier = verifiers.DefaultVerifier()

    def generate(self, params={}, **kwargs):
        """ Generate a model, and handle user-inputted and default parameters.

        Args:
            params (_type_): _description_
        """
        
        # Set model parameters
        for key in self.defaults:
            if key not in params:
                params[key] = self.defaults[key]
        # Override params directly specified by the user in keyword arguments
        for key in kwargs:
            params[key] = kwargs[key]
        
        if not self.verifier.verify(params):
            raise ValueError("Invalid parameters")

        return(self._generate(**params))
    
    def required_keys(self):
        return self.defaults.keys()
    

    def _generate(self, **params):
        # Subclass-specific implementation of the model generation
        raise NotImplementedError("Abstract class")
    

    def generate_param_text(self, exclude=[]):
        """ Generate a text representation of the model parameters.

        Args:
            exclude (list, optional): The parameters to exclude from the text. Defaults to [].

        Returns:
            str: A string representation of the model parameters
        """
        text = ""
        for key, val in self.params.items():
            if key in exclude:
                continue
            try:
                text += f"{key}: {val:.4f}\n"
            except Exception as e:
                print("Error formatting parameter:", key, val)
                continue
        return text.removesuffix("\n")


    def plot_model_params(self, axis, exclude=[], loc="upper right", pad=0.02, fontsize=10):
        """ Plot the model parameters on the given axis.

        Args:
            axis (matplotlib.axis): The axis to plot the parameters on.
            exclude (list, optional): The parameters to exclude from the plot. Defaults to [].
            loc (str, optional): The location of the text. Defaults to "upper right".
            pad (float, optional): The padding from the axis limits. Defaults to 0.02.
            fontsize (int, optional): The font size of the text. Defaults to 10.
        """
        
        ul, lr = loc.lower().split(" ")

        ha = "left" if lr == "left" else "right"
        va = "top" if ul == "upper" else "bottom"

        x = 0.02 if lr == "left" else 0.98
        y = 0.98 if ul == "upper" else 0.02

        text = self.generate_param_text(exclude=exclude)
        axis.text(x, y, text,
                    ha = ha, va = va,
                    fontsize=fontsize, 
                    transform=axis.transAxes,
                    bbox=dict(facecolor='white', alpha=0.5))


class SingleSersicModel(GalaxyModel):
    """
    SingleSersicModel is a class that generates a single Sersic galaxy model.
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
        self.verifier = verifiers.DefaultVerifier()
    
    def _generate(self, **params):
        mod, mod_params = gen_single_sersic(**params)
        self.params.update(mod_params)
        return mod, mod_params
    
    
class ExponentialDiskModel(GalaxyModel):
    """
    A model representing an exponential disk galaxy.
    Attributes:
        params (dict): A dictionary to store model parameters.
        defaults (dict): A dictionary containing default values for model parameters.
            - "MAG" (int): Default magnitude value.
            - "REFF" (int): Default effective radius value.
            - "ELLIP" (float): Default ellipticity value.
    Methods:
        __init__():
            Initializes the ExponentialDiskModel with default parameters.
        _generate(**params):
            Generates a single Sersic model with the given parameters.
            Updates the model parameters with the generated values.
            Args:
                **params: Arbitrary keyword arguments for model parameters.
            Returns:
                tuple: A tuple containing the generated model and the updated parameters.
    """
    
    def __init__(self):
        self.params = {}
        self.defaults = {
            "MAG": 21,
            "REFF": 1,
            "ELLIP": 0.3,
            "N": 1,
            "REFF_UNIT": "pix",
            "ARCCONV": 1,
        }
        self.verifier = verifiers.DefaultVerifier()

    def _generate(self, **params):
        params["N"] = 1
        mod, mod_params = gen_single_sersic(**params)
        self.params.update(mod_params)
        return mod, params
    

class EllipticalModel(GalaxyModel):
    def __init__(self):
        self.params = {}
        self.defaults = {
            "MAG": 21,
            "REFF": 1,
            "ELLIP": 0.3,
            "N": 4,
        }
        self.verifier = verifiers.DefaultVerifier()

    def _generate(self, **params):
        mod, mod_params = gen_single_sersic(**params)
        params["N"] = 4
        self.params.update(mod_params)
        return mod, params


class BulgeDiskSersicModel(GalaxyModel):
    """ 
    A model representing a galaxy with a bulge and disk component, each described by a Sersic profile.
    Attributes:
        params (dict): Dictionary to store model parameters.
        defaults (dict): Dictionary containing default values for model parameters.
        verifier (BulgeDiskVerifier): An instance of the BulgeDiskVerifier class to verify model parameters.
    Methods:
        __init__():
            Initializes the BulgeDiskSersicModel with default parameters and a verifier.
        get_bulge_disk_mags(**params):
            Calculate the bulge and disk magnitudes given the total magnitude and the bulge fraction.
            Args:
                **params: Arbitrary keyword arguments containing model parameters.
            Returns:
                tuple: A tuple containing the bulge magnitude and disk magnitude.
        _generate(**params):
            Generate the bulge and disk models based on the provided parameters.
            Args:
                **params: Arbitrary keyword arguments containing model parameters.
            Returns:
                tuple: A tuple containing the combined model and a dictionary of parameters 
                with bulge and disk suffixes.
    """

    def __init__(self):
        self.params = {}
        self.defaults = {
            "MAG": 16,          # Total Magnitude
            "FBULGE": 0.5,      # Bulge Fraction
            "REFF_BULGE": 1,    # Bulge Effective Radius in arcseconds
            "N_BULGE": 4,       # Bulge Sersic Index
            "REFF_DISK": 1,     # Disk Effective Radius in arcseconds
            "ELLIP_DISK": 0.2,       # Ellipticity
            "ELLIP_BULGE": 0.2,      # Ellipticity
        }
        self.verifier = verifiers.BulgeDiskVerifier()
    
    def get_bulge_disk_mags(self, **params):
        """
        Calculate the bulge and disk magnitudes given the total magnitude
        and the bulge fraction.
        """
        mag = params["MAG"]
        fb = params["FBULGE"]
        m_bulge = mag - 2.5 * np.log10(fb)
        m_disk = mag - 2.5 * np.log10(1 - fb)
        return m_bulge, m_disk
    

    def _generate(self, **params):
        bulge_mag, disk_mag = self.get_bulge_disk_mags(**params)

        params["PA"] = params.get("PA", np.random.uniform(0, np.pi))

        # Generate the bulge model
        bulge_params = params.copy()
        bulge_params["MAG"] = bulge_mag
        bulge_params["REFF"] = bulge_params["REFF_BULGE"]
        bulge_params["N"] = bulge_params["N_BULGE"]
        bulge_params["ELLIP"] = bulge_params["ELLIP_BULGE"]
        bulge, bulge_params = gen_single_sersic(**bulge_params)

        # Generate the disk model
        disk_params = params.copy()
        disk_params["MAG"] = disk_mag
        disk_params["REFF"] = disk_params["REFF_DISK"]
        disk_params["N"] = 1
        disk_params["ELLIP"] = disk_params["ELLIP_DISK"]
        disk, disk_params = gen_single_sersic(**disk_params)

        # Combine bulge and disk parameters with suffixes
        bulge_params = {f"BULGE_{k}": v for k, v in bulge_params.items()}
        disk_params = {f"DISK_{k}": v for k, v in disk_params.items()}

        # Combine the two models
        model = bulge + disk
        return model, {**params, **bulge_params, **disk_params}


def gen_single_sersic(**kwargs):
    """
    Generate a single Sersic model. This method is configured to more easily work with GalPRIME model classes,
    instead of just calling the Sersic2D model directly.
    Args:
        **kwargs: Arbitrary keyword arguments representing model parameters.
            - "MAG" (int): Magnitude of the galaxy, default is 22.
            - "REFF" (int): Effective radius of the galaxy, default is 1.
            - "N" (int): Sersic index, default is 1.
            - "ELLIP" (float): Ellipticity of the galaxy, default is 0.3.
            - "PA" (float): Position angle of the galaxy, default is a random value between 0 and pi.
            - "SHAPE" (tuple): Shape of the output array, default is (101, 101).
            - "x_0" (float): X-coordinate of the galaxy center, default is half of the shape's width.
            - "y_0" (float): Y-coordinate of the galaxy center, default is half of the shape's height.
            - "M0" (int): Zero-point magnitude, default is 27.

    Returns:
        tuple: A tuple containing the generated model array and a dictionary of parameters.
    """
    shape = kwargs.get("SHAPE", (101, 101))
    if not isinstance(shape, tuple):
        shape = (shape, shape)
    x_0 = kwargs.get("x_0", shape[0] / 2)
    y_0 = kwargs.get("y_0", shape[1] / 2)
    
    mag, m0 = kwargs.get("MAG", 22), kwargs.get("M0", 27)

    ELLIP = kwargs.get("ELLIP", 0.3)

    REFF = kwargs.get("REFF", 1)
    # REFF_CIRC = utils.r_circ(REFF, ELLIP)

    mod = Sersic2D(amplitude=1, r_eff=REFF, 
                   n=kwargs.get("N", 1), 
                   x_0=x_0, 
                   y_0=y_0, 
                   ellip=ELLIP, theta=kwargs.get("PA", np.random.uniform(0, np.pi)))
    ys, xs = np.mgrid[:shape[0], :shape[1]]
    
    z = mod(xs, ys) 

    # TODO : Modularize this so it can be more easily modified by the user.
    z_within_10Re = np.copy(z)
    ys, xs = np.indices(z_within_10Re.shape)
    rs = np.sqrt((xs - x_0)**2 + (ys - y_0)**2)
    z_within_10Re[rs > 10 * REFF] = np.nan

    multiplier = utils.Ltot(mag, m0=m0) / np.nansum(z_within_10Re)

    z *= multiplier


    params = {
        "MAG": mag, "M0": m0,
        "REFF": mod.r_eff.value, "N": mod.n.value,
        "ELLIP": mod.ellip.value, "PA":  mod.theta.value,
        "X0": x_0,  "Y0": y_0,
        "SHAPE": shape,
    }
    
    return z, params


def gen_gaussian(**kwargs):
    """
    Generate a single 2D Gaussian model.
    Args:
        **kwargs: Arbitrary keyword arguments representing model parameters.
            - "MAG" (int): Magnitude of the galaxy, default is 22.
            - "STDDEV" (int): Effective radius of the galaxy, default is 1.
            - "ELLIP" (float): Ellipticity of the galaxy, default is 0.3.
            - "PA" (float): Position angle of the galaxy, default is a random value between 0 and pi.
            - "SHAPE" (tuple): Shape of the output array, default is (101, 101).
            - "X_0" (float): X-coordinate of the galaxy center, default is half of the shape's width.
            - "Y_0" (float): Y-coordinate of the galaxy center, default is half of the shape's height.
            - "M0" (int): Zero-point magnitude, default is 27.

    Returns:
        tuple: A tuple containing the generated model array and a dictionary of parameters.
    """
    shape = kwargs.get("SHAPE", (101, 101))
    if not isinstance(shape, tuple):
        shape = (shape, shape)
    x_0 = kwargs.get("x_0", shape[0] / 2)
    y_0 = kwargs.get("y_0", shape[1] / 2)
    
    ellip = kwargs.get("ELLIP", 0.1)
    x_stddev = kwargs.get("STDDEV", 5)
    y_stddev = x_stddev * (1 - ellip)
    pa = kwargs.get("PA", 0)


    mod = Gaussian2D(amplitude=1, x_mean=x_0, y_mean=y_0, 
                     x_stddev=x_stddev, y_stddev=y_stddev, 
                     theta=pa)
    ys, xs = np.mgrid[:shape[0], :shape[1]]
    z = mod(xs, ys) 
    
    mag, m0 = kwargs.get("MAG", 22), kwargs.get("M0", 27)

    z *= utils.Ltot(mag, m0=m0) / np.sum(z)

    params = {
        "MAG": mag, "M0": m0,
        "REFF": mod.x_stddev.value,
        "ELLIP": ellip, "PA":  mod.theta.value,
        "X0": x_0,  "Y0": y_0,
        "SHAPE": shape,
    }
    
    return z, params


galaxy_models = {1: SingleSersicModel, 
                 2: BulgeDiskSersicModel, 
                 3: ExponentialDiskModel, 
                 4: EllipticalModel}
