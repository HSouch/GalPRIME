from ..utils import object_kde

import numpy as np


def setup_kde(model, config, table):
    required_keys = model.required_keys()
    config_keys = config.keys()
    
    # Check if all required keys are present in config
    for key in required_keys:
        if key not in config["KEYS"]:
            raise ValueError(f"Missing key {key} in config")


    columns = [table[config["KEYS"][key]] for key in required_keys]
    kde = object_kde(columns)

    return list(required_keys), kde


def default_verifier(params):
    good = params["REFF"] > 0               # Check that the effective radius is positive
    good = good and 0 < params["N"] < 10    # Check that the Sersic index is within a reasonable range
    good = good and params["MAG"] > 0       # Check that the magnitude is positive (abs. magnitude)

    return good


def sample_kde(config, keys, kde, verifier=default_verifier, attempts=100):
    verified = False

    while not verified and attempts > 0:

        sample = kde.resample(size=1).transpose()[0]
        params = dict(zip(keys, sample))
        verified = verifier(params)
        attempts -= 1
    
    if not verified:
        raise ValueError(f"Failed to sample valid parameters after {attempts} attempts")

    params = update_required(params, config)

    return params


def update_required(params, config):
    """
    Updates the 'params' dictionary with values from the 'config' dictionary.
    Args:
        params (dict): The dictionary to be updated.
        config (dict): The configuration dictionary containing the new values.
    Returns:
        dict: The updated 'params' dictionary.
    """

    params["SHAPE"] = config["MODEL"]["SIZE"]   # Update the size of the model
    params["M0"] = config["MODEL"]["ZPM"]       # Update the zero-point magnitude
    if config["MODEL"]["REFF_UNIT"].lower() == "arcsec":
        params["REFF"] /= config["MODEL"]["ARCCONV"]    # Update the effective radius to pixels

    return params


class SyntheticKDE:
    def __init__(self, kde=None, **kwargs):
        self.kde = kde


class UniformSSKDE(SyntheticKDE):
    """
    A class representing a uniform synthetic KDE (Kernel Density Estimation) model.
    
    This class inherits from the `SyntheticKDE` class and provides a method to generate a synthetic KDE
    based on uniform distributions of magnitude, r50 (half-light radius), Sersic index, and ellipticity.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    
    def generate(self, mag_lims=[20, 24], r50_lims=[0.5, 3], n_lims=[0.5, 4.5], ellip_lims=[0.05, 0.95], 
                 size=100):
        """
        Generate a synthetic KDE based on uniform distributions of magnitude, r50 (half-light radius),
        Sersic index, and ellipticity.
        
        Parameters:
        - mag_lims (list): The lower and upper limits of the magnitude distribution. Default is [20, 24].
        - r50_lims (list): The lower and upper limits of the r50 distribution. Default is [2, 10].
        - n_lims (list): The lower and upper limits of the Sersic index distribution. Default is [0.5, 4.5].
        - ellip_lims (list): The lower and upper limits of the ellipticity distribution. Default is [0.05, 0.95].
        - size (int): The number of samples to generate. Default is 100.
        
        Returns:
        - kde (object): The synthetic KDE object generated based on the specified distributions.
        """
        mags = np.random.uniform(mag_lims[0], mag_lims[1], size)
        r50s = np.random.uniform(r50_lims[0], r50_lims[1], size)
        ns = np.random.uniform(n_lims[0], n_lims[1], size)
        ellips = np.random.uniform(ellip_lims[0], ellip_lims[1], size)
        self.kde = object_kde([mags, r50s, ns, ellips])

        return self.kde
    
    
class GaussianSSKDE(SyntheticKDE):
    """
    A class representing a Gaussian synthetic KDE (Kernel Density Estimation) model.
    
    This class inherits from the `SyntheticKDE` class and provides a method to generate a synthetic KDE
    based on Gaussian distributions of magnitude, r50 (half-light radius), Sersic index, and ellipticity.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    
    def generate(self, mag_mean=22, mag_std=1, r50_mean=6, r50_std=2, n_mean=2.5, n_std=1, ellip_mean=0.5, 
                 ellip_std=0.2, size=100):
        """
        Generate a synthetic KDE based on Gaussian distributions of magnitude, r50 (half-light radius),
        Sersic index, and ellipticity.
        
        Parameters:
        - mag_mean (float): The mean of the magnitude distribution. Default is 22.
        - mag_std (float): The standard deviation of the magnitude distribution. Default is 1.
        - r50_mean (float): The mean of the r50 distribution. Default is 6.
        - r50_std (float): The standard deviation of the r50 distribution. Default is 2.
        - n_mean (float): The mean of the Sersic index distribution. Default is 2.5.
        - n_std (float): The standard deviation of the Sersic index distribution. Default is 1.
        - ellip_mean (float): The mean of the ellipticity distribution. Default is 0.5.
        - ellip_std (float): The standard deviation of the ellipticity distribution. Default is 0.2.
        - size (int): The number of samples to generate. Default is 100.
        
        Returns:
        - kde (object): The synthetic KDE object generated based on the specified distributions.
        """
        mags = np.random.normal(mag_mean, mag_std, size)
        r50s = np.random.normal(r50_mean, r50_std, size)
        ns = np.random.normal(n_mean, n_std, size)
        ellips = np.random.normal(ellip_mean, ellip_std, size)
        self.kde = object_kde([mags, r50s, ns, ellips])

        return self.kde
