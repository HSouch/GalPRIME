from ..utils import object_kde

from ..models import galaxies


def setup_kde(model, config, table):
    required_keys = model.required_keys()
    
    # Check if all required keys are present in config
    for key in required_keys:
        if key not in config["KEYS"]:
            raise ValueError(f"Missing key {key} in config")


    columns = [table[config["KEYS"][key]] for key in required_keys]
    kde = object_kde(columns)

    return list(required_keys), kde



def sample_kde(config, keys, kde, model = galaxies.SingleSersicModel(), 
               attempts=100):
    verified = False

    while not verified and attempts > 0:

        sample = kde.resample(size=1).transpose()[0]
        params = dict(zip(keys, sample))
        verified = model.verifier.verify(params)
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
