""" Configuration settings methods.

This module contains all methods required to properly set up or read a configuration file for use with TBriDGE.

"""

import urllib.request
from numpy import arange


def load_config_file(filename, verbose_test=False):
    """ Loads in a config file for TBRIDGE to run

    Args:
        filename: Filename (can be absolute or relative path, or a URL) to read config file from.
        verbose_test: If true, will create a helpful printout of params to ensure that the parameters are
            loading in properly.

    Returns:
        dict: Configuration file as a dict object.
    """
    config_values = {}

    # First try to open things locally. If that doesn't work try it as a URL
    try:
        config_lines = open(filename, "r").readlines()
    except FileNotFoundError:
        try:
            r = urllib.request.urlopen(filename)
            config_lines = []
            for line in r:
                config_lines.append(line.decode("utf-8"))
        except:
            print("Failed to get any file")
            return None

    for line in config_lines:
        line = line.strip()
        if len(line) == 0 or line[0] == "#":
            continue
        splits = line.split("=")
        config_values[splits[0].strip()] = splits[1].strip()

    for n in config_values:
        value = config_values[n]
        if value.lower() == "true":
            config_values[n] = True
            continue
        elif value.lower() == "false":
            config_values[n] = False
            continue
    config_values["SIZE"] = int(config_values["SIZE"])
    config_values["EXTRACTION_SIZE"] = int(config_values["EXTRACTION_SIZE"])
    config_values["CORES"] = int(config_values["CORES"])
    config_values["ZEROPOINT"] = int(config_values["ZEROPOINT"])
    config_values["ARC_CONV"] = float(config_values["ARC_CONV"])
    config_values["N_MODELS"] = int(config_values["N_MODELS"])
    config_values["N_BGS"] = int(config_values["N_BGS"])
    config_values["LINEAR_STEP"] = float(config_values["LINEAR_STEP"])
    config_values["ALARM_TIME"] = int(config_values["ALARM_TIME"])
    config_values["CUTOUT_FRACTION"] = float(config_values["CUTOUT_FRACTION"])
    config_values["BOX_SIZE"] = int(config_values["BOX_SIZE"])
    config_values["FILTER_SIZE"] = int(config_values["FILTER_SIZE"])

    config_values["PROVIDED_BG_TYPE"] = "None"

    for n in ("MASS_BINS", "REDSHIFT_BINS", "SFPROB_BINS"):
        """ Turn all bins in numpy aranges (just to simplify the process). Will also add a x_step parameter"""
        value_string = config_values[n].split(",")
        bin_output = arange(float(value_string[0]), float(value_string[1]), float(value_string[2]))

        config_values[n] = bin_output
        config_values[n.split("_")[0] + "_STEP"] = float(value_string[2])

    value_string = config_values["MASK_PARAMS"].split(",")
    config_values["MASK_PARAMS"] = [float(value_string[0]), float(value_string[1]), int(value_string[2])]

    if verbose_test:
        for n in config_values:
            print(n, config_values[n], type(config_values[n]))

    return config_values


def print_config(config):
    """ Print out the configuration file values"""
    for n in config:
        print(n, config[n], type(config[n]))


def default_config_params():
    """
    Dumps a dict object containing all default parameters in proper type format.
    """
    default_params = {
        "VERBOSE": True,
        "TEST_VERBOSE": False,
        "CATALOG": "cat.fits",
        "IMAGE_DIRECTORY": "images/",
        "PSF_FILENAME": "i_psfs.fits",
        "OUT_DIR": "out/",
        "SAVE_CUTOUTS": 'none',
        "CUTOUT_FRACTION": 0.2,

        "MASS_KEY": "MASSES",
        "Z_KEY": "REDSHIFTS",
        "SFPROB_KEY": "SFPROBS",
        "MAG_KEY": "i",
        "R50_KEY": "R50S",
        "N_KEY": "SERSIC_NS",
        "ELLIP_KEY": "ELLIPS",

        "SIZE": 100,
        "EXTRACTION_SIZE": 100,
        "BAND": "i",
        "ZEROPOINT": 27,
        "CORES": 4,
        "ARC_CONV": 0.2,
        "N_MODELS": 100,
        "SAME_BGS": True,
        "N_BGS": 50,

        "MASS_BINS": arange(10., 12., 0.4),
        "REDSHIFT_BINS": arange(0.1, 0.9, 0.2),
        "SFPROB_BINS": arange(0.0, 1, 0.5),
        "MASS_STEP": 0.4,
        "REDSHIFT_STEP": 0.2,
        "SFPROB_STEP": 0.5,

        "LINEAR": True,
        "LINEAR_STEP": 1,
        "USE_ALARM": True,
        "ALARM_TIME": 60,
        "MASK_PARAMS": [1, 2.0, 11],
        "BG_PARAMS": "ellipse",
        "BOX_SIZE": 40,
        "FILTER_SIZE": 3
    }

    return default_params


def dump_default_config_file(directory=""):
    """
    Dumps a default configuration file with all necessary parameters in the directory
    Args:
        directory: Directory to write file to. OPTIONAL.
    """
    lines = ["# Set verbosity printouts. VERBOSE: General printouts. TEST_VERBOSE: Additional printouts.",
             "VERBOSE             = True",
             "TEST_VERBOSE        = False",
             "",
             "# Directories and filenames -- Input and output",
             "# For SAVE CUTOUTS, options are 'none', 'mosaic', 'stitch', and 'fits'",
             "CATALOG             = cat.fits",
             "IMAGE_DIRECTORY     = images/",
             "PSF_FILENAME        = i_psfs.fits",
             "OUT_DIR             = out/",
             "SAVE_CUTOUTS        = none",
             "CUTOUT_FRACTION     = 0.2",

             "",
             "# Keys for masses, redshifts, and star-formation probability.",
             "MASS_KEY            = MASSES",
             "Z_KEY               = REDSHIFTS",
             "SFPROB_KEY          = SFPROBS",
             "",
             "# Keys for structural parameters. Magnitudes, half-light radii, Sersic index, ellipticity",
             "MAG_KEY             = i",
             "R50_KEY             = R50S",
             "N_KEY               = SERSIC_NS",
             "ELLIP_KEY           = ELLIPS",
             "",
             "# Cutout size, band, num-cores, arcsecs per pix, N models to gen, number of bgs to gen",
             "SIZE                = 100",
             "EXTRACTION_SIZE     = 100",
             "BAND                = i",
             "ZEROPOINT           = 27",
             "CORES               = 4",
             "ARC_CONV            = 0.2",
             "N_MODELS            = 100",
             "SAME_BGS            = True",
             "N_BGS               = 50",
             "",
             "# Bins to run through. (LOWER BOUND, UPPER BOUND, BIN WIDTH)",
             "# Note that the bins are defined by the LOWER BOUND to LOWER BOUND + BIN WIDTH",
             "MASS_BINS           = 10., 12., 0.4",
             "REDSHIFT_BINS       = 0.1, 0.9, 0.2",
             "SFPROB_BINS         = 0.0, 1, 0.5",
             "",
             "# Parameters for profile extraction.",
             "LINEAR              = True",
             "LINEAR_STEP         = 1",
             "USE_ALARM           = True",
             "ALARM_TIME          = 60",
             "",
             "# Parameters for Masking and background estimation ... NSIGMA, GAUSS_WIDTH, NPIX",
             "# Options for BG estimation: ellipse, circle, sigmaclip",
             "MASK_PARAMS         = 1, 2.0, 11",
             "BG_PARAMS           = ellipse",
             "BOX_SIZE            = 41",
             "FILTER_SIZE         = 3"
             ]

    with open(directory + "config.tbridge", "w+") as f:
        for n in lines:
            f.write(n + "\n")


def config_to_file(config, filename="config_out.txt"):
    """ Write a config dict to a file.

    Args:
        config (dict): The configuration parameters.
        filename (str): The location where the parameters will be written.
    """
    with open(filename, mode="w+") as f:
        for n in config:
            line = n + "\t" + str(config[n]) + "\t" + str(type(config[n])) + "\n"
            f.write(line)


def dump_default_config_file_koe(directory=""):
    """
    Dumps a default configuration file with all necessary parameters for using the KOE extraction pipeline.
    Args:
        directory: Directory to write file to. OPTIONAL.
    """
    lines = ["# Set verbosity printouts. VERBOSE: General printouts. TEST_MODE: Only extract 10 profiles as a check .",
             "VERBOSE             = True",
             "TEST_MODE        = False",
             "",
             "# Directories and filenames -- Input and output",
             "CATALOG             = cat.fits",
             "IMAGE_DIRECTORY     = images/",
             "OUT_DIR             = out/",
             "",
             "# Keys for RA, DEC, and Z",
             "RA                  = RA:",
             "DEC                 = DEC",
             "",
             "# Cutout size, band, num-cores, arcseconds per pix",
             "SIZE                = 100",
             "BAND                = i",
             "ZEROPOINT           = 27",
             "CORES               = 4",
             "ARC_CONV            = 0.2",
             "",
             "# Parameters for profile extraction.",
             "LINEAR              = True",
             "LINEAR_STEP         = 1",
             "USE_ALARM           = True",
             "ALARM_TIME          = 60",
             "",
             "# Parameters for Masking and background estimation ... NSIGMA, GAUSS_WIDTH, NPIX",
             "# Options for BG estimation: ellipse, circle, sigmaclip, 2D",
             "MASK_PARAMS         = 1, 2.0, 11",
             "BG_PARAMS           = sigmaclip",
             "BOX_SIZE            = 41",
             "FILTER_SIZE         = 6"
             ]

    with open(directory + "koe_config.tbridge", "w+") as f:
        for n in lines:
            f.write(n + "\n")


def load_config_file_koe(filename):
    """ Loads in a config file for KOE to run

    Args:
        filename: Filename (can be absolute or relative path, or a URL) to read config file from.
    Returns:
        dict: Configuration file as a dict object.
    """
    config_values = {}

    # First try to open things locally. If that doesn't work try it as a URL
    try:
        config_lines = open(filename, "r").readlines()
    except FileNotFoundError:
        try:
            r = urllib.request.urlopen(filename)
            config_lines = []
            for line in r:
                config_lines.append(line.decode("utf-8"))
        except:
            print("Failed to get any file")
            return None

    for line in config_lines:
        line = line.strip()
        if len(line) == 0 or line[0] == "#":
            continue
        splits = line.split("=")
        config_values[splits[0].strip()] = splits[1].strip()

    for n in config_values:
        value = config_values[n]
        if value.lower() == "true":
            config_values[n] = True
            continue
        elif value.lower() == "false":
            config_values[n] = False
            continue
    config_values["SIZE"] = int(config_values["SIZE"])
    config_values["CORES"] = int(config_values["CORES"])
    config_values["ZEROPOINT"] = int(config_values["ZEROPOINT"])
    config_values["ARC_CONV"] = float(config_values["ARC_CONV"])
    config_values["LINEAR_STEP"] = float(config_values["LINEAR_STEP"])
    config_values["ALARM_TIME"] = int(config_values["ALARM_TIME"])
    config_values["BOX_SIZE"] = int(config_values["BOX_SIZE"])
    config_values["FILTER_SIZE"] = int(config_values["FILTER_SIZE"])

    value_string = config_values["MASK_PARAMS"].split(",")
    config_values["MASK_PARAMS"] = [float(value_string[0]), float(value_string[1]), int(value_string[2])]

    return config_values
