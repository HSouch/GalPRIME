from configobj import ConfigObj



def dump_default_config_file():


    config = ConfigObj()
    config.filename = "default.gprime"

    config["FILES"] = {}
    config["FILES"]["CATALOGUE"] = "cat.fits"
    config["FILES"]["PSFS"] = "psfs.fits"
    config["FILES"]["BACKGROUNDS"] = "backgrounds.fits"

    config["KEYS"] = {}
    config["KEYS"]["MAG"] = "i"
    config["KEYS"]["R50"] = "R_GIM2D"
    config["KEYS"]["N"] = "SERSIC_N_GIM2D"
    config["KEYS"]["ELLIP"] = "ELL_GIM2D"

    config["BINS"] = {}
    config["Z_BEST"] = [0.1, 0.3, 0.5, 0.7, 0.9]
    config["MASS_MED"] = [10, 10.5, 11, 11.5]
    config["sfProb"] = [0, 0.5, 1.]

    config["DIRS"] = {}
    config["DIRS"]["OUTDIR"] = "gprime_out/"

    config.write()

