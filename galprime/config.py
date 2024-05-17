from configobj import ConfigObj



def dump_default_config_file(outname="default.gprime"):


    config = ConfigObj()
    config.filename = outname
    config["FILE_DIR"] = ""

    config["FILES"] = {}
    config["FILES"]["CATALOGUE"] = "cat.fits"
    config["FILES"]["PSFS"] = "psfs.fits"
    config["FILES"]["BACKGROUNDS"] = "backgrounds.fits"
    config["FILES"]["MAG_CATALOGUE"] = None

    config["KEYS"] = {}
    config["KEYS"]["RA"] = "RA_1"
    config["KEYS"]["DEC"] = "DEC_1"
    config["KEYS"]["MAG"] = "i"
    config["KEYS"]["R50"] = "R_GIM2D"
    config["KEYS"]["N"] = "SERSIC_N_GIM2D"
    config["KEYS"]["ELLIP"] = "ELL_GIM2D"

    config["PSF_INFO"] = {}
    config["PSF_INFO"]["PSF_RA"] = "RA"
    config["PSF_INFO"]["PSF_DEC"] = "DEC"

    config["BINS"] = {}
    config["BINS"]["Z_BEST"] = [0.1, 0.3, 0.5, 0.7, 0.9]
    config["BINS"]["MASS_MED"] = [10, 10.5, 11, 11.5]
    config["BINS"]["sfProb"] = [0, 0.5, 1.]

    config["DIRS"] = {}
    config["DIRS"]["OUTDIR"] = "gprime_out/"

    config["MASKING"] = {}
    config["MASKING"]["NSIGMA"] = 1
    config["MASKING"]["GAUSS_WIDTH"] = 2
    config["MASKING"]["NPIX"] = 5

    config["EXTRACTION"] = {}
    config["EXTRACTION"]["LINEAR"] = False
    config["EXTRACTION"]["STEP"] = 0.1
    config["EXTRACTION"]["NITER"] = 100

    config["BGSUB"] = {}
    

    config.write()


def read_config_file(filename):
    config = ConfigObj(filename)
    return config