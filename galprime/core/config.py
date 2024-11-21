from configobj import ConfigObj, validate
import numpy as np

def default_config():

    config = ConfigObj()
    config.filename = None
    config["FILE_DIR"] = ""
    config["NCORES"] = 2

    config["FILES"] = {}
    config["FILES"]["CATALOGUE"] = "cat.fits"
    config["FILES"]["PSFS"] = "psfs.fits"
    config["FILES"]["BACKGROUNDS"] = "backgrounds.fits"
    config["FILES"]["MAG_CATALOGUE"] = None

    config["DIRS"] = {}
    config["DIRS"]["OUTDIR"] = "gprime_out/"

    # Keys to pull information from the catalogue
    config["KEYS"] = {}
    config["KEYS"]["MAG"] = "i"
    config["KEYS"]["REFF"] = "R_GIM2D"
    config["KEYS"]["N"] = "SERSIC_N_GIM2D"
    config["KEYS"]["ELLIP"] = "ELL_GIM2D"
    config["KEYS"]["RA"] = "RA_1"
    config["KEYS"]["DEC"] = "DEC_1"

    config["PSFS"] = {}
    config["PSFS"]["PSF_RA"] = "RA"
    config["PSFS"]["PSF_DEC"] = "DEC"

    config["BINS"] = {}
    config["BINS"]["Z_BEST"] = [0.1, 0.3, 0.5, 0.7, 0.9]
    config["BINS"]["MASS_MED"] = [10, 10.5, 11, 11.5]
    config["BINS"]["sfProb"] = [0, 0.5, 1.]

    config["MODEL"] = {}
    config["MODEL"]["N_MODELS"] = 50
    config["MODEL"]["SIZE"] = 451
    config["MODEL"]["ARCCONV"] = 0.168
    config["MODEL"]["ZPM"] = 27.0
    config["MODEL"]["REFF_UNIT"] = "pixel"

    config["MASKING"] = {}
    config["MASKING"]["NSIGMA"] = 1
    config["MASKING"]["GAUSS_WIDTH"] = 2
    config["MASKING"]["NPIX"] = 5
    config["MASKING"]["BG_BOXSIZE"] = 50

    config["EXTRACTION"] = {}
    config["EXTRACTION"]["LINEAR"] = False
    config["EXTRACTION"]["STEP"] = 0.1
    config["EXTRACTION"]["MAXIT"] = 100
    config["EXTRACTION"]["FIX_CENTER"] = False
    config["EXTRACTION"]["MINSMA"] = 1
    config["EXTRACTION"]["MAXSMA"] = config["MODEL"]["SIZE"] // 2
    config["EXTRACTION"]["CONVER"] = 0.05


    config["BGSUB"] = {}
    config["BGSUB"]["BOX_SIZE"] = 42
    config["BGSUB"]["FILTER_SIZE"] = 7
    config["BGSUB"]["NSIGMA"] = 3
    config["BGSUB"]["NPIXELS"] = 10

    return config


def galprime_configspec():
    cspec = ConfigObj()

    cspec["FILE_DIR"] = ""
    cspec["OUTDIR"] = "string(default='gprime_out/')"
    cspec["NCORES"] = "integer(default=1)"
    cspec["TIME_LIMIT"] = "integer(default=10)"

    cspec["FILES"] = {}
    cspec["FILES"]["CATALOGUE"] = "string(default='cat.fits')"
    cspec["FILES"]["PSFS"] = "string(default='psfs.fits')"
    cspec["FILES"]["BACKGROUNDS"] = "string(default='backgrounds.fits')"
    cspec["FILES"]["MAG_CATALOGUE"] = "string(default=None)"

    cspec["KEYS"] = {}
    cspec["KEYS"]["RA"] = "string(default='RA_1')"
    cspec["KEYS"]["DEC"] = "string(default='DEC_1')"


    cspec["BINS"] = {}

    cspec["MODEL"] = {}
    cspec["MODEL"]["MODEL_TYPE"] = "integer(default=1)"
    cspec["MODEL"]["N_MODELS"] = "integer(default=50)"
    cspec["MODEL"]["SIZE"] = "integer(default=151)"
    cspec["MODEL"]["ARCCONV"] = "float(default=0.168)"
    cspec["MODEL"]["ZPM"] = "float(default=27.0)"
    cspec["MODEL"]["REFF_UNIT"] = "string(default='pixel')"
    
    cspec["DIRS"] = {}
    cspec["DIRS"]["OUTDIR"] = "string(default='gprime_out/')"

    cspec["MASKING"] = {}
    cspec["MASKING"]["NSIGMA"] = "float(default=1)"
    cspec["MASKING"]["GAUSS_WIDTH"] = "float(default=2)"
    cspec["MASKING"]["NPIX"] = "integer(default=5)"
    cspec["MASKING"]["BG_BOXSIZE"] = "integer(default=50)"

    cspec["EXTRACTION"] = {}
    cspec["EXTRACTION"]["LINEAR"] = "boolean(default=False)"
    cspec["EXTRACTION"]["STEP"] = "float(default=0.1)"
    cspec["EXTRACTION"]["MAXIT"] = "integer(default=100)"
    cspec["EXTRACTION"]["FIX_CENTER"] = "boolean(default=False)"
    cspec["EXTRACTION"]["MAXRIT"] = "integer(default=50)"
    cspec["EXTRACTION"]["MINSMA"] = "integer(default=1)"
    cspec["EXTRACTION"]["CONVER"] = "float(default=0.05)"
    cspec["EXTRACTION"]["INTEGRMODE"] = "string(default='bilinear')"

    cspec["BGSUB"] = {}
    cspec["BGSUB"]["BOX_SIZE"] = "integer(default=42)"
    cspec["BGSUB"]["FILTER_SIZE"] = "integer(default=7)"
    cspec["BGSUB"]["NSIGMA"] = "float(default=3)"
    cspec["BGSUB"]["NPIXELS"] = "integer(default=10)"

    return cspec

def dump_default_config_file(outname="default.gprime"):

    config = default_config()
    config.filename = outname
    config.indent_type = '    '
    config.write()


def read_config_file(filename):
    config = ConfigObj(filename, interpolation=True, configspec=galprime_configspec(), indent_type = '    ')
    
    vtor = validate.Validator()
    test = config.validate(vtor)

    for key in config["BINS"]:
        config["BINS"][key] = np.array(config["BINS"][key], dtype=float)

    return config
