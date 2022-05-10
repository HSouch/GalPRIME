"""
A command-line wrapper for GalPRIME with full arguments included.
"""

import galprime
import argparse
import warnings
from astropy.table import Table


print("GalPRIME, running version", galprime.__version__)

# ########################################################################
# ADD USER PARAMETERS  ###################################################
# ########################################################################

parser = argparse.ArgumentParser(description="Command-line wrapper for TBriDGE. " +
                                             "Please see documentation for more info.")

parser.add_argument("--CONFIG", type=str, help="Config fie filename.")

parser.add_argument("--CATALOG", type=str, help="Input Catalog")
parser.add_argument("--IMAGE_DIRECTORY", type=str, help="Image Directory (Top Level)")
parser.add_argument("--PSF_FILENAME", type=str, help="PSF Filename")
parser.add_argument("--OUT_DIR", type=str, help="Output directory")

parser.add_argument("--MASS_KEY", type=str, help="Stellar mass column name.")
parser.add_argument("--Z_KEY", type=str, help="Redshift column name.")
parser.add_argument("--SFPROB_KEY", type=str, help="Star formation probability column name.")
parser.add_argument("--MAG_KEY", type=str, help="Magnitude column name.")
parser.add_argument("--R50_KEY", type=str, help="Half light radius column name.")
parser.add_argument("--N_KEY", type=str, help="Sersic index column name.")
parser.add_argument("--ELLIP_KEY", type=str, help="Ellipticity column name.")

parser.add_argument("--BAND", type=str, help="Image band.")
parser.add_argument("--SIZE", type=int, help="Cutout Width.")
parser.add_argument("--CORES", type=int, help="Number of cores to run sims on.")
parser.add_argument("--ARC_CONV", type=float, help="Arcseconds per pixel.")
parser.add_argument("--N_MODELS", type=int, help="Number of models to simulate")

parser.add_argument("--BGS", type=str, help="Provided bg filename.")
parser.add_argument("--LINEAR", type=str, help="Extract profiles with a linear step in SMA.")
parser.add_argument("-v", action="store_true", help="Run in Verbose mode")
args = parser.parse_args()
args_dict = args.__dict__


config = tbridge.default_config_params() if args.CONFIG is None else tbridge.load_config_file(args.CONFIG)

for n in args_dict:
    if args_dict[n] is not None:
        config[n] = args_dict[n]

if args.BGS is None:
    provided_bgs, provided_bg_type = None, "None"
else:
    provided_bg_type = "None"
    try:
        provided_bgs, provided_bg_type = tbridge.load_cutouts(args.BGS), "Cutouts"
    except:
        pass
    try:
        if len(provided_bgs) == 0:
            provided_bgs, provided_bg_type = Table.read(args.BGS), "Table"
    except:
        pass

# ########################################################################
# PROCESSING #############################################################
# ########################################################################

magtable = Table.read("/home/harrison/Desktop/Research/COSMOS_V9_CLEAN.fits")["i", "ZPHOT", "MASS_MED"]

warnings.filterwarnings("ignore")

tbridge.pipeline(config, max_bins=2, mag_table=magtable, provided_bgs=provided_bgs, multiprocess_level='obj')
