import numpy as np
from astropy.cosmology import WMAP9 as cosmo
    
def kpc_per_arcsecond(z):
    return 1 / cosmo.arcsec_per_kpc_proper(z).value
