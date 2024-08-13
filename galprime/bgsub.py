from photutils.background import Background2D
from astropy.stats import sigma_clipped_stats
from astropy.convolution import convolve, Tophat2DKernel, Gaussian2DKernel

from.masking import gen_mask

import numpy as np

def dilate_mask(mask, tophat_size):
    area = np.pi * tophat_size ** 2
    kernel = Tophat2DKernel(tophat_size)
    dilated_mask = convolve(mask, kernel) >= 1. / area
    return dilated_mask


def bgsub_source_mask(data, config, tophat_sizes=[3, 5, 7]):
    # Generate very aggressive source mask and dilate it thrice
    mask = gen_mask(data, config, omit_central=False)
    dilated_mask = np.copy(mask)
    for tophat_size in tophat_sizes:
       dilated_mask = dilate_mask(dilated_mask, tophat_size)

    return dilated_mask
    


def subtract_background(data, config={}, tophat_sizes=[3, 5, 7], plot_test=False):
    box_size = config.get("BGSUB", {}).get("BOX_SIZE", 42)
    filter_size = config.get("BGSUB", {}).get("FILTER_SIZE", 7)

    source_mask = bgsub_source_mask(data, config, tophat_sizes=tophat_sizes)
    

    # Generate background object
    bkg = Background2D(data, box_size=box_size, filter_size=(filter_size, filter_size), mask=source_mask)



def estimate_background_sigclip(cutout, config=None):
    nsigma = config.get("MASKING", {}).get("NSIGMA", 1)