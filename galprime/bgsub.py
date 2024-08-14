from photutils.background import Background2D, MedianBackground, BkgZoomInterpolator
from photutils.segmentation import detect_sources

from astropy.stats import sigma_clipped_stats, SigmaClip, gaussian_fwhm_to_sigma
from astropy.convolution import convolve, Tophat2DKernel, Gaussian2DKernel

from.masking import gen_mask

import numpy as np

def dilate_mask(mask, tophat_size):
    area = np.pi * tophat_size ** 2
    kernel = Tophat2DKernel(tophat_size)
    dilated_mask = convolve(mask, kernel) >= 1. / area
    return dilated_mask


# def bgsub_source_mask(data, config, tophat_sizes=[1,3]):
def bgsub_source_mask(data, config, mask=None, filter_fwhm=None,
                     filter_size=3, kernel=None, sigclip_sigma=3.0,
                     sigclip_iters=5, dilate_size=11):
    """
        Source mask generation (from photutils 1.4) 
    """ 
    from scipy import ndimage

    nsigma = config.get("BGSUB", {}).get("NSIGMA", 1)
    npixels = config.get("BGSUB", {}).get("NPIXELS", 10)
    threshold = nsigma * sigma_clipped_stats(data)[2]

    if kernel is None and filter_fwhm is not None:
        kernel_sigma = filter_fwhm * gaussian_fwhm_to_sigma
        kernel = Gaussian2DKernel(kernel_sigma, x_size=filter_size,
                                  y_size=filter_size)
    if kernel is not None:
        kernel.normalize()

    segm = detect_sources(data, threshold, npixels)
    if segm is None:
        return np.zeros(data.shape, dtype=bool)

    selem = np.ones((dilate_size, dilate_size))
    return ndimage.binary_dilation(segm.data.astype(bool), selem)
    


def estimate_background_2D(data, config={}, tophat_sizes=[3, 5, 7],
                           exclude_percentile=90, interp=None,
                           plot_test=False):

    data = np.copy(data)
    data[np.isinf(data)] = np.nan


    box_size = config.get("BGSUB", {}).get("BOX_SIZE", 42)
    filter_size = config.get("BGSUB", {}).get("FILTER_SIZE", 7)
    nsigma = config.get("BGSUB", {}).get("NSIGMA", 3)
    npixels = config.get("BGSUB", {}).get("NPIXELS", 10)
    
    source_mask = bgsub_source_mask(data, config)
    
    if interp is None:
        interp = BkgZoomInterpolator()

    # Generate background object
    bkg = Background2D(data, box_size,
                        sigma_clip=SigmaClip(sigma=3.),
                        filter_size=filter_size,
                        bkg_estimator=MedianBackground(),
                        exclude_percentile=exclude_percentile,
                        mask=source_mask,
                        interpolator=interp)

    return source_mask, bkg


def estimate_background_sigclip(cutout, config=None):
    nsigma = config.get("MASKING", {}).get("NSIGMA", 1)
