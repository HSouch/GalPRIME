import numpy as np
from astropy.convolution import convolve
from astropy.stats import sigma_clipped_stats
from photutils import segmentation


def gen_mask(data, config=None, omit=[], omit_central=True):
    """
    Generate a mask based on the input data.

    Parameters:
    - data: numpy.ndarray
        The input data array.
    - config: dict, optional
        Configuration parameters for masking. If not provided, default values will be used.
    - omit: list, optional
       List of values to be omitted from the mask. Mainly to be used for manual masking.

    Returns:
    - mask: numpy.ndarray
        The generated mask array.

    """
    metadata = {}
    data_masked = np.copy(data)
    if config is not None:
        nsigma = float(config["MASKING"]["NSIGMA"])
        gauss_width = float(config["MASKING"]["GAUSS_WIDTH"])
        npix = int(config["MASKING"]["NPIX"])
    else:
        nsigma, gauss_width, npix = 2., 5., 5
        bg_boxsize = 50
    
    bg_mean, bg_median, bg_std = sigma_clipped_stats(data)
    threshold = bg_median + nsigma * bg_std 

    kernel = segmentation.make_2dgaussian_kernel(gauss_width, size=5)  # FWHM = 3.0
    convolved_data = convolve(data_masked, kernel)

    segment_map = segmentation.detect_sources(convolved_data, threshold, npixels=10)
    segm_deblend = segmentation.deblend_sources(convolved_data, segment_map,
                                                npixels=npix, nlevels=32, contrast=0.001,
                                                progress_bar=False)

    central_pix = (data.shape[0] // 2, data.shape[1] // 2)
    central_value = segm_deblend.data[central_pix]
    
    mask = np.zeros_like(data, dtype=bool)
    if omit_central:
        mask[np.logical_and(segm_deblend.data != central_value, segm_deblend.data != 0)] = True
    else:
        mask[segm_deblend.data != 0] = True
    for n in omit:
        mask[segm_deblend.data == n] = False

    return mask, metadata


def mask_image(data, config):
    """
    Masks the input data array based on the provided parameters.

    Parameters:
    - data (numpy.ndarray): The input data array to be masked.
    - config (dict): The configuration parameters for generating the mask.

    Returns:
    - numpy.ndarray: The masked data array with NaN values in the masked regions.
    """
    mask, metadata = gen_mask(data, config)
    data_masked = np.copy(data)
    data_masked[mask] = np.nan
    return data_masked

