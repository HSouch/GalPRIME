"""
All data-processing code methods for images, including profile extraction and masking procedures.
"""

from tqdm import tqdm
from astropy.table import Table
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma, sigma_clipped_stats
from numpy import copy, ndarray, floor, nan, sum, ndarray, nanmean, nanmedian, nanstd
from photutils import detect_threshold, detect_sources, deblend_sources, make_source_mask
from photutils.aperture import CircularAnnulus, EllipticalAnnulus


def mask_cutout(cutout, config=None, nsigma=1., gauss_width=2.0, npixels=11, omit_centre=True,
                clip_negatives=True):
    """ Masks a cutout. Users can specify parameters to adjust the severity of the mask. Default
        parameters strikes a decent balance.

    Args:
        cutout: Input cutout to mask.
        config Input config file
        nsigma: The brightness requirement for objects.
        gauss_width: The width of the gaussian kernel.
        npixels: The minimum number of pixels that an object must be comprised of to be considered a source.
        omit_centre: Set as true to leave the central object unmasked.
        clip_negatives: Remove negative pixels that are 3 sigma below the median BG level
                           (This is to get rid of artifacts).
    Returns:
        The masked cutout and the associated mask data in tuple format.
    """

    # If provided with a config file, set the parameters to what is given by the user.
    # Note that the user can just use the regular specified parameters by keeping config as None
    if config is not None:
        params = config["MASK_PARAMS"]
        nsigma = params[0]
        gauss_width = params[1]
        npixels = params[2]

    mask_data = {}
    c_x, c_y = int(floor(cutout.shape[0] / 2)), int(floor(cutout.shape[1] / 2))

    # Generate background mask and statistics

    bg_mask = make_source_mask(cutout, nsigma=2, npixels=3, dilate_size=7)


    bg_mean, bg_median, bg_std = sigma_clipped_stats(cutout, sigma=3.0, mask=bg_mask)

    # Generate source mask
    source_mask = generate_mask(cutout, nsigma=nsigma, gauss_width=gauss_width, npixels=npixels)
    source_mask = boolean_mask(source_mask, omit=[source_mask[c_x][c_y]] if omit_centre else None)
    n_masked = sum(source_mask)

    masked_cutout = copy(cutout)
    masked_cutout[source_mask] = nan

    if clip_negatives:
        limit = bg_median - (3 * bg_std)
        masked_cutout[masked_cutout <= limit] = nan

    mask_data["BG_MEAN"] = bg_mean
    mask_data["BG_MEDIAN"] = bg_median
    mask_data["BG_STD"] = bg_std
    mask_data["N_MASKED"] = n_masked
    mask_data["P_MASKED"] = n_masked / (cutout.shape[0] * cutout.shape[1])

    return masked_cutout, mask_data


def generate_mask(cutout, nsigma=1., gauss_width=2.0, npixels=5):
    """ Generates a given mask based on the input parameters


    Args:
        cutout:     input cutout (2D-Array)
        nsigma:     Number of sigma levels above the measured background pixels are required to be
                    at to be flagged as belonging to an object.
        gauss_width: The gaussian width of the smoothing kernal for the mask.
        npixels:    The number of connected pixels required in order for something to be considered an object in the
                    segmentation map.

    Returns:
        The segment array as a 2D array
    """

    sigma = gauss_width * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma).normalize()

    # Find threshold for cutout, and make segmentation map
    threshold = detect_threshold(cutout, nsigma=nsigma)
    segments = detect_sources(cutout, threshold, npixels=npixels, filter_kernel=kernel)

    # Attempt to de-blend. Return original segments upon failure.
    try:
        deb_segments = deblend_sources(cutout, segments, npixels=npixels, filter_kernel=kernel)
    except ImportError:
        print("Skimage not working!")
        deb_segments = segments
    except:
        # Don't do anything if it doesn't work
        deb_segments = segments

    segment_array = deb_segments.data

    return segment_array


def boolean_mask(mask, omit=None):
    """ Turns a given mask (photutils segment array) into a boolean array)

    Args:
        omit: If None, mask everything. Otherwise if an int or a list, do not mask those specified.
    """
    if omit is None:
        omit = []
    elif type(omit) == int:
        omit = [omit]

    bool_mask = ndarray(mask.shape, dtype="bool")

    bool_mask[:] = False
    bool_mask[mask > 0] = True
    for val in omit:
        bool_mask[mask == val] = False

    return bool_mask


def mask_cutouts(cutouts, config=None, method='standard', progress_bar=False):
    """ Mask a set of cutouts according to a certain method.

    Args:
        cutouts:
        config: User configuration file.
        method: Which method to use.
            standard : normal method. regular mask parameters. omits mask on central object
            no_central: regular mask parameters. masks central object
            background: background method: more severe mask parmeters. masks central object
        progress_bar: Use a TQDM progress bar.

    Returns:
        list of masked cutouts
    """
    masked_cutouts = []
    bg_means, bg_medians, bg_stds = [], [], []
    iterable = cutouts if not progress_bar else tqdm(cutouts)

    for cutout in iterable:
        try:
            if method == 'standard':
                masked, mask_data = mask_cutout(cutout, config=config, omit_centre=True)
            elif method == 'no_central':
                masked, mask_data = mask_cutout(cutout, config=config, omit_centre=False)
            elif method == 'background':
                masked, mask_data = mask_cutout(cutout, config=None,
                                                nsigma=0.5, gauss_width=2.0, npixels=5, omit_centre=False)
            else:
                continue

            masked_cutouts.append(masked)
            bg_means.append(mask_data["BG_MEAN"])
            bg_medians.append(mask_data["BG_MEDIAN"])
            bg_stds.append(mask_data["BG_STD"])
        except AttributeError:
            print("Cutout might be a NoneType")

    return masked_cutouts, [bg_means, bg_medians, bg_stds]
