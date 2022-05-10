""" Background estimation and subtraction module.

This module contains all methods available to estimate the backgrounds of cutouts,
and to subtract these backgrounds from cutouts or profiles as a means to correct them.

"""

import numpy as np
from astropy.table import Table
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.convolution import convolve, Tophat2DKernel

from photutils.aperture import CircularAnnulus, EllipticalAnnulus
from photutils import Background2D, MedianBackground, BkgZoomInterpolator, make_source_mask


def background_2D(img, mask, box_size, interp=None, filter_size=1,
                  exclude_percentile=90):
    ''' Run photutils background with SigmaClip and MedianBackground'''
    if interp is None:
        interp = BkgZoomInterpolator()
    return Background2D(img, box_size,
                        sigma_clip=SigmaClip(sigma=3.),
                        filter_size=filter_size,
                        bkg_estimator=MedianBackground(),
                        exclude_percentile=exclude_percentile,
                        mask=mask,
                        interpolator=interp)


class SourceMask:
    def __init__(self, img, nsigma=3., npixels=3, mask=None):
        ''' Helper for making & dilating a source mask.
             See Photutils docs for make_source_mask.'''
        self.img = img
        self.nsigma = nsigma
        self.npixels = npixels
        if mask is None:
            self.mask = np.zeros(self.img.shape, dtype=np.bool)
        else:
            self.mask = mask

    def single(self, filter_fwhm=3., tophat_size=5., mask=None):
        '''Mask on a single scale'''
        if mask is None:
            image = self.img
        else:
            image = self.img * (1 - mask)
        mask = make_source_mask(image, nsigma=self.nsigma,
                                npixels=self.npixels,
                                dilate_size=1, filter_fwhm=filter_fwhm)
        return dilate_mask(mask, tophat_size)

    def multiple(self, filter_fwhm=[3.], tophat_size=[3.], mask=None):
        """ Mask repeatedly on different scales """
        if mask is None:
            self.mask = np.zeros(self.img.shape, dtype=np.bool)
        for fwhm, tophat in zip(filter_fwhm, tophat_size):
            smask = self.single(filter_fwhm=fwhm, tophat_size=tophat)
            self.mask = self.mask | smask  # Or the masks at each iteration

        return self.mask

    def dilated(self, tophat_size=[3.], mask=None):
        """Mask using simple dilation"""
        if mask is None:
            self.mask = self.single()
        for tophat in tophat_size:
            smask = dilate_mask(self.mask, tophat)
            self.mask = self.mask | smask  # Or the masks at each iteration

        return self.mask


def dilate_mask(mask, tophat_size):
    """ Take a mask and make the masked regions bigger."""
    area = np.pi * tophat_size ** 2.
    kernel = Tophat2DKernel(tophat_size)
    dilated_mask = convolve(mask, kernel) >= 1. / area
    return dilated_mask


def estimate_background(cutout, config, model_params=None):
    """
    Estimate the background for a cutout using some of the various available methods based on what is set in the
    configuration file.

    OPTIONS:
        ellipse: uses elliptical annuli to estimate the cutout
        circle: uses a circular annuli
        sigclip (DEFAULT): measures the background using sigma-clipping.

    """
    if config["BG_PARAMS"] == "ellipse" and model_params is not None:
        bg_mean, bg_median, bg_std = estimate_bg_elliptical_annulus(cutout,
                                                                    ellipticity=model_params["ELLIP"],
                                                                    r_50=model_params["R50"],
                                                                    pa=model_params["PA"],
                                                                    width=50,
                                                                    factor=20)
    elif config["BG_PARAMS"] == 'circle':
        bg_mean, bg_median, bg_std = estimate_bg_annulus(cutout, dynamic=True, annulus_width=50)
    else:
        bg_mean, bg_median, bg_std = estimate_background_sigclip(cutout)

    return bg_mean, bg_median, bg_std


def estimate_background_sigclip(cutout, nsigma=2, npixels=3, dilate_size=7):
    """ Estimate the background mean, median, and standard deviation of a cutout using sigma-clipped-stats """

    bg_mask = make_source_mask(cutout, nsigma=nsigma, npixels=npixels, dilate_size=dilate_size)

    bg_mean, bg_median, bg_std = sigma_clipped_stats(cutout, sigma=3.0, mask=bg_mask)

    return bg_mean, bg_median, bg_std


def estimate_bg_annulus(cutout, annulus_radius=50, annulus_width=10, dynamic=True, **kwargs):
    """
    Measure the background of a cutout using a circular annulus
    """

    args = {"nsigma": 2, "npixels": 5, "dilate_size": 11, "sigclip_iters": 5}
    for arg in kwargs:
        args[arg] = kwargs[arg]

    if dynamic:
        annulus_radius = (cutout.shape[0] / 2) - annulus_width
    # First generate a background source mask for the input cutout
    bg_mask = make_source_mask(cutout,
                               nsigma=args["nsigma"],
                               npixels=args["npixels"],
                               dilate_size=args["dilate_size"],
                               sigclip_iters=args["sigclip_iters"])

    # Generate the annulus
    annulus = CircularAnnulus([cutout.shape[0] / 2, cutout.shape[1] / 2],
                              r_in=annulus_radius,
                              r_out=annulus_radius + annulus_width)

    # Generate a mask from the annulus and ensure it can be properly added to the cutout
    annulus_mask = annulus.to_mask(method='center')
    annulus_mask = annulus_mask.to_image(shape=cutout.shape)
    annulus_mask[annulus_mask == 0] = np.nan

    # Get the background pixels and remove any masked pixels
    bg_pixels = cutout * annulus_mask
    bg_pixels[bg_mask] = np.nan

    # Run our basic statistics and return them
    mean, median, std = np.nanmean(bg_pixels), np.nanmedian(bg_pixels), np.nanstd(bg_pixels)
    return mean, median, std


def estimate_bg_elliptical_annulus(cutout, ellipticity=0, r_50=50, pa=0, width=20, factor=10, return_mask=False):
    """
    Measure the background of a cutout using an elliptical annulus.
    """
    a_in = r_50 * factor
    b_in = a_in * (1 - ellipticity)

    a_out, b_out = a_in + width, b_in + width

    bg_mask = make_source_mask(cutout, nsigma=2, npixels=2, dilate_size=11)

    annulus = EllipticalAnnulus([cutout.shape[0] / 2, cutout.shape[1] / 2],
                                a_in=a_in, a_out=a_out,
                                b_in=b_in, b_out=b_out,
                                theta=pa)

    annulus_mask = annulus.to_mask(method='center')
    annulus_mask = annulus_mask.to_image(shape=cutout.shape)
    annulus_mask[annulus_mask == 0] = np.nan

    bg_pixels = cutout * annulus_mask
    bg_pixels[bg_mask] = np.nan

    mean, median, std = np.nanmean(bg_pixels), np.nanmedian(bg_pixels), np.nanstd(bg_pixels)

    if return_mask:
        return mean, median, std, bg_pixels
    else:
        return mean, median, std


def estimate_background_set(cutouts):
    """
    Estimates the background values for a set of cutouts.
    :param cutouts:
    :return:
    """
    bg_means, bg_medians, bg_stds = [], [], []

    for cutout in cutouts:
        bg_mean, bg_median, bg_std = estimate_background_sigclip(cutout)

        bg_means.append(bg_mean)
        bg_medians.append(bg_median)
        bg_stds.append(bg_std)

    return bg_means, bg_medians, bg_stds


def subtract_backgrounds(profile_set, background_array):
    """
    Generate an array of tables identical to the input except the respective backgrounds
    are subtracted from the intensity array for each table.
    :param profile_set: Set of profiles (in the photutiuls isolist format)
    :param background_array: Array of background values to subtract from each profile table.
    :return: List of profiles of length len(profile_set)
    """
    bg_subtracted_tables = []

    for i in range(0, len(profile_set)):
        this_table = profile_set[i]
        isotable_localsub = Table()

        for col in this_table.colnames:
            isotable_localsub[col] = np.copy(this_table[col])
        isotable_localsub["intens"] = ((isotable_localsub["intens"]) - background_array[i])

        bg_subtracted_tables.append(isotable_localsub)

    return bg_subtracted_tables
