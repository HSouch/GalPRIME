import sys
import numpy as np
from .isophote_l import Ellipse, EllipseGeometry

from photutils.morphology import data_properties


def isophote_fitting(data, config={}, centre_method='standard'):
    """ Wrapper for photutils.isophote methods

    Generates a table of results from isophote fitting analysis. This uses photutils Isophote procedure, which is
    effectively IRAF's Ellipse() method.
    Iterates over many possible input ellipses to force a higher success rate.

    Args:
        data: The input cutout
        config: The configuration parameters, which contains details like STEP and LINEAR for profile extraction.
        centre_method: Which method to use to determine where to place the centre of the first ellipse.
            'standard': place at the exact centre of the cutout
            'max': place at the maximum pixel in the cutout

    Returns:
        The table of results, or an empty list if not fitted successfully.

    """
    # Set-up failsafe in case of strange infinte loops in photutils
    # warnings.filterwarnings("error")

    fail_count, max_fails = 0, 1000
    linear = config.get("EXTRACTION", {}).get("LINEAR", False)
    step = config.get("EXTRACTION", {}).get("STEP", 0.1)

    # Get centre of image and cutout halfwidth
    if centre_method == 'standard':
        centre = (data.shape[0]/2, data.shape[1]/2)
    elif centre_method == 'max':
        centre = np.unravel_index(np.argmax(data), data.shape)
    else:
        centre = (data.shape[0] / 2, data.shape[1] / 2)

    cutout_halfwidth = max((np.ceil(data.shape[0] / 2), np.ceil(data.shape[1] / 2)))

    fitting_list = []

    # First, try obtaining morphological properties from the data and fit using that starting ellipse
    try:
        morph_cat = data_properties(data)
        r = 2.0
        pos = (morph_cat.xcentroid, morph_cat.ycentroid)
        a = morph_cat.semimajor_sigma.value * r
        b = morph_cat.semiminor_sigma.value * r
        theta = morph_cat.orientation.value

        geometry = EllipseGeometry(pos[0], pos[1], sma=a, eps=(1 - (b / a)), pa=theta)
        flux = Ellipse(data, geometry)
        fitting_list = flux.fit_image(maxit=100, maxsma=cutout_halfwidth, step=step, linear=linear,
                                      maxrit=cutout_halfwidth / 3)
        if len(fitting_list) > 0:
            return fitting_list

    except KeyboardInterrupt:
        sys.exit(1)
    except (RuntimeWarning, RuntimeError, ValueError, OverflowError, IndexError):
        fail_count += 1
        if fail_count >= max_fails:
            return []

    # If that fails, test a parameter space of starting ellipses
    try:
        for angle in range(0, 180, 45):
            for sma in range(2, 26, 5):
                for eps in (0.3, 0.5, 0.9):
                    geometry = EllipseGeometry(float(centre[0]), float(centre[1]), eps=eps,
                                               sma=sma, pa=angle * np.pi / 180.)
                    flux = Ellipse(data, geometry)
                    fitting_list = flux.fit_image(maxsma=cutout_halfwidth, step=step, linear=linear,
                                                  maxrit=cutout_halfwidth / 3)
                    if len(fitting_list) > 0:
                        return fitting_list

    except KeyboardInterrupt:
        sys.exit(1)
    except (RuntimeWarning, RuntimeError, ValueError, OverflowError, IndexError):

        # print("RuntimeError or ValueError")
        fail_count += 1
        if fail_count >= max_fails:
            return []

    return fitting_list


class IsophoteFitter:
    def __init__(self, data, config={}, mask=None, centre=None, **kwargs):
        self.data = data
        self.config = config
        self.mask = mask
        self.centre = centre

        if self.mask is not None:
            self.data = np.ma.masked_array(self.data, mask=self.mask)

        self.fail_count, self.max_fails = 0, 100

        self.linear = config.get("MODEL", {}).get("LINEAR", False)
        self.step = config.get("MODEL", {}).get("STEP", 0.1)

        if self.centre == None:
            self.centre = (self.data.shape[0] / 2, self.data.shape[1] / 2)
    

    def _fit_single(self, geometry, **kwargs):
        try:
            flux = Ellipse(self.data, geometry)
            fitting_list = flux.fit_image()
        except KeyboardInterrupt:
            sys.exit(1)
        except (RuntimeError, ValueError, OverflowError, IndexError):
            self.fail_count += 1
            if self.fail_count > self.max_fails:
                raise RuntimeError("Too many failed fit attempts")

    def fit(self):
        # First try fit using geometry guess
        try:
            morph_cat = data_properties(self.data, mask=self.mask, background=None)
            r = 2.0
            pos = (morph_cat.xcentroid, morph_cat.ycentroid)
            a = morph_cat.semimajor_sigma.value * r
            b = morph_cat.semiminor_sigma.value * r
            theta = morph_cat.orientation.value

            geometry = EllipseGeometry(x0=pos[0], y0=pos[1], sma=a, 
                                       eps=1 - b / a, pa=morph_cat.orientation.value)
            
            isolist = self._fit_single(geometry)
            
        except Exception as e:
            self.fail_count += 1

        # If that fails, try a brute force search
        try:
            for angle in range(0, 180, 45):
                for sma in range(2, 26, 5):
                    for eps in (0.3, 0.5, 0.9):
                        geometry = EllipseGeometry(x0=self.centre[0], y0=self.centre[1], 
                                                   sma=sma, eps=eps, pa=angle)
                        try:
                            isolist = self._fit_single(geometry)
                            
                            if len(isolist) > 0:
                                break
                            
                        except RuntimeError as e:
                            return []
                        except Exception as e:
                            continue

        except Exception as e:
            self.fail_count += 1
        