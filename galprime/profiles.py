""" Surface brightness profile methods.

This module contains all methods pertaining directly to surface brightness
profiles measured using isophote fitting.

"""


from numpy import trapz, copy, isnan, isinf, log10, arange, min, max, gradient, median, std
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from scipy.stats import norm


def get_half_light_radius(x, y):
    """ Obtains the index of the ellipse whose semimajor axis is at the half-light radius.

    Args:
        x (arr): Array of radii
        y (arr): Array of intensities

    Returns:
        best_index: the index containing the array value closest to the half-light radius
    """
    # Estimates the half-light radius for a given profile.

    half_light = trapz(y, x) / 2
    best_sec = 0     # Light contained in given section
    best_index = 0
    for index in range(0, len(x)):
        sec = trapz(y[:index], x[:index])
        if abs(sec - half_light) < abs(best_sec - half_light):
            best_sec, best_index = sec, index
    return best_index


def sb_to_i(sb, m_0=27, a_pix=1):
    """ Convert a given surface brightness to an intensity.

    Args:
        sb (float): Surface brightness
        m_0 (float): Zero-point magnitude
        a_pix (float): Area of a pixel (in arc seconds squared)

    Returns:
        Intensity value for the given surface brightness.
    """
    exponential = (sb - m_0) / -2.5
    return a_pix * (10 ** exponential)


def adjust_profile(x, y, shift=0.00, clean=False, log_x=False):
    """ Adjusts a profile to a common plotting format (surface brightness and optional shift vs regular pixels).

    Args:
        x (arr): Array of radii
        y (arr): Array of intensities
        shift (float): Introduce a shift to the intensity values (sometimes helpful if many negative numbers present)
        clean (bool): If true, clean nan and inf values from array (optional)
        log_x (bool): If True, log the radial intensities.
    Returns:
        x_new, y_new : The adjusted radial and intensity arrays
    """

    x_new, y_new = copy(x), copy(y)
    if log_x:
        x_new = log10(x_new)

    y_new = -2.5 * log10((y + shift) / (0.168 ** 2)) + 27

    if clean:
        x_clean, y_clean = [], []
        for n in range(0, len(x_new)):
            if isnan(y_new[n]) or isinf(y_new[n]):
                continue
            else:
                x_clean.append(x_new[n])
                y_clean.append(y_new[n])
        profile_adj = interp1d(x_clean, y_clean, bounds_error=False, fill_value=max(y_clean))

        return x_new, profile_adj(x_new)

    else:
        return x_new, y_new


def overlap(sma, low_1, high_1, low_2, high_2, threshold=1):
    """ Tests for non_overlap between array sets 1 and 2

    Args:
        sma (arr): Test parameters, preferably a numpy arange object
        low_1 (scipy.interpolate.interp1d): interp object for the lower bounds of set 1
        high_1 (scipy.interpolate.interp1d): interp object for the higher bounds of set 1
        low_2 (scipy.interpolate.interp1d): interp object for the lower bounds of set 2
        high_2 (scipy.interpolate.interp1d): interp object for the higher bounds of set 2
        threshold (int): How many consecutive times you see non-overlap between the two profiles.

    Returns:
        The index where the sufficient overlap is broken.
    """
    count, i = 0, 0
    for i in sma:
        overlapping = not ((high_2(i) < low_1(i)) or (low_2(i) > high_1(i)))
        if overlapping:
            count = 0
        if count == threshold:
            break
        count += 1

    return i


def smooth(arr, window=9, poly=1, iterations=2):
    """ Smooth an array with some number of Savitzky-Golay filters

    Args:
        arr (arr): Array to be be smoothed
        window (int): Filter window length
        poly (int): Order of polynomail
        iterations (int): The number of times to smooth the array

    Returns:
        Smoothed array of same length as arr

    """
    arr = copy(arr)
    for i in range(iterations):
        arr = savgol_filter(arr, window_length=window, polyorder=poly)

    return arr


def profile_quality(profile, bg_mean, bg_std, alpha=0.05, seeing=0, intens_cut=0.01):
    """ Determines the quality of a given profile compared to measured background levels.

    Args:
        profile (isotable): The input profile (still in the photutils isotable format).
        bg_mean (float): The background mean
        bg_std (float): The background standard deviation
        alpha (float): The confidence interval to define a "good" measurement
        seeing (float): The inner seeing distance within which we will not consider measurements to count.
        intens_cut (float): A minimum intensity beyond which we will stop counting measurements, inputted
                            as a fraction of the central profile intensity.

    Returns:
        good (int): The number of good measurements in the profile.
        bad (int): The number of bad measurements in the profile.
    """
    # Retrieve the background distribution as a norm dist. centered on bg_mean with std=bg_std
    bg_distribution = norm(loc=bg_mean, scale=bg_std)

    good, bad, samples = 0, 0, len(profile)
    sma_table, intens_table = profile["sma"], profile["intens"]

    max_intens = max(intens_table)
    intens_cut = max_intens * intens_cut

    # Iterate over all intensity measurements
    # If the cumulative dist. function at that point is greater than 1-alpha, then we
    # have a good measurement. Otherwise, we have a bad measurement.

    for i in range(samples):
        # Don't count measurements that are made within the seeing HWHM
        if sma_table[i] < seeing:
            continue

        intens = intens_table[i]
        if intens <= intens_cut:
            break
        if bg_distribution.sf(intens) < alpha:
            good += 1
        else:
            bad += 1

    return good, bad


def profile_asymptote(x, y, step=0.5, smooth=True, min_cut=10, threshold=1e-4):
    """ Determine the profile asymptotic value from a profile's x and y parameters.

    This method smooths out the valid data with a Savitzky-Golay filter then determines statistics on the
    values where the gradient is sufficiently low (below the threshold).

    Args:
        x (arr): Array of radii
        y (arr): Array of intensities
        step (float): When resampling the arrays, this is the step between samples.
        smooth (bool): If True, smoothing will be applied. (Recommended).
        min_cut (int): The number of measurements at the start of the profile to ignore. This is to remove any issues
                       with low gradients at the start of the profile.
        threshold (int): The gradient threshold below which measurements will be counted as sufficiently asymptotic.

    Returns:
        med (float): The median asymptote value
        std (float): The standard deviation of the asymptote measurements
        range (tuple): The range of locations (radial distances) where measurements are taken

    """
    interp = interp1d(x, y, fill_value=min(y), bounds_error=False)
    x_new = arange(min(x), max(x) + step, step)
    y_new = interp(x_new)

    # Just remove any nan values (we'll assume they their removal will not cause major problems)
    x_new, y_new = x_new[min_cut:], y_new[min_cut:]
    x_new, y_new = x_new[~isnan(y_new)], y_new[~isnan(y_new)]

    if smooth:
        # Twice smooth the data
        smoothed = savgol_filter(y_new, window_length=9, polyorder=2)
        smoothed = savgol_filter(smoothed, window_length=9, polyorder=2)
        y_new = smoothed

    grad = abs(gradient(y_new, x_new))
    y_good = y_new[grad < threshold]
    x_good = x_new[grad < threshold]

    if len(y_good) == 0:
        print("No values available to calculate gradient.")
        return None

    # Return the asymptote, the spread of values in asymptotic numbers, and the range of "good" measurements
    return median(y_good), std(y_good), (min(x_good), max(x_good))
