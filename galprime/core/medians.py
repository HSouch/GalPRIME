import numpy as np
from scipy.interpolate import interp1d


def common_sma(profiles, step=1, method='exact', dtype='isolist'):
    """ Put all profiles onto the same SMA samples
        Likely not super impactful if you are extracting profiles using the same configuration,
        but still worthwhile from a consistency standpoint.

    Args:
        profiles (array-like): List of profile objects or dictionaries containing SMA and intensity data.
        step (int, optional): Step size for SMA interpolation. Defaults to 1.
        method (str, optional): Method for aligning SMA samples. Defaults to 'exact'.
            Alternatives include 'interpolate' for linear interpolation of SMA samples.
        dtype (str, optional): Data type of profiles ('isolist' or 'table'). Defaults to 'isolist'.

    Raises:
        ValueError: If the specified method is not recognized.

    Returns:
        np.ndarray: Array of common SMA samples.
    """
    

    if dtype == 'isolist':
        profile_mins = [profiles[i].sma[0] for i in range(len(profiles))]
        profile_maxes = [profiles[i].sma[-1] for i in range(len(profiles))]
    elif dtype == 'table':
        profile_mins = [profiles[i]["sma"][0] for i in range(len(profiles))]
        profile_maxes = [profiles[i]['sma'][-1] for i in range(len(profiles))]

    if method == 'exact':
        min_sma = np.argmin(profile_mins)
        max_sma = np.argmax(profile_maxes)

        if dtype == 'isolist':
            sma_mins = profiles[min_sma].sma
            sma_maxes = profiles[max_sma].sma
        elif dtype == 'table':
            sma_mins = np.array(profiles[min_sma]["sma"])
            sma_maxes = np.array(profiles[max_sma]["sma"])
        
        sma_mins = sma_mins[sma_mins < np.nanmin(sma_maxes)]
        return np.concatenate((sma_mins, sma_maxes))

    elif method == 'interpolate':
        max_sma = np.nanmax(profile_maxes)
        return np.arange(1, max_sma, step)
    else:
        raise ValueError('Method not recognized')


def gen_profile_image(profiles, dtype='isolist', fill_value='extrapolate'):
    """ Generate a 2D image of profiles aligned on a common SMA grid.

    Args:
        profiles (array-like): List or array of profile objects or tables.
        dtype (str, optional): Type of input profiles. Defaults to 'isolist'.
        fill_value (str, optional): Fill value for interpolation. Defaults to 'extrapolate'.

    Returns:
        tuple: A tuple containing the common SMA samples and the 2D profile image.

    """

    sma = common_sma(profiles, dtype=dtype)

    prof_img = np.zeros((len(profiles), len(sma)))
    for i in range(len(profiles)):

        f = interp1d(profiles[i]['sma'], profiles[i]['intens'], 
                     kind='linear', fill_value=fill_value, bounds_error=False)
        prof_img[i] = f(sma)

    return sma, prof_img


def profile_median(profiles, dtype='isolist', fill_value='extrapolate'):
    """ Compute the median profile from a set of profiles.

    Args:
        profiles (array-like): List or array of profile objects or tables.
        dtype (str, optional): Type of input profiles. Defaults to 'isolist'.
        fill_value (str, optional): Fill value for interpolation. Defaults to 'extrapolate'.

    Returns:
        tuple: A tuple containing the median semi-major axis and the median profile.

    """
    med_sma, prof_img = gen_profile_image(profiles, dtype=dtype, fill_value=fill_value)

    noise_image = gen_profile_image(profiles, column='intens_err', dtype=dtype, fill_value=fill_value)[1]
    noise_sample = np.random.normal(loc=0.0, scale=noise_image)

    median = np.nanmedian(prof_img + noise_sample, axis=0)

    return med_sma, median


def bootstrap_median(profiles, samples=100, dtype='isolist', fill_value='extrapolate', 
                     add_individual_noise=True, stddev_factor=3.0):
    """Bootstrap the median and scatter for a set of profiles.

    Args:
        profiles (array-like): List or array of profile objects or tables.
        samples (int, optional): Number of bootstrap samples. Defaults to 100.
        dtype (str, optional): Type of input profiles. Defaults to 'isolist'.
        fill_value (str, optional): Fill value for interpolation. Defaults to 'extrapolate'.
        add_individual_noise (bool, optional): Whether to add individual noise. Defaults to True.
        stddev_factor (float, optional): Factor to scale the noise. Defaults to 3.0.

    Returns:
        tuple: A tuple containing the median semi-major axis, the bootstrapped median,
               the lower bounds, and the upper bounds.
    """
    med_sma, prof_img = gen_profile_image(profiles, dtype=dtype, fill_value=fill_value)
    noise_image = gen_profile_image(profiles, column='intens_err', dtype=dtype, fill_value=fill_value)[1]

    medians = np.zeros((samples, len(med_sma)))
    for i in range(samples):
        sample_indices = np.random.choice(len(prof_img), len(prof_img), replace=True)
        sample = prof_img[sample_indices] 
        if add_individual_noise:
            sample += np.random.normal(loc=0.0, scale=noise_image[sample_indices]) * stddev_factor
        medians[i] = np.nanmedian(sample, axis=0)
    
    medians = np.sort(medians, axis=0)

    bootstrapped_median = np.nanmedian(medians, axis=0)

    lower_index_1sig, upper_index_1sig = int(samples * 0.159), int(samples * 0.841)
    lower_index_2sig, upper_index_2sig = int(samples * 0.023), int(samples * 0.977)
    lower_index_3sig, upper_index_3sig = int(samples * 0.002), int(samples * 0.998)

    lower_1sig, upper_1sig = medians[lower_index_1sig], medians[upper_index_1sig]
    lower_2sig, upper_2sig = medians[lower_index_2sig], medians[upper_index_2sig]
    lower_3sig, upper_3sig = medians[lower_index_3sig], medians[upper_index_3sig]
    
    upper = np.vstack([upper_1sig, upper_2sig, upper_3sig])
    lower = np.vstack([lower_1sig, lower_2sig, lower_3sig])

    return med_sma, bootstrapped_median, lower, upper
