import numpy as np
from scipy.interpolate import interp1d


def common_sma_old(profiles, step=1, method='exact', dtype='isolist'):

    if dtype == 'isolist':
        profile_mins = [profiles[i].sma[0] for i in range(len(profiles))]
        profile_maxes = [profiles[i].sma[-1] for i in range(len(profiles))]
    elif dtype == 'table':
        profile_mins = [profiles[i]["sma"][0] for i in range(len(profiles))]
        profile_maxes = [profiles[i]['sma'][-1] for i in range(len(profiles))]

    if method == 'exact':
        min_sma = np.nanmax(profile_mins)
        max_sma = np.argmax(profile_maxes)
        if dtype == 'isolist':
            return profiles[max_sma].sma
        elif dtype == 'table':
            return profiles[max_sma]["sma"]
    
    elif method == 'interpolate':
        max_sma = np.nanmax(profile_maxes)
        return np.arange(1, max_sma, step)
    else:
        raise ValueError('Method not recognized')


def common_sma(profiles, step=1, method='exact', dtype='isolist'):

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
    sma = common_sma(profiles, dtype=dtype)

    prof_img = np.zeros((len(profiles), len(sma)))
    for i in range(len(profiles)):
        f = interp1d(profiles[i]['sma'], profiles[i]['intens'], kind='linear', fill_value=fill_value, bounds_error=False)
        prof_img[i] = f(sma)

    return sma, prof_img


def profile_median(profiles, dtype='isolist', fill_value='extrapolate'):
    med_sma, prof_img = gen_profile_image(profiles, dtype=dtype, fill_value=fill_value)

    median = np.nanmedian(prof_img, axis=0)

    return med_sma, median


def bootstrap_median(profiles, samples=100, dtype='isolist', fill_value='extrapolate'):
    med_sma, prof_img = gen_profile_image(profiles, dtype=dtype, fill_value=fill_value)

    medians = np.zeros((samples, len(med_sma)))
    for i in range(samples):
        sample_indices = np.random.choice(len(prof_img), len(prof_img), replace=True)
        sample = prof_img[sample_indices]
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
