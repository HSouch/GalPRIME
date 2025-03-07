import numpy as np
from scipy.interpolate import interp1d


def common_sma(profiles, step=1, method='exact', dtype='isolist'):

    if dtype == 'isolist':
        profile_set = [profiles[i].sma[-1] for i in range(len(profiles))]
    elif dtype == 'table':
        profile_set = [profiles[i]['sma'][-1] for i in range(len(profiles))]

    if method == 'exact':
        max_sma = np.argmax(profile_set)
        if dtype == 'isolist':
            return profiles[max_sma].sma
        elif dtype == 'table':
            return profiles[max_sma]["sma"]
    
    elif method == 'interpolate':
        max_sma = np.nanmax(profile_set)
        return np.arange(1, max_sma, step)
    else:
        raise ValueError('Method not recognized')
    

def gen_profile_image(profiles, dtype='isolist'):
    sma = common_sma(profiles, dtype=dtype)

    prof_img = np.zeros((len(profiles), len(sma)))
    for i in range(len(profiles)):
        f = interp1d(profiles[i]['sma'], profiles[i]['intens'], kind='linear', fill_value='extrapolate')
        prof_img[i] = f(sma)

    return sma, prof_img


def profile_median(profiles, dtype='isolist'):
    med_sma, prof_img = gen_profile_image(profiles, dtype=dtype)

    median = np.median(prof_img, axis=0)

    return med_sma, median


def bootstrap_median(profiles, samples=100, dtype='isolist'):
    med_sma, prof_img = gen_profile_image(profiles, dtype=dtype)

    medians = np.zeros((samples, len(med_sma)))
    for i in range(samples):
        sample_indices = np.random.choice(len(prof_img), len(prof_img), replace=True)
        sample = prof_img[sample_indices]
        medians[i] = np.median(sample, axis=0)
    
    medians = np.sort(medians, axis=0)

    bootstrapped_median = np.median(medians, axis=0)

    lower_index_1sig, upper_index_1sig = int(samples * 0.159), int(samples * 0.841)
    lower_index_2sig, upper_index_2sig = int(samples * 0.023), int(samples * 0.977)
    lower_index_3sig, upper_index_3sig = int(samples * 0.002), int(samples * 0.998)

    lower_1sig, upper_1sig = medians[lower_index_1sig], medians[upper_index_1sig]
    lower_2sig, upper_2sig = medians[lower_index_2sig], medians[upper_index_2sig]
    lower_3sig, upper_3sig = medians[lower_index_3sig], medians[upper_index_3sig]
    
    upper = np.vstack([upper_1sig, upper_2sig, upper_3sig])
    lower = np.vstack([lower_1sig, lower_2sig, lower_3sig])

    return med_sma, bootstrapped_median, lower, upper
