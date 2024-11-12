import numpy as np
from scipy.interpolate import interp1d


def highest_sma(profile_list):
    return np.min([max(profile.sma) for profile in profile_list])


def profiles_norm_sma(profile_list):
    sma_max = highest_sma(profile_list)
    sma_interp = np.arange(1, sma_max, 2)
    profile_interps = [interp1d(profile.sma, profile.intens, kind="linear", 
                                fill_value=0, bounds_error=False,
                                ) for profile in profile_list]
    profiles_stacked = np.vstack([interp(sma_interp) for interp in profile_interps])
    return sma_interp, profiles_stacked


def _median(profile_stack):
    return np.median(profile_stack, axis=0)


def gen_median(profile_list):
    profile_stack = profiles_norm_sma(profile_list)[1]
    return _median(profile_stack)


def bootstrap_median(profile_list, n_bootstraps=10000):
    smas, profile_stack = profiles_norm_sma(profile_list)

    # plt.imshow(profile_stack, aspect="auto", vmin=-0.005, vmax=0.005)

    sample_indices = np.random.choice(profile_stack.shape[0], (n_bootstraps, profile_stack.shape[0]), replace=True)

    medians = np.vstack([_median(profile_stack[sample]) for sample in sample_indices])

    sorted = np.sort(medians, axis=0)

    lower_index_1sig, upper_index_1sig = int(n_bootstraps * 0.159), int(n_bootstraps * 0.841)
    lower_index_2sig, upper_index_2sig = int(n_bootstraps * 0.023), int(n_bootstraps * 0.977)
    lower_index_3sig, upper_index_3sig = int(n_bootstraps * 0.002), int(n_bootstraps * 0.998)

    lower_1sig, upper_1sig = sorted[lower_index_1sig], sorted[upper_index_1sig]
    lower_2sig, upper_2sig = sorted[lower_index_2sig], sorted[upper_index_2sig]
    lower_3sig, upper_3sig = sorted[lower_index_3sig], sorted[upper_index_3sig]

    upper = np.vstack([upper_1sig, upper_2sig, upper_3sig])
    lower = np.vstack([lower_1sig, lower_2sig, lower_3sig])

    return smas, _median(profile_stack), lower, upper
