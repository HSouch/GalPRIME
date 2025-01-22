
from .. import bgsub, config
from numpy import isclose
from numpy.random import normal

import warnings

def test_estimate_background_sigclip():
    test_background = normal(size=(100, 100))

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        bg_stats = bgsub.estimate_background_sigclip(test_background, config.default_config())

        assert len(bg_stats) == 3


def test_bg_sigclip_avg():
    averages = []
    for avg in [0, 0.5, 2]:
        test_background = normal(loc=avg, size=(100, 100))

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')
            bg_stats = bgsub.estimate_background_sigclip(test_background, config.default_config())

            assert isclose(bg_stats[0], avg, atol=1e-1, rtol=1e-1)
            assert isclose(bg_stats[1], avg, atol=1e-1, rtol=1e-1)

            averages.append(avg)
    
    for i in range(1, len(averages)):
        assert averages[i] > averages[i-1]


def test_bg_sigclip_std():
    for std in [0.5, 1, 2]:
        test_background = normal(scale=std, size=(100, 100))

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')
            bg_stats = bgsub.estimate_background_sigclip(test_background, config.default_config())

            assert isclose(bg_stats[2], std, atol=1e-1, rtol=1e-1)

