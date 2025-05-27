
from scipy.stats import norm
from photutils.isophote import IsophoteList

from .fluxes import integrated_rX

from numpy import argmin, clip

def ensure_profile_format(profile):
    if isinstance(profile, IsophoteList):
        return profile.to_table()
    else:
        return profile

def profile_quality(profile, bg_info, r50_mult=5, alpha=0.05):
    """
    Evaluate the quality of a given profile based on background information and specified parameters.
    Parameters:
        profile (dict): The extracted profile in table or isolist format.
        bg_info (tuple): A tuple containing background information (mean, standard deviation).
        r50_mult (float, optional): A multiplier for the r50 value to determine the threshold. Default is 5.
            This is how far we would ideally like to extract the profile to, all things considered.
        alpha (float, optional): Significance level for the threshold calculation. Default is 0.05.
            1 - alpha is the percentage confidence that a given value does not belong to the background.
    Returns:
        float: A quality metric between 0 and 1, indicating the quality of the profile.
    """

    profile = ensure_profile_format(profile)
    sma, intens = profile['sma'], profile['intens']
    
    r50_val = integrated_rX(profile['sma'], profile['intens'], profile["ellipticity"], x=0.5)[1]

    sf = 1 - norm.cdf(intens, loc=bg_info[1], scale=bg_info[2])
    
    threshold_index = argmin(sf < alpha)
    threshold = sma[threshold_index]
    
    try:
        quality = clip(threshold / (r50_val * r50_mult), 0, 1)
    except ZeroDivisionError:
        quality = 0
    
    return quality
