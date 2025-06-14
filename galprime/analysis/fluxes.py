import numpy as np


def rX(sma, intens, X=0.5):
    """
    Calculate the radius at which the cumulative intensity reaches a fraction X of the total intensity.
    Parameters:
    sma (array-like): Semi-major axis values.
    intens (array-like): Intensity values corresponding to the semi-major axis values.
    X (float, optional): Fraction of the total intensity to find the radius for. 
    Must be between 0 and 1. Default is 0.5.
    Returns:
    tuple: A tuple containing the index and the semi-major axis value at which the 
    cumulative intensity reaches the fraction X of the total intensity.
    Raises:
    ValueError: If X is not between 0 and 1.
    """
    
    if X < 0 or X > 1:
        raise ValueError("X must be between 0 and 1")
    
    Ltot = np.trapz(intens, sma)
    for i in range(len(sma)):
        L = np.trapz(intens[:i], sma[:i])
        if L > X * Ltot:
            return i, sma[i]
    return len(sma) - 1, sma[-1]


def r50(sma, intens):
    """
    Calculate the half-light radius.
    Parameters:
    sma (array-like): Semi-major axis values.
    intens (array-like): Intensity values corresponding to the semi-major axis values.
    Returns:
    tuple: The index of the half-light radius and the semi-major axis value at that index.
    """

    return rX(sma, intens, 0.5)


def integrated_rX(sma, intens, ellip, width=1, x=0.5):
    """ Returns the rX value taking into account the INTEGRATED light profile, in annuli of width 'width'.
        So ellipticity is also needed, as it calculates the circularized profile.

        This is a slightly more accurate estimate of rX.

    Args:
        sma (arr): Semi-major axis values.
        intens (arr): Average intensity values.
        ellip (arr): Ellipticity values
        width (int, optional): Annulus width. Defaults to 1. 
            NOTE: If using linear scaling in the isophotes, it is recommended to use a width of the isophote step.
        x (float, optional): Fraction of light contained in rX. Defaults to 0.5. (50%)

    Returns:
        tuple: The index of rX and the value of rX.
    """
    sma_circ = np.sqrt(sma ** 2 * (1 - ellip))

    annuli = np.pi * ((sma_circ + width /2) ** 2 - (sma_circ - width / 2) ** 2)
    intens_integrated = annuli * intens

    return rX(sma_circ, intens_integrated, x)



def sb(intens, A_pix=1, mu_0=27):
    """
    Calculate the surface brightness.
    Parameters:
    intens (float): The intensity value.
    A_pix (float, optional): The pixel area. Default is 1.
    mu_0 (float, optional): The zero-point magnitude. Default is 27.
    Returns:
    float: The calculated surface brightness.
    """
    
    return mu_0 - 2.5 * np.log10(intens / A_pix)


def closest_to(sma, val):
    """
    Find the index of the element in the array that is closest to the given value.
    Parameters:
    sma (array-like): Semi-major axis values.
    val (float): The value to find the closest element to.
    Returns:
    int: The index of the element in the array that is closest to the given value.
    """
    
    return np.argmin(np.abs(sma - val))
