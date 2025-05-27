import numpy as np

from scipy.special import gamma, gammaincinv


def r_circ(a, ellip):
    """ Calculates the circularized radius (to normalize by ellipticity).

        R_circ = sqrt(a * b) = sqrt(a**2 * (1 - e))

    Args:
        a (float): The semimajor axis.  
        ellip (float): The ellipticity (e = 1 - b/a)

    Returns:
        float: The circularized radius
    """
    return np.sqrt(a**2 * (1 - ellip))


def to_mag(f, m_0=27):
    """ Convert a flux to a given magnitude. """
    return -2.5 * np.log10(f) + m_0


def to_sb(f, m_0=27, arcconv=0.168):
    """ Convert a flux to a surface brightness. Defaults are for the Hyper-Suprime Cam."""
    sb = -2.5 * np.log10(f / (arcconv ** 2)) + m_0
    return sb


def b(n):
    return gammaincinv(2 * n, 0.5)


def Ltot(mag, m0=27):
    return 10**((m0-mag)/2.5)


def I_e(mag, r_e, n, m0=27):
    return Ltot(mag, m0=m0) * (b(n) ** (2 * n)) / (r_e ** 2 * 2 * np.pi * n * gamma(2 * n))
