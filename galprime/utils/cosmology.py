import numpy as np


def get_angular_size_dist(z, H0=71, WM=0.27):
    """
    Return the angular size distance in Megaparsecs.
    (Stripped down version of Cosmocalc by Ned Wright and Tom Aldcroft (aldcroft@head.cfa.harvard.edu))

    :param z: The redshift.
    :type z: float
    :param H0: The Hubble constant, defaults to 71 km/s/Mpc
    :type H0: float, optional
    :param WM: matter density parameter, defaults to 0.27
    :type WM: float, optional
    """
    try:
        c = 299792.458  # velocity of light in km/sec

        if z > 100:
            z /= 299792.458  # Values over 100 are in km/s

        WV = 1.0 - WM - 0.4165 / (H0 * H0)  # Omega(vacuum) or lambda
        age = 0.0  # age of Universe in units of 1/H0

        h = H0 / 100.
        WR = 4.165E-5 / (h * h)  # includes 3 massless neutrino species, T0 = 2.72528
        WK = 1 - WM - WR - WV
        az = 1.0 / (1 + 1.0 * z)
        n = 1000  # number of points in integrals
        for i in range(n):
            a = az * (i + 0.5) / n
            adot = np.sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
            age += 1. / adot

        DCMR = 0.0

        # do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
        for i in range(n):
            a = az + (1 - az) * (i + 0.5) / n
            adot = np.sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
            DCMR = DCMR + 1. / (a * adot)

        DCMR = (1. - az) * DCMR / n

        # tangential comoving distance
        ratio = 1.0
        x = np.sqrt(abs(WK)) * DCMR
        if x > 0.1:
            if WK > 0:
                ratio = 0.5 * (np.exp(x) - np.exp(-x)) / x
            else:
                ratio = np.math.sin(x) / x
        else:
            y = x * x
            if WK < 0:
                y = -y
            ratio = 1. + y / 6. + y * y / 120.
        DCMT = ratio * DCMR
        DA = az * DCMT
        Mpc = lambda x: c / H0 * x
        DA_Mpc = Mpc(DA)

        return DA_Mpc
    except:
        raise ValueError
    

def kpc_per_arcsec(z=0.02):
    return get_angular_size_dist(z=z) * 1e3 / 206265