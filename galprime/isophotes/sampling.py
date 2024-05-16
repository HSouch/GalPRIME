import numpy as np


def eccentric_anomaly(phi, ellip=0):
    """
        Calculate the eccentric anomaly
        
        The np.cos() * 2 - 1 is a step function setting the quadrants
    """
    return np.sort(-np.arctan2(np.tan(phi) / (1-ellip), (np.cos(phi) > 0) * 2 - 1))


def ellipse_coords(isophote, angles, asint = True):
    
    x0, y0, a, pa, eps = isophote.x0, isophote.y0, isophote.sma, isophote.pa, isophote.eps
    b = a * (1 - eps)

    xs = x0 + a * np.cos(angles) * np.cos(pa) - b * np.sin(angles) * np.sin(pa)
    ys = y0 + a * np.cos(angles) * np.sin(pa) + b * np.sin(angles) * np.cos(pa)

    if asint:
        return (xs.astype(int), ys.astype(int))
    else:
        return (xs, ys)


def sample(image, isophote, thetas, min_samples=10):

    x, y = ellipse_coords(isophote, thetas, min_samples=min_samples, asint=True)
    # Remove points outside the image
    mask = (x >= 0) & (x < image.shape[1]) & (y >= 0) & (y < image.shape[0])
    x, y = x[mask], y[mask]

    return thetas[mask], image[x, y]


def intensity(image, isophote, thetas, **kwargs):
    return np.nanmean(sample(image, isophote, thetas, **kwargs)[1])


def gradient(image, isophote, thetas):
    iso_inner, iso_outer = isophote.copy(), isophote.copy()
    iso_inner.sma *= 0.95
    iso_outer.sma *= 1.05

    return intensity(image, iso_outer, thetas) - intensity(image, iso_inner, thetas)


