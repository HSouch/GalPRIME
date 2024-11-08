from scipy.signal import convolve2d


def convolve_model(model, psf):
    return convolve2d(model, psf, mode='same')