from scipy.stats import gaussian_kde

def object_kde(columns):
    """ Generate a gaussian kernel density estimate for the columns of a bin. """
    return gaussian_kde(columns)