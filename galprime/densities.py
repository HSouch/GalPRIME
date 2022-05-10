import galprime

from numpy import mgrid, copy
from skimage.filters import gaussian
from numpy import sqrt, sort


def cutout_gradient(cutout, depth_value=5):
    """
    Estimate the gradient of a first-smoothed cutout according to
    the following formula:

        G(x,y) = ( I(x,y) - I_0 )^2 / d
            where d is the distance from (x,y) to the cutout centre.
        The returned gradient is the Nth largest gradient value, where N is
        defined by the user as depth_value.

    Args:
        cutout: The cutout to measure
        depth_value: Nth largest gradient value to return (default is 5).

    Returns:
        The gradient of the submitted cutout.
    """
    x, y = mgrid[:cutout.shape[0], :cutout.shape[1]]
    cutout_centre = (cutout.shape[0] / 2, cutout.shape[1] / 2)
    c_x, c_y = int(cutout_centre[0]), int(cutout_centre[1])

    cutout = gaussian(cutout, sigma=5)

    dists = sqrt((cutout_centre[0] - x) ** 2 + (cutout_centre[1] - y) ** 2)

    central_value = cutout[c_x][c_y]

    grads = (cutout - central_value) ** 2 / dists
    grads = sort(grads, axis=None)

    return grads[len(grads) - 1 - depth_value]


def cutout_gradients(cutouts, mask_cutouts=False):
    """
    Return the gradients for a set of cutouts. (Convenience function).

    Args:
        cutouts: List of cutouts.
        mask_cutouts: Set to True if the input cutouts need to have light from
            nearby objects masked first.

    Returns:
        list: the gradients for the supplied cutouts.
    """

    if mask_cutouts:
        cutouts, bg_info = galprime.mask_cutouts(cutouts, method='no_central')

    gradients = []
    for cutout in cutouts:
        gradients.append(cutout_gradient(cutout))
    return gradients
