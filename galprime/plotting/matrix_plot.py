from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.table import Table

from scipy.interpolate import interp1d
from numpy import copy, log10, isnan, isinf, round


def matrix_plot(directory, x_bins, y_bins,
                stubs=["bare", "bgadded", "noisy", "bgsub"],
                colours=["red", "blue", "green", "violet"], extra_params="_0.0-0.5"):

    print(x_bins, y_bins)

    n_x, n_y = len(x_bins), len(y_bins)

    fig, ax = plt.subplots(n_x, n_y, sharex=True, sharey=True)
    fig.set_figheight(4)
    fig.set_figwidth(5)
    for stub in stubs:
        filename_prefix = directory + stub + "_profiles/"

        for x in range(0, n_x):
            for y in range(0, n_y):
                # ax[x][y].text(0.5, 0.5, str(round(x_bins[n_x - x - 1], 2)) + "," + str(round(y_bins[y], 2)),
                #               color="red")
                pass
    plt.show()
    return None


def load_median_info(filename):
    bare_hdul = fits.open(filename)
    bare_median, l, u = Table.read(bare_hdul[1]), Table.read(bare_hdul[2]), Table.read(bare_hdul[3])
    bare_median_adj = adjust_profile(bare_median["SMA"], bare_median["INTENS"])
    bare_l_1sig_adj = adjust_profile(l["SMA"], l["INTENS_1SIG"])
    bare_u_1sig_adj = adjust_profile(u["SMA"], u["INTENS_1SIG"])
    bare_l_2sig_adj = adjust_profile(l["SMA"], l["INTENS_2SIG"])
    bare_u_2sig_adj = adjust_profile(u["SMA"], u["INTENS_2SIG"])
    bare_l_3sig_adj = adjust_profile(l["SMA"], l["INTENS_3SIG"])
    bare_u_3sig_adj = adjust_profile(u["SMA"], u["INTENS_3SIG"])
    bare_l_5sig_adj = adjust_profile(l["SMA"], l["INTENS_5SIG"])
    bare_u_5sig_adj = adjust_profile(u["SMA"], u["INTENS_5SIG"])


def adjust_profile(x, y, clean=False, shift=0.00):
    # hlr_index, hlr_value = get_half_light_radius(x, y)
    x_new, y_new = copy(x), copy(y)
    y_new = -2.5 * log10((y + shift) / (0.168 ** 2)) + 27

    if clean:
        x_clean, y_clean = [], []
        for n in range(0, len(x_new)):
            if isnan(y_new[n]) or isinf(y_new[n]):
                continue
            else:
                x_clean.append(x_new[n])
                y_clean.append(y_new[n])
        profile_adj = interp1d(x_clean, y_clean, bounds_error=False, fill_value=max(y_clean))

        return x_new, profile_adj(x_new)

    else:
        return x_new, y_new