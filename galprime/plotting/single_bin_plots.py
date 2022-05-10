from matplotlib import pyplot as plt
from astropy.table import Table
from numpy import copy, log10
import tbridge


def adjust_profile(x, y, m_0=27, arcconv=0.168):

    x_new, y_new = copy(x), copy(y)
    y_new = -2.5 * log10(y / (arcconv ** 2)) + m_0

    return x_new, y_new


def single_bin_plot(table_arrays, colours=None, ind_profile_alpha=1, medians=False, to_sb=True,
                    bootstraps=False, ylim=(30, 20),
                    xlabel=None, ylabel=None, output="", dpi=150):

    plt.figure(figsize=(8, 6))

    for i in range(0, len(table_arrays)):
        for j in range(0, len(table_arrays[i])):
            this_prof = table_arrays[i][j]

            sma, intens = this_prof["sma"], this_prof["intens"]
            sma, intens = adjust_profile(sma, intens)

            plt.plot(sma, intens, c=colours[i], alpha=ind_profile_alpha, lw=1)

        if medians:
            interps = tbridge.as_interpolations(table_arrays[i])
            median_sma, median_interp = tbridge.get_median(interps, tbridge.bin_max(table_arrays[i]))
            median_intens = median_interp(median_sma)

            if to_sb:
                median_sma, median_intens = adjust_profile(median_sma, median_intens)

            plt.plot(median_sma, median_intens, lw=3, c=colours[i])

    if ylim is None:
        pass
    else:
        plt.ylim(ylim[0], ylim[1])
    plt.xlim(1,)
    plt.xscale("log")

    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)

    plt.tight_layout()

    if output == "":
        plt.show()
    else:
        plt.savefig(output, dpi=dpi)


def single_prof(x, y, error=None):
    plt.figure(figsize=(8, 6))

    if error is None:
        plt.plot(x, y, color="orange", lw=2)
    else:
        plt.errorbar(x, y, yerr=error)

    plt.tight_layout()
    plt.show()
