import numpy as np
import matplotlib as mpl


def pyplot_style():
    """ If you just call pyplot style, return the default. """
    default_pyplot_style()


def default_pyplot_style():
    # xticks
    mpl.rcParams["xtick.major.size"] = 3
    mpl.rcParams["xtick.major.width"] = 2
    mpl.rcParams["xtick.minor.size"] = 2
    mpl.rcParams["xtick.minor.width"] = 1

    # yticks
    mpl.rcParams["ytick.major.size"] = 3
    mpl.rcParams["ytick.major.width"] = 2
    mpl.rcParams["ytick.minor.size"] = 2
    mpl.rcParams["ytick.minor.width"] = 1

    mpl.rcParams["axes.linewidth"] = 1.5

    mpl.rc("xtick", labelsize=12)
    mpl.rc("ytick", labelsize=12)

    font = {
        "family": "serif",
        #'weight': 'bold',
        "size": 12,
    }

    mpl.rc("font", **font)

    mpl.rc("lines", linewidth=1, linestyle="solid", color="black")


def lavender_cmap(step_1=175):
    x, y = np.mgrid[:251, :251]

    lavender = [86, 82, 100]
    indigo = [29, 0, 51]
    yellow = [100, 100, 0]

    step_2 = 256 - step_1

    r, g, b = [], [], []

    r = np.array(np.append(np.linspace(lavender[0], indigo[0], step_1), 
                     np.linspace(indigo[0], yellow[0], step_2))) / 100
    g = np.array(np.append(np.linspace(lavender[1], indigo[1], step_1), 
                     np.linspace(indigo[1], yellow[1], step_2))) / 100
    b = np.array(np.append(np.linspace(lavender[2], indigo[2], step_1), 
                     np.linspace(indigo[2], yellow[2], step_2))) / 100

    a = np.array([1 for i in range(len(r))])

    full_arr = np.asarray([[r[i], g[i], b[i], a[i]] for i in range(0, len(r))])

    new_cmap = mpl.colors.ListedColormap(full_arr)
    return new_cmap