import numpy as np
import matplotlib as mpl


def pyplot_style(**kwargs):
    """ If you just call pyplot style, return the default. """
    default_pyplot_style(**kwargs)


def default_pyplot_style(**kwargs):
    # xticks
    mpl.rcParams["xtick.major.size"] = kwargs.get("xtick_major_size", 2)
    mpl.rcParams["xtick.major.width"] = kwargs.get("xtick_major_width", 1)
    mpl.rcParams["xtick.minor.size"] = kwargs.get("xtick_minor_size", 1.5)
    mpl.rcParams["xtick.minor.width"] = kwargs.get("xtick_minor_width", 0.75)

    # yticks
    mpl.rcParams["ytick.major.size"] = kwargs.get("ytick_major_size", 2)
    mpl.rcParams["ytick.major.width"] = kwargs.get("ytick_major_width", 1)
    mpl.rcParams["ytick.minor.size"] = kwargs.get("ytick_minor_size", 1.5)
    mpl.rcParams["ytick.minor.width"] = kwargs.get("ytick_minor_width", 0.75)

    mpl.rcParams["axes.linewidth"] = kwargs.get("axes_linewidth", 1.25)

    mpl.rc("xtick", labelsize=kwargs.get("xtick_labelsize", 12))
    mpl.rc("ytick", labelsize=kwargs.get("ytick_labelsize", 12))

    font = {
        "family": "serif",
        #'weight': 'bold',
        "size": kwargs.get("font_size", 12),
    }

    mpl.rc("font", **font)

    mpl.rc("lines", linewidth=kwargs.get("linewidth", 1.5), 
           linestyle="solid", color="black")


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
