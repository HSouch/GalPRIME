from astropy.visualization import ZScaleInterval
from matplotlib import pyplot as plt

import numpy as np

def show_cutouts(cutouts, nrows=5, ncols=5, method="zscale", cmap="gray_r", **kwargs):
    print(len(cutouts.cutouts))

    outname = kwargs.get("outname", None)
    dpi = kwargs.get("dpi", 150)

    vmin, vmax = kwargs.get("vmin", -3), kwargs.get("vmax", 1)

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 2.5, nrows * 2.5), facecolor="white", dpi=dpi)
    for i in range(nrows):
        for j in range(ncols):
            
            cutout = cutouts.cutouts[i * ncols + j]
            
            if method == "zscale":
                cutout = cutout
                interval = ZScaleInterval()
                vmin, vmax = interval.get_limits(cutout)

            if method == "linear":
                cutout = cutout
                vmin, vmax = vmin, vmax
            
            if method == "log":
                cutout = np.log10(cutout)
                vmin, vmax = vmin, vmax

            axes[i, j].imshow(cutout, cmap=cmap, vmin=vmin, vmax=vmax)
            axes[i, j].set_xticks([])
            axes[i, j].set_yticks([])


    plt.tight_layout()
    if outname is not None:
        plt.savefig(outname)
    else:
        plt.show()