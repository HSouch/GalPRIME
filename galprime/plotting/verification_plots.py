from matplotlib import pyplot as plt


def plot_binlist(binlist, cols=6, key="R50", **kwargs):
    """ Check that the binlist rebinning is working as expected.

    Args:
        binlist (BinList): A BinList object.
        cols (int): Number of columns. The number of rows is calculated automatically from this
        key (str): The key to plot. Default is "R50".

    Returns:
        _type_: _description_
    """
    
    color = kwargs.get("color", "black")
    outname = kwargs.get("outname", None)

    rows = int(len(binlist.bins)/cols) + (len(binlist.bins) % cols > 0)
    fig, ax = plt.subplots(rows, cols, figsize=(cols*2.5, rows*2.5))

    for i in range(rows):
        for j in range(cols):
            index = i*cols + j
            if index >= len(binlist.bins):
                break
            bin = binlist.bins[index]
            ax[i, j].hist(bin.objects[key], bins=30, histtype='step', stacked=True, fill=False, color=color)
            
            xmin, xmax = ax[i, j].get_xlim()
            ymin, ymax = ax[i, j].get_ylim()
            dx, dy = xmax - xmin, ymax - ymin
            ax[i, j].text(xmin + 0.1*dx, ymin + 0.9*dy, f'{i}:{key}: {bin.bin_info[key]}', fontsize=8)

    for j in range(cols):
        ax[-1, j].set_xlabel(key)

    plt.tight_layout()
    if outname is not None:
        plt.savefig(outname)
    else:
        plt.show()

    return fig, ax
