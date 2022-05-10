from numpy import ceil, sqrt, log
from matplotlib import pyplot as plt


def view_cutouts(cutouts, output="", dpi=150, show_index=False, log_scale=True, limits=(1, -3), cmap="magma_r"):
    """ Generate cutouts plot for input filename"""

    width = int(ceil(sqrt(len(cutouts))))

    print(len(cutouts), "cutouts to visualize. Generating", str(width) + "x" + str(width), "figure.")

    fig, ax = plt.subplots(width, width)

    fig.set_figheight(width)
    fig.set_figwidth(width)

    index = 0
    for x in range(0, width):
        for y in range(0, width):
            ax[x][y].set_xticks([])
            ax[x][y].set_yticks([])
            try:
                image = cutouts[index]
                if log_scale:
                    ax[x][y].imshow(log(image), cmap=cmap, vmax=limits[0], vmin=limits[1])
                else:
                    ax[x][y].imshow(image, cmap=cmap, vmax=limits[0], vmin=limits[1])

                if show_index:
                    ax[x][y].text(5, 0, str(index), color="red", fontweight="bold",
                                  fontsize=10, alpha=0.7, **{'fontname': 'Helvetica'})
            except:
                pass

            index += 1
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.2, wspace=0.2)

    if output == "":
        plt.show()
    else:
        print("Saving figure to " + output)
        plt.savefig(output, dpi=dpi)

