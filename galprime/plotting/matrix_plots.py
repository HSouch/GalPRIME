from matplotlib import pyplot as plt
import numpy as np
import galprime as gp


class MatrixPlot:
    def __init__(self, config, **kwargs):
        self.config = config

        self.nrows = kwargs.get("nrows", 1)
        self.ncols = kwargs.get("ncols", 1)

        self.outname = kwargs.get("outname", "plot.pdf")
    
    def auto_figsize(self, width=12):
        return (width, width/self.ncols*self.nrows)
    
    
    def _populate(self, i, j, axis, **kwargs):
        pass

    def plot(self, **kwargs):
        figsize = kwargs.get("figsize", self.auto_figsize())
        fig, ax = plt.subplots(self.nrows, self.ncols, facecolor="white", figsize=figsize, 
                               sharey=kwargs.get("sharey", True), sharex=kwargs.get("sharex", True))
        
        ax = np.asanyarray(ax)
        for i in range(self.nrows):
            for j in range(self.ncols):
                i_plot = self.nrows - i - 1 if kwargs.get("reverse_y", False) else i
                j_plot = j if kwargs.get("reverse_x", False) else j

                ind = i_plot*self.ncols + j_plot
                axis = ax.flatten()[ind]
                
                self._populate(i, j, axis, **kwargs)

        if kwargs.get("title", None):
            fig.suptitle(kwargs.get("title"), fontsize=kwargs.get("title_size", 16))
        
        return fig, ax
    
    def cleanup(self):
        plt.tight_layout()
    
    def save(self, **kwargs):
        plt.savefig(self.outname, dpi=200, bbox_inches="tight")
        plt.close()



class KDEPlot(MatrixPlot):
    def __init__(self, config, cat=None, **kwargs):
        super().__init__(config, **kwargs)
        self.config = gp.read_config_file(config)

        self.cols = self.config["KEYS"]

        bins, shape = {}, []
        for key in self.config["BINS"].keys()[:2]:
            bins[key] = np.array(self.config["BINS"][key])
            shape.append(len(bins[key]) - 1)
        
        try:
            self.binlist = gp.bin_catalogue(cat, bin_params=config["BINS"], params=config["KEYS"]) 
        except Exception as e:
            print(f'Failed to bin catalogue: {e}')

        self.nrows, self.ncols = shape


        self.val = kwargs.get("val", "REFF")

    def plot(self, **kwargs):
        fig, ax = super().plot(**kwargs)
        
        self.cleanup()
        self.save(**kwargs)

        return fig, ax
    
    def _populate(self, i, j, axis, **kwargs):
        try:
            indices = [i, j]
            indices.extend(kwargs.get("indices", []))
            b = self.binlist.get_bin_from_indices(indices)
        except:
            return 

        key = self.cols[self.val]
        kde = gp.object_kde(b.objects[key]).resample(10000)[0]
        kde_min, kde_max = np.nanmin(kde), np.nanmax(kde)
        d_kde = (kde_max - kde_min) / 10
        output = axis.hist(kde, bins=20, histtype="step", color=kwargs.get("color", "black"), lw=3,
                            label=f"{key} KDE", density=True)
        
        if j == self.ncols - 1:
            axis.set_ylabel(b.bin_id_stub(0), fontsize=16)
            axis.yaxis.set_label_position("right")
        
        if i == 0:
            axis.set_xlabel(b.bin_id_stub(1), fontsize=16)
            axis.xaxis.set_label_position("top")
        
        if i == self.nrows - 1:
            axis.set_xlabel(key, fontsize=16)
            axis.xaxis.set_label_position("bottom")


    @staticmethod
    def plot_all_kdes(config, cat=None, outname_prefix="KDEPLOT", **kwargs):
        
        for key, val in config["KEYS"].items():
            outname = f"{outname_prefix}_{key}.pdf"
            
            kdeplot = KDEPlot(config, cat=cat, val=key, outname=outname, **kwargs)
            kdeplot.plot(title=f"{key} KDE", **kwargs)
