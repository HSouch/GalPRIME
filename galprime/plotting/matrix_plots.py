import os

from matplotlib import pyplot as plt
import numpy as np
import galprime as gp

from astropy.io import fits
from astropy.table import Table


class MatrixPlot:
    def __init__(self, config, **kwargs):
        self.config = config

        self.nrows = kwargs.get("nrows", 1)
        self.ncols = kwargs.get("ncols", 1)

        self.use_run_id = kwargs.get("use_run_id", True)
        self.outname = kwargs.get("outname", "plot.pdf")
        self.ylims = kwargs.get("ylims", None)

    def auto_figsize(self, width=12):
        return (width, width/self.ncols*self.nrows)
    
    def _populate(self, i, j, axis, **kwargs):
        pass

    def plot(self, **kwargs):
        figsize = kwargs.get("figsize", self.auto_figsize(width=kwargs.get("width", 12)))
        fig, ax = plt.subplots(self.nrows, self.ncols, facecolor="white", figsize=figsize, 
                               sharey=kwargs.get("sharey", True), sharex=kwargs.get("sharex", True))

        ax = np.asanyarray(ax)
        for i in range(self.nrows):
            for j in range(self.ncols):
                i_plot = self.nrows - i - 1 if kwargs.get("reverse_y", False) else i
                j_plot = j if kwargs.get("reverse_x", False) else -j

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
        if len(shape) == 1:
            shape.append(1)
        
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


class ProfilePlot(MatrixPlot):
    def __init__(self, outdir=None, config=None, **kwargs):
        self.outdir = outdir
        
        if outdir is None and config is None:
            raise ValueError("Outdir must be specified")

        if config is None:
            self.files = gp.gen_filestructure(outdir, generate=False)
            self.config_filename = None
            for f in os.listdir(self.files["ADDL_DATA"]):
                if f.startswith("config_"):
                    self.config_filename = os.path.join(self.files["ADDL_DATA"], f)
            
            if self.config_filename is None:
                raise ValueError("Could not find config file in ADDL_DATA directory")
            
            self.run_id = kwargs.get("run_id", gp.get_run_id(self.outdir))
            self.config = gp.read_config_file(self.config_filename)
        else:
            self.config = gp.read_config_file(config)

        self.x_index = kwargs.get("x_index", 0)
        self.y_index = kwargs.get("y_index", 1)
        self.bin_keys = self.config["BINS"].keys()

        self.xbins = self.config["BINS"][self.bin_keys[self.x_index]]
        self.ybins = self.config["BINS"][self.bin_keys[self.y_index]]

        super().__init__(self.config, 
                         ncols=len(self.ybins) - 1,
                         nrows=len(self.xbins) - 1,
                         **kwargs)

        self.ylims = kwargs.get("ylim", [32, 25])


    def plot(self, **kwargs):
        fig, ax = super().plot(**kwargs)

        self.draw_bin_labels(ax, **kwargs)
        self.cleanup()

        xmin, xmax = ax[0, 0].get_xlim()
        ymin, ymax = ax[0, 0].get_ylim()
        dx, dy = xmax - xmin, ymax - ymin
        outtext = ""
        if kwargs.get("plot_run_id", True):
            outtext += f'Run ID = {self.run_id}\n'
        if kwargs.get("plot_nmodels", True):
            outtext += r'$N_{models}$' + f' = {kwargs.get("nmodels", self.config["MODEL"]["N_MODELS"])}\n'

        ax[0, -1].text(xmax - 0.05 * dx, ymin + 0.98 * dy, outtext, 
                                    fontsize=kwargs.get("fontsize", 10), ha="right", va="top")

    def draw_bin_labels(self, ax, **kwargs):
        ylabel = kwargs.get("ylabel", r'$\log_{10}\left(\frac{F}{A_{pix}}\right)$ [mag / arcsec $^2$]')
        xlabel = kwargs.get("xlabel", r'$\log_{10}(R)$ [pix]')

        rows_param = kwargs.get("row_label", "z")
        cols_param = kwargs.get("col_label", r"$\log_{10}(M_*)$")

        for i in range(self.nrows):
            ax[i, 0].set_ylabel(ylabel, fontsize=kwargs.get("fontsize", 12))

            ax[i, -1].set_ylabel(f'{self.xbins[i]} < {rows_param} < {self.xbins[i + 1]}',
                                fontsize=kwargs.get("fontsize", 12))
            ax[i, -1].yaxis.set_label_position("right")

        for j in range(self.ncols):
            ax[-1, j].set_xlabel(xlabel, fontsize=kwargs.get("fontsize", 12))
            ax[0, j].set_title(f'{self.ybins[j]} < {cols_param} < {self.ybins[j + 1]}')
        
        ax[-1, -1].legend(fontsize=kwargs.get("legend_fontsize", 10), frameon=False)


    def plot_profs(self, tabs, axis, colors=None, labels=None, to_sb=True):
        for i in range(len(tabs)):
            tab = tabs[i]
            color = colors[i]
            label = labels[i]

            sma, median = tab["SMA"], tab["MEDIAN"]
            low_1sig, up_1sig = tab["LOW_1SIG"], tab["UP_1SIG"]
            low_2sig, up_2sig = tab["LOW_2SIG"], tab["UP_2SIG"]
            low_3sig, up_3sig = tab["LOW_3SIG"], tab["UP_3SIG"]
            
            if to_sb:
                median = gp.to_sb(median, m_0=self.config["MODEL"]["ZPM"], arcconv=self.config["MODEL"]["ARCCONV"])
            axis.plot(sma, median, color=color, label=label)
            
            for low, up in zip([low_1sig, low_2sig, low_3sig], [up_1sig, up_2sig, up_3sig]):
                if to_sb:
                    low = gp.to_sb(low, m_0=self.config["MODEL"]["ZPM"])
                    up = gp.to_sb(up, m_0=self.config["MODEL"]["ZPM"])
                axis.fill_between(sma, low, up, color=color, alpha=0.2)


    def _populate(self, i, j, axis, **kwargs):
        bin_indices = [kwargs.get("third_index", 0)] * len(self.bin_keys)
        bin_indices[self.x_index] = i
        bin_indices[self.y_index] = j
        
        bin_suffix = ""
        for suffix in bin_indices:
            bin_suffix += "_{}".format(suffix)
        
        run_id = self.run_id if self.use_run_id else ""

        median_set = f'{self.files["MEDIANS"]}{run_id}{bin_suffix}.fits'

        with fits.open(median_set) as hdul:
            bare = Table.read(hdul[1])[1:]
            coadd = Table.read(hdul[2])[1:]
            bgsub = Table.read(hdul[3])[1:]

        bare_color = kwargs.get("bare_color", "red")
        coadd_color = kwargs.get("coadd_color", "blue")
        bgsub_color = kwargs.get("bgsub_color", "gold")

        if kwargs.get("xlim", None) is not None:
            axis.set_xlim(kwargs.get("xlim"))

        axis.set_ylim(kwargs.get("ylim", self.ylims))
        axis.set_xscale("log")

        self.plot_profs([bare, coadd, bgsub], axis, 
                        colors=[bare_color, coadd_color, bgsub_color], 
                        labels=["Model", "Coadd", "BG Sub"],
                        to_sb=kwargs.get("to_sb", True))
        
        if kwargs.get("grid", True):
            axis.grid(zorder=0, which="both", linestyle="--", alpha=0.3)
    

class DiffPlot(ProfilePlot):
    def __init__(self, outdir=None, **kwargs):
        super().__init__(outdir, **kwargs)
        self.ylims = [-1, 1]

    def plot_profs(self, tabs, axis, colors=None, labels=None):
        bare_tab = tabs[0]
        bare_sma, bare_median = bare_tab["SMA"], bare_tab["MEDIAN"]
        bare_low_1sig, bare_up_1sig = bare_tab["LOW_1SIG"], bare_tab["UP_1SIG"]
        bare_low_2sig, bare_up_2sig = bare_tab["LOW_2SIG"], bare_tab["UP_2SIG"]
        bare_low_3sig, bare_up_3sig = bare_tab["LOW_3SIG"], bare_tab["UP_3SIG"]

        for i in range(1, len(tabs)):
            tab = tabs[i]
            sma, median = tab["SMA"], tab["MEDIAN"]
            low_1sig, up_1sig = tab["LOW_1SIG"], tab["UP_1SIG"]
            low_2sig, up_2sig = tab["LOW_2SIG"], tab["UP_2SIG"]
            low_3sig, up_3sig = tab["LOW_3SIG"], tab["UP_3SIG"]

            if len(bare_median) != len(median):
                continue

            diff = (median - bare_median) / bare_median
            diff_low_1sig = (low_1sig - bare_low_1sig) / bare_median
            diff_up_1sig = (up_1sig - bare_up_1sig) / bare_median
            diff_low_2sig = (low_2sig - bare_low_2sig) / bare_median
            diff_up_2sig = (up_2sig - bare_up_2sig) / bare_median
            diff_low_3sig = (low_3sig - bare_low_3sig) / bare_median
            diff_up_3sig = (up_3sig - bare_up_3sig) / bare_median

            axis.plot(sma, diff, color=colors[i ], label=labels[i])
            axis.fill_between(sma, diff_low_1sig, diff_up_1sig, color=colors[i], alpha=0.2)
            axis.fill_between(sma, diff_low_2sig, diff_up_2sig, color=colors[i], alpha=0.2)
            axis.fill_between(sma, diff_low_3sig, diff_up_3sig, color=colors[i], alpha=0.2)

        axis.axhline(0, color="black", linestyle="--", alpha=0.5, zorder=0)
