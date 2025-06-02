import galprime as gp

import numpy as np
from astropy.table import Table
from astropy.visualization import ZScaleInterval
import matplotlib.pyplot as plt

import warnings


def load_mag_kde(config):
    """Load the magnitude KDE if specified in the config."""
    if (str(config["FILES"]["MAG_CATALOGUE"]).lower() != "none" and 
        config["FILES"]["MAG_CATALOGUE"] is not None):
        mag_table = Table.read(f'{config["FILE_DIR"]}{config["FILES"]["MAG_CATALOGUE"]}')
        mags = mag_table[config["KEYS"]["MAG"]]
        return gp.object_kde([mags])
    else:
        return None


def load_and_trim_table(config, logger):
    """ Load the catalogue and trim it according to the config. """
    table = Table.read(f'{config["FILE_DIR"]}{config["FILES"]["CATALOGUE"]}')
    table = gp.trim_table(table, config)
    logger.info(f'Loaded catalogue with {len(table)} entries')
    return table


class GPrimeSingle:
    """ A single instance of a GalPRIME iteration. """
    
    def __init__(self, config, model, params, bg=None, psf=None, 
                 logger=None, id=None, save_output=False, metadata={}):
        self.config = config
        self.model = model
        self.params = params

        self.save_output = save_output

        self.id = id if id is not None else np.random.randint(1e9, 1e10)

        self.bg = bg
        self.psf = psf

        self.logger = logger

        self.stop_code = 0
        self.isophote_lists = []

        self.metadata = metadata
        self.metadata["ID"] = self.id


    def process(self):
        """
        Executes the full simulation processing pipeline for a single object.
        The processing steps include:
            1. Model generation and PSF convolution.
            2. Addition of model to background, background estimation, and subtraction.
            3. Mask generation for both background-added and background-subtracted images.
            4. Extraction of isophotal profiles from the convolved model, background-added, and 
                background-subtracted images.
        At each stage, updates internal state and handles errors by raising RuntimeError with 
                informative messages.
        The following attributes are updated during processing:
            - self.model_image: Generated model image.
            - self.model_params: Parameters used for model generation.
            - self.convolved_model: Model image after PSF convolution.
            - self.bg_added_model: Model image with background added.
            - self.bg_params: Estimated background statistics (mean, median, std).
            - self.params["BG_MEAN"], self.params["BG_MED"], self.params["BG_STD"]: Background statistics.
            - self.source_mask: Source mask from background estimation.
            - self.background: Background model object.
            - self.bg_model: 2D background model array.
            - self.bgsub: Background-subtracted image.
            - self.mask_bgadded, self.mask_data_bgadded: Masks for background-added image.
            - self.mask_bgsub, self.mask_data_bgsub: Masks for background-subtracted image.
            - self.isophote_lists: List of isophotal profile results for each processed image.
            - self.stop_code: Integer code indicating the current processing stage or completion.
        Raises:
            RuntimeError: If any processing step fails, with a message indicating the stage and error details.
        """

        # Generate model and convolve with PSF
        try:
            self.stop_code = 1
            self.model_image, self.model_params = self.model.generate(self.params)
            self.convolved_model = gp.convolve_model(self.model_image, self.psf)
        except Exception as e:
            raise RuntimeError(f'{self.id} failed convolution: {e}')

        # Add model to background, estimate background, and subtract
        try:    
            self.stop_code = 2
            self.bg_added_model = self.convolved_model + self.bg

            self.bg_params = gp.estimate_background_sigclip(self.bg_added_model, self.config)
            self.params["BG_MEAN"], self.params["BG_MED"], self.params["BG_STD"] =  self.bg_params

            self.source_mask, self.background = gp.estimate_background_2D(self.bg_added_model, self.config)
            
            self.bg_model = self.background.background
            self.bgsub = self.bg_added_model - self.bg_model
        except Exception as e:
            raise RuntimeError(f'{self.id} failed bg-modeling: {e}')

        # Mask image(s)
        try:    
            self.stop_code = 3
            self.mask_bgadded, self.mask_data_bgadded = gp.gen_mask(self.bg_added_model, config=self.config)
            self.mask_bgsub, self.mask_data_bgsub = gp.gen_mask(self.bgsub, config=self.config)
        except Exception as e:
            raise RuntimeError(f'{self.id} failed masking: {e}')
        
        # Extract profiles
        try:
            self.stop_code = 4
            # Extract profiles
            for dataset in [self.convolved_model, 
                            np.ma.array(self.bg_added_model, mask=self.mask_bgadded), 
                            np.ma.array(self.bgsub, mask=self.mask_bgsub)]:
                with warnings.catch_warnings():     # Suppress warnings from astropy fitting
                    warnings.simplefilter("ignore")
                    isolist = gp.isophote_fitting(dataset, self.config)
                    self.isophote_lists.append(isolist)

        except Exception as e:
            raise RuntimeError(f'{self.id} failed extraction: {e}')
        
        self.stop_code = 10

    def condensed_output(self):
        """
        Generate a condensed dictionary output containing key simulation results and metadata.
        Returns:
            dict: A dictionary with the following keys:
                - "ISOLISTS": List of ISOLIST values extracted from each item in self.isophote_lists.
                - "PARAMS": The parameters used for the simulation (self.params).
                - "ITERATION": The iteration number from metadata, or 999 if not present.
                - "BG_INDEX": The background index from metadata, or 999 if not present.
                - "PSF_INDEX": The PSF index from metadata, or 999 if not present.
                - "METADATA": The full metadata dictionary.
        Raises:
            ValueError: If an error occurs during the creation of the condensed output.
        """

        try:
            output = {"ISOLISTS": [n["ISOLIST"] for n in self.isophote_lists],
                "PARAMS": self.params,
                "ITERATION": self.metadata.get("ITERATION", 999),
                "BG_INDEX": self.metadata.get("BG_INDEX", 999),
                "PSF_INDEX": self.metadata.get("PSF_INDEX", 999),
                "METADATA": self.metadata}
        except Exception as e:
            raise ValueError(f"Error when creating condensed output: {e}")
        return output
    

    def plot_results(self, **kwargs):
        """
        Plot various stages and components of the simulation results in a 2x4 grid of subplots.
        This method visualizes the convolved model, background, coadded model, source mask, background model,
        background-subtracted image, masked background-subtracted image, and surface brightness profiles.
        It provides options to customize the appearance and save the resulting figure.
        Parameters
        ----------
        fontsize : int, optional
            Font size for plot titles and annotations. Default is 14.
        figsize : tuple, optional
            Size of the figure in inches (width, height). Default is (12, 6.5).
        cmap : str or Colormap, optional
            Colormap to use for image displays. Default is "Greys".
        prof_ylims : tuple, optional
            Y-axis limits for the profile plot (min, max). Default is (31, 24).
        filename : str, optional
            If provided, the figure will be saved to this file path instead of being displayed.
        dpi : int, optional
            Dots per inch for saving the figure. Default is 75.
        Notes
        -----
        - The method uses ZScaleInterval for image scaling.
        - Surface brightness profiles are plotted with error bands.
        - If an error occurs during profile plotting, a warning is issued.
        Raises
        ------
        Warning
            If there is an error plotting the surface brightness profiles, a warning is issued.
        """


        fontsize = kwargs.get("fontsize", 14)

        fig, ax = plt.subplots(2, 4, figsize=kwargs.get("figsize", (12, 6.5)), facecolor="white")

        lims = ZScaleInterval(contrast=0.1).get_limits(self.bg)

        cmap = kwargs.get("cmap", "Greys")

        ax[0, 0].imshow(self.convolved_model, origin="lower", cmap=cmap, vmin=lims[0], vmax=lims[1])
        self.model.plot_model_params(axis=ax[0, 0], fontsize=fontsize - 4, exclude=["X0", "Y0", "SHAPE"])
        ax[0, 0].set_title("Convolved Model", fontsize=fontsize)

        img = ax[0, 1].imshow(self.bg, origin="lower", cmap=cmap, vmin=lims[0], vmax=lims[1])
        ax[0, 1].text(0.05, 0.95, f'BG_INDEX = {self.metadata["BG_INDEX"]}', transform=ax[0, 1].transAxes,
                    fontsize=fontsize - 4, color='black', va='top', ha='left', bbox=dict(facecolor='white', alpha=0.5))
        ax[0, 1].set_title("Background", fontsize=fontsize)

        cbar_ax = ax[0, 1].inset_axes([0.05, 0.05, 0.5, 0.05])  # [left, bottom, width, height]

        coadd = self.bg + self.convolved_model
        coadd[self.mask_bgadded] = np.nan

        ax[0, 2].imshow(coadd, origin="lower", cmap=cmap, vmin=lims[0], vmax=lims[1])
        ax[0, 2].set_title("Coadded Model (Masked)", fontsize=fontsize)

        ax[0, 3].imshow(self.source_mask, origin="lower", cmap="gray", vmin=0, vmax=1)
        ax[0, 3].set_title("Source Mask", fontsize=fontsize)

        bg_model_lims = ZScaleInterval(contrast=0.1).get_limits(self.bg_model)
        ax[1, 0].imshow(self.bg_model, origin="lower", cmap=cmap, vmin=bg_model_lims[0], vmax=bg_model_lims[1])
        ax[1, 0].set_title("Background Model", fontsize=10)


        ax[1, 1].imshow(self.bgsub, origin="lower", cmap=cmap, vmin=lims[0], vmax=lims[1])
        ax[1, 1].set_title("Background Subtracted", fontsize=fontsize)
        # print(obj.metadata)
        plt.suptitle(f"Object ID: {self.metadata['ID']}", fontsize=fontsize)
        
        bgsub = self.bgsub.copy()
        bgsub[self.mask_bgsub] = np.nan
        ax[1, 2].imshow(bgsub, origin="lower", cmap=cmap, vmin=lims[0], vmax=lims[1])
        ax[1, 2].set_title("BG Subtracted (Masked)", fontsize=fontsize)

        try:
            for n, color, label in zip(self.isophote_lists, ["red", "blue", "gold"], ["MODEL", "COADD", "BGSUB"]):
                n = n["ISOLIST"]
                sb = gp.to_sb(n.intens)
                low = gp.to_sb(n.intens - 3 * n.int_err)
                high = gp.to_sb(n.intens + 3 * n.int_err)
                ax[1, 3].plot(n.sma, sb, color=color)
                ax[1, 3].fill_between(n.sma, low, high, color=color, alpha=0.2)
        except Exception as e:
            warnings.warn(f"Error plotting profiles: {e}")
            
        prof_ylims = kwargs.get("prof_ylims", (31, 24))
        ax[1, 3].set_ylim(prof_ylims[0], prof_ylims[1])
        ax[1, 3].set_title("Profiles", fontsize=fontsize)
        ax[1, 3].set_xscale("log")

        for axis in ax.flatten()[:-1]:
            axis.set(xticks=[], yticks=[])
        
        plt.tight_layout()

        filename = kwargs.get("filename", None)
        if filename:
            plt.savefig(filename, dpi=kwargs.get("dpi", 75), bbox_inches="tight")
            plt.close()
        else:
            plt.show()


gprime_stop_codes = {
    0: "NOT RUN",
    1: "Started Model Gen",
    2: "Started BG handling",
    3: "Started masking",
    4: "Started profile extraction",

    10: "Finished process"
}