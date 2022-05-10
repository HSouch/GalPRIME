import galprime
import galprime.plotting as plotter
from astropy.table import Table
from astropy.io import fits

import time
import traceback

from pebble import ProcessPool
from concurrent.futures import TimeoutError

from numpy import round, array
from numpy.random import choice, randint, seed


def pipeline(config, max_bins=None, mag_table=None,
             provided_bgs=None, provided_psfs=None,
             progress_bar=False, multiprocess_level='obj'):
    """
    Runs the entire simulation pipeline assuming certain data exists.
    :param config: Values from properly loaded configuration file.
    :param max_bins: The number of bins to process (useful if running tests).
    :param mag_table: Optional array of magnitudes.
    :param mag_table_keys Keys for optional mag table binning.
    :param provided_bgs: A set of provided background cutouts [OPTIONAL].
    :param provided_psfs: A set of provided PSFs related to the provided backgrounds [OPTIONAL].
    :param progress_bar: Have a TQDM progress bar.
    :param multiprocess_level: Where in the simulations to divide into cores
        'obj'  - Divide at the object level, where each core handles a single object in each bin.
        'bin'  - Divide at the bin level, so each core is responsible for a single bin
        'none' - Do not do any multiprocessing at all (SLOW).

        The object level is less memory intensive, but bins are processed one by one instead of simultaneously.

    EXAMPLE USAGE:

    config_values = tbridge.load_config_file("path/to/config/file.tbridge")
    tbridge.pipeline(config_values, max_bins=10)

    """

    binned_objects = galprime.bin_catalog(config)
    max_bins = len(binned_objects) if max_bins is None else max_bins

    verbose = config["VERBOSE"]
    if verbose:
        print(max_bins, "bins to process.")

    if config["SAME_BGS"] and provided_bgs is None:
        provided_bgs, provided_psfs = galprime.get_backgrounds(config, n=int(config["N_BGS"]),
                                                              return_psfs=True, return_bg_info=False)
        if verbose:
            print(len(provided_bgs), "obtained with associated PSFs.")

    for b in binned_objects[:max_bins]:

        if mag_table is not None:
            separate_mags = galprime.bin_mag_catalog(mag_table, b,
                                                    mag_table_keys=[config["MASS_KEY"], config["Z_KEY"]],
                                                    bin_keys=["MASSES", "REDSHIFTS"])[config["MAG_KEY"]]
        else:
            separate_mags = None

        _process_bin(b, config, separate_mags=separate_mags,
                     provided_bgs=provided_bgs, provided_psfs=provided_psfs,
                     progress_bar=progress_bar, multiprocess=True)

    galprime.config_to_file(config, filename=config["OUT_DIR"] + "tbridge_config.txt")


def _process_bin(b, config_values, separate_mags=None, provided_bgs=None, provided_psfs=None,
                        progress_bar=False, multiprocess=False, profiles_per_row=3):
    """
        Process a single bin of galaxies. (Tuned for pipeline usage but can be used on an individual basis.
        :param b: Bin to obtain full profiles from.
        :param config_values: Values from properly loaded configuration file.
        :param separate_mags: Optional array of magnitudes.
        :param provided_bgs: Array of provided backgrounds
        :param progress_bar: Use a TQDM progress bar (note with multithreading this might get funky).
        :param multiprocess: Run in multiprocessing mode.
            This means both the model creation AND the
        :return:
        """
    t_start = time.time()
    verbose = config_values["VERBOSE"]

    # Load in bin information, and prepare all necessary structural parameters.
    keys, columns = b.return_columns()
    mags, r50s, ns, ellips = galprime.pdf_resample(galprime.structural_parameters(keys, columns), resample_size=1000)
    if separate_mags is not None:
        mags = galprime.pdf_resample(separate_mags, resample_size=len(r50s))[0]

    # Simulate the model rows, using multiprocessing to speed things up#
    if verbose:
        print("Processing", config_values["N_MODELS"], "models for: ", b.bin_params)

    # Prepare containers for simulations
    job_list = []
    full_profile_list, bg_infolist, cutout_infolist, bg_2D_profiles = [[] for i in range(profiles_per_row)], [], [], []
    masked_cutouts, unmasked_cutouts = [], []
    param_results = Table(data=None, names=["MAG", "R50", "I_R50", "N", "ELLIP", "PA"])

    # Get the original models
    models, model_parameters = galprime.simulate_sersic_models(mags, r50s, ns, ellips,
                                                              config_values, n_models=config_values["N_MODELS"])

    # Run multithreaded simulation code
    with ProcessPool(max_workers=config_values["CORES"]) as pool:
        for i in range(len(models)):
            job_list.append(pool.schedule(_process_model,
                                          args=(models[i], config_values, model_parameters[i],
                                                provided_bgs, provided_psfs),
                                          timeout=config_values["ALARM_TIME"]))
    # Collect the results
    for i in range(len(job_list)):
        try:
            result = job_list[i].result()
            profiles = result["PROFILES"]
            if len(profiles) != profiles_per_row:
                continue
            # If we got enough profiles, append everything to the appropriate arrays
            for i in range(len(full_profile_list)):
                full_profile_list[i].append(profiles[i])

            parameters = result["INFO"]
            row = [parameters[key] for key in param_results.colnames]
            param_results.add_row(row)

            bg_infolist.append(result["BG_DATA"])
            cutout_infolist.append(result["CUTOUT_DATA"])

            # masked_cutouts.append(result["MASKED_CUTOUT"])
            # unmasked_cutouts.append(result["UNMASKED_CUTOUT"])
            bg_2D_profiles.append(result["BGSUB_PROFILE"])

        except Exception as error:
            if config_values["TEST_VERBOSE"]:
                print(error.args, i, traceback.print_exc())
            else:
                print(error.args, i)

    bg_info, cutout_info = [[], [], []], [[], [], []]

    for i in range(0, len(bg_infolist)):
        bg_info[0].append(bg_infolist[i]["BG_MEAN"])
        bg_info[1].append(bg_infolist[i]["BG_MEDIAN"])
        bg_info[2].append(bg_infolist[i]["BG_STD"])

    for i in range(0, len(cutout_infolist)):
        cutout_info[0].append(cutout_infolist[i]["BG_MEAN"])
        cutout_info[1].append(cutout_infolist[i]["BG_MEDIAN"])
        cutout_info[2].append(cutout_infolist[i]["BG_STD"])

    # Subtract the median values from the bgadded profiles
    if config_values["BG_PARAMS"] == "2D":
        print("Using 2D profiles")
        bg_sub_profiles = bg_2D_profiles
    else:
        bg_sub_profiles = galprime.subtract_backgrounds(full_profile_list[2], bg_info[1])

    full_profile_list.append(bg_sub_profiles)

    if verbose:
        print("Finished extraction. Saving info for", len(full_profile_list[0]),
              "successful extractions out of", len(models), ".")

    # Save the profiles to the required places
    galprime.save_profiles(full_profile_list,
                          bin_info=b.bin_params,
                          out_dir=config_values["OUT_DIR"],
                          keys=["bare", "noisy", "bgadded", "bgsub"],
                          bg_info=bg_info, cutout_info=cutout_info,
                          structural_params=param_results)

    if verbose:
        print("Finished", b.bin_params, "-- Time Taken:", round((time.time() - t_start) / 60, 2), "minutes.")
        print()


def _process_model(sersic_model, config, model_params, provided_bgs=None, provided_psfs=None):
    """
    Run processing for a single model.

    :param sersic_model: The input model to process
    :param config: Configuration values loaded from config file
    :param provided_bgs: OPTIONAL --- set of provided backgrounds
    """
    # First make the input models

    if provided_bgs is None:
        background, psf, bg_data = galprime.get_background(config)
        convolved_model = galprime.convolve_models(sersic_model, psf=psf)
        bg_added_model = convolved_model + background

    else:
        if config["PROVIDED_BG_TYPE"] == "Cutouts":
            # We will now randomly select a background and psf, and use that as our background
            seed()
            bg_index = randint(0, len(provided_bgs),)
            convolved_model = galprime.convolve_models(sersic_model, psf=provided_psfs[bg_index])
            bg_added_model = convolved_model + provided_bgs[bg_index]

        else:
            seed()
            cutout, row_data = galprime.cutout_from_table(provided_bgs)

            with fits.open(config["PSF_FILENAME"]) as psfs:
                psf = galprime.get_closest_psf(psfs, row_data["RA"], row_data["DEC"])
            convolved_model = galprime.convolve_models(sersic_model, config_values=config, psf=psf)
            bg_added_model = convolved_model + cutout

    # Estimate the background using the sigclip method before we do anything to it with 2D Background subtraction
    bg_mean, bg_median, bg_std = galprime.estimate_background_sigclip(bg_added_model)
    bg_data = {"BG_MEAN": bg_mean, "BG_MEDIAN": bg_median, "BG_STD": bg_std}

    # If we're doing 2D background subtraction we'll need to estimate and then subtract the background BEFORE
    # profile extraction.
    bgsub_model = None
    if config["BG_PARAMS"] == "2D":
        sm = galprime.SourceMask(bg_added_model)
        mask = sm.multiple(filter_fwhm=[1, 3, 5], tophat_size=[4, 2, 1])
        bg = galprime.background_2D(bg_added_model, mask,
                                   box_size=config["BOX_SIZE"], filter_size=config["FILTER_SIZE"])
        bgsub_model = bg_added_model - bg.background

    # Now we will make the cutout smaller for extraction purposes
    convolved_model = galprime.central_cutout(convolved_model, width=config["EXTRACTION_SIZE"])
    bg_added_model = galprime.central_cutout(bg_added_model, width=config["EXTRACTION_SIZE"])
    bgsub_model = galprime.central_cutout(bgsub_model, width=config["EXTRACTION_SIZE"])

    # Generate the noise model and mask the background-added model
    noisy_model = galprime.add_to_noise(convolved_model)
    masked_model, mask_data = galprime.mask_cutout(bg_added_model, config=config)
    model_row = [convolved_model, noisy_model, masked_model]

    # Then extract everything we can
    profile_extractions = []
    for i in range(0, len(model_row)):
        try:
            extraction = galprime.isophote_fitting(model_row[i], config=config)
            profile_extractions.append(extraction.to_table())
        except Exception as error:
            print(error.args)
            continue

    # We'll extract the BG-subtracted profile is necessary
    bgsub_profile = None
    if config["BG_PARAMS"] == "2D":
        try:
            bgsub_profile = galprime.isophote_fitting(galprime.mask_cutout(bgsub_model, config=config)[0],
                                                     config=config).to_table()
        except Exception as error:
            print(error.args)

    # put the results into a dictionary format and return everything
    return {"PROFILES": profile_extractions,
            "BGSUB_PROFILE": bgsub_profile,
            "MASK_DATA": mask_data,
            "BG_DATA": bg_data,
            "CUTOUT_DATA": {"BG_MEAN": bg_mean, "BG_MEDIAN": bg_median, "BG_STD": bg_std},
            # "UNMASKED_CUTOUT": bg_added_model,
            # "MASKED_CUTOUT": masked_model,
            "INFO": model_params}
