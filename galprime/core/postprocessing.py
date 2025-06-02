import galprime as gp

from astropy.table import Table
from astropy.io import fits

from tqdm import tqdm
import os, shutil


good_colnames = ['sma', 'intens', 'intens_err', 'ellipticity', 'ellipticity_err', 'pa', 'pa_err',
                 'x0', 'x0_err', 'y0', 'y0_err', 'ndata', 'nflag', 'niter', 'stop_code']


def gen_median_table(smas, median, low, up):
    """
    Generate a table containing median and bootstrapped uncertainty values.

    Parameters
    ----------
    smas : list or array-like
        Semi-major axes values.
    median : list or array-like
        Median values.
    low : list of lists or array-like
        Lower bootstrapped uncertainty (1-sigma, 2-sigma, 3-sigma).
    up : list of lists or array-like
        Upper bootstrapped uncertainty (1-sigma, 2-sigma, 3-sigma).

    Returns
    -------
    astropy.table.Table
    """

    colnames = ["SMA", "MEDIAN", "LOW_1SIG", "LOW_2SIG", "LOW_3SIG", "UP_1SIG", "UP_2SIG", "UP_3SIG"]
    return Table([smas, median, *low, *up], names=colnames)


def gen_median_hdul(bare_table, coadd_table, bgsub_table, outname):
    """
    Generate a FITS HDUList containing median tables and write to a file.

    Parameters
    ----------
    bare_table : astropy.table.Table
        The bare table data to be included in the HDUList.
    coadd_table : astropy.table.Table
        The coadded table data to be included in the HDUList.
    bgsub_table : astropy.table.Table
        The background-subtracted table data to be included in the HDUList.
    outname : str
        The output filename for the FITS file.

    Returns
    -------
    None
    """

    hdul = fits.HDUList()
    for table, name in zip([bare_table, coadd_table, bgsub_table], ["BARE", "COADD", "BGSUB"]):
        hdul.append(fits.BinTableHDU(data=table, name=name))
    hdul.writeto(outname, overwrite=True)



def dict_to_header(test):
    # For converting model params to a FITS header file
    header = fits.Header()
    for key, val in test.items():
        if isinstance(val, tuple):
            header[key] = val[0]
            header.comments[key] = val[1]
        else:
            header[key] = val
    return header


def handle_output(results, outdirs, config, bin_id="0"):
    """
    Handles the output of the processing results by saving individual profiles and generating median tables.
    Parameters
    ----------
    results : list of dict
        List of dictionaries containing the results with "ISOLISTS" key.
    outdirs : dict
        Dictionary containing output directories with keys "MODEL_PROFS", "COADD_PROFS", "BGSUB_PROFS", and "MEDIANS".
    config : dict
        Configuration dictionary containing the "RUN_ID".
    bin_id : str, optional
        Identifier for the bin, by default "0".
    Returns
    -------
    None
    """
    
    # Grab all of the tables and sort them into the appropriate locations
    bare_profiles  = [n["ISOLISTS"][0].to_table() for n in results]
    coadd_profiles = [n["ISOLISTS"][1].to_table() for n in results]
    bgsub_profiles = [n["ISOLISTS"][2].to_table() for n in results]

    # mod_params = [n["PARAMS"] for n in results]
    mod_params = []
    for n in results:
        row = n["PARAMS"]
        row.update(n["METADATA"])
        mod_params.append(row)

    # Save individual profiles
    model_hdul, coadd_hdul, bgsub_hdul = fits.HDUList(), fits.HDUList(), fits.HDUList()
    for isolists, hdul in zip([bare_profiles, coadd_profiles, bgsub_profiles], [model_hdul, coadd_hdul, bgsub_hdul]):
        for i, t in enumerate(isolists):
            for col in t.colnames:
                if col not in good_colnames:
                    t.remove_column(col)
            header = dict_to_header(mod_params[i])
            hdul.append(fits.BinTableHDU(data = t, name=f"ISOLIST_{i}", header=header))
            
    model_hdul.writeto(f'{outdirs["MODEL_PROFS"]}{config["RUN_ID"]}_{bin_id}.fits', overwrite=True)
    coadd_hdul.writeto(f'{outdirs["COADD_PROFS"]}{config["RUN_ID"]}_{bin_id}.fits', overwrite=True)
    bgsub_hdul.writeto(f'{outdirs["BGSUB_PROFS"]}{config["RUN_ID"]}_{bin_id}.fits', overwrite=True)

    bare_table = gen_median_table(*gp.bootstrap_median(bare_profiles, dtype="table"))
    coadd_table = gen_median_table(*gp.bootstrap_median(coadd_profiles, dtype="table"))
    bgsub_table = gen_median_table(*gp.bootstrap_median(bare_profiles, dtype="table"))

    gen_median_hdul(bare_table, coadd_table, bgsub_table, 
                    f'{outdirs["MEDIANS"]}{config["RUN_ID"]}_{bin_id}.fits')
    

    # Generate metadata file
    meta_file = gen_model_data_file(mod_params)
    meta_file.write(f'{outdirs["MODEL_PARAMS"]}{config["RUN_ID"]}_{bin_id}_MODEL_PARAMS.fits', 
                    format='fits', overwrite=True)


def gen_model_data_file(param_list):
    params = list(param_list[0].keys())

    rows = []
    for p in param_list:
        rows.append([p[param] for param in params])

    return Table(rows=rows, names=params)


def hdul_to_table(hdul):
    # Convert HDUList to Astropy Table
    return Table.read(hdul, format='fits')




def combine_profile_sets(loc_1, loc_2, run_id_1, run_id_2, outdir):
    """
    Combine profile sets from two directories and compute the median profiles.
    This function reads FITS files from two specified directories, combines the profiles,
    and writes the combined profiles to an output directory. It also computes the median
    profiles using bootstrap resampling.
    Parameters
    ----------
    loc_1 : str
        The path to the first output GalPRIME directory.
    loc_2 : str
        The path to the second output GalPRIME directory.
    run_id_1 : str
        The run identifier prefix for first GalPRIME run.
    run_id_2 : str
        The run identifier prefix for second GalPRIME run.
    outdir : str
        The output directory where combined FITS files will be saved.
    Returns
    -------
    medians : dict
        A dictionary where keys are the FITS file names and values are lists containing
        the semi-major axis (smas), median profile, lower bound, and upper bound.
    """

    files_1 = [f.removeprefix(f'{run_id_1}') for f in os.listdir(loc_1) if f.endswith('.fits')]
    files_2 = [f.removeprefix(f'{run_id_2}') for f in os.listdir(loc_2) if f.endswith('.fits')]
    
    medians = {}

    for f in tqdm(files_1, desc=f"Combining profiles in {loc_1.split('/')[-2]}"):
        if f not in files_2:
           continue
        
        with fits.open(f"{loc_1}/{run_id_1}{f}") as hdul:
            data_1 = hdul[1:]

            with fits.open(f"{loc_2}/{run_id_2}{f}") as hdul:
                data_2 = hdul[1:]
        
                hdul_out_new = fits.HDUList()
                hdul_out_new.append(fits.PrimaryHDU())

                hdul_out_new.extend(data_1)
                hdul_out_new.extend(data_2)
                hdul_out_new.writeto(f"{outdir}/{run_id_1}{f}", overwrite=True)

                profiles = [hdul_to_table(hdul) for hdul in data_1]
                profiles.extend([hdul_to_table(hdul) for hdul in data_2])

        smas, median, low, up = gp.bootstrap_median(profiles, dtype="table")

        medians[f] = [smas, median, low, up]

    return medians


def get_run_id(outdir):
    """ Automatically determine the run IDs for a given directory

    Args:
        outdir (str): GalPRIME output directory

    Raises:
        ValueError: If no run ID is found, the method will raise a ValueError./

    Returns:
        int: The largest run ID found for the given output directory.
    """

    directory = gp.gen_filestructure(outdir, generate=False)
    
    config_dir = directory["ADDL_DATA"]
    run_ids = []
    try:
        for fn in os.scandir(config_dir):
            fn = fn.name
            if not fn.startswith("config"):
                continue
            
            run_id = int(fn.split(".")[0].split("_")[1])
            run_ids.append(run_id)
    except Exception as e:
        print(f"Error getting run IDs: {e}")
    
    if len(run_ids) == 0:
        raise ValueError("Could not find any index in this output directory.")
    else:
        return max(run_ids)


def combine_outputs(filedir_1, filedir_2, outdir, run_id_1=None, run_id_2=None, suffix=""):
    """
    Combine output profiles from two different runs and generate median tables and HDUs.
    Parameters
    ----------
    filedir_1 : str
        Directory path for the first set of files.
    filedir_2 : str
        Directory path for the second set of files.
    run_id_1 : str
        Identifier for the first run.
    run_id_2 : str
        Identifier for the second run.
    outdir : str
        Output directory where the combined results will be saved.
    Returns
    -------
    None
    """

    output_files = gp.gen_filestructure(outdir)

    

    files_1 = gp.gen_filestructure(filedir_1, generate=False)
    if run_id_1 is None:
        run_id_1 = get_run_id(filedir_1)
    
    files_2 = gp.gen_filestructure(filedir_2, generate=False)
    if run_id_2 is None:
        run_id_2 = get_run_id(filedir_2)

    
    for fset in [files_1, files_2]:
        for f in os.listdir(fset["ADDL_DATA"]):
            shutil.copy(f"{fset['ADDL_DATA']}{f}", f"{output_files['ADDL_DATA']}{f}")

    to_combine = ["MODEL_PROFS", "COADD_PROFS", "BGSUB_PROFS"]
    median_sets = []
    for key in to_combine:
        median_sets.append(combine_profile_sets(files_1[key], files_2[key], 
                                                run_id_1, run_id_2, 
                                                output_files[key]))
    
    # Generate median tables and HDUs
    bare_sets, coadd_sets, bgsub_sets = median_sets
    for f in bare_sets.keys():
        bare_table = gp.gen_median_table(*bare_sets[f])
        coadd_table = gp.gen_median_table(*coadd_sets[f])
        bgsub_table = gp.gen_median_table(*bgsub_sets[f])

        gp.gen_median_hdul(bare_table, coadd_table, bgsub_table, 
                           outname=f'{output_files["MEDIANS"]}{f}')
