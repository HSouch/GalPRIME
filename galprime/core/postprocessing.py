import galprime as gp

from astropy.table import Table
from astropy.io import fits

from tqdm import tqdm
import os


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
    
    bare_profiles = [n["ISOLISTS"][0] for n in results]
    coadd_profiles = [n["ISOLISTS"][1] for n in results]
    bgsub_profiles = [n["ISOLISTS"][2] for n in results]

    # Save individual profiles
    model_hdul, coadd_hdul, bgsub_hdul = fits.HDUList(), fits.HDUList(), fits.HDUList()
    for isolists, hdul in zip([bare_profiles, coadd_profiles, bgsub_profiles], [model_hdul, coadd_hdul, bgsub_hdul]):
        for i in range(len(isolists)):
            t = isolists[i].to_table()
            for col in t.colnames:
                if col not in good_colnames:
                    t.remove_column(col)
            hdul.append(fits.BinTableHDU(data = t, name=f"ISOLIST_{i}"))
    model_hdul.writeto(f'{outdirs["MODEL_PROFS"]}{config["RUN_ID"]}_{bin_id}.fits', overwrite=True)
    coadd_hdul.writeto(f'{outdirs["COADD_PROFS"]}{config["RUN_ID"]}_{bin_id}.fits', overwrite=True)
    bgsub_hdul.writeto(f'{outdirs["BGSUB_PROFS"]}{config["RUN_ID"]}_{bin_id}.fits', overwrite=True)

    bare_table = gen_median_table(*gp.bootstrap_median(bare_profiles))
    coadd_table = gen_median_table(*gp.bootstrap_median(coadd_profiles))
    bgsub_table = gen_median_table(*gp.bootstrap_median(bare_profiles))

    gen_median_hdul(bare_table, coadd_table, bgsub_table, f'{outdirs["MEDIANS"]}{config["RUN_ID"]}_{bin_id}.fits')


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
                hdul_out_new.writeto(f"{outdir}/{f}", overwrite=True)

                profiles = [hdul_to_table(hdul) for hdul in data_1]
                profiles.extend([hdul_to_table(hdul) for hdul in data_2])

        smas, median, low, up = gp.bootstrap_median(profiles, dtype="table")

        medians[f] = [smas, median, low, up]

    return medians


def combine_outputs(filedir_1, filedir_2, run_id_1, run_id_2, outdir):
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
    files_2 = gp.gen_filestructure(filedir_2, generate=False)
    
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
