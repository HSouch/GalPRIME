import galprime as gp

from astropy.table import Table
from astropy.io import fits


good_colnames = ['sma', 'intens', 'intens_err', 'ellipticity', 'ellipticity_err', 'pa', 'pa_err',
                 'x0', 'x0_err', 'y0', 'y0_err', 'ndata', 'nflag', 'niter', 'stop_code']


def handle_output(results, outdirs, config, bin_id="0"):
    
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

    # Save medians
    bare_smas, bare_median, bare_low, bare_up = gp.bootstrap_median(bare_profiles)
    coadd_smas, coadd_median, coadd_low, coadd_up = gp.bootstrap_median(coadd_profiles)
    bgsub_smas, bgsub_median, bgsub_low, bgsub_up = gp.bootstrap_median(bgsub_profiles)

    colnames = ["SMA", "MEDIAN", "LOW_1SIG", "LOW_2SIG", "LOW_3SIG", "UP_1SIG", "UP_2SIG", "UP_3SIG"]
    bare_table = Table([bare_smas, bare_median, *bare_low, *bare_up], names=colnames)
    coadd_table = Table([coadd_smas, coadd_median, *coadd_low, *coadd_up], names=colnames)
    bgsub_table = Table([bgsub_smas, bgsub_median, *bgsub_low, *bgsub_up], names=colnames)

    median_hdul = fits.HDUList()
    for table, name in zip([bare_table, coadd_table, bgsub_table], ["BARE", "COADD", "BGSUB"]):
        median_hdul.append(fits.BinTableHDU(data=table, name=name))
    median_hdul.writeto(f'{outdirs["MEDIANS"]}{config["RUN_ID"]}_{bin_id}.fits', overwrite=True)
    