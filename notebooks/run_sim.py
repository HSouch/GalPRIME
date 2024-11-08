import numpy as np
import galprime as gp

from astropy.table import Table

import time
import argparse
import multiprocessing as mp

import os

from joblib import Parallel, delayed

import pebble


global logger

parser = argparse.ArgumentParser(description="Run GalPRIME simulation")
parser.add_argument("config_filename", type=str, help="Path to config file")
parser.add_argument("--max_bins", type=int, default=None, help="Maximum number of bins to process")
parser.add_argument("--run_id", type=int, default=gp.get_dt_intlabel(), help="Run ID")
parser.add_argument("--log_level", type=int, default=20, help="Logging level. 10=DEBUG, 20=INFO, 30=WARNING")
parser.add_argument("--verbose", action="store_true", help="Verbose output")
parser.add_argument("--keep_temp", action="store_true", help="Delete temporary files")
args = parser.parse_args()

start_time = time.perf_counter()

config_filename = "myconfig.gprime"


class GPrimeSingle:
    def __init__(self, config, model, params, bg=None, psf=None, logger=None):
        self.config = config
        self.model = model
        self.params = params

        self.bg = bg
        self.psf = psf

        self.logger = logger

    def process(self):

        model_image, model_params = self.model.generate(self.params)
        # Convolve model with psf
        convolved_model = gp.convolve_model(model_image, self.psf)
        
        # Add to background
        bg_added_model = convolved_model + self.bg
        
        # Create bg-subtracted image
        source_mask, background = gp.estimate_background_2D(bg_added_model, self.config)
        
        background = background.background
        bgsub = bg_added_model - background

        # Mask image(s)
        mask_bgadded, mask_data_bgadded = gp.gen_mask(bg_added_model, config=self.config)
        mask_bgsub, mask_data_bgsub = gp.gen_mask(bgsub, config=self.config)

        # Extract profiles
        for dataset in [convolved_model, np.ma.array(bg_added_model, mask=mask_bgadded), np.ma.array(bgsub, mask=mask_bgsub)]:
            gp.isophote_fitting(dataset, self.config)
        
        # Save outputs

        


def process_single(fn):
    logger.debug(f"Processing {fn}")
    gprime_single = gp.load_object(fn)
    gprime_single.logger = logger
    gprime_single.process()

    return gprime_single
    



if __name__ == '__main__':
    config = gp.read_config_file(config_filename)
    outfiles = gp.gen_filestructure(config["DIRS"]["OUTDIR"])

    run_id = args.run_id

    logger = gp.setup_logging(run_id, args.log_level, 
                            log_filename=f'{config["DIRS"]["OUTDIR"]}output_{run_id}.log')

    logger.info(f"Starting run ID:{run_id}, GalPRIME Version: {gp.__version__}", )

    # Load in backgrounds, PSFS, and object catalogue
    bgs = gp.Cutouts.from_file(f'{config["FILE_DIR"]}{config["FILES"]["BACKGROUNDS"]}', logger=logger)
    psfs = gp.Cutouts.from_file(f'{config["FILE_DIR"]}{config["FILES"]["PSFS"]}', logger=logger)
    psfs.get_ra_dec(ra_key=config["PSFS"]["PSF_RA"], dec_key=config["PSFS"]["PSF_DEC"])

    table = Table.read(f'{config["FILE_DIR"]}{config["FILES"]["CATALOGUE"]}')
    table = gp.trim_table(table, config)
    logger.info(f'Loaded catalogue with {len(table)} entries')

    print(f"GalPRIME v{gp.__version__} -- Starting run ID:{run_id}")
    print(f'Logfile saved to: {config["DIRS"]["OUTDIR"]}output_{run_id}.log')

    binlist = gp.bin_catalogue(table, bin_params=config["BINS"], params=config["KEYS"], logger=logger)
    max_bins = len(binlist.bins) if args.max_bins is None else min(args.max_bins, len(binlist.bins))


    model = gp.galaxy_models[config["MODEL"]["MODEL_TYPE"]]
    logger.info(f"Using model {model.__name__}")


    def process_bin(b):
        cores = config["NCORES"]
        n_objects = config["MODEL"]["N_MODELS"]
        bg_indices = np.random.randint(0, len(bgs.cutouts), n_objects)
        psf_indices = np.random.randint(0, len(psfs.cutouts), n_objects)  
        
        model_template = model()

        keys, kde = gp.setup_kde(model_template, config, b.objects)

        # Create and pickle the gprime single objects
        to_process = []
        for i in range(n_objects):
            filename = f'{outfiles["TEMP"]}{run_id}_{b.bin_id()}_{i}.pkl'

            bg = bgs.cutouts[bg_indices[i]]
            psf = psfs.cutouts[psf_indices[i]]

            params = gp.sample_kde(config, keys, kde)
            params = gp.update_required(params, config)

            gprime_single = GPrimeSingle(config, model(), params, bg=bg, psf=psf)
            gp.save_object(gprime_single, filename)
 
            to_process.append(filename)
        
        iterator = None
        # Process the gprime single objects (mutiprocessed)
        with pebble.ProcessPool(max_workers=cores) as pool:
            future = pool.map(process_single, to_process, timeout=15)
            try:
                iterator = future.result()
            except pebble.TimeoutError as error:
                logger.warning(f"Timeout error: {error}")
            except pebble.ProcessExpired as error:
                logger.warning(f"Process expired: {error}")

        # Remove temporary files if specified
        if not args.keep_temp:
            logger.info(f"Removing temporary files for bin {b.bin_id()}")
            for fn in to_process:
                os.remove(fn)

        

    # Go through the bins and process them
    for i in range(max_bins):
        b = binlist.bins[i]
        logger.info(f'Processing bin {b.bin_id()}')
        # Process the bin
        try:
            process_bin(b)
        except Exception as e:
            logger.error(f"Critical error processing bin {b.bin_id()}: {e}")
            continue

        logger.info(f"Finished processing bin {b.bin_id()}")

    time_elapsed = time.perf_counter() - start_time
    unit = "seconds" if time_elapsed < 60 else "minutes"
    if time_elapsed > 60:
        time_elapsed /= 60

    


    logger.info("Finished processing all bins")
    logger.info(f"Run complete: Time Elapsed: {time.perf_counter() - start_time:.2f} {unit}")

    print(f'Run complete: Time Elapsed: {time.perf_counter() - start_time:.2f} seconds')
    print("Results saved to: ", config["DIRS"]["OUTDIR"])

