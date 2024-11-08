import numpy as np
import galprime as gp

from astropy.table import Table

import time
import argparse
import multiprocessing as mp

import os



parser = argparse.ArgumentParser(description="Run GalPRIME simulation")
parser.add_argument("config_filename", type=str, help="Path to config file")
parser.add_argument("--max_bins", type=int, default=None, help="Maximum number of bins to process")
parser.add_argument("--run_id", type=int, default=gp.get_dt_intlabel(), help="Run ID")
parser.add_argument("--log_level", type=int, default=20, help="Logging level")
parser.add_argument("--verbose", action="store_true", help="Verbose output")
parser.add_argument("--keep_temp", action="store_true", help="Delete temporary files")
args = parser.parse_args()

start_time = time.perf_counter()

config_filename = "myconfig.gprime"


class GPrimeSingle:
    def __init__(self, config, model, params, keys, bg=None, psf=None, logger=None):
        self.config = config
        self.model = model
        self.params = params
        self.keys = keys

        self.bg = bg
        self.psf = psf

        self.logger = logger

    def process(self):

        model_image = self.model.generate()
        
        # Convolve model with psf
        convolved_model = gp.convolve_model(self.model, self.psf)
        print(convolved_model.shape)
        # Add to background

        # Create bg-subtracted image

        # Mask image(s)

        # Extract profiles

        # Save outputs

        pass


def process_single(fn, logger):
    logger.debug(f"Processing {fn}")
    gprime_single = gp.load_object(fn)
    gprime_single.logger = logger
    gprime_single.process()
    


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

    print(f"Starting run ID:{run_id}")
    print(f'Logfile saved to: {config["DIRS"]["OUTDIR"]}output_{run_id}.log')

    binlist = gp.bin_catalogue(table, bin_params=config["BINS"], params=config["KEYS"], logger=logger)
    max_bins = len(binlist.bins) if args.max_bins is None else min(args.max_bins, len(binlist.bins))


    model = gp.galaxy_models[config["MODEL"]["MODEL_TYPE"]]
    logger.info(f"Using model {model.__name__}")


    def process_bin(b):
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

            gprime_single = GPrimeSingle(config, model(), b, config["KEYS"], bg, psf)
            gp.save_object(gprime_single, filename)

            to_process.append(filename)
        
        # Process the gprime single objects (mutiprocessed)
        with mp.Pool(processes=config["NCORES"]) as pool:
            pool.map(process_single,  to_process)

        # Remove temporary files if specified
        if not args.keep_temp:
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

