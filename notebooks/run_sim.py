import galprime as gp

import argparse

parser = argparse.ArgumentParser(description="Run GalPRIME simulation")
parser.add_argument("config_filename", type=str, help="Path to config file")
parser.add_argument("--max_bins", type=int, default=None, help="Maximum number of bins to process")
parser.add_argument("--run_id", type=int, default=gp.get_dt_intlabel(), help="Run ID")
args = parser.parse_args()


config_filename = "myconfig.gprime"


sim = gp.GPrime(config_filename, log_level=20, run_id=args.run_id)
run = sim.run(max_bins=args.max_bins)
