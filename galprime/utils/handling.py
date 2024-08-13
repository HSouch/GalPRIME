import os
import logging


def setup_logging(run_id, log_level, log_filename=None):
    logger = logging.getLogger(f"RUN_{run_id}")

    os.remove(log_filename)
    logging.basicConfig(filename=log_filename, level=log_level,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    return logger
