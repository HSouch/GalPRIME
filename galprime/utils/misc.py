from scipy.stats import gaussian_kde

import os

def object_kde(columns):
    """ Generate a gaussian kernel density estimate for the columns of a bin. """
    return gaussian_kde(columns)


def gen_filestructure(outdir):
    os.makedirs(outdir, exist_ok=True)

    file_dict = {"MODEL_PROFS": f"{outdir}/model_profiles/",
                 "COADD_PROFS": f"{outdir}/coadd_profiles/",
                 "BGSUB_PROFS": f"{outdir}/bgsub_profiles/",
                 "MODEL_MEDIANS": f"{outdir}/model_medians/",
                 "COADD_MEDIANS": f"{outdir}/coadd_medians/",
                 "BGSUB_MEDIANS": f"{outdir}/bgsub_medians/",
                 "ADDL_DATA": f"{outdir}/additional_data/",
                 "TEMP": f"{outdir}/tempfiles/"}
    
    for key, value in file_dict.items():
        os.makedirs(value, exist_ok=True)

    return file_dict
