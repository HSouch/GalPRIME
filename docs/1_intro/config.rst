Config
======

GalPRIME primarily runs using a config file, which handles all key steps of the automated extraction process.
You can automatically generate 

.. code-block:: python

    import galprime as gp
    gp.dump_default_config_file(outname="config.gprime")



.. note::
    If parameters are missing, they are automatically handled by the CONFIGSPEC. The full config is dumped into the
    output directory as a pickled object under **additional_data/**


Config Parameters
-----------------

In the following subsection we describe some of the paramters you will find in the default config.
Please note that the format in the actual ConfigObj is a nested dict, so FILES.CATALOGUE in the
configobj *config* is actually ``config["FILES"]["CATALOGUE"]``.

* ``FILES.CATALOGUE``: The FITS catalogue containing necessary structural parameters for binning / modeling.
* ``FILES.PSFS``: Filename of a FITS datacube containing PSFS for convolving input models.
* ``FILES.BACKGROUNDS``: Filename of a FITS datacube containing the input backgrounds to add models to.
* ``FILES.MAG_CATALOGUE``: Filename for a separate catalogue containing a set of magnitudes. Can be None.

* ``DIRS.OUTDIR``: The output directory that everything is saved to (default is ``gprime_out/``).

KEYS
^^^^

The keys section contains all the necessary keys for model gen. ``GalPRIME`` cuts the input catalogue 
down to just these keys (as well as the columns used in binning) while running to save on memory. The
format for this section is ``["KEYS"]["KEY"] = "colname"`` where ``KEY`` is the structural parameter
needed by GalPRIME, and ``colname`` is the corresponding column in the input catalogue. For example,
the Sersic index might be passed into ``GalPRIME`` by including:

.. code-block:: python 

    config["KEYS"]["N"] = "SERSIC_N_GIM2D"

Where in this case, Sersic indices were measured using GIM2D.


BINS 
^^^^

The bins section is where the user specifies how to bin the input catalogue. The format is as follows:

.. code-block:: python

    config["BINS"]["COLNAME"] = [...]

Where ``COLNAME`` is the column in the input catalogue to bin by, and the array contains a list of
bin edges set by the user. For example, to bin by log stellar mass (``MASS_MED``) in steps of 0.5 
from 10 to 11.5, and then by redshift (``Z_BEST``) in steps of 0.2 from 0.1 to 0.9,
the user would include the following:

.. code-block:: python

    config["BINS"]["MASS_MED"] = [10, 10.5, 11, 11.5]
    config["BINS"]["Z_BEST"] = [0.1, 0.3, 0.5, 0.7, 0.9]

which will create 3x4=12 bins of mass and redshift.


MODEL
^^^^^

The MODEL section is where the user sets the size of the model, the number of models to generate per 
bin, as well as the type of model being generated. Other important structural parameters, including 
the zero-point magnitude and the pixel scale (arcseconds per pixel) are specified here.
For more information on the models available in ``GalPRIME``, please see :doc:`../2_processing/models`.

