Running GalPRIME
================


Running from main binary script
-------------------------------

GalPRIME automatically comes with a binary script that can be used to run simulations with a given configuration file.
Once the code is installed, the ``run_galprime`` script should be included in your path and can be called directly from
anywhere in a terminal. 
To execute the script, use the following command in your terminal:

.. code-block:: bash

    run_galprime config.gprime

You can see the available arguments included in the binary script by adding ``-h`` to the command.

Number of objects and multiple sims
"""""""""""""""""""""""""""""""""""

It is recommended to have a maximum of ~500 objects per bin for a given GalPRIME simulation. For users who wish to have
more, they should run multiple individual runs and them combine the outputs together 
(See  :doc:`../2_processing/postprocessing` for info on combining sim outputs).

Extraction of Real Galaxy Light Profiles 
----------------------------------------

GalPRIME also has a pipeline script for the extraction of real profiles. Note that while the configuration file is 
mostly the same for both pipelines, there are some minor differences (please refer to the documentation section on 
config files for more information).

.. code-block:: bash

    gp_extract extract_config.gprime

