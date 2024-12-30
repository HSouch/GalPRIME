Running GalPRIME
================


Running from main binary script
-------------------------------

GalPRIME automatically comes with a binary script that can be used to run simulations with a given configuration file. 
To execute the script, use the following command in your terminal:

.. code-block:: bash

    run_galprime config.gprime




Extraction of Real Galaxy Light Profiles 
----------------------------------------

GalPRIME also has a pipeline script for the extraction of real profiles. Note that while the configuration file is 
mostly the same for both pipelines, there are some minor differences (please refer to the documentation section on 
config files for more information).

.. code-block:: bash

    gp_extract extract_config.gprime

