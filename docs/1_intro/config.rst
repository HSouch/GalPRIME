Configuration
=============

GalPRIME primarily runs using a config file, which handles all key steps of the automated extraction process.
You can automatically generate 

.. code-block:: python

    import galprime as gp
    gp.dump_default_config_file(outname="config.gprime")




.. note::
    If parameters are missing, they are automatically handled by the CONFIGSPEC. The full config is dumped into the
    output directory as a pickled object under **additional_data/**