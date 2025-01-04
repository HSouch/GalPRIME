Postprocessing
--------------

GalPRIME has a healthy collection of built-in processes to handle outputs. Please refer to the ``notebooks`` section on 
the Github for other uses of the code. 

.. toctree::
    :maxdepth: 1

    medians
    plotting
    

Combining Outputs
^^^^^^^^^^^^^^^^^

Sometimes it helps to run multiple instances of ``run_galprime`` with a smaller number of objects per bin. In this case,
the outputs need to be combined together for plotting and analysis. To do this, you can use the GalPRIME function
``combine_outputs``. Below is an example of combining together two seperate outputs with their respective run IDs:

.. code-block:: python

    import galprime as gp

    output_1 = "test/gprime_1/"
    output_2 = "test/gprime_2/"

    run_id_1, run_id_2 = 1, 2

    gp.combine_outputs(output_1, output_2, run_id_1, run_id_2, "output_combined/"

This will generate unified sets of profiles (in FITS file format) as well as new sets of medians based on the combined
data.

