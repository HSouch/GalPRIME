.. GalPRIME documentation master file, created by
   sphinx-quickstart on Thu Apr  4 11:36:50 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

GalPRIME
========

.. raw:: html

   <img src="_static/gprime_logo.png" width="90%"  style="margin-bottom: 32px;"/>

GalPRIME is a code designed to automate the process of 1D Surface Brightness Profile extraction. It has tools to
automatically extract real profiles from a set of images (given an input catalogue). Its primary purpose, however,
is to determine the efficiency of profile extraction by injecting mock images into either real or simulated 
background images.

GalPRIME operates by generating N-dimensional KDEs from a set of input parameters (typically a galaxy catalogue with 
requisite data columns to construct a morphological model. These KDEs are often generated from binned data, which is
fully customizable by the user. 

For example, you can use GalPRIME to:

#. Test profile extraction performance for typical galaxies by creating mass and redshift bins.
#. Test background subtraction performance by running GalPRIME with bins of size, ellipticity, and Sersic index.
#. Extract profiles of real galaxies from a set of images, and an input catalogue.
#. Test profile extraction performance for disk (n ~ 1) vs elliptical (n ~ 4) galaxies.

.. toctree::
   :maxdepth: 1

   1_intro/installation
   1_intro/config 

.. toctree::
   :maxdepth: 2

   2_processing/running
   2_processing/models


.. toctree::
   :maxdepth: 2

   2_processing/postprocessing

.. toctree::
   :maxdepth: 1

   3_misc/papers


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
