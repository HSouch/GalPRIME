Models
======

**GalPRIME** is primarily designed to work with Sersic models, of the form 

.. math:: 
    I(r) = I_{R_e}\exp\left\{b_n\left[1 - \left(\frac{r}{R_e}\right)^{1/n}\right]\right\}

Where **n** is the *Sersic index* and sets the overall profile shape. For most galaxies, the Sersic index will lie 
between 0.5 and 5. Disk galaxies typically have incides of **n=1**, which returns the typical exponential profile.
Elliptical galaxies typically have indices around **n=4**, which gives higher light buildup in the center and in the 
outskirts of the galaxy.

**GalPRIME** has a set of builtin models that will cover most typical usages of the code. However, a user can define 
their own model, but it will require some additional work to modify the simulation pipeline.

**GalPRIME** models (``galprime.models.galaxies.GalaxyModel``) are essentially wrappers for 2D
`Astropy <https://docs.astropy.org/en/stable/modeling/predef_models2D.html>`_ models. 
These models are designed to be automatically generated using KDEs and therefore require a few additional bells and 
whistles. 

Each model has a directory of **defaults** set by the user. If the user does not modify any paramters for the model, 
each model creation will use the default parameters instead. You can see what defaults exist by:

.. code-block:: python

    import galprime as gp 
    mod = gp.SingleSersicModel()
    mod.defaults

Furthermore, in order for a given ``GalPRIME`` sim to work, the input catalogue is required to have all of the 
parameters specified under the defaults. Setting up a config file properly involves checking what parameters are
requested by the model, and customizing the ``config["KEYS"]`` section of the config file accordingly. For example,
with a Single Sersic model, the user might have the following config section:

.. code-block:: bash

    [KEYS]
        MAG = i
        REFF = R_GIM2D
        N = SERSIC_N_GIM2D
        ELLIP = ELL_GIM2D

Built-In Models
---------------

To change which model ``GalPRIME`` uses, the user must modify the ``config["MODEL"]["MODEL_TYPE"] = index`` parameter in the
configuration file. These are the current builtin models and their respective indices:

#. SingleSersicModel: 1
#. BulgeDiskSersicModel: 2
#. ExponentialDiskModel: 3
#. EllipticalModel: 4


Single-Sersic Model
^^^^^^^^^^^^^^^^^^^

The main model used in ``GalPRIME``. It is a single instance of a Sersic model, requiring keys of magnitude, effecitve
radius in pixels, arcseconds or kpc, Sersic index, and ellipticity, given by

.. math::
    \epsilon = 1 - \frac{b}{a}

where b and a are the semiminor and semimajor axies of a given isophotal ellipse.

.. note::

    If the effective radius unit is given as kpc by specifying ``config["MODEL"]["REFF_UNIT"] = "kpc"``, note that a
    redshift must always be specified to convert the effective radius into units of pixels. It is preferable to 
    make a version of the table with effective radii given in pixels.

Bulge-Disk Model
^^^^^^^^^^^^^^^^

The bulge-disk model is a combination of two Sersic models. It requires a slighlty different set of inputs which are
shown below.

.. code-block:: bash
    
    [KEYS]
        MAG = gg2d
        FBULGE = __B_T_g
        REFF_BULGE = Rd
        N_BULGE = nb
        REFF_DISK = Re
        ELLIP_DISK = ellip_d
        ELLIP_BULGE = ellip_b

These include the TOTAL magnitude and the bulge fraction (i.e. how much of the light belongs to to the bulge between 0 
and 1). All other paramters should be explanatory. 
This setup is designed to work with catalogues such as the bulge-disk 
decomposition work found in `Simard+2011 <https://iopscience.iop.org/article/10.1088/0067-0049/196/1/11/pdf>`_.
**GalPRIME** automatically determines the magnitudes of each component, and will add together two separate Sersic 
models with a shared position angle. 



Exponential Disk and Elliptical Models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
These are Sersic models with a fixed Sersic index. N=1 for exponential disks, and N=4 for elliptical models. Thus,
these parameters are not needed when using these models.