Models and Generators
=====================


GalPRIME Models
---------------

**GalPRIME** is primarily designed to work with `Sersic models <https://en.wikipedia.org/wiki/S%C3%A9rsic_profile>`_, 
of the form 

.. math:: 
    I(r) = I_{R_e}\exp\left\{b_n\left[1 - \left(\frac{r}{R_e}\right)^{1/n}\right]\right\}

Where **n** is the *Sersic index* and sets the overall profile shape. For most galaxies, the Sersic index will lie 
between 0.5 and 5. Disk galaxies typically have incides of **n=1**, which returns the usual exponential profile.
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
^^^^^^^^^^^^^^^

To change which model ``GalPRIME`` uses, the user must modify the ``config["MODEL"]["MODEL_TYPE"] = index`` parameter in the
configuration file. These are the current builtin models and their respective indices:

#. SingleSersicModel: 1
#. BulgeDiskSersicModel: 2
#. ExponentialDiskModel: 3
#. EllipticalModel: 4


Single-Sersic Model
"""""""""""""""""""

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
""""""""""""""""

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
""""""""""""""""""""""""""""""""""""""
These are Sersic models with a fixed Sersic index. N=1 for exponential disks, and N=4 for elliptical models. Thus,
these parameters are not needed when using these models.


Building Your Own Model 
-----------------------

If the user wishes to build their own model, they are required to create a subclass of the 
``galprime.models.galaxies.GalaxyModel`` class. The subclass then needs to modify the following details. For example,
if we wanted to instantiate a Gaussian model, it could be done using the following setup:

.. code-block:: python
    
    import galprime as gp

    def gen_gaussian(**kwargs):

        shape = kwargs.get("SHAPE", (101, 101))
        if not isinstance(shape, tuple):
            shape = (shape, shape)
        x_0 = kwargs.get("x_0", shape[0] / 2)
        y_0 = kwargs.get("y_0", shape[1] / 2)
        
        ellip = kwargs.get("ELLIP", 0.1)
        x_stddev = kwargs.get("STDDEV", 5)
        y_stddev = x_stddev * (1 - ellip)
        pa = kwargs.get("PA", 0)

        mod = Gaussian2D(amplitude=1, x_mean=kwargs.get("x_0", shape[0] / 2), y_mean=kwargs.get("y_0", shape[1] / 2),
                            x_stddev=x_stddev, y_stddev=y_stddev, theta=pa)
        ys, xs = np.mgrid[:shape[0], :shape[1]]
        z = mod(xs, ys) 
        
        mag, m0 = kwargs.get("MAG", 22), kwargs.get("M0", 27)

        z *= gp.Ltot(mag, m0=m0) / np.sum(z)

        params = {
            "MAG": mag, "M0": m0,
            "STDDEV": mod.x_stddev.value,
            "ELLIP": ellip, "PA":  mod.theta.value,
            "X0": x_0,  "Y0": y_0,
            "SHAPE": shape,
        }
        
        return z, params
    
    class GaussianVerifier(gp.ParamVerifier):
    
        def __init__(self):
            super().__init__()
            self.conditions = [
                self.mag_condition,
                self.stddev_condition,
                self.ellip_condition,
            ]

        def mag_condition(self, p):
            return p["MAG"] > 0
        
        def stddev_condition(self, p):
            return p["STDDEV"] > 0

        def ellip_condition(self, p):
            return 0 <= p["ELLIP"] <= 1

    class GaussianModel(gp.GalaxyModel):

        def __init__(self):       
            self.params = {}
            self.defaults = {
                "MAG": 20.0,
                "STDDEV": 5.0,
                "ELLIP": 0.0,
            }
            self.verifier = GaussianVerifier()
        
        def _generate(self, **params):
            mod, mod_params = gen_gaussian(**params)
            self.params.update(mod_params)
            return mod, mod_params

Here there are three things to consider. The first is the method ``gen_gaussian`` which actually invokes the Astropy
model and handles all keyword arguments. The second class is a ``verifier`` class which will check if the input
parameters are valid for model generation. This is especially helpful for parameters generated from a KDE, which can 
sometimes return unphysical parameter sets (ellipticity < 0, for example). The last class is the actual GaussianModel
class, which defines the default parameters and invokes a unique method called ``_generate``.



Generating Synthetic Images
---------------------------

.. toctree::
    :maxdepth: 1

    generators