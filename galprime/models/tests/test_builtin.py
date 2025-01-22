# Test all built-in GalPRIME models

import numpy as np
from .helpers import ModelTestBase
from .. import galaxies


class TestSingleSersicModel(ModelTestBase):
    name = "SingleSersicModel"
    model = galaxies.SingleSersicModel


class TestBulgeDiskModel(ModelTestBase):
    name = "BulgeDiskModel"
    model = galaxies.BulgeDiskSersicModel

    # Have to override this test because the model is a composite model
    def test_verify_bad_params(self):
        mod = self.model()
        params = {**mod.defaults, "MAG": -20}
        assert not mod.verifier.verify(params)
        
        params = {**mod.defaults, "REFF_BULGE": -1}
        assert not mod.verifier.verify(params)

        params = {**mod.defaults, "REFF_DISK": -1}
        assert not mod.verifier.verify(params)


class TestExponentialDiskModel(ModelTestBase):
    name = "ExponentialDiskModel"
    model = galaxies.ExponentialDiskModel

    def test_sersic_index(self):
        mod, mod_params = self.model()._generate()

        assert mod_params["N"] == 1


class TestEllipticalModel(ModelTestBase):
    name = "EllipticalGalaxyModel"
    model = galaxies.EllipticalModel

    def test_sersic_index(self):
        mod, mod_params = self.model()._generate()

        assert mod_params["N"] == 4



def test_gen_single_sersic():

    param_set_1 = {"MAG": 20, "REFF": 1, "ELLIP": 0.3, "N": 1}
    param_set_2 = {"MAG": 22, "REFF": 1.5, "ELLIP": 0.3, "N": 4}
    param_set_3 = {"MAG": 20, "REFF": 1, "ELLIP": 1.2, "N": 2, "PA": np.pi / 4}


    for pset in [param_set_1, param_set_2, param_set_3]:
        mod, mod_params = galaxies.gen_single_sersic(**pset)

        for key, val in pset.items():
            assert mod_params[key] == val
        