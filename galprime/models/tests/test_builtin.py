# Test all built-in GalPRIME models

from.helpers import ModelTestBase
from .. import galaxies


class TestSingleSersicModel(ModelTestBase):
    name = "SingleSersicModel"
    model = galaxies.SingleSersicModel


class TestBulgeDiskModel(ModelTestBase):
    name = "BulgeDiskModel"
    model = galaxies.BulgeDiskSersicModel


class TestExponentialDiskModel(ModelTestBase):
    name = "ExponentialDiskModel"
    model = galaxies.ExponentialDiskModel

    def test_sersic_index(self):
        mod, mod_params = self.model()._generate()

        assert mod_params["N"] == 1


class TestEllipticalGalaxyModel(ModelTestBase):
    name = "EllipticalGalaxyModel"
    model = galaxies.EllipticalGalaxyModel

    def test_sersic_index(self):
        mod, mod_params = self.model()._generate()

        assert mod_params["N"] == 4
