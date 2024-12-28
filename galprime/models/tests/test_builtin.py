# Test all built-in GalPRIME models

from.helpers import ModelTestBase
from .. import galaxies


class TestSingleSersicModel(ModelTestBase):
    name = "SingleSersicModel"
    model = galaxies.SingleSersicModel


class TestBulgeDiskModel(ModelTestBase):
    name = "BulgeDiskModel"
    model = galaxies.BulgeDiskSersicModel