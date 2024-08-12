
from .. import galaxies
from ...config import default_config

class ModelTestBase:
    name = "GalaxyModel"
    model = None
    test_params = {}
    config = default_config()

    def test_generate(self):
        mod = self.model
        mod.config = self.config
        mod.params = self.test_params

        model = mod.generate()    
        

class TestSingleSersicModel(ModelTestBase):
    name = "SingleSersicModel"
    model = galaxies.SingleSersicModel()
    
    test_params = {
        "mag": 20,
        "r50": 5,
        "n": 2,
        "ellip": 0.5
    }

    def test_generate(self):
        assert True

    def test_verify_params(self):
        assert True

    def test_gen_multiple(self):
        assert True


