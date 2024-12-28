import numpy as np

class ModelTestBase:
    """
    ModelTestBase is a base class for testing builtin GalPRIme models.
    Attributes:
        name (str): The name of the model.
        model (class): The model class to be tested.
    Methods:
        test_defaults():
            Tests that the model's defaults attribute is not None and is a dictionary.
        test_default_gen():
            Tests that the model can generate data and parameters, and that the generated data is a 2D numpy array.
        test_different_size():
            Tests that the model can be generated with different sizes and that the generated data has the correct shape.
        test_diff_mags():
            Tests that the model can be generated with different magnitudes, that the generated data is 
            different for different magnitudes, and that each successive model has lower total brightness.
    """

    name = None
    model = None


    def test_defaults(self):
        assert self.model().defaults is not None
        assert isinstance(self.model().defaults, dict)
    
    def test_default_gen(self):
        model_gen, model_params = self.model().generate()

        assert model_gen is not None
        assert isinstance(model_gen, np.ndarray)
        assert len(model_gen.shape) == 2

        assert model_params is not None
        assert isinstance(model_params, dict)
    
    def test_different_size(self):
        # Test that the model can be generated with different sizes
        for size in [151, 251, 301]:
            mod = self.model()
            params = {**mod.defaults, "SHAPE": size}
            model_gen, _ = mod.generate(params=params)
            
            assert model_gen is not None
            assert len(model_gen.shape) == 2
            assert model_gen.shape == (int(size), int(size))

    def test_diff_mags(self):
        # Test that the model can be generated with different magnitudes
        model_mags = []
        for mag in [20, 22, 24]:
            mod = self.model()
            params = {**mod.defaults, "MAG": mag}
            model_gen, _ = mod.generate(params=params)
            
            assert model_gen is not None
            assert len(model_gen.shape) == 2

            model_mags.append(model_gen)
        
        assert not np.allclose(model_mags[0], model_mags[1])
        
        # Check that each successive model is of lower total brightness
        for i in range(1, len(model_mags)):
            assert np.sum(model_mags[i]) < np.sum(model_mags[i-1])
