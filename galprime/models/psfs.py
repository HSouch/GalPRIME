from astropy.modeling import models


class PSFModel:
    def __init__(self, params, **kwargs):
        self.params = params

    def generate(self):
        raise NotImplementedError("Abstract class")
    

class SingleGaussianPSF(PSFModel):
    
    def __init__(self, params={}, **kwargs):
        super().__init__(params, **kwargs)

    def generate(self):
        return models.Gaussian2D(**self.params)

