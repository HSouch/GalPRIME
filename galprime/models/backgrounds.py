from .. import cutouts


class BackgroundModels(cutouts.Cutouts):

    def __init__(self, **kwargs):
        super().__init__()

    @staticmethod
    def populate():
        raise NotImplementedError("Abstract class")
    

class XGradientBackgrounds(BackgroundModels):
    pass


class YGradientBackgrounds(BackgroundModels):
    pass


class DualGradientBackgrounds(BackgroundModels):
    pass


class RadialGradientBackgrounds(BackgroundModels):
    pass
