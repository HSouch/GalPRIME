

class ParamVerifier:
    """
    A verifier class that checks if given parameters meet all specified conditions.
    Attributes:
        conditions (list): A list of conditions (functions) that the parameters must satisfy.
    Methods:
        verify(params):
            Checks if all conditions are met for the given parameters.
            Args:
                params: The parameters to be verified.
            Returns:
                bool: True if all conditions are met, False otherwise.
    """
    
    def __init__(self):
        self.conditions = []

    def verify(self, params):
        is_valid = all(condition(params) for condition in self.conditions)
        return is_valid
    
    
class NullVerifier(ParamVerifier):
    def __init__(self):
        self.conditions = []


class DefaultVerifier(ParamVerifier):

    def __init__(self):
        super().__init__()
        self.conditions = [
            lambda p: p["MAG"] > 0,
            lambda p: p["REFF"] > 0,
            lambda p: 0 < p["N"] < 10,
            lambda p: 0 < p["ELLIP"] < 1,
        ]

class BulgeDiskVerifier(ParamVerifier):
    
        def __init__(self):
            super().__init__()
            self.conditions = [
                lambda p: p["MAG"] > 0,
                lambda p: 0 < p["FBULGE"] < 1,
                lambda p: p["REFF_BULGE"] > 0,
                lambda p: p["REFF_DISK"] > 0,
                lambda p: 0 < p["ELLIP_DISK"] < 1,
                lambda p: 0 < p["ELLIP_BULGE"] < 1,
            ]
