

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
            self.mag_condition,
            self.reff_condition,
            self.n_condition,
            self.ellip_condition,
        ]
    
    def mag_condition(self, p):
        return p["MAG"] > 0

    def reff_condition(self, p):
        return p["REFF"] > 0

    def n_condition(self, p):
        return 0 < p["N"] < 10

    def ellip_condition(self, p):
        return 0 < p["ELLIP"] < 1

class BulgeDiskVerifier(ParamVerifier):

    def __init__(self):
        super().__init__()
        self.conditions = [
            self.mag_condition,
            self.fbulge_condition,
            self.reff_bulge_condition,
            self.reff_disk_condition,
            self.ellip_disk_condition,
            self.ellip_bulge_condition,
        ]

    def mag_condition(self, p):
        return p["MAG"] > 0

    def fbulge_condition(self, p):
        return 0 < p["FBULGE"] < 1

    def reff_bulge_condition(self, p):
        return p["REFF_BULGE"] > 0

    def reff_disk_condition(self, p):
        return p["REFF_DISK"] > 0

    def ellip_disk_condition(self, p):
        return 0 < p["ELLIP_DISK"] < 1

    def ellip_bulge_condition(self, p):
        return 0 < p["ELLIP_BULGE"] < 1
