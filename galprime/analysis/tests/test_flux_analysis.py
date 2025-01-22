
from .. import fluxes

from numpy import arange
from astropy.modeling.models import Sersic1D


mock_xs = arange(1, 50, 1)
mock_sersic = Sersic1D(amplitude=1, r_eff=1, n=1)(mock_xs)


def test_rX_valid():
    # Check that rX returns a valid index and value

    r10_idx, r10_val = fluxes.rX(mock_xs, mock_sersic, 0.1)
   
    assert r10_idx >= 0
    assert r10_val >= 0

    
def test_rX_proper_math():
    # Check that rX monotonically increases with X

    r10_idx, r10_val = fluxes.rX(mock_xs, mock_sersic, 0.1)
    r50_idx, r50_val = fluxes.rX(mock_xs, mock_sersic, 0.5)
    r90_idx, r90_val = fluxes.rX(mock_xs, mock_sersic, 0.9)

    assert r10_val <= r50_val
    assert r50_val <= r90_val


def test_rX_limits():
    # Check that rX returns the last value for x=1

    r100_idx, r100_val = fluxes.rX(mock_xs, mock_sersic, 1)
    assert r100_idx == len(mock_xs) - 1


def test_r50():
    # Check that r50 is the same as rX with x=0.5

    r50_idx, r50_val = fluxes.rX(mock_xs, mock_sersic, 0.5)
    r50_builtin_idx, r50_builtin_val = fluxes.r50(mock_xs, mock_sersic)

    assert r50_idx == r50_builtin_idx
    assert r50_val == r50_builtin_val