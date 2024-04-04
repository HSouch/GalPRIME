import numpy as np


def to_sb(y, m_0=27, arcconv=0.168):
    y_new = -2.5 * np.log10(y / (arcconv ** 2)) + m_0
    return y_new