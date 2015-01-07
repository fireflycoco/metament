import numpy as np

def q_ext(m, alpha):
    """Calculate Rayleigh extinction efficiency.
    m: complex refractive index
    alpha: dimensionless size parameter, alpha = pi * Dp / lambda
    """
    m2 = m ** 2
    m4 = m ** 4
    A = (m2 - 1.0) / (m2 + 2.0)
    B = (m4 + 27.0 * m2 + 38.0) / (2.0 * m2 + 3.0)
    res = 4.0 * alpha * np.imag(A * (1.0 + (alpha**2)/ 15.0 * A*B)) + \
            8.0 / 3.0 * (alpha ** 4) * np.real(A ** 2)
    return res

def q_scat(m , alpha):
    """Calculate Rayleigh scattering efficiency.
    m: complex refractive index
    alpha: dimensionless size parameter, alpha = pi * Dp / lambda
    """
    m2 = m ** 2
    A = (m2 - 1.0) / (m2 + 2.0)
    res = 8.0 / 3.0 * (alpha ** 4) * (np.abs(A) ** 2)
    return res
