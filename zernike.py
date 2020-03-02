"""
Created on Mon May 19 13:30:30 2014
@author: WenLv
"""

from scipy.misc import factorial
from numpy import meshgrid, sqrt, arctan2, zeros, cos, sin, array, where, linspace


def polar_array(T, num):
    """
    Generate polar coordinates.
    """
    x = linspace(-T, T, num)
    X, Y = meshgrid(x, x)
    r = sqrt(X ** 2 + Y ** 2)
    th = arctan2(Y, X)
    return r, th


def rnm(n, m, rho):
    """
    Return an array with the zernike Rnm polynomial calculated at rho points.
    """
    Rnm = zeros(rho.shape)
    S = (n - abs(m)) / 2
    for s in range(0, S + 1):
        CR = (
            pow(-1, s)
            * factorial(n - s)
            / (
                factorial(s)
                * factorial(-s + (n + abs(m)) / 2)
                * factorial(-s + (n - abs(m)) / 2)
            )
        )
        p = CR * pow(rho, n - 2 * s)
        Rnm = Rnm + p
    Rnm[rho > 1.0] = 0
    return Rnm


def zernike(n, m, rho, theta):
    """
    **ARGUMENTS:** 
    
    ===== ==========================================     
    n     n order of the Zernike polynomial
    m     m order of the spatial frequency
    rho   Matrix containing the radial coordinates. 
    theta Matrix containing the angular coordinates.
    ===== ==========================================
 
    """
    Rnm = rnm(n, m, rho)
    Nnm = sqrt(2 * (n + 1) / (1 + where(m == 0, 1, 0)))

    if m >= 0:
        Zmn = Nnm * Rnm * cos(m * theta)
    elif m < 0:
        Zmn = -Nnm * Rnm * sin(m * theta)

    return Zmn


def i2nm(i):
    """
    Return the n and m orders of the i'th zernike polynomial
       
    ========= == == == == == == == == == == == == == == == ===
    i          0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 ...
    n-order    0  1  1  2  2  2  3  3  3  3  4  4  4  4  4 ...
    m-order    0 -1  1 -2  0  2 -3 -1  1  3 -4 -2  0  2  4 ...
    ========= == == == == == == == == == == == == == == == ===
    """
    ia = array(i)
    n = (1 + (sqrt(8 * (ia) + 1) - 3) / 2).astype(int)
    ni = n * (n + 1) / 2
    m = -n + 2 * (i - ni)
    return n, m


def zerniken(i, rho, theta):
    """
    Return the normalized zernike function by the order i
    """
    n, m = i2nm(i)
    return zernike(n, m, rho, theta)
