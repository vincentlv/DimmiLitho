import pathlib

import numpy as np
import pp
import scipy.signal as sg


def smooth(ndarray, pixels_kernel=21):
    """ returns a 2D convolution with a square

    .. plot::
        :include-source:

        import pp
        import litho

        c = pp.c.verniers()
        c.flatten()
        gdspath = pp.write_gds(c, gdspath)

        ndarray = litho.gds_to_np(gdspath)
        nf = smooth(ndarray)
        pp.plotgds(gdspathout)
    """
    xx = np.linspace(-1, 1, pixels_kernel)
    X, Y = np.meshgrid(xx, xx)
    R = X ** 2 + Y ** 2
    G = np.exp(-10 * R)
    return sg.convolve2d(0.9 * ndarray + 0.05, G, "same") / np.sum(G)


def overlay(gdspathin, gdspathout):
    c = pp.Component()

    ci = pp.import_gds(gdspathin)
    co = pp.import_gds(gdspathout)

    cir = c << ci
    cor = c << co

    cir.x = 0
    cir.y = 0
    cor.x = 0
    cor.y = 0
    return c


def _demo():
    """ overlay nominal and smooth designs """
    import litho

    gdspath1 = "gds/gds1.gds"
    gdspath2 = "gds/gds2.gds"

    layer = (1, 0)

    c = pp.c.verniers()
    c = pp.c.ring(layer=layer)
    c = pp.c.dbr(n=5)
    c.flatten()
    pp.write_gds(c, gdspath1)

    layer = 1
    threshold = 0.5
    pixels_per_um = 100
    pixels_kernel = int(np.ceil(pixels_per_um * 100e-3))  # 50nm = 50e-3

    ndarray = litho.to_np(gdspath1, layer=layer, pixels_per_um=pixels_per_um)
    ns = smooth(ndarray=ndarray, pixels_kernel=pixels_kernel)

    # litho.plot(ns)
    litho.to_gds(ns, gdspath=gdspath2, threshold=threshold, pixels_per_um=pixels_per_um)
    c = overlay(gdspath1, gdspath2)
    pp.show(c)


if __name__ == "__main__":
    _demo()
