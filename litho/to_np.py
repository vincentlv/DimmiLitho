from litho import plot
from litho.mask import Mask


def to_np(
    gdspath, pixels_per_um=1, layer=1, boundary=0.16,
):
    """ from GDS returns np.ndarray

    Args:
        gdspath:
        pixels_per_um: needs to be 10X bigger than the minimum feature we want to capture
        layer_in: integer
        boundary: image margin percentage

    Returns:
        np.ndarray

    .. plot::
      :include-source:

      import pp
      import litho

      c = pp.c.verniers()
      c.flatten()
      pp.write_gds(c, gdspath)

      d = gds_to_np(gdspath, pixels_per_um=20, layer_in=1, boundary=0.16)
      litho.plot(d)
      litho.plt.show()

    """
    m = Mask()
    m.openGDS(str(gdspath), layer, boundary=boundary, pixels_per_um=pixels_per_um)
    return m.data.T


def _demo_to_np():
    """ from GDS -> np.ndarray """
    import pp

    gdspath = "gds/verniers.gds"
    c = pp.c.verniers()

    # w = 1000
    # c = pp.c.waveguide(width=w, length=w)
    c = pp.c.circle(radius=10, layer=(1, 0))

    c.flatten()
    pp.write_gds(c, gdspath)

    d = to_np(gdspath, pixels_per_um=20, layer=1, boundary=0.16)
    plot(d)


if __name__ == "__main__":
    _demo_to_np()
