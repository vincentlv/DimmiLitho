"""dimmilitho - litho simulation"""

import pp
from litho.config import CONFIG
from litho.mask import Mask
from pp.rotate import rotate
from skimage.measure import find_contours

__version__ = "0.0.2"
__author__ = "vincentlv <vicentlv>"
__all__ = ["CONFIG", "Mask"]


def smooth(
    gdspathin,
    gdspathout,
    threshold=0.3,
    xmax=300,
    ymax=300,
    gridsize=10,
    layer_in=1,
    layer_out=1,
    scale_factor=1,
    angle=-90,
):
    """ returns a litho simulated version of gdspathin

    Args:
        gdspathin: input GDS
        gdspathout: output GDS
        threshold: from analog to binary, lower threshold makes bigger features
        xmax:
        ymax:
        gridsize: 10
        layer: 1

    Mask Default Parameters
    1. x/y_range means the computing area. Different value are supported
    2. x/y_gridsize the simulated size of the area. Different value are supported
    3. CD infomation is usable for method poly2mask
    """
    m = Mask()
    m.x_range = [-xmax, xmax]
    m.y_range = [-ymax, ymax]
    m.x_gridsize = gridsize
    m.y_gridsize = gridsize
    m.openGDS(gdspathin, layer_in)
    m.maskfft()

    m.smooth()
    mth = m.sdata > threshold
    contours = find_contours(mth, level=threshold)

    # ci = pp.import_gds(gdspathin)
    # xscale = ci.xsize / gridsize
    # scalerate = 0.6428571428571429

    c = pp.Component(f"{gdspathin.stem}_smooth")
    for ci in contours:
        c.add_polygon(ci / gridsize ** 2 / scale_factor, layer=layer_out)

    cr = rotate(c, angle=angle)
    pp.write_gds(cr, gdspathout)
    return cr


def overlay(gdspathin, gdspathout):
    c = pp.Component()
    ci = pp.import_gds(gdspathin)
    co = pp.import_gds(gdspathout)
    cir = c << ci
    cor = c << co
    cir.xmin = 0
    cir.ymin = 0
    cor.xmin = 0
    cor.ymin = 0

    cir.x = 0
    cir.y = 0
    cor.x = 0
    cor.y = 0
    return c


if __name__ == "__main__":
    gdspathin = CONFIG["samples"] / "verniers.gds"
    gdspathout = CONFIG["samples"] / "verniers_out.gds"
    layer = 1
    threshold = 0.35
    ymax = xmax = 600
    gridsize = 1
    scale_factor = 20 * 23 / 7.1
    angle = 90

    # gdspathin = CONFIG["gdslib"] / "AND2_X4.gds"
    # gdspathout = CONFIG["gdslib"] / "AND2_X4_smooth.gds"
    # layer = 10
    # threshold = 0.4
    # ymax = xmax = 600
    # gridsize = 5  # smaller features
    # gridsize = 10
    # scale_factor = 1
    # angle = -90

    c = smooth(
        gdspathin,
        gdspathout,
        threshold=threshold,
        xmax=xmax,
        ymax=ymax,
        gridsize=gridsize,
        scale_factor=scale_factor,
        layer_in=1,
        layer_out=2,
        angle=angle,
    )
    c = overlay(gdspathin, gdspathout)

    # pp.show(gdspathin)
    pp.show(c)

    # import pp
    # pp.show(gdspathout)
