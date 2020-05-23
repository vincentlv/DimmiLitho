"""dimmilitho - litho simulation"""
import pathlib

import pp
from pp.rotate import rotate
from skimage.measure import find_contours

from litho.config import CONFIG
from litho.mask import Mask

__version__ = "0.0.2"
__author__ = "vincentlv <vicentlv>"
__all__ = ["CONFIG", "Mask", "smooth"]


def smooth(
    gdspathin,
    gdspathout=None,
    threshold=0.3,
    xmax=300,
    ymax=300,
    gridsize=10,
    layer_in=1,
    layer_out=1,
    scale_factor=65.5,
    angle=-90,
):
    """ returns a litho simulated Component of gdspathin

    Args:
        gdspathin: input GDS
        gdspathout: output GDS
        threshold: from analog to binary, lower threshold makes bigger features
        xmax:
        ymax:
        gridsize: 10
        layer_in: 1
        layer_out: 1
        scale_factor: 1
        angle: in degrees -90

    .. plot::
        :include-source:

        import pp
        from litho.samples.verniers import verniers

        c = verniers()
        c.x = 0
        c.y = 0
        pp.plotgds(c)


    .. plot::
        :include-source:

        import pp
        from litho import CONFIG
        from litho import smooth
        from litho import overlay

        gdspathin = CONFIG["samples"] / "verniers.gds"
        gdspathout = CONFIG["samples"] / "verniers_out.gds"

        layer = 1
        threshold = 0.35
        ymax = xmax = 600
        gridsize = 1
        scale_factor = 20 * 23 / 7.1
        angle = 90

        c = smooth(
            gdspathin,
            gdspathout,
            threshold=threshold,
            xmax=xmax,
            ymax=ymax,
            gridsize=gridsize,
            scale_factor=scale_factor,
            layer_in=layer,
            layer_out=2,
            angle=angle,
        )

        pp.plotgds(gdspathout)

    """
    gdspathin = pathlib.Path(gdspathin)
    gdspathout = gdspathout or gdspathin.with_suffix(".litho.gds")
    m = Mask()
    m.x_range = [-xmax, xmax]
    m.y_range = [-ymax, ymax]
    m.x_gridsize = gridsize
    m.y_gridsize = gridsize
    m.openGDS(str(gdspathin), layer_in)
    m.maskfft()

    m.smooth()
    mth = m.sdata > threshold
    contours = find_contours(mth, level=threshold)

    # ci = pp.import_gds(gdspathin)
    # xscale = ci.xsize / gridsize
    # scalerate = 0.6428571428571429

    c = pp.Component(f"{gdspathin.stem}_smooth")
    for ci in contours:
        c.add_polygon(ci * gridsize / scale_factor, layer=layer_out)

    cr = rotate(c, angle=angle)
    cr.x = 0
    cr.y = 0
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


def _demo_gate():
    gdspathin = CONFIG["gdslib"] / "AND2_X4.gds"
    gdspathout = CONFIG["gdslib"] / "AND2_X4_smooth.gds"
    layer = 10
    threshold = 0.4
    ymax = xmax = 600
    gridsize = 5  # smaller features
    gridsize = 10
    angle = -90

    c = smooth(
        gdspathin,
        gdspathout,
        threshold=threshold,
        xmax=xmax,
        ymax=ymax,
        gridsize=gridsize,
        layer_in=layer,
        layer_out=2,
        angle=angle,
    )
    c = overlay(gdspathin, gdspathout)

    # pp.show(gdspathin)
    pp.show(c)


def _remap_gate():
    gdspathin = CONFIG["gdslib"] / "AND2_X4.gds"
    ci = pp.import_gds(gdspathin)
    ci = ci.remove_layers([1, 2, 3, 4, 5, 9, 11, (63, 63), 235])
    pp.write_gds(ci, "and2.gds")


if __name__ == "__main__":
    gdspathin = CONFIG["samples"] / "verniers.gds"
    gdspathout = CONFIG["samples"] / "verniers_out.gds"
    layer = 1
    threshold = 0.35
    ymax = xmax = 600
    gridsize = 1
    angle = 90

    # gdspathin = CONFIG["samples"] / "and2.gds"
    # gdspathout = CONFIG["gdslib"] / "and2_smooth.gds"
    # layer = 10
    # threshold = 0.4
    # ymax = xmax = 600
    # gridsize = 5  # smaller features
    # gridsize = 10
    # angle = -90

    c = smooth(
        gdspathin,
        gdspathout,
        threshold=threshold,
        xmax=xmax,
        ymax=ymax,
        gridsize=gridsize,
        layer_in=layer,
        layer_out=2,
        angle=angle,
    )
    c = overlay(gdspathin, gdspathout)

    # pp.show(gdspathin)
    pp.show(c)
    pp.plotgds(c)

    # import pp
    # pp.show(gdspathout)
