"""dimmilitho - litho simulation"""
import pathlib

import gdspy
import pp
from pp.rotate import rotate
from skimage.measure import find_contours

from litho.config import CONFIG
from litho.mask import Mask

__version__ = "0.0.2"
__author__ = "vincentlv <vicentlv>"
__all__ = ["CONFIG", "Mask", "smooth"]


def carve_holes(gdspathin, gdspathout=None, layer=(1, 0)):
    """ takes a GDSpath where holes are polygons defined twice
    returns GDSpath where GDS has carved the holes
    it basically does an XOR

    Args:
        gdspathin:
        gdspathout: if None appends _litho to gdspathin
        layer (tuple): where polygons are
    """
    gdspathin = pathlib.Path(gdspathin)
    gdspathout = gdspathout or gdspathin.with_suffix(".fix.gds")

    g = pp.import_gds(gdspathin)

    # default scale is pm, scale it back up to microns
    si_layer = {"layer": layer[0], "datatype": layer[1]}

    # clear gds cell library for conflicts
    [
        gdspy.current_library.remove(each_cell)
        for each_cell in list(gdspy.current_library.cells)
    ]
    lib = gdspy.GdsLibrary()

    cell = lib.new_cell(gdspathin.stem + "_litho")

    # declare empty lithography polygon instance
    litho_polyinst = gdspy.Polygon([], **si_layer)
    litho_polypts = g.get_polygons()
    for i in range(len(litho_polypts)):
        # scale points back to microns
        temp_pts = litho_polypts[i]
        temp_pts = temp_pts.tolist()
        if len(litho_polypts[i]) > 2:
            litho_polyinst = gdspy.boolean(
                litho_polyinst,
                gdspy.Polygon(temp_pts, layer=layer[0], datatype=layer[1]),
                operation="xor",
                layer=layer[0],
                datatype=layer[1],
            )

    cell.add(litho_polyinst)

    lib.write_gds(str(gdspathout))
    return gdspathout


def smooth(
    gdspathin,
    gdspathout=None,
    threshold=0.3,
    xmax=300,
    ymax=300,
    gridsize=1,
    layer_in=1,
    layer_out=1,
    scalerate=45 / 70,
    angle=0,
    remove_layers=None,
):
    """ returns a litho simulated Component of gdspathin

    Make sure your gridsize is in the order of the feature that you want to calculate.
    If you choose a gridsize too small it will take long time
    If your gridsize is too big for your features the litho simulated will not make any sense

    Args:
        gdspathin: input GDS
        gdspathout: output GDS
        threshold: from analog to binary, lower threshold makes bigger features
        xmax:
        ymax:
        gridsize (nm): resolution.
        layer_in: 1
        layer_out: 1
        scalerate: 45/70
        angle (deg): rotation
        remove_layers: list of tuples of layers to remove before the smooth operation

    .. plot::
        :include-source:

        import pp
        from litho import CONFIG
        from litho import smooth
        from litho import overlay
        from litho.samples.verniers import verniers

        gdspathin = CONFIG["samples"] / "verniers.gds"
        gdspathout = CONFIG["samples"] / "verniers_out.gds"

        layer = 1
        threshold = 0.35
        ymax = xmax = 600
        gridsize = 1
        scalerate = 45/70
        angle = 90
        remove_layers = [(111,0)]

        c = verniers()
        c.flatten()
        pp.write_gds(c, gdspathin)

        c = smooth(
            gdspathin,
            gdspathout,
            threshold=threshold,
            xmax=xmax,
            ymax=ymax,
            gridsize=gridsize,
            scalerate=scalerate,
            layer_in=layer,
            layer_out=2,
            angle=angle,
            remove_layers=remove_layers
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
    m.openGDS(str(gdspathin), layer_in, boundary=0.16, scalerate=scalerate)
    m.maskfft()

    m.smooth()
    mth = m.sdata > threshold
    contours = find_contours(mth.T, level=threshold)

    c = pp.Component(f"{gdspathin.stem}_smooth")
    for ci in contours:
        c.add_polygon(ci / scalerate / 100 * gridsize, layer=layer_out)

    cr = rotate(c, angle=angle)
    cr.x = 0
    cr.y = 0
    pp.write_gds(cr, gdspathout)
    return cr


def np_to_gds(ndarray, gdspath, **kwargs):
    """ returns a gdspath from a numpy array """
    c = np_to_component(ndarray, **kwargs)
    pp.write_gds(c, gdspath)


def np_to_component(
    ndarray,
    threshold=0.99,
    scalerate=45 / 70,
    gridsize=1,
    layer_out=(2, 0),
    name="litho",
):
    """ returns a gdsfactory component from a numpy array """
    c = pp.Component(name)

    for contour in find_contours(ndarray, level=threshold):
        c.add_polygon(contour / scalerate / 100 * gridsize, layer=layer_out)
    return c


def gds_to_numpy_array(
    gdspathin,
    gridsize=1,
    scalerate=45 / 70,
    layer_in=1,
    xmax=300,
    ymax=300,
    boundary=0.16,
):
    """ rasterizes a GDS and returns a numpy array representation of the image

    Args:
        gdspathin
        gridsize
        scalerate
        boundary: how much extra space 16% adds to the image margin

    Returns:
        numpy array
    """
    m = Mask()
    m.x_range = [-xmax, xmax]
    m.y_range = [-ymax, ymax]
    m.x_gridsize = gridsize
    m.y_gridsize = gridsize
    m.openGDS(str(gdspathin), layer_in, boundary=boundary, scalerate=scalerate)
    return m.data.T


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
    """ the gate has small features
    therefore a gridsize of 1nm
    """
    gdspathin = CONFIG["gdslib"] / "AND2_X4.gds"
    gdspathin2 = CONFIG["gdslib"] / "AND2_X4_simpler.gds"
    gdspathout = CONFIG["gdslib"] / "AND2_X4_smooth.gds"

    c = pp.import_gds(gdspathin)
    c = c.remove_layers([1, 2, 3, 4, 5, 9, 11, (63, 63), 235])
    pp.write_gds(c, gdspathin2)

    layer = 10
    threshold = 0.4

    c = smooth(
        gdspathin2, gdspathout, threshold=threshold, layer_in=layer, layer_out=2,
    )
    c = overlay(gdspathin2, gdspathout)
    pp.show(c)


def _demo_verniers():
    from litho.samples.verniers import verniers

    gdspathin = CONFIG["samples"] / "verniers.gds"
    gdspathout = CONFIG["samples"] / "verniers_out.gds"
    layer = 1
    threshold = 0.35
    ymax = xmax = 600
    gridsize = 1
    c = verniers()
    c.flatten()
    pp.write_gds(c, gdspathin)

    c = smooth(
        gdspathin,
        gdspathout,
        threshold=threshold,
        xmax=xmax,
        ymax=ymax,
        gridsize=gridsize,
        layer_in=layer,
        layer_out=2,
    )
    c = overlay(gdspathin, gdspathout)
    pp.show(c)


def _demo_gds_to_numpy_array():
    """ from GDS -> np.ndarray """
    import matplotlib.pyplot as plt

    gdspathin = CONFIG["samples"] / "verniers.gds"
    d = gds_to_numpy_array(gdspathin)
    plt.imshow(d.T)
    plt.show()


def _demo_gds_np_gds():
    """ from GDS -> np.ndarray -> GDS """
    import matplotlib.pyplot as plt

    gdspathin = CONFIG["samples"] / "verniers.gds"
    gdspathout = CONFIG["samples"] / "verniers2.gds"
    scalerate = 0.5

    d = gds_to_numpy_array(gdspathin, scalerate=scalerate, boundary=0.1)
    np_to_gds(d, gdspathout, scalerate=scalerate, threshold=0.99)
    c = overlay(gdspathin, gdspathout)
    pp.show(c)


if __name__ == "__main__":
    # _demo_gate()
    # _demo_verniers()
    # _demo_gds_to_numpy_array()
    _demo_gds_np_gds()
