import pathlib

import gdspy
import pp
from skimage.measure import find_contours

from litho.plot import plot
from litho.to_np import to_np


def carve_holes(gdspathin, gdspathout=None, layer=(1, 0)):
    """ takes a GDSpath where holes are polygons defined twice
    returns GDSpath where GDS has carved the holes
    it basically does an XOR

    Args:
        gdspathin:
        gdspathout: if None uses gdspathin
        layer (tuple): where polygons are
    """
    gdspathin = pathlib.Path(gdspathin)
    gdspathout = gdspathout or gdspathin

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


def to_gds(
    ndarray,
    gdspath,
    threshold=0.99,
    pixels_per_um=1,
    layer=(2, 0),
    component_name="litho",
):
    """ saves a numpy array into GDS

    Args:
        ndarray:
        gdspath:
        threshold: 0.99
        pixels_per_um: 1
        layer_out: (2, 0)
        name: "litho"
    """
    c = pp.Component(component_name)

    for contour in find_contours(ndarray, level=threshold):
        c.add_polygon(contour / pixels_per_um, layer=layer)

    pp.write_gds(c, gdspath)
    carve_holes(gdspathin=gdspath, layer=layer)
    return c


def component_to_np(c, pixels_per_um=1, layer_in=1, boundary=0.16):
    """ returns a numpy array from a component
    """
    layers = c.get_layers()
    c.remove_layers(layers - {layer_in})
    c._bb_valid = False
    pp.show(c)
    xsize, ysize = c.size
    gdspath = pp.write_gds(c)

    return to_np(
        gdspath, pixels_per_um=pixels_per_um, layer_in=layer_in, boundary=boundary,
    )


def _demo_component_to_np():
    """ from gdsfactory component -> np.ndarray
    """
    c = pp.c.waveguide(width=1, length=10)
    c = pp.c.dbr()
    c.flatten()
    d = component_to_np(c, pixels_per_um=1, layer_in=(1, 0))
    plot(d)


def _demo_to_gds():
    """ from GDS -> np.ndarray -> GDS """
    """
    for verniers you need at least 20 pixels_per_um
    circle is hard to render with no loss
    """
    from litho.to_np import to_np
    import pp

    layer = (1, 0)
    gdspath = "gds/verniers.gds"
    c = pp.c.verniers()

    # w = 1000
    # c = pp.c.waveguide(width=w, length=w)
    c = pp.c.circle(radius=10, layer=layer)
    c = pp.c.ring(radius=10, layer=layer)

    c.flatten()
    pp.write_gds(c, gdspath)

    pixels_per_um = 20

    ndarray = to_np(gdspath, pixels_per_um=pixels_per_um)
    c2 = to_gds(
        ndarray=ndarray, gdspath=gdspath, layer=layer, pixels_per_um=pixels_per_um
    )
    pp.show(c2)


if __name__ == "__main__":
    _demo_to_gds()
    # _demo_gds_np_gds()
    # _demo_component_to_np()
