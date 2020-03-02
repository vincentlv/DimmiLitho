import numpy as np
import pp
from pp.rotate import rotate


@pp.autoname
def verniers_width(width_min=0.1, width_max=0.5, gap=0.1, size_max=11):
    c = pp.Component()
    y = 0

    widths = np.linspace(width_min, width_max, int(size_max / (width_max + gap)))

    for width in widths:
        w = c << pp.c.waveguide(width=width, length=size_max)
        y += width / 2
        w.y = y
        y += width / 2 + gap

    return c


if __name__ == "__main__":
    c = verniers_width()
    c.flatten()
    pp.write_gds(c, 'verniers.gds')
    pp.show(c)
