""" litho sim

calculates litho effects
"""

import matplotlib.pyplot as plt
from litho.config import CONFIG
from litho.mask import Mask

"""polygon 2 mask"""
# mp = [ [[-1,6],[-1, 2],[1, 2],[1, 1],[6, 1],[6, 0],[0, 0],[0, 1],[-2, 1],[-2, 6],[-1, 6]], \
#   [[6, -1],[6, -2],[1, -2],[1, -3],[4, -3],[4, -6],[3, -6],[3, -4],[0, -4],[0, -1],[6, -1]] ]
# m = Mask()
# m.x_range = [-300.0,300.0]
# m.y_range = [-300.0,300.0]
# m.x_gridsize = 1.5
# m.y_gridsize = 1.5
# m.CD = 45
# m.polygons = mp
# m.poly2mask()

"""from GDS"""
m = Mask()
m.x_range = [-300.0, 300.0]
m.y_range = [-300.0, 300.0]
m.x_gridsize = 10
m.y_gridsize = 10
m.openGDS(CONFIG["gdslib"] / "AND2_X4.gds", 10)
m.maskfft()

plt.imshow(
    m.data,
    extent=(m.x_range[0], m.x_range[1], m.y_range[0], m.y_range[1]),
    cmap="hot",
    interpolation="none",
)

m.smooth()
plt.imshow(
    m.sdata,
    extent=(m.x_range[0], m.x_range[1], m.y_range[0], m.y_range[1]),
    cmap="hot",
    interpolation="none",
)
plt.show()
