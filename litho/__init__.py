"""dimmilitho - litho simulation"""

from litho.config import CONFIG
from litho.mask import Mask

__version__ = "0.0.1"
__author__ = "vincentlv <vicentlv>"
__all__ = ["CONFIG", "Mask"]


def smooth(gdspathin, gdspathout, threshold=0.3):
    """ returns a litho simulated version of gdspathin
    """
    m = Mask()
    m.x_range = [-300.0, 300.0]
    m.y_range = [-300.0, 300.0]
    m.x_gridsize = 10
    m.y_gridsize = 10
    m.openGDS(CONFIG["gdslib"] / "AND2_X4.gds", 10)
    m.maskfft()

    m.smooth()
    mth = m.sdata > threshold


if __name__ == "__main__":
    gdspathin = CONFIG["samples"] / "verniers.gds"
    gdspathout = CONFIG["samples"] / "verniers_out.gds"
    smooth(gdspathin, gdspathout)
