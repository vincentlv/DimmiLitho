""" inverse litho

calculate needed mask to account for litho effects
"""

import time

import matplotlib.pyplot as plt
from litho.config import CONFIG
from litho.ilt import RobustILT
from litho.lens import LensList
from litho.mask import Mask
from litho.source import Source
from litho.tcc import TCCList

m = Mask()
m.x_gridsize = 2.5
m.y_gridsize = 2.5
m.openGDS(CONFIG["samples"] / "verniers.gds", layername=1, boundary=0.3)
m.maskfft()

s = Source()
s.na = 1.35
s.maskxpitch = m.x_range[1] - m.x_range[0]
s.maskypitch = m.y_range[1] - m.y_range[0]
s.type = "annular"
s.sigma_in = 0.7
s.sigma_out = 0.9
s.smooth_deta = 0.00
s.shiftAngle = 0
s.update()
s.ifft()

o = LensList()
o.na = s.na
o.maskxpitch = s.maskxpitch
o.maskypitch = s.maskypitch
o.focusList = [0]
o.focusCoef = [1]
o.calculate()


tic = time.time()
print("Calculating TCC and SVD kernels")
t = TCCList(s, o)
t.calculate()
print("###taking %1.3f seconds" % (time.time() - tic))


print("Calculating ILT")
iterations = 20
i = RobustILT(m, t)
i.image.resist_a = 100
i.image.resist_tRef = 0.9
i.stepSize = 0.4
i.image.doseList = [0.9, 1, 1.1]
i.image.doseCoef = [0.3, 1, 0.3]
i.run(iterations)

plt.figure()
plt.imshow(i.maskdata, origin="lower")

plt.figure()
plt.imshow(i.maskdata > 0.9, origin="lower")
plt.show()
