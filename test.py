"""
Created on Wen Apr 27 2016
@author: WenLv (wenlv@hust.edu.cn)
"""

from lens import LensList
from tcc import TCCList
from mask import Mask
from source import Source 
from ilt import RobustILT

import numpy as np    
import time


m = Mask()
m.x_gridsize = 2.5
m.y_gridsize = 2.5
m.openGDS('./NanGateLibGDS/NOR2_X2.gds',layername=11,boundary=0.3)
m.maskfft()


s = Source()
s.na = 1.35
s.maskxpitch = m.x_range[1] - m.x_range[0]
s.maskypitch = m.y_range[1] - m.y_range[0]
s.type = 'annular'
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
print "Calculating TCC and SVD kernels"
t = TCCList(s,o)
t.calculate()
print "###taking %1.3f seconds" %(time.time()-tic)


print "Calculating ILT"
i = RobustILT(m,t)
i.image.resist_a = 100
i.image.resist_tRef = 0.9
i.stepSize = 0.4
i.image.doseList = [0.9,1,1.1]
i.image.doseCoef = [0.3,1, 0.3]
i.run(100)




