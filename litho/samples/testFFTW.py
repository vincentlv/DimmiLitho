"""
Created on Fri Apr 29 2016
@author: WenLv (wenlv@hust.edu.cn)

We need PRE-PROCESS the matrix BEFORE and AFTER fftw, because unless
1) The zero frequency is un-centered
2) The frequency pattern is with oscillation

Here are some tricks:
1) One way is using fftshift both before and after. Adopted here!
2) Maybe, put a(i,j) = a(i,j)*(-1)^(j+k) will center the zero freqency
   without need of fftshift, faster but with sacrificing accuracy slightly.

"""

import matplotlib.pyplot as plt
import numpy as np
import pyfftw

# FFTW syntax
a = pyfftw.empty_aligned((101, 101), dtype="complex128")
b = pyfftw.empty_aligned((101, 101), dtype="complex128")
fft_o = pyfftw.FFTW(a, b, axes=(0, 1))

a[:] = np.zeros((101, 101))
a[40:60, 40:60] = 1

# kk = np.ones((101,101))*(-1)
# kk[::2,::2] = 1
# kk[1::2,1::2] = 1
# a[:] = a*kk

a[:] = np.fft.ifftshift(a)
fft_o()
b[:] = np.fft.fftshift(b)


plt.imshow(b.real)
plt.colorbar()
plt.show()
