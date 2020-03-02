"""
Created on Mon Dec 29 2014
@author: WenLv (wenlv@hust.edu.cn)

NOTE: This is for Scalar Pupil Assumption, not Jones Pupil
      LensList is a Lens container, used for, e.g., robust mask synthesis
"""

import numpy as np
import copy
from zernike import zerniken


class Lens:
    def __init__(self):
        self.na = 1.35
        self.nLiquid = 1.414
        self.wavelength = 193.0  # nm
        self.defocus = 0.0  # nm
        self.maskxpitch = 1000  # nm
        self.maskypitch = 1000
        self.Zn = [9]
        self.Cn = [0.0]

    def update(self):
        self.detaf = self.wavelength / (self.maskxpitch * self.na)
        self.detag = self.wavelength / (self.maskypitch * self.na)
        self.fnum = int(np.ceil(2 / self.detaf))
        self.gnum = int(np.ceil(2 / self.detag))

    def calPupil(self, shiftx=0, shifty=0):
        fx = np.linspace(
            -self.fnum * self.detaf, self.fnum * self.detaf, 2 * self.fnum + 1
        )
        fy = np.linspace(
            -self.gnum * self.detag, self.gnum * self.detag, 2 * self.gnum + 1
        )
        FX, FY = np.meshgrid(fx - shiftx, fy - shifty, indexing="xy")

        R = np.sqrt(FX ** 2 + FY ** 2)
        TH = np.arctan2(FY, FX)
        H = copy.deepcopy(R)
        H = np.where(H > 1.0, 0.0, 1.0)
        R[R > 1.0] = 0.0

        W = np.zeros((2 * self.gnum + 1, 2 * self.fnum + 1), dtype=np.complex)

        for ii in range(len(self.Zn)):
            W = W + zerniken(self.Zn[ii], R, TH) * self.Cn[ii]

        if self.na < 1:
            W = W + self.defocus / self.wavelength * (
                np.sqrt(1 - (self.na ** 2) * (R ** 2)) - 1
            )
        elif self.na >= 1:
            # W = W + self.defocus/self.wavelength*\
            #         (np.sqrt(self.nLiquid**2-(self.na**2)*(R**2))-self.nLiquid)
            W = W + (self.na ** 2) / (2 * self.wavelength) * self.defocus * (R ** 2)
        self.fdata = H * np.exp(-1j * 2 * (np.pi) * W)

    def calPSF(self):
        normlize = 1  # self.detaf * self.detag
        self.data = (
            np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(self.fdata))) * normlize
        )


class LensList(Lens):
    def __init__(self):
        Lens.__init__(self)
        self.focusList = [0.0]
        self.focusCoef = [1.0]
        self.fDataList = []
        self.sDataList = []

    def calculate(self):
        self.update()
        for ii in self.focusList:
            self.defocus = ii
            self.calPupil()
            self.fDataList.append(self.fdata)
            self.calPSF()
            self.sDataList.append(self.data)


if __name__ == "__main__":
    o = LensList()
    o.na = 0.85
    o.focusList = [-50, 0, 50]
    o.focusCoef = [0.5, 1, 0.5]
    o.calculate()

    O = Lens()
    O.na = 0.85
    O.update()
    O.calPupil()
    O.calPSF()
