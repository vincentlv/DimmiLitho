"""

"""

import numpy as np
import pyfftw
import scipy.signal as sg


class ImageHopkins:
    """ImageHopkinsList is a container, used for, e.g., robust mask synthesis
    This is for Scalar Assumption, not Vector
    """

    def __init__(self, mask, tcc):
        self.tcc = tcc  # TCC
        self.mask = mask  # Mask
        self.order = tcc.order  # TCC Order
        self.resist_a = 80  # Resist model parameter: sharpness
        self.resist_t = 0.6  # Resist model parameter: threshold
        self.kernels = tcc.kernels  # Kernels
        self.coefs = tcc.coefs  # Coefs

        self.norm = self.mask.y_gridnum * self.mask.x_gridnum
        self.x1 = int(np.floor(self.mask.x_gridnum / 2) - self.tcc.s.fnum)
        self.x2 = int(np.floor(self.mask.x_gridnum / 2) + self.tcc.s.fnum + 1)
        self.y1 = int(np.floor(self.mask.y_gridnum / 2) - self.tcc.s.gnum)
        self.y2 = int(np.floor(self.mask.y_gridnum / 2) + self.tcc.s.gnum + 1)

        self.spat_part = pyfftw.empty_aligned(
            (self.mask.y_gridnum, self.mask.x_gridnum), dtype="complex128"
        )
        self.freq_part = pyfftw.empty_aligned(
            (self.mask.y_gridnum, self.mask.x_gridnum), dtype="complex128"
        )
        self.ifft_image = pyfftw.FFTW(
            self.freq_part, self.spat_part, axes=(0, 1), direction="FFTW_BACKWARD"
        )

    def calAI(self):  # much faster than calAIold()
        AI_freq_dense = np.zeros(
            (self.mask.y_gridnum, self.mask.x_gridnum), dtype=np.complex128
        )
        AI_freq_sparse = np.zeros(
            (int(self.y2 - self.y1), int(self.x2 - self.x1)), dtype=np.complex128
        )
        for ii in range(self.order):
            self.x1 = int(self.x1)
            self.x2 = int(self.x2)
            self.y1 = int(self.y1)
            self.y2 = int(self.y2)
            e_field = (
                self.kernels[:, :, ii]
                * self.mask.fdata[self.y1 : self.y2, self.x1 : self.x2]
            )
            e_field_conj = (
                np.conj(np.rot90(self.kernels[:, :, ii], 2))
                * self.mask.fdata[self.y1 : self.y2, self.x1 : self.x2]
            )
            AA = sg.convolve2d(e_field, e_field_conj, "same", "wrap")
            AI_freq_sparse += self.coefs[ii] * AA
        AI_freq_dense[self.y1 : self.y2, self.x1 : self.x2] = AI_freq_sparse

        self.freq_part[:] = np.fft.ifftshift(AI_freq_dense)
        self.ifft_image()
        self.AI = np.real(np.fft.fftshift(self.spat_part)) / self.norm

    def calAIold(self):
        AI = np.zeros((self.mask.y_gridnum, self.mask.x_gridnum))
        for ii in range(self.order):
            e_field = np.zeros(
                (self.mask.y_gridnum, self.mask.x_gridnum), dtype=np.complex128
            )
            e_field[self.y1 : self.y2, self.x1 : self.x2] = (
                self.kernels[:, :, ii]
                * self.mask.fdata[self.y1 : self.y2, self.x1 : self.x2]
            )
            AA = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(e_field)))
            AI += self.coefs[ii] * np.abs(AA * np.conj(AA))
        self.AI = AI

    def calRI(self):
        self.RI = 1 / (1 + np.exp(-self.resist_a * (self.AI - self.resist_t)))


class ImageHopkinsList(ImageHopkins):
    def __init__(self, mask, tccList):
        self.mask = mask
        self.tcc = tccList
        self.order = tccList.order
        self.kernelList = tccList.kernelList
        self.coefList = tccList.coefList
        self.focusList = tccList.focusList
        self.focusCoef = tccList.focusCoef
        self.doseList = [1.0]
        self.doseCoef = [1.0]
        self.AIList = []
        self.RIList = []
        self.resist_a = 80
        self.resist_tRef = 0.5

        self.norm = self.mask.y_gridnum * self.mask.x_gridnum
        self.x1 = np.floor(self.mask.x_gridnum / 2) - self.tcc.s.fnum
        self.x2 = np.floor(self.mask.x_gridnum / 2) + self.tcc.s.fnum + 1
        self.y1 = np.floor(self.mask.y_gridnum / 2) - self.tcc.s.gnum
        self.y2 = np.floor(self.mask.y_gridnum / 2) + self.tcc.s.gnum + 1

        self.spat_part = pyfftw.empty_aligned(
            (self.mask.y_gridnum, self.mask.x_gridnum), dtype="complex128"
        )
        self.freq_part = pyfftw.empty_aligned(
            (self.mask.y_gridnum, self.mask.x_gridnum), dtype="complex128"
        )
        self.ifft_image = pyfftw.FFTW(
            self.freq_part, self.spat_part, axes=(0, 1), direction="FFTW_BACKWARD"
        )

    def calculate(self):
        length = len(self.focusList)
        for ii in range(length):
            self.kernels = self.kernelList[ii]
            self.coefs = self.coefList[ii]
            self.calAI()
            self.AIList.append(self.AI)
            self.RIList.append([])
            for jj in self.doseList:
                self.resist_t = self.resist_tRef * jj
                self.calRI()
                self.RIList[ii].append(self.RI)


if __name__ == "__main__":
    from litho.lens import Lens
    from litho.mask import Mask
    from litho.source import Source
    from litho.tcc import TCC

    mp = [
        [
            [-1, 6],
            [-1, 2],
            [1, 2],
            [1, 1],
            [6, 1],
            [6, 0],
            [0, 0],
            [0, 1],
            [-2, 1],
            [-2, 6],
            [-1, 6],
        ],
        [
            [6, -1],
            [6, -2],
            [1, -2],
            [1, -3],
            [4, -3],
            [4, -6],
            [3, -6],
            [3, -4],
            [0, -4],
            [0, -1],
            [6, -1],
        ],
    ]
    m = Mask()
    m.x_range = [-300.0, 300.0]
    m.y_range = [-400.0, 300.0]
    m.x_gridsize = 1.0
    m.y_gridsize = 1.0
    m.CD = 40
    m.polygons = mp
    m.poly2mask()
    m.smooth()
    m.maskfft()

    """nominal ILT setting"""
    s = Source()
    s.na = 1.35
    s.maskxpitch = 600.0
    s.maskypitch = 800.0
    s.type = "annular"
    s.sigma_in = 0.6
    s.sigma_out = 0.8
    s.update()
    s.ifft()

    o = Lens()
    o.na = s.na
    o.maskxpitch = s.maskxpitch
    o.maskypitch = s.maskypitch
    o.update()
    o.calPupil()
    o.calPSF()

    t = TCC(s, o)
    t.calMutualIntensity()
    t.calSpatTCC()
    t.svd()

    i = ImageHopkins(m, t)
    i.calAI()

    """robust ILT setting"""
    # from litho.lens import LensList
    # from litho.tcc import TCCList

    # s = Source()
    # s.na = 1.25
    # s.maskxpitch = 600.0
    # s.maskypitch = 1000
    # s.type = 'annular'
    # s.sigma_in = 0.5
    # s.sigma_out = 0.8
    # s.update()
    # s.ifft()
    #
    # o = LensList()
    # o.na = s.na
    # o.maskxpitch = s.maskxpitch
    # o.maskypitch = s.maskypitch
    # o.focusList = [-50, 0, 50]
    # o.focusCoef = [0.5, 1, 0.5]
    # o.calculate()
    #
    # t = TCCList(s,o)
    # t.calculate()
    #
    # i = ImageHopkinsList(m,t)
    # i.doseList = [0.95, 1, 1.05]
    # i.doseCoef = [0.5, 1, 0.5]
    # i.calculate()
