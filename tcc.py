"""
Created on Sta Jun 28 2015
@author: WenLv (wenlv@hust.edu.cn)
"""
import numpy as np
from source import Source
from lens import LensList
import scipy as sci
import pyfftw

class TCC:
    def __init__(self,source,lens):
        self.s = source
        self.s.update()
        self.s.ifft()

        self.l = lens
        self.l.update()
        self.l.calPupil()
        self.l.calPSF()
        
        self.order = 7
        self.psf = lens.data
        
    def calMutualIntensity(self):
        self.gnum,self.fnum = self.s.data.shape
        J = np.zeros((self.gnum,self.fnum,self.gnum,self.fnum),dtype = np.complex)       
        for ii in range(self.gnum):
            for jj in range(self.fnum):
                J[:,:,ii,jj]=self.s.spatMutualData.real[(self.gnum-ii-1):(2*self.gnum-ii-1),(self.fnum-jj-1):(2*self.fnum-jj-1)]
        self.jsource = np.reshape(J,(self.gnum*self.fnum,self.gnum*self.fnum))
        
    def calSpatTCC(self):
        H = np.reshape(self.psf,(self.psf.size,1))
        self.tcc2d = self.jsource * np.dot(H,H.transpose())/self.s.detaf/self.s.detag
        
    def svd(self):
        self.spat_part = pyfftw.empty_aligned((self.gnum,self.fnum,self.gnum,self.fnum),\
                                               dtype='complex128')
        self.freq_part = pyfftw.empty_aligned((self.gnum,self.fnum,self.gnum,self.fnum),\
                                               dtype='complex128')
        self.fft_svd = pyfftw.FFTW(self.spat_part,self.freq_part,axes=(0,1,2,3))    
        

        tcc4d = self.tcc2d.reshape((self.gnum,self.fnum,self.gnum,self.fnum))
        self.spat_part[:] = np.fft.ifftshift(tcc4d)
        self.fft_svd()
        tcc4df = np.fft.fftshift(self.freq_part)
        # tcc4df = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(tcc4d)))

        tcc2df = tcc4df.reshape((self.gnum*self.fnum,self.gnum*self.fnum))
        
        # U,S,V = np.linalg.svd(tcc2df)
        U,S,V = sci.sparse.linalg.svds(tcc2df, self.order) #faster than svd
        self.coefs = S[0:self.order]
        self.kernels = np.zeros((self.gnum,self.fnum,self.order),dtype = np.complex)
        for ii in range(self.order):
            self.kernels[:,:,ii] = np.reshape(U[:,ii],(self.gnum,self.fnum))
        

class TCCList(TCC):
    def __init__(self, source, lensList):
        self.s = source
        self.PSFList = lensList.sDataList
        self.order = 7
        self.focusList = lensList.focusList
        self.focusCoef = lensList.focusCoef
        self.kernelList = []
        self.coefList = []
    
    def calculate(self):
        self.calMutualIntensity()
        for ii in self.PSFList:
            self.psf = ii
            self.calSpatTCC()
            self.svd()
            self.coefList.append(self.coefs)
            self.kernelList.append(self.kernels)
        
        
if __name__ == "__main__":
    s = Source()
    s.type = 'annular'
    s.sigma_in = 0.5
    s.sigma_out = 0.8
    s.na = 1.35
    s.update()
    s.ifft()
    
    o = LensList()
    o.maskxpitch = s.maskxpitch
    o.maskypitch = s.maskypitch
    o.na =s.na
    o.focusList = [-50, 0,50]
    o.focusCoef = [0.5, 1, 0.5]
    o.calculate()
    
    tcc = TCCList(s,o)
    tcc.calculate()
    
    
    