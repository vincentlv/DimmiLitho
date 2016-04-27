"""
Created on Sta Jun 28 2015
@author: WenLv (wenlv@hust.edu.cn)

NOTE: This is for Scalar Assumption, not Vector
      This is for Hopkins
      ImageHopkinsList is a container, used for, e.g., robust mask synthesis
"""

import numpy as np

class ImageHopkins:
    def __init__(self, mask, tcc):
        self.tcc = tcc               # TCC
        self.mask = mask             # Mask
        self.order = tcc.order       # TCC Order
        self.resist_a = 80           # Resist model parameter: shapness
        self.resist_t = 0.6          # Resist model parameter: threshold
        self.kernels = tcc.kernels   # Kerneks
        self.coefs = tcc.coefs       # Coefs
        
    def calAI(self):        
        x1 = np.floor(self.mask.x_gridnum/2) - self.tcc.s.fnum
        x2 = np.floor(self.mask.x_gridnum/2) + self.tcc.s.fnum  + 1
        y1 = np.floor(self.mask.y_gridnum/2) - self.tcc.s.gnum 
        y2 = np.floor(self.mask.y_gridnum/2) + self.tcc.s.gnum  + 1
        
        AI = np.zeros((self.mask.y_gridnum,self.mask.x_gridnum))       
        for ii in range(self.order):
            e_field = np.zeros((self.mask.y_gridnum,self.mask.x_gridnum),dtype=np.complex) 
            e_field[y1:y2,x1:x2] = self.kernels[:,:,ii]*self.mask.fdata[y1:y2,x1:x2]
            AA = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(e_field)))
            AI += self.coefs[ii]*np.abs(AA*np.conj(AA))     
        self.AI = AI
    
    def calRI(self):
        self.RI = 1/(1+np.exp(-self.resist_a*(self.AI - self.resist_t)))
 
class ImageHopkinsList(ImageHopkins):
    def __init__(self,mask, tccList):
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
        
    def calculate(self):
        length = len(self.focusList)
        for ii in xrange(length):
            self.kernels = self.kernelList[ii]
            self.coefs = self.coefList[ii]
            self.calAI()
            self.AIList.append(self.AI)
            self.RIList.append([])
            for jj in self.doseList:
                self.resist_t = self.resist_tRef*jj
                self.calRI()
                self.RIList[ii].append(self.RI)
                        
if __name__ == "__main__":
    from tcc import TCCList,TCC
    from mask import Mask 
    from source import Source
    from lens import LensList, Lens
    
    mp = [ [[-1,6],[-1, 2],[1, 2],[1, 1],[6, 1],[6, 0],[0, 0],[0, 1],[-2, 1],[-2, 6],[-1, 6]], \
       [[6, -1],[6, -2],[1, -2],[1, -3],[4, -3],[4, -6],[3, -6],[3, -4],[0, -4],[0, -1],[6, -1]] ]
    m = Mask()
    m.x_range = [-300.0,300.0]
    m.y_range = [-300.0,300.0]
    m.x_gridsize = 2.0
    m.y_gridsize = 2.0
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
    s.type = 'annular'
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
    
    t = TCC(s,o)
    t.calMutualIntensity()
    t.calSpatTCC()
    t.svd()
    
    i = ImageHopkins(m,t)
    i.calAI()
    
    
    """robust ILT setting"""
    #s = Source()
    #s.na = 1.25
    #s.maskxpitch = 600.0
    #s.maskypitch = 1000
    #s.type = 'annular'
    #s.sigma_in = 0.5
    #s.sigma_out = 0.8
    #s.update()
    #s.ifft()
    #
    #o = LensList()
    #o.na = s.na
    #o.maskxpitch = s.maskxpitch
    #o.maskypitch = s.maskypitch
    #o.focusList = [-50, 0, 50]
    #o.focusCoef = [0.5, 1, 0.5]
    #o.calculate()
    #
    #t = TCCList(s,o)
    #t.calculate()
    #
    #i = ImageHopkinsList(m,t)
    #i.doseList = [0.95, 1, 1.05]
    #i.doseCoef = [0.5, 1, 0.5]
    #i.calculate()