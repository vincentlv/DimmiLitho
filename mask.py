"""
Created on Sta Jun 28 2015
@author: WenLv (wenlv@hust.edu.cn)

Note: Binary Mask
"""

import numpy as np
import math
import scipy.signal as sg
import pyfftw

"""This Libarary is Used Only for Polygon Data Processing
"""
from PIL import Image, ImageDraw


class Mask:
    def __init__(self):
        """Default Parameters
           1. x/y_range means the computing area. Different value are supported
           2. x/y_gridsize the simulated size of the area. Different value are supported
           3. CD infomation is usable for method poly2mask
        """
        self.x_range = [-500,500] # nm
        self.y_range = [-500,500]
        self.x_gridsize = 2 # nm
        self.y_gridsize = 2
        self.CD = 45 #nm 
       
        
    def poly2mask(self):
        """Get Pixel-based Mask Image from Polygon Data
           The Poylgon Data Form are sensitive
           Similar to poly2mask in Matlab
        """
        self.x_gridnum = int((self.x_range[1] - self.x_range[0])/self.x_gridsize)
        self.y_gridnum = int((self.y_range[1] - self.y_range[0])/self.y_gridsize)
        img = Image.new('L',(self.x_gridnum,self.y_gridnum),0)
        
        self.perimeter = 0.0;
        for ii in self.polygons:
            pp = np.array(ii)*self.CD  #polygon 
            polygonlen = len(pp)
            self.perimeter += np.sum(np.abs(pp[0:-1] - pp[1:polygonlen]))     
            pp[:,0] = (pp[:,0] - self.x_range[0])/self.x_gridsize
            pp[:,1] = (pp[:,1] - self.y_range[0])/self.y_gridsize
            vetex_list = list(pp)
            polygon = [tuple(y) for y in vetex_list]            
            ImageDraw.Draw(img).polygon(polygon,outline=1,fill=1)
                 
        self.data = np.array(img)
        self.data = np.float64(self.data)

        self.spat_part = pyfftw.empty_aligned((self.y_gridnum,self.x_gridnum),\
                                               dtype='complex128')
        self.freq_part = pyfftw.empty_aligned((self.y_gridnum,self.x_gridnum),\
                                               dtype='complex128')
        self.fft_mask = pyfftw.FFTW(self.spat_part,self.freq_part,axes=(0,1)) 

    def openGDS(self, gdsdir, layername, boundary = 0.16, scalerate = 45/70.0):
        
        from gdsii.library import Library
        with open(gdsdir, 'rb') as stream:
            lib = Library.load(stream)
            
        a = lib.pop(0)
        b = []
        xmin = []
        xmax = []
        ymin = []
        ymax = []
        for ii in xrange(0,len(a)):
            if a[ii].layer == layername:
                #if hasattr(a[ii],'data_type'):
                if len(a[ii].xy)>1:
                    aa = np.array(a[ii].xy)/10.0*scalerate
                    b.append(aa)                    
                    xmin.append(min([k for k,v in aa]))
                    xmax.append(max([k for k,v in aa]))
                    ymin.append(min([v for k,v in aa]))
                    ymax.append(max([v for k,v in aa]))                     
        self.polylist = b
               
        xmin = min(xmin)
        xmax = max(xmax)
        ymin = min(ymin)
        ymax = max(ymax)
        self.xmin = xmin - boundary*(xmax - xmin)
        self.xmax = xmax + boundary*(xmax - xmin)
        self.ymin = ymin - boundary*(ymax - ymin)
        self.ymax = ymax + boundary*(ymax - ymin)
        self.x_range = [self.xmin, self.xmax]
        self.y_range = [self.ymin, self.ymax]
        
        self.x_gridnum = int((self.xmax - self.xmin)/self.x_gridsize)
        self.y_gridnum = int((self.ymax - self.ymin)/self.y_gridsize)
        img = Image.new('L',(self.x_gridnum,self.y_gridnum),0)
        
        self.perimeter = 0.0;
        for ii in self.polylist:
            pp = np.array(ii)  #polygon 
            polygonlen = len(pp)
            self.perimeter += np.sum(np.abs(pp[0:-1] - pp[1:polygonlen]))
                 
            pp[:,0] = (pp[:,0] - self.xmin)/self.x_gridsize
            pp[:,1] = (pp[:,1] - self.ymin)/self.y_gridsize
            vetex_list = list(pp)
            polygon = [tuple(y) for y in vetex_list]            
            ImageDraw.Draw(img).polygon(polygon,outline=1,fill=1)
        
            self.perimeter += np.sum(np.abs(pp[0:-1] - pp[1:polygonlen]))
                 
        self.data = np.array(img)
        # Creat fourier transform pair, pyfftw syntax
        self.spat_part = pyfftw.empty_aligned((self.y_gridnum,self.x_gridnum),\
                                               dtype='complex128')
        self.freq_part = pyfftw.empty_aligned((self.y_gridnum,self.x_gridnum),\
                                               dtype='complex128')
        self.fft_mask = pyfftw.FFTW(self.spat_part,self.freq_part,axes=(0,1))            
     
    # use the fftw packages            
    def maskfft(self):
        self.spat_part[:] = np.fft.ifftshift(self.data)
        self.fft_mask()
        self.fdata = np.fft.fftshift(self.freq_part)
            
    def maskfftold(self):
        self.fdata = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(self.data)))

    def smooth(self):
        xx = np.linspace(-1,1,21)
        X, Y = np.meshgrid(xx,xx)
        R = X**2 + Y**2
        G =  np.exp(-10*R)
        D = sg.convolve2d(0.9*self.data+0.05, G,'same')/np.sum(G)
        self.sdata = D   

        
if __name__ == '__main__': 
    """polygon 2 mask"""                    
    #mp = [ [[-1,6],[-1, 2],[1, 2],[1, 1],[6, 1],[6, 0],[0, 0],[0, 1],[-2, 1],[-2, 6],[-1, 6]], \
    #   [[6, -1],[6, -2],[1, -2],[1, -3],[4, -3],[4, -6],[3, -6],[3, -4],[0, -4],[0, -1],[6, -1]] ]
    #m = Mask()
    #m.x_range = [-300.0,300.0]
    #m.y_range = [-300.0,300.0]
    #m.x_gridsize = 1.5
    #m.y_gridsize = 1.5
    #m.CD = 45
    #m.polygons = mp
    #m.poly2mask()
    
    """from GDS"""
    m = Mask()
    m.x_range = [-300.0,300.0]
    m.y_range = [-300.0,300.0]
    m.x_gridsize = 10
    m.y_gridsize = 10
    m.openGDS('./NanGateLibGDS/AND2_X4.gds',10)
    m.maskfft()

    import matplotlib.pyplot as plt
    plt.imshow(m.data,\
               extent=(m.x_range[0],m.x_range[1],m.y_range[0],m.y_range[1]),\
               cmap='hot',\
               interpolation='none')
    plt.show()
    #m.smooth()






