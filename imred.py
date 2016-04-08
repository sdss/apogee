import numpy as np
from astropy.io import fits
from astropy.io import ascii
from astropy.modeling import models, fitting
from astropy.convolution import convolve, Box1DKernel
import matplotlib.pyplot as plt
import glob
import bz2
import os
import pdb
try: 
    import pyds9
except:
    print 'pyds9 is not available, proceeding'

def setup(name='*',dir='./',idet=1) :
    """ 
    Sets up global root, indir, and detector iD

    Args:
      name (str) : root file name
    
    Keyword args:
      dir (str) : input directory name [default='./']

    Returns:
      nothing
    """
    global root, indir, det
    root = name
    indir = dir+'/'
    det = getdet(idet)

def rd(num,ext=0,bias=True,verbose=False) :
    """ 
    Reads image by name or number, by default subtracts bias using overscan

    Args:
      num (int or str) : number or name of image to read

    Keyword args:
      ext (int) : extension to read (default=0)
      bias (bool) : subtract overscan? (default=True)
      verbose (bool) : use verbose mode? (default=False)

    Returns:
      HDU : HDU of image
    """

    if type(num) == str :
        file=num
    else :
        file=glob.glob(indir+'/'+root+'.'+det.formstr.format(num)+'.fits*')
        if len(file) > 0 : file=file[0]
    if verbose: print 'Reading: ', file
    if '.bz2' in file :
        hdu=fits.open(bz2.BZ2File(file),ignore=True,ignore_missing_end=True)
    else :
        hdu=fits.open(file,ignore=True,ignore_missing_end=True)

    if 'object' in hdu[0].header.cards and hdu[0].header['object'] == '' : 
        hdu[0].header['object'] = os.path.basename(file)
    if bias :
        if det.biastype == 0 :
            b=det.biasbox.mean(hdu[ext].data)
            if verbose: print 'subtracting overscan: ', b
            hdu[ext].data -= b
        elif det.biastype == 1 :
            over=np.median(hdu[ext].data[:,det.biasbox.xmin:det.biasbox.xmax],axis=1)
            over=stretch(over,ncol=hdu[ext].data.shape[1])
            hdu[ext].data -= over
    return hdu[ext]

def reduce(num,bias=None,flat=None,trim=False,verbose=False) :
    """ 
    Reads with overscan subtraction, subtracts bias and divides by flat as specified
 
    Args:
      num (int or str) : number or name of image to reduce

    Keyword args:
      bias (numpy array) : bias frame to subtract
      flat (numpy array) : flat frame to divide
      trim (bool) : trim frame? (default=False)
      verbose (bool) : be verbose? (default=False)

    Returns:
      HDU : hdu of reduced image
    """
    if verbose: print 'reading: ', num
    hdu=rd(num,verbose=verbose)
    if bias is not None :
        if verbose: print 'subtracting bias'
        hdu.data -= bias
    if flat is not None :
        if verbose: print 'dividing by flat'
        hdu.data /= flat
    if trim :
        out= window(hdu,det.trimbox)
    else :
        out= hdu  
    return out

def combine(ims,norm=False,bias=None,flat=None,trim=False,verbose=False) :
    """ 
    Combines input list of images (names or numbers) by median, 
    optionally normalizes before combination

    Args:
      ims (list of int or str): list of images to combine

    Keyword args:
      norm (bool) : normalize images before combining? (default=False)
      bias (numpy array) : bias frame to subtract? (default=None)
      flat (numpy array) : flat frame to divide? (default=None)

    Returns:
      median of input data arrays, after reduction as requested
    """
    cube=[]
    for im in ims :
        print 'Reading image: ', im
        h=reduce(im,bias=bias,flat=flat,trim=trim,verbose=verbose) 
        if norm :
            b=det.normbox
            norm=np.median(h.data[b.ymin:b.ymax,b.xmin:b.xmax])
            print 'Normalizing image by : ', norm
            cube.append(h.data/norm)
        else :
            cube.append(h.data)
    print 'Combining: ', ims
    return np.median(cube,axis=0)
    
def specflat(flat,rows=True,indiv=False,wid=100) :
    """
    Removes spectral signature from a flat by dividing by smoothed version

    Args:
      flat: input flat fields
   
    Keyword args:
      rows (bool) : specifies if smoothing is along rows (default), otherwise columns
      indiv (bool) : specifies if smoothing is done row-by-row, or single for whole image (default)
      wid (int) : boxcar kernel width (default=100)

    Returns:
      flat with spectral signature removed
    """     
    boxcar = Box1DKernel(wid)
    smooth=flat
    if rows :
        if indiv :
            for row in range(flat.shape[0]) :
                smooth[row,:] /= convolve(flat[row,:],boxcar,boundary='extend')
        else : 
            c=convolve(flat.sum(axis=0),boxcar,boundary='extend')
            for row in range(flat.shape[0]) :
                smooth[row,:] /= c

    else :
        print 'smoothing by columns not yet implemented!'
        pdb.set_trace()
      
    return smooth 

def abx(im,box) :
    """ 
    Returns dictionary with image statistics in box 
    """
    return {'mean': box.mean(im),
            'stdev': box.stdev(im),
            'max': box.max(im),
            'min': box.min(im),
            'peakx': np.unravel_index(im[box.ymin:box.ymax,box.xmin:box.xmax].argmax(),(box.nrow(),box.ncol()))[1]+box.xmin,
            'peaky': np.unravel_index(im[box.ymin:box.ymax,box.xmin:box.xmax].argmax(),(box.nrow(),box.ncol()))[0]+box.ymin}

class BOX() :
    """ Defines BOX class"""
    def __init__(self) :
        self.xmin = -1
        self.xmax = -1
        self.ymin = -1
        self.ymax = -1

    def set(self,xmin,xmax,ymin,ymax):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

    def nrow(self):
        return(self.ymax-self.ymin)

    def ncol(self):
        return(self.xmax-self.xmin)

    def show(self):
        print '[{:4d}:{:4d},{:4d}:{:4d}]'.format(self.ymin,self.ymax,self.xmin,self.xmax)

    def mean(self,data):
        return data[self.ymin:self.ymax,self.xmin:self.xmax].mean() 

    def stdev(self,data):
        return data[self.ymin:self.ymax,self.xmin:self.xmax].std() 

    def max(self,data):
        return data[self.ymin:self.ymax,self.xmin:self.xmax].max() 

    def min(self,data):
        return data[self.ymin:self.ymax,self.xmin:self.xmax].min() 

    def median(self,data):
        return np.median(data[self.ymin:self.ymax,self.xmin:self.xmax])

class DET() :
    """ Defines detector class """
    def __init__(self) :
        self.gain = 0.
        self.rn = 0.
        self.biastype = 0
        self.biasbox = BOX()
        self.normbox = BOX()
        self.trimbox = BOX()
        self.formstr = "{:04d}"

def getdet(idet) :
    """ returns detector object given input detector index"""
    d=DET()
    if idet == 11 :
       # APO SPICAM
       d.gain=3.8
       d.rn=6
       d.biasbox.set(1040,1070,10,1000)
       d.normbox.set(400,600,400,600)
    elif idet == 16 :
       # APO ARCES
       d.gain=3.8
       d.rn=7
       d.biasbox.set(2052,2057,20,2028)
       d.trimbox.set(200,1850,0,2047)
    elif idet == 17 :
       # APO 1m APOGEE
       d.gain=3.8
       d.rn=6
       d.biastype=0
       d.biasbox.set(520,540,10,500)
       d.normbox.set(400,600,400,600)
       d.formstr = "{:03d}"
    elif idet == 36 :
       # APO ARCTIC
       d.gain=3.8
       d.rn=6
       d.biastype=1
       d.biasbox.set(2060,2090,10,2040)
       d.normbox.set(900,1100,900,1100)
    elif idet == 44 :
       # DIS blue
       d.gain=1.71
       d.rn=3.9
       d.biasbox.set(1030,1050,0,2047)
       d.trimbox.set(0,2047,0,1023)
       d.formstr='{:04d}b'
    elif idet == 45 :
       # DIS red
       d.gain=1.71
       d.rn=3.9
       d.biasbox.set(1030,1050,0,2047)
       d.trimbox.set(0,2047,0,1023)
       d.formstr='{:04d}r'
    return d

def gfit(data,xcen,ycen,size=5,sub=True) :
    """ Does gaussian fit to input data given initial xcen,ycen"""
    g_init=models.Gaussian2D(x_mean=xcen,y_mean=ycen,x_stddev=1,y_stddev=1,amplitude=data[ycen,xcen])+models.Const2D(0.)
    fit=fitting.LevMarLSQFitter()
    #fit=fitting.SLSQPLSQFitter()
    y,x=np.mgrid[ycen-size:ycen+size,xcen-size:xcen+size]
    z=data[ycen-size:ycen+size,xcen-size:xcen+size]
    g=fit(g_init,x,y,z)
    if sub :
        out=data
        out[ycen-size:ycen+size,xcen-size:xcen+size]-=g[0](x,y)
        return out
    return g[0](x,y)+g[1](x,y)

def look(disp,pause=True,ds9=False,i1=None,i2=None,min=None, max=None) :
    """ Displays series of files """
    if i1 is not None and i2 is not None :
        files=range(i1,i2+1)
    else :
        files=glob.glob(indir+'/*.fits*')
    for file in files :
       hd=rd(file)
       tv(disp,hd,min=min,max=max)
       plt.plot(hd.data[:,2060:2080].mean(axis=1))
       print hd.data[:,2060:2080].mean(),hd.data[:,2060:2080].std()
       plt.draw()
       if pause :
           pdb.set_trace()

def getfiles(type,listfile=None,filter=None,verbose=False) :
    """ Get all files of desired type from specified directory.
        If file is specified, read numbers from that file, else use IMAGETYP card"""
    if listfile is None :
        list=[]
        if verbose: print 'directory: ', indir
        for file in glob.glob(indir+'/'+root+'*.fits*') :
            head=fits.open(file)[0].header
            #if verbose :
            #   print 'file: ', file
            #   print 'IMAGETYPE: ', head['IMAGETYP']

            try :
                if head['IMAGETYP'] == type :
                    if filter is None or head['FILTER'] == filter :
                        list.append(file)
            except :
                pass
    else :
        list=ascii.read(indir+listfile,Reader=ascii.NoHeader)['col1']
    return list

def window(hdu,box) :
    """
    Reduce size of image and header accordingly
    """
    hdu.header['CRVAL1'] = box.xmin
    hdu.header['CRVAL2'] = box.ymin
    hdu.header['NAXIS1'] = box.ncol()
    hdu.header['NAXIS2'] = box.nrow()
    box.show()
    new = fits.PrimaryHDU(hdu.data[box.ymin:box.ymax+1,box.xmin:box.xmax+1],hdu.header)
    return new

def tv(disp,hd,min=None,max=None) :
    """ Displays HDU or data array on specified ds9/tv device"""
    if isinstance(hd, (np.ndarray)) :
        data=hd
    else :
        data=hd.data

    if type(disp) is pyds9.DS9 :
        disp.set_np2arr(data)
        if min is not None and max is not None :
            disp.set("scale limits {:5d} {:5d}".format(min,max))
    else :
        disp.tv(data,min=min,max=max)

def stretch(a,ncol=None,nrow=None) :
    """ Stretches a 1D image into a 2D image along rows or columns """
    if nrow is None and ncol is None :
        print 'Must specify either nrow= or ncol='
        return
    if nrow is not None and ncol is not None :
        print 'Must specify only one of nrow= or ncol='
        return
    if ncol is not None :
        out=np.zeros([a.shape[0],ncol])
        for i in range(ncol) :
            out[:,i]=a
    if nrow is not None :
        out=np.zeros([nrow,a.shape[0]])
        for i in range(nrow) :
            out[i,:]=a
    return out
