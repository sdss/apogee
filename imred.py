import numpy as np
import astropy
from astropy.io import fits
from astropy.io import ascii
from astropy.modeling import models, fitting
from astropy.convolution import convolve, Box1DKernel
import matplotlib.pyplot as plt
import glob
import bz2
import os
import pdb
import image
try: 
    import pyds9
except:
    print('pyds9 is not available, proceeding')

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

def read(num,ext=0,bias=True,verbose=False) :
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
    if verbose: 
        print('root: ', root)
        print('Reading: ', file)
    if '.bz2' in file :
        hdu=fits.open(bz2.BZ2File(file),ignore=True,ignore_missing_end=True)
    else :
        hdu=fits.open(file,ignore=True,ignore_missing_end=True)

    if 'object' in hdu[0].header.cards and hdu[0].header['object'] == '' : 
        hdu[0].header['object'] = os.path.basename(file)
    if bias :
        if det.biastype == 0 :
            b=det.biasbox.mean(hdu[ext].data)
            if verbose: print('subtracting overscan: ', b)
            hdu[ext].data = hdu[ext].data.astype(float)-b
        elif det.biastype == 1 :
            over=np.median(hdu[ext].data[:,det.biasbox.xmin:det.biasbox.xmax],axis=1)
            boxcar = Box1DKernel(10)
            over=convolve(over,boxcar,boundary='extend')
            over=image.stretch(over,ncol=hdu[ext].data.shape[1])
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
    if verbose: print('reading: ', num)
    hdu=read(num,verbose=verbose)
    if bias is not None :
        if verbose: print('subtracting bias')
        hdu.data -= bias
    if flat is not None :
        if verbose: print('dividing by flat')
        hdu.data /= flat
    if trim :
        out= window(hdu,det.trimbox)
    else :
        out= hdu  
    return out

def combine(ims,norm=False,bias=None,flat=None,trim=False,verbose=False,
            disp=None,min=None,max=None,div=False) :
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
        print('Reading image: ', im)
        h=reduce(im,bias=bias,flat=flat,trim=trim,verbose=verbose) 
        if norm :
            b=det.normbox
            norm=np.median(h.data[b.ymin:b.ymax,b.xmin:b.xmax])
            print('Normalizing image by : ', norm)
            cube.append(h.data/norm)
        else :
            cube.append(h.data)
    print('Combining: ', ims)
    comb = np.median(cube,axis=0)
    if disp is not None :
        for im in ims :
            print(im)
            h=reduce(im,bias=bias,flat=flat,trim=trim,verbose=verbose) 
            if norm :
                b=det.normbox
                norm=np.median(h.data[b.ymin:b.ymax,b.xmin:b.xmax])
                h.data /= norm
            print('Normalizing image by : ', norm)
            if div :
                disp.tv(h.data/comb,min=min,max=max)
            else :
                disp.tv(h.data-comb,min=min,max=max)
            pdb.set_trace()
        disp.tv(comb,min=min,max=max)
        pdb.set_trace()
    return comb
    
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
        print('smoothing by columns not yet implemented!')
        pdb.set_trace()
      
    return smooth 

class DET() :
    """ 
    Defines detector class 
    """
    def __init__(self) :
        self.gain = 0.
        self.rn = 0.
        self.biastype = 0
        self.biasbox = image.BOX()
        self.normbox = image.BOX()
        self.trimbox = image.BOX()
        self.formstr = "{:04d}"

def getdet(idet) :
    """ 
    returns detector object given input detector index
    """
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
    elif idet == 32 :
       # APO 1m Leach
       d.gain=1
       d.rn=7
       d.biastype=0
       d.biasbox.set(2060,2090,10,2000)
       d.normbox.set(900,1200,900,1200)
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

def look(tv,pause=True,files=None,list=None,min=None, max=None) :
    """ 
    Displays series of files 
    """
    if files is None and list is None :
        files=glob.glob(indir+'/*.fits*')
    if list is not None :
        f=open(indir+list,'r')
        files=[]
        for line in f :
            files.extend(np.arange(int(line.split()[0]),int(line.split()[1])))
        f.close()
        pdb.set_trace()
 
    for file in files :
       hd=read(file,verbose=True,bias=False)
       disp(tv,hd,min=min,max=max)
       if pause :
           pdb.set_trace()

def getfiles(type,listfile=None,filter=None,verbose=False) :
    """ 
    Get all files of desired type from specified directory. If file is specified, read numbers from that file, else use IMAGETYP card
    """
    if listfile is None :
        list=[]
        if verbose: print('directory: ', indir)
        for file in glob.glob(indir+'/'+root+'*.fits*') :
            head=fits.open(file)[0].header
            #if verbose :
            #   print('file: ', file)
            #   print('IMAGETYPE: ', head['IMAGETYP'])

            try :
                if head['IMAGETYP'] == type :
                    if filter is None or head['FILTER'] == filter :
                        list.append(file)
            except :
                pass
    else :
        list=ascii.read(indir+listfile,Reader=ascii.NoHeader)['col1']
    return list

def disp(tv,hd,min=None,max=None,sky=False) :
    """ 
    Displays HDU or data array on specified ds9/tv device
    """

    if type(tv) is pyds9.DS9 :
        if isinstance(hd, (np.ndarray)) :
            tv.set_np2arr(hd)
            data=hd
        elif isinstance(hd, (astropy.io.fits.hdu.hdulist.HDUList)) :
            tv.set_pyfits(hd)
            data=hd[0].data
        elif isinstance(hd, (astropy.io.fits.hdu.image.PrimaryHDU)) :
            tv.set_np2arr(hd.data)
            data=hd.data
        else :
            print('Unrecognized data type for display: ',type(hd)) 
        if sky :
           skyval = mmm.mmm(data)
           min = skyval[0]-5*skyval[1]
           max = skyval[0]+20*skyval[1]
        if min is not None and max is not None :
            tv.set("scale limits {:5d} {:5d}".format(min,max))
        else :
            tv.set("scale histequ")

    else :
        if isinstance(hd, (astropy.io.fits.hdu.hdulist.HDUList)) :
            tv.tv(hd[0],min=min,max=max)
        else :
            tv.tv(hd,min=min,max=max)

