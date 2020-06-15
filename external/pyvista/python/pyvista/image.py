from __future__ import print_function
import numpy as np
import copy
from tools import plots
from astropy.io import fits
from astropy.io import ascii
from astropy.modeling import models, fitting
from astropy.convolution import convolve, Box1DKernel
from astropy.nddata import StdDevUncertainty, support_nddata
import scipy.signal
import scipy.ndimage
import matplotlib.pyplot as plt
import glob
import bz2
import os
import pdb
try: 
    import pyds9
except:
    print('pyds9 is not available, proceeding')

class BOX() :
    """ 
    Defines BOX class
    """
    def __init__(self,n=None,nr=None,nc=None,sr=1,sc=1,cr=None,cc=None,xr=None,yr=None) :
        """ Define a BOX

            Args :
               n (int) : size of box (if square)
               nr (int) : number of rows 
               nc (int) : number of cols
               sr (int) : start row
               sc (int) : start column
               cr (int) : central row (supercedes sr)
               cc (int) : central column (supercedes sc)
               xr       : [xmin,xmax]  (supercedes cc and sc)
               yr       : [ymin,ymax]  (supercedes cr and sr)
        """
        if nr is None and nc is None and n is None and xr is None and yr is None:
            print('You must specify either n=, or nr= and nc=')
            return
        elif nr is None and nc is None :
            nr=n
            nc=n
        elif nr is None :
            print('You much specify nr= with nc=')
        elif nc is None :
            print('You much specify nc= with nr=')

        if cr is not None and cc is not None :
            sr=cr-nr/2
            sc=cc-nr/2

        if xr is not None :
            self.xmin=xr[0]
            self.xmax=xr[1]
        else :
            self.xmin = sc
            self.xmax = sc+nc-1
        if yr is not None :
            self.ymin=yr[0]
            self.ymax=yr[1]
        else :
            self.ymin = sr
            self.ymax = sr+nr-1

    def set(self,xmin,xmax,ymin,ymax):
        """ Resets limits of a box
        
            Args:
                xmin : lower x value
                xmax : higher x value
                ymin : lower y value
                ymax : higher xyvalue
        """
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

    def nrow(self):
        """ Returns number of rows in a box

            Returns :
                number of rows 
        """
        return(self.ymax-self.ymin+1)

    def ncol(self):
        """ Returns number of columns in a box

            Returns :
                number of columns
        """
        return(self.xmax-self.xmin+1)

    def show(self):
        """ Prints box limits
        """
        print('    SC    NC    SR    NR  Exp       Date     Name')
        print('{:6d}{:6d}{:6d}{:6d} '.format(self.xmin,self.ncol(),self.ymin,self.nrow()))

    def mean(self,data):
        """ Returns mean of data in box

            Args :
                data : input data (CCDData or np.array)
 
            Returns:
                mean of data in box
        """
        if self.nrow() <= 0 or self.ncol() <= 0 : return 0.
        return data[self.ymin:self.ymax+1,self.xmin:self.xmax+1].mean() 

    def stdev(self,data):
        """ Returns standard deviation of data in box

            Args :
                data : input data (CCDData or np.array)

            Returns:
                standard deviation of data in box
        """
        if self.nrow() == 0 or self.ncol() == 0 : return 0.
        return data[self.ymin:self.ymax+1,self.xmin:self.xmax+1].std() 

    def max(self,data):
        """ Returns maximum of data in box

            Args :
                data : input data (CCDData or np.array)

            Returns:
                maximum of data in box
        """
        if self.nrow() == 0 or self.ncol() == 0 : return 0.
        return data[self.ymin:self.ymax+1,self.xmin:self.xmax+1].max() 

    def min(self,data):
        """ Returns minimum of data in box

            Args :
                data : input data (CCDData or np.array)

            Returns:
                minimum of data in box
        """
        if self.nrow() == 0 or self.ncol() == 0 : return 0.
        return data[self.ymin:self.ymax+1,self.xmin:self.xmax+1].min() 

    def median(self,data):
        """ Returns median of data in box

            Args :
                data : input data (CCDData or np.array)

            Returns:
                median of data in box
        """
        if self.nrow() == 0 or self.ncol() == 0 : return 0.
        return np.median(data[self.ymin:self.ymax+1,self.xmin:self.xmax+1])

    def setval(self,data,val):
        """ Sets data in box to specified value
        """
        if self.nrow() == 0 or self.ncol() == 0 : return 0.
        data[self.ymin:self.ymax+1,self.xmin:self.xmax+1] = val

@support_nddata
def abx(data,box) :
    """
    Returns dictionary with image statistics in box.

    Args :
        data  : input data (CCDData or np.array)
        box   : pyvista BOX

    Returns :
        dictionary with image statistics : 'mean', 'stdev', 'min', 'max', 'peakx', 'peaky'
    """
    return {'mean': box.mean(data),
            'stdev': box.stdev(data),
            'max': box.max(data),
            'min': box.min(data),
            'peakx': np.unravel_index(
                        data[box.ymin:box.ymax,box.xmin:box.xmax].argmax(),
                        (box.nrow(),box.ncol()) )[1]+box.xmin,
            'peaky': np.unravel_index(
                        data[box.ymin:box.ymax,box.xmin:box.xmax].argmax(),
                        (box.nrow(),box.ncol()) )[0]+box.ymin}

def gfit(data,x0,y0,size=5,fwhm=3,sub=True,plot=None,fig=1,scale=1,pafixed=False) :
    """ 
    Does gaussian fit to input data given initial xcen,ycen
    """
    fit=fitting.LevMarLSQFitter()
    #fit=fitting.SLSQPLSQFitter()
    z=data[int(y0)-size:int(y0)+size,int(x0)-size:int(x0)+size]
    xcen,ycen=np.unravel_index(np.argmax(z),z.shape)
    xcen+=(int(x0)-size)
    ycen+=(int(y0)-size)
    y,x=np.mgrid[ycen-size:ycen+size,xcen-size:xcen+size]
    z=data[ycen-size:ycen+size,xcen-size:xcen+size]
    g_init=models.Gaussian2D(x_mean=xcen,y_mean=ycen,
                             x_stddev=fwhm/2.354,y_stddev=fwhm/2.354,
                             amplitude=data[ycen,xcen],theta=0.,
                             fixed={'theta':pafixed})+models.Const2D(0.)
    g=fit(g_init,x,y,z)
    xfwhm=g[0].x_stddev*2.354*scale
    yfwhm=g[0].y_stddev*2.354*scale
    fwhm=np.sqrt(xfwhm*yfwhm)
    xcen=g[0].x_mean.value
    ycen=g[0].y_mean.value
    theta=(g[0].theta.value % (2*np.pi)) * 180./np.pi
    print('xFWHM:{:8.2f}   yFWHM:{:8.2f}   FWHM:{:8.2f}  SCALE:{:8.2f}  PA:{:8.2f}'.format(xfwhm,yfwhm,fwhm,scale,theta))
    if plot is not None:
        r = np.sqrt((y-ycen)**2 + (x-xcen)**2)
        plots.plotp(plot,r,z,xt='R(pixels)',yt='Intensity')
        r = np.arange(0.,5*fwhm/2.354/scale,0.1)
        peak=g[0].amplitude
        plot.plot(r,peak*np.exp(-np.power(r, 2.) / (2 * np.power(g[0].x_stddev, 2.)))+g[1].amplitude)
        plot.plot(r,peak*np.exp(-np.power(r, 2.) / (2 * np.power(g[0].y_stddev, 2.)))+g[1].amplitude)
        plot.text(0.9,0.9,'x: {:7.1f} y: {:7.1f} fw: {:8.2f}'.format(xcen,ycen,fwhm),transform=plot.transAxes,ha='right')
        plt.draw()
       
    if sub :
        out=data
        out[ycen-size:ycen+size,xcen-size:xcen+size]-=g[0](x,y)
        return out
    return g

def tvstar(tv,plot=None,size=11,fwhm=5,scale=1,pafixed=False) :
    """ Fit gaussian and show radial profile of stars marked interactively
    """
    key=''
    print('Hit key near star center, "q" to quit')
    while key != 'q' :
        key,x,y=tv.tvmark()
        tv.plotax1.cla()
        gfit(tv.img,x,y,size=size,fwhm=fwhm,scale=scale,plot=tv.plotax1,sub=False,pafixed=pafixed)

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

def stretch(a,ncol=None,nrow=None) :
    """ 
    Stretches a 1D image into a 2D image along rows or columns 
    """
    if nrow is None and ncol is None :
        print('Must specify either nrow= or ncol=')
        return
    if nrow is not None and ncol is not None :
        print('Must specify only one of nrow= or ncol=')
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

def __get_cnpix(a) :
    """
    Gets CNPIX cards from HDU a, sets them to 1 if they don't exist
    """
    try:
        cnpix1=a.header['CNPIX1']
    except:
        a.header['CNPIX1']=1
    try:
        cnpix2=a.header['CNPIX2']
    except:
        a.header['CNPIX2']=1

    return a.header['CNPIX1'],a.header['CNPIX2']

def __get_overlap(a,b,dc=0,dr=0,box=None) :
    """
    Returns overlap coordinates from two input HDUs
    """
    
    if box is not None :
        print('need to implement box=!')
        return
 
    a_cnpix1,a_cnpix2 = __get_cnpix(a)
    b_cnpix1,b_cnpix2 = __get_cnpix(b)

    ixmin = max(a_cnpix1,b_cnpix1+dc)
    iymin = max(a_cnpix2,b_cnpix2+dr)
    ixmax = min(a_cnpix1+a.header['NAXIS1'],b_cnpix1+dc+b.header['NAXIS1'])
    iymax = min(a_cnpix2+a.header['NAXIS2'],b_cnpix2+dr+b.header['NAXIS2'])

    return (iymin-a_cnpix2,iymax-a_cnpix2,ixmin-a_cnpix1,ixmax-a_cnpix1,
           iymin-b_cnpix2-dr,iymax-b_cnpix2-dr,ixmin-b_cnpix1-dc,ixmax-b_cnpix1-dc)

def __check_hdu(a) :
    """
    Checks if input variable is an HDU
    """

    if type(a) is fits.hdu.image.PrimaryHDU :
        return True
    else:
        print('Input must be HDU type, with header and data!')
        return False

def add(a,b,dc=0,dr=0,box=None) :
    """ 
    Adds b to a, paying attention to CNPIX 
    """
    if __check_hdu(a) is False or __check_hdu(b) is False : return
    ay1,ay2,ax1,ax2,by1,by2,bx1,bx2 = __get_overlap(a,b,dr=dr,dc=dc,box=box)
    a.data[ay1:ay2,ax1:ax2] += b.data[by1:by2,bx1:bx2]

def sub(a,b,dc=0,dr=0,box=None) :
    """ 
    Subracts b from a, paying attention to CNPIX 
    """
    if __check_hdu(a) is False or __check_hdu(b) is False : return
    ay1,ay2,ax1,ax2,by1,by2,bx1,bx2 = __get_overlap(a,b,dr=dr,dc=dc,box=box)
    a.data[ay1:ay2,ax1:ax2] -= b.data[by1:by2,bx1:bx2]

def mul(a,b,dc=0,dr=0,box=None) :
    """ 
    Multiplies b by a, paying attention to CNPIX 
    """
    if __check_hdu(a) is False or __check_hdu(b) is False : return
    ay1,ay2,ax1,ax2,by1,by2,bx1,bx2 = __get_overlap(a,b,dr=dr,dc=dc,box=box)
    a.data[ay1:ay2,ax1:ax2] *= b.data[by1:by2,bx1:bx2]

def div(a,b,dc=0,dr=0,box=None) :
    """ 
    Divides a by b, paying attention to CNPIX 
    """
    if __check_hdu(a) is False or __check_hdu(b) is False : return
    ay1,ay2,ax1,ax2,by1,by2,bx1,bx2 = __get_overlap(a,b,dr=dr,dc=dc,box=box)
    a.data[ay1:ay2,ax1:ax2] /= b.data[by1:by2,bx1:bx2]

def clip(hd,min=None,max=None,vmin=None,vmax=None,box=None) :
    """
    Clipping tasks: sets all values above or below input values to specified values

    Args:
         hd : input HDU

    Keyword args:
         min=  (float) : clip values below min
         vmin= (float) : values to clip min values to. If min= is not given clips values <vmin to vmin
         max=  (float) : clip values above max
         vmax= (float) : values to clip max values to. If max= is not given clips values >vmax to vmax
    """
    if __check_hdu(hd) is False : return
   
    if box is not None :
        print('need to implement box=!')
        return
 
    if min is not None or vmin is not None :
        if vmin is None: 
            clipval=0
        else :
            clipval=vmin
        if min is None:
            min=vmin
        iy,ix=np.where(hd.data > min)
        hd.data[iy,ix]=clipval

    if max is not None or vmax is not None :
        if vmax is None: 
            clipval=0
        else :
            clipval=vmax
        if max is None:
            max=vmax
        iy,ix=np.where(hd.data > max)
        hd.data[iy,ix]=clipval

def buf(hd) :
    """
    Display information about HDU
    """ 
    if __check_hdu(hd) is False : return

    print('    SC    NC    SR    NR  Exp       Date     Name')
    cnpix1,cnpix2 = __get_cnpix(hd)
    npix1 = hd.header['NAXIS1']
    npix2 = hd.header['NAXIS2']
    print('{:6d}{:6d}{:6d}{:6d}'.format(cnpix1,npix1,cnpix2,npix2))

    #dict=globals()
    #for key in dict :
    #    if type(dict[key]) is fits.hdu.image.PrimaryHDU : print(key)


def rd(file,ext=0) :
    """
    Read file into HDU
    """
    try:
        return fits.open(file)[ext]
    except :
        print('cannot open file: ', file, ' extension: ', ext)

def create(box=None,n=None,nr=None,nc=None,sr=1,sc=1,cr=None,cc=None,const=None) :
    """
    Creates a new HDU
    """
    if box is not None:
        nr=box.nrow()
        nc=box.ncol()
        sc=box.xmin
        sr=box.ymin
    else :
        if nr is None and nc is None :
            try :
                nr=n
                nc=n
            except:
                print('You must specify either box=, n=, or nr= and nc=')
                return
        if cr is not None and cc is not None :
            sr=cr-nr/2
            sc=cc-nr/2
    try :
        im=np.zeros([nr,nc])
    except :
        print('must specify image size ')
        return
    hd=fits.PrimaryHDU(im)
    hd.header['CNPIX1'] = sc
    hd.header['CNPIX2'] = sr
    if const is not None :
        hd.data += const
    return hd

def sky(im,box=None,max=None,min=None,plot=None):
    """
    Estimate sky value in an image by fitting parabola to peak of histogram

    Args:
        im (HDU or numpy array): input image data 

    Keyword args:
        box=   : only use values within specified box (default=None)
        min=   : ignore values below min in sky computation (default=None)
        max=   : ignore values above max in sky computation (default=None)
        plot=  : matplotlib axes to view histogram and fit (default=None)
    """

    if type(im) is fits.hdu.image.PrimaryHDU :
        data = im.data
    else :
        data = im
    if box is not None :
        reg = data[box.ymin:box.ymax+1,box.xmin:box.xmax+1]
    else :
        reg = data

    if min is None: min = reg.min()
    if max is None: max = reg.max()
    if min > max :
        raise ValueError("min must be less than max")

    gd = np.where((reg >min) & (reg<max))
    if len(gd[0]) < 1 :
        raise ValueError("no pixels between min and max")

    # get median and stdev in desired region
    med = np.median(reg[gd])
    sig = reg[gd].std()
    print('initial median, sigma: ', med, sig)

    # create histogram around median and find peak
    gd = np.where((reg.flatten() > med-2*sig) & (reg.flatten() < med+2*sig))[0]
    hist,bins = np.histogram(reg.flatten()[gd],bins=np.arange(med-2*sig,med+2*sig))
    max = np.max(hist)
    imax = np.argmax(hist)

    # find half power points on either side of peak
    i1=imax
    while hist[i1] > max/2. and i1 > 0 :
        i1-=1
    i2=imax
    while hist[i2] > max/2. and i2 < len(hist) :
        i2+=1

    # fit parabola to peak, and determine location of fit max
    binwidth=bins[1]-bins[0]
    p_init=models.Polynomial1D(degree=2)
    fit=fitting.LinearLSQFitter()
    p=fit(p_init,bins[i1:i2+1]+binwidth,hist[i1:i2+1])
    sky=-p.parameters[1]/(2.*p.parameters[2])
    if plot is not None:
        plot.plot(bins[i1:i2+1]+binwidth,hist[i1:i2+1])
        plot.plot(bins[i1:i2+1]+binwidth,p(bins[i1:i2+1]+binwidth))
        plt.draw()

    return sky

def getdata(hd) :

    if isinstance(hd, (np.ndarray)) :
        data=hd
    elif isinstance(hd, (astropy.io.fits.hdu.hdulist.HDUList)) :
        data=hd[0].data
    elif isinstance(hd, (astropy.io.fits.hdu.image.PrimaryHDU)) :
        data=hd.data
    else :
        print('Unrecognized data type: ',type(hd))
    return(data)

def xcorr(a,b,lags,medfilt=0) :
    """ Cross correlation function between two arrays, calculated at lags

        Args:
            a, b : input 1D arrays
            lags : array (1D) of x-corrlation lags
            medfilt : size of median filter for arrays (default=0)

        Returns :
            fit peak of cross-correlation (quadratic fit)
            1D cross-correlation function
    """

    # compute xcorr with starting and ending position to allow full range of lags
    xs = -lags[0]
    xe = np.min([a.shape[-1],b.shape[-1]])-lags[-1]
    atmp=np.atleast_2d(a)
    btmp=np.atleast_2d(b)
    if medfilt>0 :
        atmp=np.atleast_2d(atmp-scipy.signal.medfilt(a,kernel_size=[1,medfilt]))
        btmp=np.atleast_2d(btmp-scipy.signal.medfilt(b,kernel_size=[1,medfilt]))

    if atmp.shape[0] == btmp.shape[0] :
        shift=np.zeros([1,len(lags)])
        for i,lag in enumerate(lags) :
            shift[0,i]=np.sum(atmp[:,xs:xe]*btmp[:,xs+lag:xe+lag])
    elif atmp.shape[0] == 1 :
        shift=np.zeros([btmp.shape[0],len(lags)])
        for row in range(btmp.shape[0]) :
            print('cross correlating row: {:d}'.format(row),end='\r')
            for i,lag in enumerate(lags) :
                shift[row,i]=np.sum(atmp[0,xs:xe]*btmp[row,xs+lag:xe+lag])
    else:
        raise ValueError('input arrays must have same nrows, or first must have 1 row')
        return
    fitpeak=np.zeros(shift.shape[0])
    for row in range(shift.shape[0]) :
        peak=shift[row,:].argmax()
        try :
            fit=np.polyfit(range(-3,4),shift[row,peak-3:peak+4],2)
            fitpeak[row]=peak+-fit[1]/(2*fit[0])
        except TypeError :
            fitpeak[row]=peak

    return fitpeak,np.squeeze(np.array(shift))

def xcorr2d(a,b,lags) :
    """ Two-dimensional cross correlation

        Args:
            a, b : input CCDData frames
            lags : array (1D) of x-corrlation lags

        Returns:
            (x,y) position of cross correlation peak from quadratic fit to x-correlation
            2D cross correlation function
    """
    # do x-corrlation over section of image that fits within input lag array
    xs = -lags[0]
    xe = np.min([a.shape[1],b.shape[1]])-lags[-1]
    ys = -lags[0]
    ye = np.min([a.shape[0],b.shape[0]])-lags[-1]

    # compute x-corrleation
    shift = np.zeros([len(lags),len(lags)])
    for i, xlag in enumerate(lags) :
        for j, ylag in enumerate(lags) :
            shift[j,i] = np.sum(a.data[ys:ye,xs:xe]*b.data[ys+ylag:ye+ylag,xs+xlag:xe+xlag])

    # quadratic fit and determine peak
    fit=fitting.LinearLSQFitter()
    mod=models.Polynomial2D(degree=2)
    y,x=np.meshgrid(lags,lags)
    yp,xp=np.unravel_index(shift.argmax(),shift.shape)
    print(yp,xp)
    p=fit(mod,x[yp-1:yp+2,xp-1:xp+2],y[yp-1:yp+2,xp-1:xp+2],shift[yp-1:yp+2,xp-1:xp+2])
    a = np.array([ [2*p.parameters[2], p.parameters[5]], [p.parameters[5],2*p.parameters[4]] ])
    b = np.array([-p.parameters[1],-p.parameters[3]])
    peak=np.linalg.solve(a,b)+(xp,yp)+(lags[0],lags[0])

    return peak,shift


def zap(hd,size,nsig=3,mask=False) : 
    """ Median filter array and replace values > nsig*uncertainty
    """
    filt=scipy.signal.medfilt(hd.data,size)
    if nsig >= 0 : bd = np.where(np.atleast_2d(hd.data)-filt > nsig*hd.uncertainty.array)
    else : bd = np.where(np.atleast_2d(hd.data)-filt < nsig*hd.uncertainty.array)
    np.atleast_2d(hd.data)[bd[0],bd[1]] = np.atleast_2d(filt)[bd[0],bd[1]]
    if mask :
        if hd.mask is None : hd.mask=np.zeros(hd.data.shape,dtype=bool)
        hd.mask[bd[0],bd[1]] = True

def smooth(hd,size) :
    """ Boxcar smooth image
    """
    hd.data=scipy.ndimage.uniform_filter(hd.data,size=size)
    npix=1
    for dim in size : npix*=dim
    hd.uncertainty=StdDevUncertainty(np.sqrt(scipy.ndimage.uniform_filter(hd.uncertainty.array**2,size=size)/npix))
