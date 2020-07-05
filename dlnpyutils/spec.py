#!/usr/env python

# Imports
import numpy as np
from scipy.special import erf
from scipy.optimize import curve_fit, least_squares
from astropy.io import fits
from scipy import ndimage
from scipy.interpolate import interp1d
from scipy.special import wofz
from . import utils


def trace(im,yestimate=None,yorder=2,sigorder=4,step=10):
    """ Trace the spectrum.  Spectral dimension is assumed to be on the horizontal axis."""
    ny,nx = im.shape
    if yestimate is None:
        ytot = np.sum(im,axis=1)
        yestimate = np.argmax(ytot)
    # Smooth in spectral dimension
    # a uniform (boxcar) filter with a width of 50
    smim = ndimage.uniform_filter1d(im, 50, 1)
    nstep = nx//step
    # Loop over the columns in steps and fit Gaussians
    tcat = np.zeros(nstep,dtype=np.dtype([('x',float),('pars',float,4)]))
    for i in range(nstep):
        pars,cov = dln.gaussfit(y[yestimate-10:yestimate+10],im[yestimate-10:yestimate+10,step*i+step//2])
        tcat['x'][i] = step*i+step//2
        tcat['pars'][i] = pars
    # Fit polynomail to y vs. x and gaussian sigma vs. x
    ypars = np.polyfit(tcat['x'],tcat['pars'][:,1],yorder)
    sigpars = np.polyfit(tcat['x'],tcat['pars'][:,2],sigorder)
    # Model
    mcat = np.zeros(nx,dtype=np.dtype([('x',float),('y',float),('sigma',float)]))
    xx = np.arange(nx)
    mcat['x'] = xx
    mcat['y'] = np.poly1d(ypars)(xx)
    mcat['sigma'] = np.poly1d(sigpars)(xx)
    return tcat, ypars, sigpars,mcat

def boxcar(im):
    """ Boxcar extract the spectrum"""
    ny,nx = im.shape
    ytot = np.sum(im,axis=1)
    yest = np.argmax(ytot)
    # Background subtract
    yblo = int(np.maximum(yest-50,0))
    ybhi = int(np.minimum(yest+50,ny))
    med = np.median(im[yblo:ybhi,:],axis=0)
    medim = np.repeat(med,ny).reshape(ny,nx)
    subim = im.copy()-medim
    # Sum up the flux
    ylo = int(np.maximum(yest-20,0))
    yhi = int(np.minimum(yest+20,ny))
    flux = np.sum(subim[ylo:yhi,:],axis=0)
    return flux

def linefit(x,y,initpar,bounds,err=None):
    # Fit Gaussian profile to data with center and sigma fixed.
    # initpar = [height, center, sigma, constant offset]
    cen = initpar[1]
    sigma = initpar[2]
    def gline(x, amp, const=0):
        """1-D gaussian: gaussian(x, amp, cen, sig)"""
        return amp * np.exp(-(x-cen)**2 / (2*sigma**2)) + const
    line_initpar = [initpar[0],initpar[3]]
    lbounds, ubounds = bounds
    line_bounds = ([lbounds[0],lbounds[3]],[ubounds[0],ubounds[3]])
    return curve_fit(gline, x, y, p0=line_initpar, bounds=line_bounds, sigma=err)


def extract(im,imerr=None,mcat=None,nobackground=False):
    """ Extract a spectrum"""
    ny,nx = im.shape
    x = np.arange(nx)
    y = np.arange(ny)
    # No trace information input, get it
    if mcat is None:
        tcat,ypars,sigpars,mcat=trace(im)
    # Loop over the columns and get the flux using the trace information
    cat = np.zeros(nx,dtype=np.dtype([('x',int),('pars',float,2),('perr',float,2),
                                      ('flux',float),('fluxerr',float)]))
    for i in range(nx):
        line = im[:,i].flatten()
        if imerr is not None:
            lineerr = imerr[:,i].flatten()
        else:
            lineerr = np.ones(len(line))   # unweighted
        # Fit the constant offset and the height of the Gaussian
        #  fix the central position and sigma
        ycen = mcat['y'][i]
        ysigma = mcat['sigma'][i]
        ht0 = np.maximum(line[int(np.round(ycen))],0.01)
        initpar = [ht0,ycen,ysigma,np.median(line)]
        if nobackground is True:
            initpar = [ht0,ycen,ysigma,0]
        # Only fit the region right around the peak
        y0 = int(np.maximum(ycen-50,0))
        y1 = int(np.minimum(ycen+50,ny))
        bnds = ([0,ycen-1e-4,ysigma-1e-4,0],[1.5*ht0,ycen,ysigma,1.5*initpar[3]])
        if nobackground is True:
            bnds = ([0,ycen-1e-4,ysigma-1e-4,0],[1.5*ht0,ycen,ysigma,0.1])
        pars,cov = linefit(y[y0:y1],line[y0:y1],initpar=initpar,bounds=bnds,err=lineerr[y0:y1])
        perr = np.sqrt(np.diag(cov))
        # Gaussian area = ht*wid*sqrt(2*pi)
        flux = pars[0]*ysigma*np.sqrt(2*np.pi)
        fluxerr = perr[0]*ysigma*np.sqrt(2*np.pi)
        cat['x'][i] = i
        cat['pars'][i] = pars
        cat['perr'][i] = perr
        cat['flux'][i] = flux
        cat['fluxerr'][i] = fluxerr
    return cat

def emissionlines(spec,thresh=None):
    """Measure the emission lines in an arc lamp spectrum. """
    nx = len(spec)
    x = np.arange(nx)
    
    # Threshold
    if thresh is None:
        thresh = np.min(spec) + (np.max(spec)-np.min(spec))*0.05
    
    # Detect the peaks
    sleft = np.hstack((0,spec[0:-1]))
    sright = np.hstack((spec[1:],0))
    peaks, = np.where((spec>sleft) & (spec>sright) & (spec>thresh))
    npeaks = len(peaks)
    print(str(npeaks)+' peaks found')
    
    # Loop over the peaks and fit them with Gaussians
    gcat = np.zeros(npeaks,dtype=np.dtype([('x0',int),('x',float),('xerr',float),('pars',float,4),('perr',float,4),
                                           ('flux',float),('fluxerr',float)]))
    resid = spec.copy()
    gmodel = np.zeros(nx)
    for i in range(npeaks):
        x0 = peaks[i]
        xlo = np.maximum(x0-6,0)
        xhi = np.minimum(x0+6,nx)
        initpar = [spec[x0],x0,1,0]
        bnds = ([0,x0-3,0.1,0],[1.5*initpar[0],x0+3,10,1e4])
        pars,cov = dln.gaussfit(x[xlo:xhi],spec[xlo:xhi],initpar,bounds=bnds,binned=True)
        perr = np.sqrt(np.diag(cov))
        gmodel1 = dln.gaussian(x[xlo:xhi],*pars)
        gmodel[xlo:xhi] += (gmodel1-pars[3])
        resid[xlo:xhi] -= (gmodel1-pars[3])
        # Gaussian area = ht*wid*sqrt(2*pi)
        flux = pars[0]*pars[2]*np.sqrt(2*np.pi)
        fluxerr = perr[0]*pars[2]*np.sqrt(2*np.pi)
        gcat['x0'][i] = x0
        gcat['x'][i] = pars[1]
        gcat['xerr'][i] = perr[1]
        gcat['pars'][i] = pars
        gcat['perr'][i] = perr
        gcat['flux'][i] = flux
        gcat['fluxerr'][i] = fluxerr
        
    return gcat, gmodel


def continuum(spec,bin=50,perc=60,norder=4):
    """ Derive the continuum of a spectrum."""
    nx = len(spec)
    x = np.arange(nx)
    # Loop over bins and find the maximum
    nbins = nx//bin
    xbin1 = np.zeros(nbins,float)
    ybin1 = np.zeros(nbins,float)
    for i in range(nbins):
        xbin1[i] = np.mean(x[i*bin:i*bin+bin])
        ybin1[i] = np.percentile(spec[i*bin:i*bin+bin],perc)
    # Fit polynomial to the binned values
    coef1 = np.polyfit(xbin1,ybin1,norder)
    cont1 = np.poly1d(coef1)(x)
    
    # Now remove large negative outliers and refit
    gdmask = np.zeros(nx,bool)
    gdmask[(spec/cont1)>0.8] = True
    xbin = np.zeros(nbins,float)
    ybin = np.zeros(nbins,float)
    for i in range(nbins):
        xbin[i] = np.mean(x[i*bin:i*bin+bin][gdmask[i*bin:i*bin+bin]])
        ybin[i] = np.percentile(spec[i*bin:i*bin+bin][gdmask[i*bin:i*bin+bin]],perc)
    # Fit polynomial to the binned values
    coef = np.polyfit(xbin,ybin,norder)
    cont = np.poly1d(coef)(x)
    
    return cont,coef
