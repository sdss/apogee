#!/usr/bin/env python
#
# PLOTTING.PY - plotting functions
#

from __future__ import print_function

__authors__ = 'David Nidever <dnidever@noao.edu>'
__version__ = '20200209'  # yyyymmdd

#import os
#import sys
import numpy as np
import warnings
#from scipy.signal import medfilt
#from scipy.ndimage.filters import median_filter,gaussian_filter1d
from astropy.utils.exceptions import AstropyWarning
#from scipy.ndimage.filters import convolve
#from . import utils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy import stats
import copy
import dlnpyutils.utils as dln

# Ignore these warnings, it's a bug
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

def zscaling(im,contrast=0.25,nsample=500000):
    """
    This finds the IRAF display z1,z2 scalings
    for an image
    
    Parameters:
    im         2D image
    =contrast  Contrast to use (contrast=0.25 by default)
    =nsample   The number of points of the image to use for the sample
                (nsample=5e4 by default).

    Returns:
     [z1,z1]         The minimum and maximum value to use for display scaling.

    Examples:
      zcals = zscale(im)


    From IRAF display help

    If  the  contrast  is  not  zero  the  sample  pixels  are ranked in
    brightness to form the function I(i) where i  is  the  rank  of  the
    pixel  and  I is its value.  Generally the midpoint of this function
    (the median) is very near the peak of the image histogram and  there
    is  a  well defined slope about the midpoint which is related to the
    width of the histogram.  At the ends of the I(i) function there  are
    a  few very bright and dark pixels due to objects and defects in the
    field.  To determine  the  slope  a  linear  function  is  fit  with
    iterative rejection;
    
            I(i) = intercept + slope * (i - midpoint)
    
    If  more  than half of the points are rejected then there is no well
    defined slope and the full range of the sample defines  z1  and  z2.
    Otherwise  the  endpoints  of the linear function are used (provided
    they are within the original range of the sample):
    
            z1 = I(midpoint) + (slope / contrast) * (1 - midpoint)
            z2 = I(midpoint) + (slope / contrast) * (npoints - midpoint)

    The actual IRAF program is called zsc_zlimits and is in the file
    /net/astro/iraf2.12/iraf/pkg/images/tv/display/zscale.x


    By D.Nidever  Oct 2007  (using IRAF display zscale algorithm) 
    """

    if im.ndim != 2:
        raise ValueError('The input must be 2 dimensiona')
    nx,ny = im.shape
    n = nx*ny

    nsample = np.minimum(nsample, n)

    xind = np.round(np.random.random(nsample)*(nx-1)).astype(int)
    yind = np.round(np.random.random(nsample)*(ny-1)).astype(int)    

    f = im[xind,yind]
    si = np.argsort(f)
    f2 = f[si]
    x = np.arange(nsample)
    midpoint = np.round(nsample*0.5)
    med = np.median(f)

    zmin = np.min(im)
    zmax = np.max(im)
    zmed = np.median(im)

    # Robust fitting program
    coef = dln.poly_fit(x,f2,1,robust=True)

    # y = m*x + b
    # I = intercept + slope * (i-midpoint)
    # I = intercept + slope * i - slope*midpoint
    # I = slope * i + (intercept - slope*midpoint)
    # b = intercept - slope*midpoint
    # intercept = b + slope*midpoint
    
    slope = coef[0]
    intercept = coef[1] + slope*midpoint

    #    z1 = I(midpoint) + (slope / contrast) * (1 - midpoint)
    #    z2 = I(midpoint) + (slope / contrast) * (npoints - midpoint)
    #z1 = f2[midpoint] + (slope/contrast) * (1L - midpoint)
    #z2 = f2[midpoint] + (slope/contrast) * (nsample - midpoint)
    z1 = zmed + (slope/contrast) * (1 - midpoint)
    z2 = zmed + (slope/contrast) * (nsample - midpoint)

    z1 = np.maximum(z1, zmin)
    z2 = np.minimum(z2, zmax)
    
    return [z1,z2]

    

def hist2d(x,y,z=None,statistic='count',xrange=None,yrange=None,dx=None,dy=None,nx=200,ny=200,
           zscale=None,log=None,noerase=False,zmin=None,zmax=None,center=True,force=True,cmap=None,
           xtitle=None,ytitle=None,title=None,origin='lower',aspect='auto'):
    """
    Make a 2D histogram
    
    Parameters:
    x
    y
    z
    statistic
    xrange
    yrange
    dx
    dy
    nx
    ny
    zscale
    log
    noerase
    zmin
    zmax
    center
    force
    cmap
    xtitle
    ytitle
    title
    origin
    aspect

    """

    # Getting the current figure, creating a new one if necessary
    fig = plt.gcf()
    
    if noerase is False:
        plt.clf()   # clear the plot

    # Input dx/dy
    if dx is None:
        dx0 = None
    else:
        dx0 = copy.copy(dx)
    if dy is None:
        dy0 = None
    else:
        dy0 = copy.copy(dy)
        
    # Min and Max's
    xmin = np.min(x)
    xmax = np.max(x)
    ymin = np.min(y)
    ymax = np.max(y)    

    # Temporary DX and DY if not input
    # X-axis
    if (dx is None):
        # XRANGE set
        if xrange is not None:
            dx = (xrange[1]-xrange[0])/(nx-1)
        # XRANGE not set
        else:
            dx = (xmax-xmin)/(nx-1)
    # Y-axis
    if (dy is None):
        # YRANGE set
        if yrange is not None:
            dy = (yrange[1]-yrange[0])/(ny-1)
        # YRANGE not set
        else:
            dy = (ymax-ymin)/(ny-1)

    # Axes reversed
    if xrange is not None:
        # X axis reversed
        if xrange[1]<xrange[0]:
            xrange = np.flip(xrange)
            xflip = 1
    if yrange is not None:
        # Y axis reversed
        if yrange[1]<yrange[0]:
            yrange = np.flip(yrange)
            yflip = 1



    # Cutting based on the xrange/yrange if set
    #   Otherwise the image size might get too large if there are
    #   points far outside the xrange/yrange
    # Xrange set
    if xrange is not None:
        mask1 = ((x>=(xrange[0]-dx)) & (x<=(xrange[1]+dx)))
        if np.sum(mask1)>0:
            x = x[mask1]
            y = y[mask1]
            if z is not None:
                z = z[mask1]
        else:
            raise ValueError('No points left')

        # Min and Max's
        ymin = np.min(y)
        ymax = np.max(y)
        xmin = np.min(x)
        xmax = np.max(x)

        # Setting DX/DY if not input
        if dx0 is None: dx=(xmax-xmin)/(nx-1)
        if dy0 is None: dy=(ymax-ymin)/(ny-1)
    # yrange set
    if yrange is not None:
        mask2 = ((y>=(yrange[0]-dy)) & (y<=(yrange[1]+dy)))
        if np.sum(mask2)>0:
            x = x[mask2]
            y = y[mask2]
            if z is not None:
                z = z[mask2]
        else:
            raise ValueError('No points left')

        # Min and Max's
        ymin = np.min(y)
        ymax = np.max(y)
        xmin = np.min(x)
        xmax = np.max(x)

        # Setting DX/DY if not input
        if dx0 is None: dx=(xmax-xmin)/(nx-1)
        if dy0 is None: dy=(ymax-ymin)/(ny-1)


    # Step sizes must be positive
    dx = np.abs(dx)
    dy = np.abs(dy)


    # Setting the range
    if xrange is not None:
        # Default is that the xrange/yrange are for the centers
        # Move xmin/xmax accordingly
        if center is True:
            off = 0.5*dx
        else:
            off = 0.0

        # xmin and xmax must be xrange+/-integer*dx
        if force is True:
            # If xmin < xrange[0] then move it down an integer number of dx's
            if (xmin < (xrange[0]-off)):
                diff = (xrange[0]-off)-xmin
                xmin = (xrange[0]-off)-np.ceil(diff/dx)*dx
            else:
                xmin = xrange[0]-off
            # If xmax > xrange[1] then move it up an integer number of dx's
            if (xmax > (xrange[1]+off)):
                diff = xmax-(xrange[1]+off)
                xmax = (xrange[1]+off)+np.ceil(diff/dx)*dx
            else:
                xmax = xrange[1]+off

        # Don't force the xrange
        else:
            if ((xrange[0]-off) < xmin) | ((xrange[1]+off) > xmax):
                xmin = np.minimum((xrange[0]-off), xmin)
                xmax = np.maximum((xrange[1]+off), xmax)

    # Setting the range
    if yrange is not None:
        # Default is that the xrange/yrange are for the centers
        # Move ymin/ymax accordingly
        if center is True:
            off = 0.5*dx
        else:
            off = 0.0

        # xmin and xmax must be xrange+/-integer*dx
        if force is True:
            # If ymin < yrange[0] then move it down an integer number of dy's
            if (ymin < (yrange[0]-off)):
                diff = (yrange[0]-off)-ymin
                ymin = (yrange[0]-off)-np.ceil(diff/dy)*dy
            else:
                ymin = yrange[0]-off
            # If ymax > yrange[1] then move it up an integer number of dy's
            if (ymax > (yrange[1]+off)):
                diff = ymax-(yrange[1]+off)
                ymax = (yrange[1]+off)+np.ceil(diff/dy)*dy
            else:
                ymax = yrange[1]+off

        # Don't force the xrange
        else:
            if ((yrange[0]-off) < ymin) | ((yrange[1]+off) > ymax):
                ymin = np.minimum((yrange[0]-off), ymin)
                ymax = np.maximum((yrange[1]+off), ymax)


    # Setting final DX/DY and NX/NY
    # X-axis
    if dx0 is None:
        dx = (xmax-xmin)/(nx-1.)
    else:
        nx = np.floor((xmax-xmin)/dx)+ 1  # only want bins fully within the range
    # Y-axis
    if dy0 is None:
        dy = (ymax-ymin)/(ny-1)
    else:
        ny = np.floor((ymax-ymin)/dy)+ 1  # only want bins fully within the range


    # Final xrange/yrange, if not already set
    if xrange is None:
        xrange = [xmin,xmax]
    if yrange is None:
        yrange = [ymin,ymax]
        
                
    # No z input
    if z is None:
        im, xedges, yedges = np.histogram2d(x,y,range=[xrange,yrange],bins=[nx,ny])
    else:
        im, xedges, yedges, binnumber = stats.binned_statistic_2d(x,y,z,statistic=statistic,range=[xrange,yrange],bins=[nx,ny])

    # Plot the image
    #fig, ax = plt.subplots()
    norm = None
    if log is True:
        norm = colors.LogNorm(vmin=zmin,vmax=zmax)
    if zscale is True:
        zmin,zmax = zscaling(im)
    extent = [xrange[0], xrange[1], yrange[0], yrange[1]]
    plt.imshow(im,cmap=cmap,norm=norm,aspect=aspect,vmin=zmin,vmax=zmax,origin=origin,extent=extent)
        
    # Axis titles
    if xtitle is not None:
        plt.xlabel(xtitle)
    if ytitle is not None:
        plt.ylabel(ytitle)        
    if title is not None:
        plt.title(title)
        
    # Add the colorbar
    plt.colorbar()

    return

def display(im,x=None,y=None,log=False,xrange=None,yrange=None,noerase=False,zscale=False,zmin=None,zmax=None,
            xtitle=None,ytitle=None,title=None,origin='lower',aspect='auto',cmap=None):
    """ Display an image."""

    # Getting the current figure, creating a new one if necessary
    fig = plt.gcf()
    
    if noerase is False:
        plt.clf()   # clear the plot

    nx,ny = im.shape        
        
    # No X/Y inputs
    if x is None:
        x = np.arange(nx)
    if y is None:
        y = np.arange(ny)

    # Min's and Max's
    xmin = np.min(x)
    xmax = np.max(x)
    ymin = np.min(y)
    ymax = np.max(y)

    # Trim based on xrange/yrange
    if xrange is not None:
        val0,xlo = dln.closest(x,xrange[0])
        val1,xhi = dln.closest(x,xrange[1])        
        im = im[xlo:xhi+1,:]
        x = x[xlo:xhi+1]
        nx,ny = im.shape
        xmin = np.min(x)
        xmax = np.max(x)
    if yrange is not None:
        val0,ylo = dln.closest(y,yrange[0])
        val1,yhi = dln.closest(y,yrange[1])        
        im = im[:,ylo:yhi+1]
        y = y[ylo:yhi+1]
        nx,ny = im.shape        
        ymin = np.min(y)
        ymax = np.max(y)
        
    # Plot the image
    #fig, ax = plt.subplots()
    norm = None
    if log is True:
        norm = colors.LogNorm(vmin=zmin,vmax=zmax)
    if zscale is True:
        zmin,zmax = zscaling(im)

    extent = [xmin, xmax, ymin, ymax]
    plt.imshow(im,cmap=cmap,norm=norm,aspect=aspect,vmin=zmin,vmax=zmax,origin=origin,extent=extent)
        
    # Axis titles
    if xtitle is not None:
        plt.xlabel(xtitle)
    if ytitle is not None:
        plt.ylabel(ytitle)        
    if title is not None:
        plt.title(title)
        
    # Add the colorbar
    plt.colorbar()
    
    return


def plot(x,y,z=None,marker=None,log=False,noerase=False,zmin=None,zmax=None,linewidth=None,
         xtitle=None,ytitle=None,title=None,cmap=None,alpha=None):
    """ Create a line or scatter plot.  like plotc.pro"""

    # xerr, yerr, symbol size
    
    # Getting the current figure, creating a new one if necessary
    fig = plt.gcf()
    
    if noerase is False:
        plt.clf()   # clear the plot

    # Make the plot
    norm = None
    if log is True:
        norm = colors.LogNorm(vmin=zmin,vmax=zmax)
    plt.scatter(x,y,c=z,marker=marker,cmap=cmap,norm=norm,vmin=zmin,vmax=zmax,alpha=alpha,linewidth=linewidth)

    # Axis titles
    if xtitle is not None:
        plt.xlabel(xtitle)
    if ytitle is not None:
        plt.ylabel(ytitle)        
    if title is not None:
        plt.title(title)
    
    # Add the colorbar
    if z is not None:
        plt.colorbar()
    
    return
