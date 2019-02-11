# encoding: utf-8
#
# @Author: Jon Holtzman
# @Date: July 2018
# @Filename: norm.py
# @License: BSD 3-Clause
# @Copyright: Jon Holtzman

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import numpy as np
import matplotlib.pyplot as plt
import copy
import pdb
from astropy.io import fits
from astropy.io import ascii
from apogee.aspcap import aspcap
from apogee.aspcap import ferre
from scipy.ndimage.filters import median_filter
from scipy import interpolate
from tools import plots

def cont(spec,specerr,chips=False,order=4,poly=True) :
    """ Returns continuum normalized spectrum
    """
    x = np.arange(0,len(spec))
   
    if chips :
        cont=np.full_like(spec,np.nan)
        pranges=aspcap.gridPix(apStar=True)
        for prange in pranges :
            s = spec[prange[0]:prange[1]]
            serr = specerr[prange[0]:prange[1]]
            xx = x[prange[0]:prange[1]]
            if poly :
                cont[prange[0]:prange[1]] = polyfit(xx,s,serr,order)
    else :
        if poly :
            cont = polyfit(x,spec,specerr,order)

    return cont

def polyfit(x,y,yerr,order) :
    """ continuum via polynomial fit
    """
    gd = np.where(np.isfinite(y))[0]
    # note unconventional definition in numpy.polyfit for weights!
    p = np.poly1d(np.polyfit(x[gd],y[gd],order,w=1./yerr[gd]))
    return p(x)

def correct(field,libfile,plot=True,write=None,width=151) :
    """ Read raw FERRE files and create correction to normalization based on ratio
    """
    lib=ferre.rdlibhead(libfile+'.hdr')
    out=ferre.read(field,libfile+'.hdr')
    nspec=out[1]['obs'].shape[0]
    norm=copy.copy(out[1]['obs'])
    if plot : fig,ax=plots.multi(1,3)
    for i in range(nspec) :
        p1=0
        for ichip in range(3) :
            npix=lib[1][ichip]['NPIX']
            obs=out[1]['obs'][i,p1:p1+npix]
            # replace bad pixels with median filtered value
            med=median_filter(obs,[5*width],mode='nearest')
            bd=np.where(obs < 0.01)[0]
            obs[bd] = med[bd]
            # get ratio of observed / model, make sure edge is reasonable number
            mdl=out[1]['mdl'][i,p1:p1+npix]
            ratio=obs/mdl
            ratio[0]=np.median(ratio[0:9])
            ratio[-1]=np.median(ratio[-9:0])
            corr=median_filter(obs/mdl,[width],mode='nearest')
            norm[i,p1:p1+npix] = corr
            if plot :
                ax[ichip].cla()
                plots.plotl(ax[ichip],np.arange(npix),obs,yr=[0.95,1.05],color='r')
                plots.plotl(ax[ichip],np.arange(npix),mdl,yr=[0.95,1.05],color='b')
                plots.plotl(ax[ichip],np.arange(npix),obs/mdl,yr=[0.95,1.05],color='g')
                plots.plotl(ax[ichip],np.arange(npix),corr,color='k',linewidth=3,yr=[0.95,1.05])
                plt.draw()
                plt.show()
            p1+=npix
        if plot : pdb.set_trace()
    if write is not None:  
        # make sure we have no zeros
        bd=np.where(norm.flatten() < 0.01)[0]
        norm.flatten()[bd]=1.
        ferre.writespec(write,norm)
    return norm

