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
from astropy.io import fits
from astropy.io import ascii
from apogee.utils import apload
from apogee.aspcap import aspcap

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
    p = np.poly1d(np.polyfit(x[gd],y[gd],order,w=1./yerr))
    return p(x)

