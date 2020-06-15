# routines for simulating images

import numpy as np
import pdb
from astropy.nddata import CCDData

def gauss2d(data, coords,fwhm=1,back=0,noise=False,rn=0) :
    """ Add 2d gaussians to data given input coords, fwhm, noise"""
    if type(coords[0]) is not list : coords=[coords]
    ypix,xpix = np.mgrid[0:data.shape[0],0:data.shape[1]]
    sig2 = (fwhm/2.354)**2
    for coord in coords: 
        amp = coord[0]/2./np.pi/sig2
        x = coord[1]
        y = coord[2]
        dist2 = (xpix-x)**2 + (ypix-y)**2
        gd = np.where(dist2 < 100*sig2)
        data[gd[0],gd[1]] += amp*np.exp(-dist2[gd[0],gd[1]]/(2.*sig2))
    data = data.astype(float) + back

    if noise:
        data = np.random.poisson(data)
        if rn>0 : data = data.astype(float) +  \
                         rn*np.random.normal(size=data.shape)

    return CCDData(data,unit='photon')
