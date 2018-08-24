# encoding: utf-8
#
# @Author: Jon Holtzman
# @Date: March 2018
# @Filename: spectra.py
# @License: BSD 3-Clause
# @Copyright: Jon Holtzman

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import numpy as np
import copy

# utility routines for working with spectra

def fits2vector(header,axis) :
    """ Routine to return vector of axis values from a FITS header CRVAL, CDELT, NAXIS for specified axis
    """
    caxis='{:1d}'.format(axis)
    return header['CRVAL'+caxis]+header['CDELT'+caxis]*np.arange(header['NAXIS'+caxis])

def vector(start,delta,n) :
    """ Routine to return vector of values given start, delta, n
    """
    return float(start)+np.arange(int(n))*float(delta)

def add_dim(header,crval,cdelt,crpix,ctype,idim) :
    """ Add a set of CRVAL/CDELT,CRPIX,CTYPE cards to header
    """
    header.append(('CRVAL{:d}'.format(idim),crval))
    header.append(('CDELT{:d}'.format(idim),cdelt))
    header.append(('CRPIX{:d}'.format(idim),crpix))
    header.append(('CTYPE{:d}'.format(idim),ctype))

def vactoair(wave_vac) :
    """ Convert vacuum wavelengths to air wavelengths

        Corrects for the index of refraction of air under standard conditions.  
        Wavelength values below 2000 A will not be altered.  Accurate to about 10 m/s.

        From IDL Astronomy Users Library, which references Ciddor 1996 Applied Optics 35, 1566
    """
    if not isinstance(wave_vac, np.ndarray) : 
        vac = np.array([wave_vac])
    else :
        vac = wave_vac

    air = copy.copy(vac)
    g = np.where(vac >= 2000)[0]     #Only modify above 2000 A
    sigma2 = (1.e4/vac[g] )**2.       #Convert to wavenumber squared

    # Compute conversion factor
    fact = 1. +  5.792105E-2/(238.0185E0 - sigma2) + 1.67917E-3/( 57.362E0 - sigma2)
    
    air[g] = vac[g]/fact
    return air 

def airtovac(wave_air) :
    """ Convert air wavelengths to vacuum wavelengths

        Corrects for the index of refraction of air under standard conditions.  
        Wavelength values below 2000 A will not be altered.  Accurate to about 10 m/s.

        From IDL Astronomy Users Library, which references Ciddor 1996 Applied Optics 35, 1566
    """
    if not isinstance(wave_air, np.ndarray) : 
        air = np.array([wave_air])
    else :
        air = wave_air

    vac = copy.copy(air)
    g = np.where(vac >= 2000)[0]     #Only modify above 2000 A

    for iter in range(2) :
        sigma2 = (1e4/vac[g])**2.     # Convert to wavenumber squared
        # Compute conversion factor
        fact = 1. +  5.792105E-2/(238.0185E0 - sigma2) + 1.67917E-3/( 57.362E0 - sigma2)

        vac[g] = air[g]*fact              #Convert Wavelength
    return vac
