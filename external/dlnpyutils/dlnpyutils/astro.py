#!/usr/bin/env python
#
# ASTRO.PY - astronomy utility functions
#

from __future__ import print_function

__authors__ = 'David Nidever <dnidever@noao.edu>'
__version__ = '20191226'  # yyyymmdd

#import re
#import logging
#import os
#import sys
import numpy as np
import warnings
#from astropy.io import fits
#from astropy.table import Table, Column
#from astropy import modeling
#from astropy.convolution import Gaussian1DKernel, convolve
#from glob import glob
#from scipy.signal import medfilt
#from scipy.ndimage.filters import median_filter,gaussian_filter1d
#from scipy.optimize import curve_fit, least_squares
#from scipy.special import erf
#from scipy.interpolate import interp1d
#from scipy.linalg import svd
from astropy.utils.exceptions import AstropyWarning
#import socket
#from scipy.signal import convolve2d
#from scipy.ndimage.filters import convolve
#import astropy.stats
from . import utils


# Ignore these warnings, it's a bug
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

# Convert wavelengths in air to vacuum
def airtovac(wave):
    """
    Convert air wavelengths to vacuum wavelengths 

    Wavelengths are corrected for the index of refraction of air under 
    standard conditions.  Wavelength values below 2000 A will not be 
    altered.  Uses relation of Ciddor (1996).

    INPUT/OUTPUT:
      WAVE_AIR - Wavelength in Angstroms, scalar or vector
              If this is the only parameter supplied, it will be updated on
              output to contain double precision vacuum wavelength(s). 

    EXAMPLE:
      If the air wavelength is  W = 6056.125 (a Krypton line), then 
      AIRTOVAC, W yields an vacuum wavelength of W = 6057.8019

    METHOD:
	Formula from Ciddor 1996, Applied Optics 62, 958

    NOTES: 
      Take care within 1 A of 2000 A.   Wavelengths below 2000 A *in air* are
      not altered.       
    REVISION HISTORY
      Written W. Landsman                November 1991
      Use Ciddor (1996) formula for better accuracy in the infrared 
          Added optional output vector, W Landsman Mar 2011
      Iterate for better precision W.L./D. Schlegel  Mar 2011
    """

    wave_air = np.atleast_1d(wave).copy()  # makes sure it's an array
    wave_vac = np.atleast_1d(wave).copy()  # initialize
    
    g,ng = utils.where(wave_vac >= 2000)     #Only modify above 2000 A
    
    if ng>0:
        for iter in range(2):
            sigma2 = (1e4/wave_vac[g] )**2     # Convert to wavenumber squared
            
            # Compute conversion factor
            fact = 1.0 +  5.792105e-2/(238.0185e0 - sigma2) + 1.67917e-3/( 57.362e0 - sigma2)
            
            wave_vac[g] = wave_air[g]*fact              # Convert Wavelength

    return wave_vac


def vactoair(wave_vac):
    """
    Convert vacuum wavelengths to air wavelengths

    Corrects for the index of refraction of air under standard conditions.  
    Wavelength values below 2000 A will not be altered.  Accurate to 
    about 10 m/s.


    INPUT/OUTPUT:
	WAVE_VAC - Vacuum Wavelength in Angstroms, scalar or vector
		If the second parameter is not supplied, then this will be
               updated on output to contain double precision air wavelengths.

    EXAMPLE:
	If the vacuum wavelength is  W = 2000, then 

	IDL> VACTOAIR, W 

	yields an air wavelength of W = 1999.353 Angstroms

    METHOD:
	Formula from Ciddor 1996  Applied Optics , 35, 1566

    REVISION HISTORY
      Written, D. Lindler 1982 
      Documentation W. Landsman  Feb. 1989
      Use Ciddor (1996) formula for better accuracy in the infrared 
           Added optional output vector, W Landsman Mar 2011
    """

  
    wave_vac = np.atleast_1d(wave_vac).copy()  # makes sure it's an array
    wave_air = np.atleast_1d(wave_vac).copy()  # initialize
    g,ng = utils.where(wave_air >= 2000)     # Only modify above 2000 A
    
    if ng>0:
        sigma2 = (1e4/wave_vac[g] )**2   # Convert to wavenumber squared

        # Compute conversion factor
        fact = 1.0 +  5.792105e-2/(238.0185e0 - sigma2) + 1.67917e-3/( 57.362e0 - sigma2)
    
        # Convert wavelengths
        wave_air[g] = wave_vac[g]/fact

    return wave_air


def xyz2lbd(x,y,z,R0=8.5):
    """ Convert galactocentric X/Y/Z coordinates to l,b,dist."""
    ll = np.atleast_1d(x).copy()*0.0
    bb = np.atleast_1d(x).copy()*0.0
    dd = np.atleast_1d(x).copy()*0.0    
    for i in range(len(ll)):
        xx = np.float64(x[i])
        yy = np.float64(y[i])
        zz = np.float64(z[i])
        rho = np.sqrt( (xx+R0)**2 + yy**2)  # distance from sun in X/Y plane

        lrad = np.arctan2(yy,xx+R0)
        
        brad = 0.5*np.pi - np.arctan2(rho,zz)      # this is more straighforward
        if brad > 0.5*np.pi:
            brad = brad-np.pi
        if brad < -0.5*np.pi:
            brad = brad+np.pi

        # This doesn't work if z=0
        #if cos(0.5*!dpi-brad) ne 0.0 then d = zz/cos(0.5*!dpi-brad)
        #if cos(0.5*!dpi-brad) eq 0.0 then d = abs(zz)
        d = np.sqrt( (xx+R0)**2 + yy**2 + zz**2 )

        b = np.rad2deg(brad)
        l = np.rad2deg(lrad)
        
        if l < 0.:
            l = l+360.
        
        ll[i] = l
        bb[i] = b
        dd[i] = d

    return ll,bb,dd

def vgsr2vhelio(vgsr,lon,lat,vcirc=240.0):
    """
    Given an input array called vhel with values of heliocentric
    radial velocities and input arrays of the same length with
    the Galactic coordinates of the object (gl,gb),
    this code calculates Vlsr and Vgsr.

    INPUTS:
    vgsr     Galactocentric Standard of Rest Radial Velocity
    lon      Galactic longitude (in degrees)
    lat      Galactic latitude (in degrees)
    =vcirc   MW circular velocity.  220 km/s by default.
   
    OUTPUTS:
    Vhelio  The heliocentric velocity in km/s.

    USAGE:
    IDL>vhelio = vgsr2vhelio(vgsr,lon,lat)

    By D.Nidever 2005
    """

    # code assumes gl & gb given in degrees.

    gl = np.deg2rad(lon)
    gb = np.deg2rad(lat)

    cgl=np.cos(gl)
    sgl=np.sin(gl)
    cgb=np.cos(gb)
    sgb=np.sin(gb)

    #  This equation takes the solar motion w.r.t. the LSR as
    #  (9,11,6) km/sec (Mihalas & Binney)
    #vlsr=vhel+((9.0d0*cgb*cgl)+(11.0d0*cgb*sgl)+(6.0d0*sgb))

    #  This equation takes the rotation velocity of the LSR to
    #  be 220 km/sec

    vlsr = vgsr-(vcirc*cgb*sgl)

    #   Using updated values from Dehnen & Binney 1998, MNRAS, 298, 387
    #   U = 10.00 +/- 0.36 km/s
    #   V = 5.25 +/- 0.62 km/s
    #   W = 7.17 +/- 0.38 km/s
    #   This is in a right-handed system
    #vlsr=vhel+((9.0d0*cgb*cgl)+(11.0d0*cgb*sgl)+(6.0d0*sgb))
    #vhel=vlsr-( (9.0d0*cgb*cgl)+(11.0d0*cgb*sgl)+(6.0d0*sgb) )
    vhel = vlsr-( (10.0*cgb*cgl)+(5.25*cgb*sgl)+(7.17*sgb) )
    
    return vhel


def vgsr2vlsr(v,l,b,dir,vcirc=240.0):
    """

    PURPOSE:
    Given an input array called vhel with values of heliocentric
    radial velocities and input arrays of the same length with
    the Galactic coordinates of the object (gl,gb),
  this code calculates Vlsr and Vgsr.

    INPUTS:
    vel     Vgsr/Vlsr velocity to convert
    l       Galacitc longitude in degrees
    b       Galactic latitutde in degrees
    dir     Direction of conversion
             dir=1    Vgsr->Vlsr
             dir=-1   Vlsr->Vgsr
    =vcirc  The rotation speed of the MW disk at the sun's
                 radius.  Vcirc=240 km/s by default.

    OUTPUTS:
    Vlsr/Vgsr velocity array

    USAGE:
    IDL>vlsr = vgsr2vlsr(vgsr,lon,lat)

    By D.Nidever Feb 2006
    """

    gl = l
    gb = b
    
    gl = np.deg2rad(gl)
    gb = np.deg2rad(gb)

    cgl = np.cos(gl)
    sgl = np.sin(gl)
    cgb = np.cos(gb)
    sgb = np.sin(gb)

    #  This equation takes the solar motion w.r.t. the LSR as
    #  (9,11,6) km/sec (Mihalas & Binney)
    #
    #vlsr=vhel+((9.0d0*cgb*cgl)+(11.0d0*cgb*sgl)+(6.0d0*sgb))


    #  This equation takes the rotation velocity of the LSR to
    #  be 220 km/sec
    #
    #vgsr=vlsr+(220.0d0*cgb*sgl)

    #vgsr=vlsr+(220.0d0*cgb*sgl)
    #vlsr=vgsr-(220.0d0*cgb*sgl)

    # input vgsr, output vlsr
    if dir == 1:
        v2 = v-(vcirc*cgb*sgl)

    # input vlsr, output vgsr
    if dir == -1:
        v2 = v+(vcirc*cgb*sgl)

    return v2


def galaxy_model(flat=False,vcirc=240.0,R0=8.5,rscale=12.5,zscale=5.0,rmin=1.5,vdisk=240.0,vdisp=0.0,rhelcut=None,nstars=100000):
    """
    Model of the proper motions and rvs of the galaxy.
    Mostly copied from ggss/ggss_model.pro

    INPUT:
    =vdisk    The circular velocity of the disk (default vdisk=vcirc)
    =vcirc    The circular velocity of the sun (default vcirc=220 km/s)
    =R0       The distance from the sun to the GC (default R0=8.5 kpc)
    =rscale   The upper radial cutoff of the disk (default rscale=12.5 kpc)
    =zscale   The z-height cutoff of the disk (default zscale=5.0 kpc)
    =vdisp    Velocity dispersion (default is vdisp=0 km/s)
    =rmin     The lower radial cutoff of the disk(default rmin=1.5 kpc)
    =rhelcut  Only keep stars within this distance of the sun.
    =nstars   The number of stars in the simulation (default nstars=100,000)
    /flat     Simulate a "flat" galaxy, i.e. no thickness, all stars at b=0.

    OUTPUT:
    x       Galactocentric X direction (in kpc)
    y       Galactocentric Y direction (in kpc)
    z       Galactocentric Z direction (in kpc)
    vx      Velocity in the galactocentric X direction (in km/s)
    vy      Velocity in the galactocentric Y direction (in km/s)
    vz      Velocity in the galactocentric Z direction (in km/s)
    l       Galactic longitude
    b       Galactic latitude
    vhel    Heliocentic radial velocity (in km/s)
    vlsr    Radial velocity wrt the local standard of rest (LSR) (in km/s)
    vgsr    Radial velocity wrt the galactic standard of rest (GSR) (in km/s)
    mul     True proper motion in galactic longitude (in mas/yr)
    mub     True proper motion in galactic latitude (in mas/yr)

    Written by D.Nidever  March 2006
    """

    #The Coordinate System (right-handed):
    # x = toward galactic center from sun
    # y = toward l=90
    #  z = toward NGP
    #  phi = 0 at x, positive toward +y (against galactic rotation)
    #  r = radial distance from GC


    # Getting the random positions
    np.random.seed(0)
    x = (np.random.rand(2*nstars)-0.5)*2.0*rscale
    y = (np.random.rand(2*nstars)-0.5)*2.0*rscale    
    # Flat galaxy
    if flat is True:
        z = np.zeros(nstars,float)
    # Galaxy extended in the Z direction
    else:
        z = np.random.randn(2*nstars)*zscale

    # Only want stars within a radius RSCALE and greater than RMIN*RSCALE 
    gd, = np.where( (np.sqrt(x**2.+y**2.) < rscale) & (np.sqrt(x**2.+y**2.) > rmin))
    x = x[gd]
    y = y[gd]
    z = z[gd]

    # RHEL CUT
    if rhelcut is not None:
        xhel = x+R0
        rhel = np.sqrt(xhel**2.+y**2.+z**2.)
        gd, = np.where(rhel <= rhelcut)
        if len(gd)==0:
            raise ValueError('No stars within '+str(rhelcut)+' kpc of the sun')
        x = x[gd]
        y = y[gd]
        z = z[gd]

    # Trim to number of stars that we want
    if len(x)>nstars:
        x = x[0:nstars]
        y = y[0:nstars]
        z = z[0:nstars]
        
    # Cylindrical and spherical coordinates
    rho = np.sqrt(x**2.+y**2.)    # distance from z axis
    gphi = np.arctan2(y,x)         # the galactic azimuthal coordinate (0-2*pi), in radians
    gtheta = np.arctan2(rho,z)     # the galactic polar coordinate (0-pi), in radians
    r = np.sqrt(x**2.+y**2.+z**2.) # radial distance, in kpc

    # Get Galactic coordinates
    l,b,d = xyz2lbd(x,y,z,R0=R0)
    
    # Use a radius-dependent rotation curve
    vcoef = np.array([0.0, 55.0, -4.0, 0.09])
    vdisk = np.poly1d(np.flip(vcoef))(rho)

    # Cartesian velocities
    #  Moving in the -gphi direction
    #  phihat = -sin(phi)*xhat + cos(phi)*yhat
    vx = (-vdisk)*(-np.sin(gphi))
    vy = (-vdisk)*np.cos(gphi)
    vz = x.copy()*0.
    # add velocity dispersion
    if vdisp > 0.0:
        np.random.seed(1)
        vx += np.random.randn(len(vx))*vdisp/np.sqrt(3.)
        vy += np.random.randn(len(vx))*vdisp/np.sqrt(3.)
        vz += np.random.randn(len(vx))*vdisp/np.sqrt(3.)
    vtot = np.sqrt(vx**2+vy**2+vz**2)     # total velocity, should be vcirc

    # *** NOW WE "OBSERVE" THE STARS ***

    # RADIAL VELOCITY
    
    #Calculating the heliocentric radial velocity
    # Vrad = Vrel dotted with unit vector from sun to the star
    # Vrad = V*cos(angle)

    ## Relative velocity of the star wrt the sun
    # THIS ASSUMES THAT THE SUN IS IN COMPLETE CIRCULAR MOTION
    vrelvec = np.zeros((len(x),3),np.float64)
    vrelvec[:,0] = vx
    vrelvec[:,1] = vy-vcirc    # the sun is moving in +y direction
    vrelvec[:,2] = vz

    # unit vector for our line-of-sight toward the star, rhat
    #  in helio cartesian coordinate, xhel=x+R0
    # rhat = (x/r)*xhat + (y/r)*yhat + (z/r)*zhat
    xhel = x+R0
    rhel = np.sqrt(xhel**2+y**2+z**2)
    rhat = np.zeros((len(x),3),np.float64)
    rhat[:,0] = xhel/rhel 
    rhat[:,1] = y/rhel
    rhat[:,2] = z/rhel

    # ****THIS IS VLSR NOT VHELIO*****

    # Take the dot product of the two vectors
    # vvec dot rhat = vrelvec[0]*rhat[0] + vrelvec[1]*rhat[1] + vrelvec[2]*rhat[2]
    vlsr = vrelvec*rhat
    vlsr = np.sum(vlsr,axis=1)

    # Getting Vgsr and Vhelio velocity
    vgsr = vgsr2vlsr(vlsr,l,b,-1)
    vhel = vgsr2vhelio(vgsr,l,b,vcirc=vcirc)


    # PROPER MOTIONS

    # Calculating the heliocentric transverse velocity and proper motion

    # heliocentric spherical coordinates
    htheta = 90-b   # 0 at NP, 180 at SP
    hphi = l

    # In the +LON direction

    # unit vector in the +LON direction at the star's position, lhat
    # lhat = -sin(hphi)*xhat + cos(hphi)*yhat
    lhat = np.zeros((len(x),3),np.float64)
    lhat[:,0] = -np.sin(np.deg2rad(hphi))
    lhat[:,1] = np.cos(np.deg2rad(hphi))
    lhat[:,2] = 0.0

    # transverse velocity in +LON direction
    #  relative velocity dot lhat
    vtanl = vrelvec*lhat
    vtanl = np.sum(vtanl,axis=1)

    # proper motion in +LON direction
    #   mu(arcsec/year) = Vtan(km/s)/(4.74*d(pc))
    mul = vtanl/(4.74*rhel*1000.0)          # in arcsec/year
    mul = mul * 1000.0                      # in mas/year

    # in +LAT direction

    # unit vector in the +LAT direction at the star's position, bhat
    # bhat = cos(htheta)*cos(hphi)*xhat + cos(htheta)*sin(hphi)*yhat - sin(htheta)*zhat
    bhat = np.zeros((len(x),3),np.float64)
    bhat[:,0] = np.cos(np.deg2rad(htheta))*np.cos(np.deg2rad(hphi))
    bhat[:,1] = np.cos(np.deg2rad(htheta))*np.sin(np.deg2rad(hphi))
    bhat[:,2] = -np.sin(np.deg2rad(htheta))

    # transverse velocity in +LAT direction
    #  relative velocity dot lhat
    vtanb = vrelvec*bhat
    vtanb = np.sum(vtanb,axis=1)

    # proper motion in +LAT direction
    #   mu(arcsec/year) = Vtan(km/s)/(4.74*d(pc))
    mub = vtanb/(4.74*rhel*1000.0)          # in arcsec/year
    mub = mub * 1000.0                      # in mas/year

    # The position angle of the proper motion, E of N
    pa = np.rad2deg(np.arctan2(mul,mub))
    bd, = np.where(pa < 0.0)
    if len(bd)>0:
        pa[bd] = pa[bd]+360.0

    # Make the output structure
    dt = np.dtype([('x',float),('y',float),('z',float),('l',float),('b',float),('d',float),('rho',float),('phi',float),
                   ('vx',float),('vy',float),('vz',float),('vtot',float),('vhelio',float),('vlsr',float),('vgsr',float),('pml',float),('pmb',float)])
    out = np.zeros(len(l),dtype=dt)
    out['x'] = x
    out['y'] = y
    out['z'] = z
    out['l'] = l
    out['b'] = b
    out['d'] = d
    out['rho'] = rho
    phi = 180-np.rad2deg(gphi)
    bd, = np.where(phi > 180)
    if len(bd)>0:
        phi[bd] -= 360       # left-handed system
    out['phi'] = phi
    out['vx'] = vx
    out['vy'] = vy
    out['vz'] = vz
    out['vtot'] = np.sqrt(vx**2+vy**2+vz**2)
    out['vhelio'] = vhel
    out['vlsr'] = vlsr
    out['vgsr'] = vgsr
    out['pml'] = mul
    out['pmb'] = mub

    return out
