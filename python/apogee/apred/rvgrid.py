import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.io import ascii
from astropy.io import fits
from apogee.aspcap import aspcap
from apogee.speclib import sample
from apogee.speclib import synth
from apogee.speclib import lsf
from apogee.utils import spectra
from scipy import interpolate, ndimage

import pdb

def mkgrid(apred='r10', telescope='apo25m', lsfid=5440020, waveid=2420038,clobber=True) :

    # create parameters of sample and read in
    if clobber : sample.sample(gridclass='rv')
    stars=ascii.read('test_rv')

    # synthesize spectra and read in. Will fail on hot stars, but add those below
    if clobber : synth.mksynth('test_rv')
    grid=fits.open('test_rv.fits')[1].data
    pars=fits.open('test_rv.fits')[0].data

    # load output parameters of successful syntheses
    n=len(stars)
    outpar = np.zeros(n,dtype=[
                           ('teff','f4'),
                           ('logg','f4'),
                           ('mh','f4'),
                           ('am','f4'),
                           ('cm','f4')])
    for i in range(pars.shape[0]) :
        outpar['teff'][i] = pars[i,0]
        outpar['logg'][i] = pars[i,1]
        outpar['mh'][i] = pars[i,2]
        outpar['am'][i] = pars[i,3]
        outpar['cm'][i] = pars[i,4]
    nout=pars.shape[0]

    # Add in hot stars
    a=fits.open(os.environ['APOGEE_SPECLIB']+'/synth/synspec/kurucz/solarisotopes/sBA/ap00cp00np00vp20.fits')
    wav=spectra.fits2vector(a[0].header,1)
    hot=np.where(stars['Teff']>8000)[0]
    init=0
    for i in hot :
        te=stars['Teff'][i]
        logg=stars['logg'][i]
        mh=stars['[M/H]'][i]
        ite=int(round((te-a[0].header['CRVAL2'])/a[0].header['CDELT2']))
        ilogg=int(round((logg-a[0].header['CRVAL3'])/a[0].header['CDELT3']))
        imh=int(round((mh-a[0].header['CRVAL4'])/a[0].header['CDELT4']))
        print(te,logg,mh,ite,ilogg,imh)
        if ite>=0 & ilogg >=0 & imh >=0 :
            spec=a[0].data[imh,ilogg,ite,:].reshape(1,38001)
            # need to try to remove continuum from hot grid
            baseline= np.polynomial.Polynomial.fit(np.log(wav),spec[0,:],1)
            gd=np.where(spec[0,:]/baseline(np.log(wav)) > 0.99)[0]
            baseline= np.polynomial.Polynomial.fit(np.log(wav[gd]),spec[0,gd],1)
 
            grid=np.append(grid,spec/baseline(np.log(wav)),axis=0)
            outpar['teff'][nout] = te
            outpar['logg'][nout] = logg
            outpar['mh'][nout] = mh
            outpar['am'][nout] = 0.
            outpar['cm'][nout] = 0.
            nout += 1
    # Trim to successful syntheses
    outpar=outpar[0:nout]

    # output FITS file
    hdu=fits.HDUList()
    h=fits.PrimaryHDU()
    h.header['COMMENT'] = 'Synthetic RV grids'
    h.header['COMMENT'] = 'HDU #1: parameter table'
    #h.header['COMMENT'] = 'HDU #2: high resolution, unsmoothed'
    h.header['COMMENT'] = 'HDU #2: LSF smoothed, apStar resolution (log, 6.e-6 dispersion)'
    h.header['COMMENT'] = 'HDU #3: LSF smoothed, apVisit type resolution (log, 4.75e-6 dispersion)'
    hdu.append(h)

    # parameter table
    hdu.append(fits.TableHDU(outpar))

    # raw spectra (currently commented out)
    h=fits.ImageHDU(grid)
    h.header['CRVAL1'] = 15100.
    h.header['CDELT1'] = 0.05
    h.header['CTYPE1'] = 'WAVELENGTH'
    h.header['LSFID'] = lsfid
    h.header['WAVEID'] = waveid
    h.header['APRED'] = apred
    #hdu.append(h)

    # put on logarithmic wavelength scale and smooth to R=30000
    nspec=grid.shape[0]
    resol=30000.
    ws=spectra.airtovac(np.linspace(15100.,17000.,38001))
    wout=10.**(np.arange(aspcap.logw0,np.log10(17000.),aspcap.dlogw/9.))
    smooth=np.zeros([nspec,len(wout)])
    sig=1./(aspcap.dlogw/9.*np.log(10))/resol/2.354
    for i in range(nspec) :
        # Interpolate the input spectrum
        ip= interpolate.InterpolatedUnivariateSpline(ws, grid[i,:], k=3)
        smooth[i,:] = ndimage.gaussian_filter1d(ip(wout),sig,mode='nearest')

    # now ressample to apstar and (approx) apvisit wavelenth scales
    wapstar=10.**(aspcap.logw0+np.arange(aspcap.nw_apStar)*aspcap.dlogw)
    apstardata=np.zeros([nspec,len(wapstar)])
    wapvisit=10.**(np.arange(aspcap.logw0,np.log10(17000.),4.75e-6))
    apvisitdata=np.zeros([nspec,len(wapvisit)])
    for i in range(nspec) :
        # Interpolate the input spectrum
        ip= interpolate.InterpolatedUnivariateSpline(wout, smooth[i,:], k=3)
        apstardata[i,:] = ip(wapstar)
        apvisitdata[i,:] = ip(wapvisit)

    # convolve and sample to apStar grid in vacuum
    #ws=spectra.airtovac(np.linspace(15100.,17000.,38001))
    #x,ls= synth.getlsf(lsfid,waveid,apred=apred,prefix='lsf_',fill=True)
    #apstardata=lsf.convolve(ws,grid,lsf=ls,xlsf=x)
    h=fits.ImageHDU(apstardata)
    h.header['CRVAL1'] = aspcap.logw0
    h.header['CDELT1'] = aspcap.dlogw
    h.header['CTYPE1'] = 'LOG(WAVELENGTH)'
    #h.header['LSFID'] = lsfid
    #h.header['WAVEID'] = waveid
    #h.header['APRED'] = apred
    h.header['RESOL'] = resol
    hdu.append(h)

    # convolve and sample to higher resolution for apVisit spectra
    #wout=10.**(np.arange(4.179,np.log10(17000.),4.75e-6))
    #smoothdata=lsf.convolve(ws,grid,lsf=ls,xlsf=x,highout=True)
    #nwav=smoothdata.shape[1]
    #nspec=smoothdata.shape[0]
    #apvisitdata=np.zeros([nspec,len(wout)])
    #wav=10.**(aspcap.logw0+np.arange(nwav)*aspcap.dlogw/9.)
    #for i in range(nspec) :
    #    # Interpolate the input spectrum
    #    bd=np.where(np.isnan(smoothdata[i,:]))[0]
    #    smoothdata[i,bd]=1.
    #    ip= interpolate.InterpolatedUnivariateSpline(wav, smoothdata[i,:], k=3)
    #    apvisitdata[i,:] = ip(wout)
    h=fits.ImageHDU(apvisitdata)
    h.header['CRVAL1'] = 4.179
    h.header['CDELT1'] = 4.75e-6
    h.header['CTYPE1'] = 'LOG(WAVELENGTH)'
    h.header['RESOL'] = resol
    hdu.append(h)

    # output
    hdu.writeto('apg_synthgrid.fits',overwrite=True)

