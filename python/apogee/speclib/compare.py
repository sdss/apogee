from astropy.io import fits
import os
import numpy as np
import matplotlib.pyplot as plt
from apogee.aspcap import aspcap
from apogee.speclib import synth
from apogee.utils import spectra

te=5750
logg=4.5
mh=0.
#old=fits.open('kurucz/solarisotopes/tsGK_150714/ap00cp00np00vp10.fits')[0].data
#new=fits.open('marcs/solarisotopes/tdGK_180901/ap00cp00np00vp12.fits')[1].data
# synthesis
old=synth.mkturbospec(te,logg,mh,0.,0.,0.,kurucz=True,linelist='20150714')[1]
oldmarcs=synth.mkturbospec(te,logg,mh,0.,0.,0.,kurucz=False,linelist='20150714')[1]
new=synth.mkturbospec(te,logg,mh,0.,0.,0.,kurucz=False,linelist='20180901')[1]
w=np.linspace(15100,17000,38001)
plt.clf()
plt.plot(spectra.airtovac(w),new/old)
plt.plot(spectra.airtovac(w),oldmarcs/old)

# LSF convolved (with lsfc)
irot=6
olds=fits.open(os.environ['APOGEE_SPECLIB']+'/synth/turbospec/kurucz/solarisotopes/tsGK_150714_lsfc/ap00cp00np00vp10.fits')
it=int((te-olds[1].header['CRVAL2'])/olds[1].header['CDELT2'])
ig=int((logg-olds[1].header['CRVAL3'])/olds[1].header['CDELT3'])
im=int((mh-olds[1].header['CRVAL4'])/olds[1].header['CDELT4'])
spec=np.append(olds[1].data[irot,0,im,ig,it,:],olds[2].data[irot,0,im,ig,it,:])
spec=np.append(spec,olds[3].data[irot,0,im,ig,it,:])
spec=aspcap.aspcap2apStar(spec)

news=fits.open(os.environ['APOGEE_SPECLIB']+'/synth/turbospec/marcs/solarisotopes/tdGK_180901_lsfc_l33/ap00cp00np00vp12.fits')
it=int((te-news[0].header['CRVAL2'])/news[0].header['CDELT2'])
ig=int((logg-news[0].header['CRVAL3'])/news[0].header['CDELT3'])
im=int((mh-news[0].header['CRVAL4'])/news[0].header['CDELT4'])
newspec=news[0].data[irot,im,ig,it,:]

plt.plot(aspcap.apStarWave(),newspec/spec)

# output to fits
hdulist=fits.HDUList()
hdu=fits.ImageHDU(old)
spectra.add_dim(hdu.header,15100.,0.05,1,'Wavelength',1)
hdulist.append(hdu)
hdu=fits.ImageHDU(oldmarcs)
spectra.add_dim(hdu.header,15100.,0.05,1,'Wavelength',1)
hdulist.append(hdu)
hdu=fits.ImageHDU(new)
spectra.add_dim(hdu.header,15100.,0.05,1,'Wavelength',1)
hdulist.append(hdu)
hdu=fits.ImageHDU(spec)
spectra.add_dim(hdu.header,4.179,6.e-6,1,'log Wavelength',1)
hdulist.append(hdu)
hdu=fits.ImageHDU(newspec)
spectra.add_dim(hdu.header,4.179,6.e-6,1,'log Wavelength',1)
hdulist.append(hdu)
hdulist.writeto('solarcomp.fits')
