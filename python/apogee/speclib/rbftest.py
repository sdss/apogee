import numpy as np
import matplotlib.pyplot as plt
import pdb
from astropy.io import fits
import os
from tools import plots

# read in holefiles
hole=fits.open(os.environ['APOGEE_SPECLIB']+'/atmos/marcs/MARCS_v3_2016/MARCS_GK_holefile.fits')[0].data
rbfhole=fits.open('rbf_MARCS_GK_holefile.fits')[0].data

neg=np.where(rbfhole < 1.e-5)
zer=np.where(np.isclose(rbfhole,0.))

# choose a [C/M] and [alpha/M]
j=np.where((neg[0] == 0) & (neg[1] == 4))[0]
a=fits.open('ap00cm10np00vp12.fits')[0].data
rbf=fits.open('rbf_ap00cm10np00vp12.fits')[0].data
icm=0
iam=4

# loop through spectra in this file
for imh in range(10,a.shape[0]) :
  for ilogg in range(a.shape[1]) :
    for iteff in range(a.shape[2]) :
        print(hole[icm,iam,imh,ilogg+1,iteff],rbfhole[icm,iam,imh,ilogg+1,iteff])
        plt.clf()
        plt.plot(rbf[imh,ilogg,iteff,:]/a[imh,ilogg,iteff,:]/(rbf[imh,ilogg,iteff,:]/a[imh,ilogg,iteff,:]).mean())
        plt.ylim(0.9,1.1)
        plt.draw()
        pdb.set_trace()
    

