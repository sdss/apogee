# routine to get dither shifts from successfully combined pairs

import glob
from astropy.io import fits
import pdb
import matplotlib.pyplot as plt
import numpy as np

files=glob.glob('*/*/5[789]??[0-9]/apPlate-a*.fits')
mjd=[]
dith=[]
print(len(files))
pdb.set_trace()
for ifile,file in enumerate(files) :
    print(ifile,file)
    try:
        a=fits.open(file)
        for shift in a[14].data['RELSHIFT'] :
            mjd.append(a[0].header['JD-MID'])
            dith.append(shift)
    except : pass
plt.scatter(mjd,dith,marker='o',color='r',s=10)
np.savetxt('dither.dat',np.vstack([mjd,dith]).T)

pdb.set_trace()
