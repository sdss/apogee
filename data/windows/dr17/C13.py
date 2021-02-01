import numpy as np
from apogee.utils import spectra

a=np.loadtxt('C13.dat')
vac=spectra.airtovac(a)
fp=open('C13.wind','w')
for wind in vac :
  fp.write('{:12.3f}{:12.3f}{:3}\n'.format(wind[0],wind[1],1))
fp.close()

wave=np.loadtxt('wave.dat')
for w in wave[:,0] :
  weight=0.
  for wind in vac :
    if w >wind[0] and w<wind[1] : weight=1.
  print('{:12.3}'.format(weight))


