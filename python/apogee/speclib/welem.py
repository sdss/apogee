#
# Routine to make list of windows for a given element from a FERRE filter file
#
from astropy.io import ascii
import pdb

el='Ce'
halfwidth=1.2

wave=ascii.read('wave.dat',Reader=ascii.NoHeader)['col1']
pix=ascii.read('wave.dat',Reader=ascii.NoHeader)['col2']
filt=ascii.read(el+'.filt',Reader=ascii.NoHeader)['col1']

fp=open(el+'.wave','w')
start=-1
cens=[]
for i in range(len(filt) ):
  if filt[i] > 0 and start < 0:
    start=wave[i]-halfwidth
  elif filt[i] == 0. and start >=0 :
    wind=[start,wave[i]+halfwidth]
    cens.append(wind)
    start=-1

for c in cens :
    fp.write('{:8.3f} {:8.3f}\n'.format(c[0],c[1]))
    print(c,c[1]-c[0])

fp.close()
pdb.set_trace()
