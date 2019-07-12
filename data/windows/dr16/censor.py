from astropy.io import ascii
import pdb

pix=ascii.read('wave.dat',Reader=ascii.NoHeader)['col2']
filt=ascii.read('Fe.filt',Reader=ascii.NoHeader)['col1']

start=0
cens=[]
for i in range(len(filt) ):
  if filt[i] > 0 and start >= 0:
    skip=[start,pix[i]]
    cens.append(skip)
    start=-1
  elif filt[i] == 0. and start <0 :
    start = pix[i]

cens.append([start,-1])

pdb.set_trace()
