#
# Routine to make list of windows for a given element from a FERRE filter file
#
from astropy.io import ascii
from apogee.aspcap import aspcap
import numpy as np
import pdb

def mkwind(el,halfwidth=1.2) :
    """ Create default window file from master filt file
    """
    chipswave=aspcap.gridWave()
    wave=[]
    for ichip in range(3) : wave.extend(chipswave[ichip])

    filt=ascii.read(el+'.filt',Reader=ascii.NoHeader)['col1']

    pdb.set_trace()
    fp=open(el+'.wind','w')
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
        fp.write('{:8.3f} {:8.3f}   1\n'.format(c[0],c[1]))
        print(c,c[1]-c[0])

    fp.close()

def mkmask(el) :
    """ Given input master filt file, and window file with desired windows, output mask file
    """
    chipswave=aspcap.gridWave()
    wave=[]
    for ichip in range(3) : wave.extend(chipswave[ichip])

    filt=ascii.read(el+'.filt',Reader=ascii.NoHeader)['col1']
    wind=ascii.read(el+'.wind',Reader=ascii.NoHeader)
    fp=open(el+'.mask','w')
    bd=np.where(wind['col3'] == 0)[0]
    for w in bd :
        j=np.where( (wave>wind['col1'][w]) & (wave<wind['col2'][w]) )[0]
        filt[j] = 0.
    np.savetxt(fp,filt,fmt='%12.8f')
    fp.close()

