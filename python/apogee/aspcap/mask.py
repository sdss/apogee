#
# Routine to make list of windows for a given element from a FERRE filter file
#
from astropy.io import ascii
from apogee.aspcap import aspcap
import numpy as np
import pdb

def mkwind(el,halfwidth=1.2,invert=False) :
    """ Create default window file with unmasked regions from master filt file
        use with halfwidth to get windows for minigrid synthesis
        use with invert to get masked regions
    """
    chipswave=aspcap.gridWave()
    wave=[]
    for ichip in range(3) : wave.extend(chipswave[ichip])

    filt=ascii.read(el+'.filt',Reader=ascii.NoHeader)['col1']

    fp=open(el+'.wind','w')
    start=-1
    cens=[]
    for i in range(len(filt) ):
        if invert :
            if filt[i] == 0 and start < 0:
              start=wave[i]-halfwidth
            elif filt[i] > 0. and start >=0 :
              wind=[start,wave[i]+halfwidth]
              cens.append(wind)
              start=-1
        else :
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

def mkmask(el,globalmask=None) :
    """ Given input master filt file, and window file with desired windows, output mask file
    """
    chipswave=aspcap.gridWave()
    wave=[]
    for ichip in range(3) : wave.extend(chipswave[ichip])

    filt=ascii.read(el+'.filt',format='no_header')['col1']
    wind=ascii.read(el+'.wind',format='no_header')
    fp=open(el+'.mask','w')
    bd=np.where(wind['col3'] == 0)[0]
    for w in bd :
        j=np.where( (wave>wind['col1'][w]) & (wave<wind['col2'][w]) )[0]
        filt[j] = 0.
    # add in bad pixelsl from global mask if we have one 
    if globalmask is not None :
        g=ascii.read(globalmask,format='no_header')['col1']
        bd=np.where(g == 0)[0]
        filt[bd] = 0.
    np.savetxt(fp,filt,fmt='%12.8f')
    fp.close()

#    elems=['C','CI','N','O','Na','Mg','Al','Si','P','S','K','Ca','Ti','TiII','V','Cr','Mn','Fe','Co','Ni','Cu','Ge','Rb','Ce','Nd','Yb']

def mkparam(out='param.mask',edgefile='global.mask', globalfile='mask_v02_aspcap.txt', els=['CI','O','Na','Al','Si','P','S','K','Ca','Ti','TiII','V','Cr','Mn','Co','Ni','Cu','Ge','Rb','Ce','Nd','Yb']) :
    """ Make param mask from edge, global and element masks
    """
    mask=ascii.read(edgefile,format='no_header')['col1']
    gmask=ascii.read(globalfile,format='no_header')['col1']
    bd=np.where(gmask == 0)[0]
    mask[bd] = 0.
    for el in els :
        elem=ascii.read(el+'.filt',format='no_header')['col1']
        bd=np.where(elem > 0.)[0]
        mask[bd] = 0.
    fp=open(out,'w')
    np.savetxt(fp,mask,fmt='%8.2f')
    fp.close()
   


