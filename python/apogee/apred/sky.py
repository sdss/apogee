from tools import plots
from apogee.utils import bitmask
from astropy.io import fits
import numpy as np

def check(frame,fiber) :

    pixmask=bitmask.PixelBitMask()
    fig,ax=plots.multi(1,1,figsize=(12,4))
    for chip in ['a','b','c'] :
        a=fits.open('apCframe-{:s}-{:08d}.fits'.format(chip,frame))
        ax.plot(a[4].data[300-fiber,:],a[5].data[300-fiber,:])
        ax.plot(a[4].data[300-fiber,:],a[1].data[300-fiber,:])
        mask=np.where((a[3].data[300-fiber,:]&pixmask.getval('SIG_SKYLINE')) > 0)[0]
        ax.plot(a[4].data[300-fiber,mask],a[1].data[300-fiber,mask],'ro')

