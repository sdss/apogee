import numpy as np
from tools import plots
import pdb

def chip(a,row=150,ax=None,pixel=False,visit=False,color=None) :
    """ Routine to plot 3 chips in 3 panels
    """
    if ax is None : fig,ax=plots.multi(1,3,hspace=0.3)
    if visit : rows=range(3)
    else : rows=[row,row,row]
    chips=['a','b','c']
    x=np.arange(a['a'][1].header['NAXIS1'])
    for ichip,chip in enumerate(chips) :
        if pixel :
            plots.plotl(ax[ichip],x,a[chip][1].data[rows[ichip],:],color=color)
        else :
            plots.plotl(ax[ichip],a[chip][4].data[rows[ichip],:],a[chip][1].data[rows[ichip],:],color=color)

    try : return fig,ax
    except : return
