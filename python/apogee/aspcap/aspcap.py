# encoding: utf-8
#
# @Author: Jon Holtzman
# @Date: March 2018
# @Filename: synth.py
# @License: BSD 3-Clause
# @Copyright: Jon Holtzman

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import numpy as np
import os
import pdb
import scipy.ndimage.filters
from astropy.io import fits
from astropy.io import ascii
#from holtz.tools import struct
#from holtz.tools import plots
from apogee.utils import apload

def params() :
    '''
    Define the order of the parameter arrays, with the associated FERRE names, tag names, and flag names
    '''
    tagnames=np.array(['TEFF','LOGG','LOGVMICRO','M_H','C_M','N_M','ALPHA_M','LGVSINI','PARAM_O'])
    flagnames=np.array(['TEFF','LOGG','VMICRO','M_H','C_M','N_M','ALPHA_M','VSINI','O'])
    params=np.array(['TEFF','LOGG','LOG10VDOP','METALS','C','N','O Mg Si S Ca Ti','LGVSINI','O'])
    return params,tagnames,flagnames

def elems(nelem=0) :
    '''
    Define the order of the element arrays
    '''

    elems=['C','CI','N','O','Na','Mg','Al','Si','P','S','K','Ca','Ti','TiII','V','Cr','Mn','Fe','Co','Ni','Cu','Ge','Ce','Rb','Y','Nd']
    #return,['C','N','O','Na','Mg','Al','Si','S','K','Ca','Ti','V','Mn','Fe','Ni','Nd']
    #return,['C','N','O','Na','Mg','Al','Si','S','K','Ca','Ti','V','Mn','Fe','Ni']
    elemtoh=[0,0,0,0,1,0,1,0,1,0,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] #0,0,0,0,0,0,0,0,0]
    tagnames=[]
    elemfitnames=[]
    for i in range(len(elems) ) :
        if elemtoh[i] :
            tagnames.append(elems[i]+'_Fe')
            elemfitnames.append('['+elems[i]+'/Fe]')
        else :
            tagnames.append(elems[i]+'_M')
            elemfitnames.append('['+elems[i]+'/M]')
    return elems,elemtoh,tagnames,elemfitnames

logw0=4.179
dlogw=6.e-6
nw_apStar=8575
def apStarWave() :
    """ Returns apStar wavelengths
    """
    return 10.**(logw0+np.arange(nw_apStar)*dlogw)

logw0_chip=np.array([4.180476,4.200510,4.217064])
nw_chip=np.array([3028,2495,1991])
def gridWave() :
    """ Returns aspcap grid wavelengths
    """
    return [10.**(logw0_chip[0]+np.arange(nw_chip[0])*dlogw),
            10.**(logw0_chip[1]+np.arange(nw_chip[1])*dlogw),
            10.**(logw0_chip[2]+np.arange(nw_chip[2])*dlogw)]

def gridPix(apStar=True) :
    """ Returns chip pixel ranges in apStar or aspcap grid
    """
    if apStar :
        w=np.log10(apStarWave())
        s1 = np.where(np.isclose(w,logw0_chip[0],rtol=0.))[0][0]
        s2 = np.where(np.isclose(w,logw0_chip[1],rtol=0.))[0][0]
        s3 = np.where(np.isclose(w,logw0_chip[2],rtol=0.))[0][0]
        e1 = np.where(np.isclose(w,logw0_chip[0]+nw_chip[0]*dlogw,rtol=0.))[0][0]
        e2 = np.where(np.isclose(w,logw0_chip[1]+nw_chip[1]*dlogw,rtol=0.))[0][0]
        e3 = np.where(np.isclose(w,logw0_chip[2]+nw_chip[2]*dlogw,rtol=0.))[0][0]
        return [[s1,e1],[s2,e2],[s3,e3]]
    else :
        return [[0,nw_chip[0]],[nw_chip[0],nw_chip[0]+nw_chip[1]],[nw_chip[0]+nw_chip[1],nw_chip[0]+nw_chip[1]+nw_chip[2]]]

def aspcap2apStar(aspcap):
    """ from input aspcap spectrum on aspcap grid, return spectrum on apStar grid 
    """
    apstar=np.zeros(nw_apStar)
    pix_out=gridPix()
    pix_in=gridPix(apStar=False)
    for pin,pout in zip(pix_in,pix_out) :
        apstar[pout[0]:pout[1]] = aspcap[pin[0]:pin[1]]
    return apstar

def apStar2aspcap(apstar):
    """ from input aspcap spectrum on aspcap grid, return spectrum on apStar grid 
    """
    aspcap=np.zeros(nw_chip.sum())
    pix_out=gridPix()
    pix_in=gridPix(apStar=False)
    for pin,pout in zip(pix_in,pix_out) :
        aspcap[pin[0]:pin[1]] = apstar[pout[0]:pout[1]] 
    return aspcap

def readstars(starlist,libpar) :
    '''
    Runs stars in starlist through FERRE using libpar
    '''
    for star in starlist :
        spec,err=readstar(star)
        cont=cont_normalize(spec)

def ferre(spec,libpar,te0=None,logg0=None,mh0=None) :
    '''
    Runs FERRE for a set of input spectra
    '''
    np.savetxt(out+'.frd',spec)
    np.savetxt(out+'.err',spec)
    np.savetxt(out+'.con',spec)

    # write ferre control file
    ferre.writenml()

    out=ferre.readferre(name)

def elemsens(els=None,plot=None,ylim=[0.1,-0.3],teff=4750,logg=2.,feh=-1.,smooth=None) :
    '''
    Returns and optionally plots wavelength sensitivity to changes in elemental abundances for specified elements from MOOG mini-elem grid
    '''
    elem=fits.open(os.environ['APOGEE_REDUX']+'/speclib/moog/elemsens.fits')
    if els is None :
        els = elems()[0]
    elif type(els) == str :
        els = [els]
    wave=[]
    out=[]
    for el in els :
        for i in range(1,25) :
            card='HDU{:02d}'.format(i)
            try :
              if elem[0].header[card].strip().upper() == el.strip().upper() :
                it=int(round((teff-elem[i].header['CRVAL2'])/elem[i].header['CDELT2']))
                ig=int(round((logg-elem[i].header['CRVAL3'])/elem[i].header['CDELT3']))
                ife=int(round((feh-elem[i].header['CRVAL4'])/elem[i].header['CDELT4']))
                diff=elem[i].data[ife,ig,it,:]
                if smooth is not None:
                    diff=scipy.ndimage.filters.gaussian_filter(diff,smooth)
                wave=elem[i].header['CRVAL1']+np.arange(elem[i].header['NAXIS1'])*elem[i].header['CDELT1']
                if plot is not None:
                    #plot.plot(wave,diff,color='g')
                    plot.plot(wave,diff)
                    plot.set_ylim(ylim[0],ylim[1])
                out.append(diff)
            except: pass
    if len(out) == 1 :
        return wave, out[0]
    else :
        return wave, out

def sensplot(ax=None,offset=0) :
    if ax is None :
        fig,ax=plots.multi(1,2,hspace=0.001,sharex=True)
    els=['O','Mg','Si','S','Ca','Ti','Na','Al','K','P']
    cols=['r','g','b','c','y','m','r','g','b','c']
    ls=['-','-','-','-','-','-',':',':',':',':']
    for i in range(len(els)) :
        w,s=elemsens(els=els[i])
        plots.plotl(ax[0],w,s+offset,label=els[i],color=cols[i],ls=ls[i])
    ax[0].legend(fontsize='small')

    elems=['C','CI','N','O','Na','Mg','Al','Si','P','S','K','Ca','Ti','TiII','V','Cr','Mn','Fe','Co','Ni','Cu','Ge','Ce','Rb','Y','Nd']
    els=['V','Cr','Mn','Co','Ni','Cu','Ge','Ce','Rb','Nd']
    cols=['r','g','b','c','y','m','r','g','b','c']
    ls=['-','-','-','-','-','-',':',':',':',':']
    for i in range(len(els)) :
        w,s=elemsens(els=els[i])
        plots.plotl(ax[1],w,s+offset,label=els[i],color=cols[i],ls=ls[i])
    ax[1].legend(fontsize='small')
    #elems=['C','CI','N','O','Na','Mg','Al','Si','P','S','K','Ca','Ti','TiII','V','Cr','Mn','Fe','Co','Ni','Cu','Ge','Ce','Rb','Y','Nd']
    

def data(str,loc=None) :       
    '''
    Add apogeeObject data to structure
    '''
    add=np.empty(1,dtype=[('RA','f4'),('DEC','f4'),('J','f4'),('H','f4'),('K','f4'),('AK_TARG','f4'),('SFD_EBV','f4')])
    new=struct.add_cols(str,add)
    for i in range(len(str)) :
        name = str['APOGEE_ID'][i]
        print(i, name)
        apogee_id = name.split('.')[0].split('-')[2]
        if loc is None :
            loc = name.split('.')[1].split('_')[0]
        ap=apload.ApLoad()
        s=ap.apStar(loc,apogee_id)
        field=s[0].header['FIELD']
        try :
            obj=fits.open(os.environ['APOGEE_TARGET']+'/apogeeObject/apogeeObject_'+field+'.fits')[1].data
        except :
            obj=fits.open(os.environ['APOGEE_TARGET']+'/apogee2object/apogee2object_'+field+'.fits')[1].data
        j=np.where(obj['APOGEE_ID'] == apogee_id)[0]
        for card in add.dtype.names :
            new[card][i]=obj[card][j]
    return new

def elemmask(el,maskdir='filters_26112015',plot=None,yr=[0,1]) :
    '''
    '''
    mask=np.loadtxt(os.getenv('SPECLIB_DIR')+'/lib/'+maskdir+'/'+el+'.filt') 
    wave=np.loadtxt(os.getenv('SPECLIB_DIR')+'/lib/'+maskdir+'/wave.dat')
    if plot is not None :
        plots.plotl(plot,wave,mask,yr=yr)
    return wave,mask
