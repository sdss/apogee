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
import glob
import os
import pdb
from shutil import copyfile
import subprocess
import scipy.ndimage.filters
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
#from holtz.tools import struct
from tools import plots
from apogee.utils import apload
from apogee.speclib import isochrones
try: from apogee.aspcap import ferre
except: pass

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

def getparams(name,lib,coarse=None,n=None) :

    objs=ascii.read(name+'.ipf')['col1']
    if n is not None : objs = objs[0:n]

    # if we have a coarse library, run FERRE with it first
    #   to get starting guesses
    if coarse is not None :
        l=ferre.rdlibhead(coarse+'.hdr')[0]
        link(name+'.obs',name+'_coarse.obs')
        link(name+'.err',name+'_coarse.err')
        ferre.writeipf(name+'_coarse',coarse+'.hdr',objs)
        ferre.writenml(name+'_coarse.nml',name+'_coarse',l,algor=3,renorm=4,obscont=1,ncpus=32,init=1)
        copyfile(name+'_coarse.nml','input.nml')
        subprocess.call(['ferre.x'],shell=False)
        out,outspec,outwave=ferre.read(name+'_coarse',coarse+'.hdr')
        ferre.writeipf(name+'_1',lib+'.hdr',objs,param=out['FPARAM']) 
    else :
        ferre.writeipf(name+'_1',lib+'.hdr',objs)
    l=ferre.rdlibhead(lib+'.hdr')[0]
    ferre.writenml(name+'_1.nml',name+'_1',l,algor=3,renorm=4,obscont=1,ncpus=32,init=1)
    copyfile(name+'_1.nml','input.nml')
    link(name+'.obs',name+'_1.obs')
    link(name+'.err',name+'_1.err')
    subprocess.call(['ferre.x'],shell=False)
    out,outspec,outwave=ferre.read(name+'_1',lib+'.hdr')


def link(src,dest) :
    try :
        os.remove(dest)
    except :
        pass
    os.symlink(src,dest)
   
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

def intplot(a=None,param='FPARAM',indir='cal',apred='r10',aspcap='t33b',verbose=False) :
    """ Given input structure, plot HR diagram, and enter event loop to mark stars to plot spectra
    """

    load=apload.ApLoad(apred=apred,aspcap=aspcap,verbose=verbose)
    if a is None : a=load.allCal()[1].data

    fig,ax = hr(a,param=param)
    plots.event(fig)
    sf,sa=plots.multi(1,1)
    nplot = 11
    hf,ha=plots.multi(1,nplot,figsize=(8.5,11),hspace=0.2)
    ha2=[]
    for i in range(nplot) : ha2.append(ha[i].twinx())
    print('hit any key near object to plot in HR diagram, q to quit')
    while (1) :
        ret=plots.mark(fig)
        if ret[2] == 'q' : break
        ind=plots._index[0]
        try :
            load.settelescope='apo25m'
            try : f=load.aspcapField(a['ALTFIELD'][ind])
            except : f=load.aspcapField(a['FIELD'][ind])
        except :
            load.settelescope='lco25m'
            try : f=load.aspcapField(a['ALTFIELD'][ind])
            except : f=load.aspcapField(a['FIELD'][ind])
        if f is None :
            f=glob.glob(indir+'/*'+a['APOGEE_ID'][plots._index[0]]+'*')
            dir=os.path.dirname(f[0])
            f=glob.glob(dir+'/*aspcapField*.fits')
            f=fits.open(f[0])
        print('f: ',f)
        data=f[1].data
        j=np.where(data['APOGEE_ID'] == a['APOGEE_ID'][plots._index[0]])[0][0]
        sa.cla()
        for i in range(11) : 
            ha[i].cla()
            ha2[i].cla()
            ha[i].set_ylabel('Flux')
            ha2[i].set_ylabel(r'$\chi^2$')
            ha[i].set_ylim(0.5,1.3)
            ha2[i].set_ylim(0.,20.)
        plot(10.**f[3].data['WAVE'][0],f[2].data['SPEC'][j,:],ax=ha,sum=True,color='k')
        plot(10.**f[3].data['WAVE'][0],f[2].data['SPEC_BESTFIT'][j,:],ax=ha,sum=True,color='b')
        chi2 = (f[2].data['SPEC'][j,:]-f[2].data['SPEC_BESTFIT'][j,:])**2/f[2].data['ERR'][j,:]**2
        plot(10.**f[3].data['WAVE'][0],chi2,ax=ha2,sum=True,alpha=0.4)
        plots.plotl(sa,10.**f[3].data['WAVE'][0],f[2].data['SPEC'][j,:])
        plots.plotl(sa,10.**f[3].data['WAVE'][0],f[2].data['SPEC_BESTFIT'][j,:])
        text1=r'ID: {:s} FIELD: {:s} SNR: {:6.1f} $\chi^2$: {:6.1f}'.format(
             data['APOGEE_ID'][j],data['FIELD'][j],data['SNR'][j],data['PARAM_CHI2'][j])
        text2=r'Teff: {:5.0f} logg: {:5.1f} [M/H]: {:5.2f} [$\alpha$/M]: {:5.2f} [C/M]: {:5.2f} [N/M]: {:5.2f}'.format(
             data[param][j,0],data[param][j,1],data[param][j,3],data[param][j,6],data[param][j,4],data[param][j,5])
        sf.suptitle(text1+'\n'+text2)
        hf.suptitle(text1+'\n'+text2)
        plt.draw()
        plt.show()
    plt.close(hf)
    plt.close(sf)
    plt.close(fig)

def hr(a,param='FPARAM',colorbar=False,zt='[M/H]',zr=None,iso=False, hard=None, grid=None) :
    """ Plot an HR diagram from input structure

        Args:
            all  : structure that includes stellar parameter array with (Teff, logg, ...)
            param : tag to use (default='FPARAM')
            colorbar : show colorbar? (default= False)
    """
    fig,ax = plots.multi(1,1)
    if grid is None :
        teff=a[param][:,0]
        logg=a[param][:,1]
    else :
        teff=a['FPARAM_CLASS'][:,grid,0]
        logg=a['FPARAM_CLASS'][:,grid,1]
    if zt == '[M/H]' : 
        z=a[param][:,3]
        if zr is None : zr=[-2,0.5]
    elif zt == 'chi2' : 
        z=a['PARAM_CHI2']
        if zr is None : zr=[0,10]
    plots.plotc(ax,teff,logg,z,xr=[8000,3000],yr=[6,-1],zr=zr,
                xt='Teff',yt='log g',zt=zt,colorbar=colorbar)
    plots._data = a
    if iso:
        colors=['b','g','k','r']
        for i,z in enumerate(['zm20','zm10','zp00','zp05']) :
            isodata=isochrones.read(os.environ['ISOCHRONE_DIR']+'/'+z+'.dat',agerange=[8.99,9.01])
            isochrones.plot(ax,isodata,'teff','logg',color=colors[i],alpha=0.3)
            isodata=isochrones.read(os.environ['ISOCHRONE_DIR']+'/'+z+'.dat',agerange=[9.99,10.01])
            isochrones.plot(ax,isodata,'teff','logg',color=colors[i],alpha=0.3)
    if hard is not None: fig.savefig(hard)
    return fig,ax

def multihr(a,param='FPARAM',colorbar=False,hard=None) :
    fig,ax = plots.multi(2,4,hspace=0.001,wspace=0.001,figsize=(8,12))

    z=a[param][:,3]
    zr=[-2,0.5]
    zt='[M/H]'
    plots.plotc(ax[0,0],a[param][:,0],a[param][:,1],z,xr=[8000,3000],yr=[6,-1],zr=zr,
                xt='Teff',yt='log g',zt=zt,colorbar=colorbar)
    ax[0,0].text(0.05,0.9,zt,transform=ax[0,0].transAxes)

    z=a['PARAM_CHI2']
    zr=[0,10]
    zt='CHI2'
    plots.plotc(ax[0,1],a[param][:,0],a[param][:,1],z,xr=[8000,3000],yr=[6,-1],zr=zr,
                xt='Teff',zt=zt,colorbar=colorbar)
    ax[0,1].text(0.05,0.9,zt,transform=ax[0,1].transAxes)

    z=a[param][:,4]
    zr=[-1,1.0]
    zt='[C/M]'
    plots.plotc(ax[1,0],a[param][:,0],a[param][:,1],z,xr=[8000,3000],yr=[6,-1],zr=zr,
                xt='Teff',yt='log g',zt=zt,colorbar=colorbar)
    ax[1,0].text(0.05,0.9,zt,transform=ax[1,0].transAxes)

    z=a[param][:,5]
    zr=[-1,1.0]
    zt='[N/M]'
    plots.plotc(ax[1,1],a[param][:,0],a[param][:,1],z,xr=[8000,3000],yr=[6,-1],zr=zr,
                xt='Teff',zt=zt,colorbar=colorbar)
    ax[1,1].text(0.05,0.9,zt,transform=ax[1,1].transAxes)

    z=a[param][:,6]
    zr=[-1,1.0]
    zt=r'[$\alpha$/M]'
    plots.plotc(ax[2,0],a[param][:,0],a[param][:,1],z,xr=[8000,3000],yr=[6,-1],zr=zr,
                xt='Teff',yt='log g',zt=zt,colorbar=colorbar)
    ax[2,0].text(0.05,0.9,zt,transform=ax[2,0].transAxes)

    z=a[param][:,4]-a[param][:,5]
    zr=[-1,1.0]
    zt='[C/N]'
    plots.plotc(ax[2,1],a[param][:,0],a[param][:,1],z,xr=[8000,3000],yr=[6,-1],zr=zr,
                xt='Teff',zt=zt,colorbar=colorbar)
    ax[2,1].text(0.05,0.9,zt,transform=ax[2,1].transAxes)

    z=10.**a[param][:,2]
    zr=[0.3,4]
    zt='vmicro'
    plots.plotc(ax[3,0],a[param][:,0],a[param][:,1],z,xr=[8000,3000],yr=[6,-1],zr=zr,
                xt='Teff',zt=zt,colorbar=colorbar)
    ax[3,0].text(0.05,0.9,zt,transform=ax[3,0].transAxes)

    z=10.**a[param][:,7]
    zr=[0,10]
    zt='vrot'
    plots.plotc(ax[3,1],a[param][:,0],a[param][:,1],z,xr=[8000,3000],yr=[6,-1],zr=zr,
                xt='Teff',zt=zt,colorbar=colorbar)
    ax[3,1].text(0.05,0.9,zt,transform=ax[3,1].transAxes)

    if hard is not None: fig.savefig(hard)

#def hrclass(a) :
#    fig,ax=plots.multi(2,3)
#    for class in ['GKg','GKg_O','Mg','Md','GKd','Fd'] :

def plot(wave,spec,color=None,ax=None,hard=None,sum=False,title=None,alpha=None) :
    """  Multipanel plots of APOGEE spectra
    """
    if sum : ny=11
    else : ny=10
    if ax is None : fig,ax=plots.multi(1,ny,figsize=(8.5,11),hspace=0.2)
    for i in range(10) :
        plots.plotl(ax[i],wave,spec,xr=[15000+i*200,15200+i*200],color=color,linewidth=0.3,alpha=alpha)
        ax[i].xaxis.label.set_size(6)
        ax[i].yaxis.label.set_size(6)
        ax[i].tick_params(axis = 'both', which = 'major', labelsize = 6)
        ax[i].xaxis.set_minor_locator(plt.MultipleLocator(10.))
    if sum :
        plots.plotl(ax[10],wave,spec,xr=[15100,17000],color=color,linewidth=0.3,alpha=alpha)
        ax[10].xaxis.label.set_size(6)
        ax[10].yaxis.label.set_size(6)
        ax[10].tick_params(axis = 'both', which = 'major', labelsize = 6)
        ax[10].xaxis.set_minor_locator(plt.MultipleLocator(100.))
    try: 
        if title is not None : fig.suptitle(title)
    except: pass
    if hard is not None : fig.savefig(hard)

    try: return fig,ax
    except: return
