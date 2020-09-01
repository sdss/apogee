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
import matplotlib
try: matplotlib.use('Agg')
except: pass

import copy
import numpy as np
import glob
import os
import shutil
import pdb
import time
import yaml
from shutil import copyfile
import subprocess
import scipy.ndimage.filters
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table, Column
#from holtz.tools import struct
from tools import plots
from tools import match
from tools import html
from apogee.utils import apload
from apogee.utils import bitmask
from apogee.utils import gaia
from apogee.utils import spectra
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

    elems=['C','CI','N','O','Na','Mg','Al','Si','P','S','K','Ca','Ti','TiII','V','Cr','Mn','Fe','Co','Ni','Cu','Ge','Rb','Ce','Nd','Yb']
    #return,['C','N','O','Na','Mg','Al','Si','S','K','Ca','Ti','V','Mn','Fe','Ni','Nd']
    #return,['C','N','O','Na','Mg','Al','Si','S','K','Ca','Ti','V','Mn','Fe','Ni']
    elemtoh=[0,0,0,0,1,0,1,0,1,0,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1] #0,0,0,0,0,0,0,0,0]
    tagnames=[]
    elemfitnames=[]
    for i in range(len(elems) ) :
        if elemtoh[i] :
            tagnames.append(elems[i]+'_Fe')
            elemfitnames.append('['+elems[i]+'/Fe]')
        else :
            tagnames.append(elems[i]+'_M')
            elemfitnames.append('['+elems[i]+'/M]')
    return np.array(elems),np.array(elemtoh),np.array(tagnames),np.array(elemfitnames)

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
    aspcap=np.zeros(nw_chip.sum(),dtype=apstar.dtype)
    pix_out=gridPix()
    pix_in=gridPix(apStar=False)
    for pin,pout in zip(pix_in,pix_out) :
        aspcap[pin[0]:pin[1]] = apstar[pout[0]:pout[1]] 
    return aspcap

def link(src,dest) :
    try : os.remove(dest)
    except : pass
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
    """ Reads element window file, optional plot
    """
    mask=np.loadtxt(os.getenv('SPECLIB_DIR')+'/lib/'+maskdir+'/'+el+'.filt') 
    wave=np.loadtxt(os.getenv('SPECLIB_DIR')+'/lib/'+maskdir+'/wave.dat')
    if plot is not None :
        plots.plotl(plot,wave,mask,yr=yr)
    return wave,mask

def apField2aspcapField(planfile,nobj=None,minerr=0.005,apstar_vers='stars') :
    """ create initial aspcapField file from apField file
    """

    # read configuration file
    plan=yaml.safe_load(open(planfile,'r'))
    apred=plan['apred_vers']
    aspcap_vers=plan['aspcap_vers']
    aspcap_config=plan['aspcap_config']
    instrument=plan['instrument']
    telescope=plan['telescope']
    field=plan['field']

    # setup reader and load apField file
    load = apload.ApLoad(apred=apred,aspcap=aspcap_vers,telescope=telescope)
    apfieldname=load.filename('Field',field=field)
    if apstar_vers != 'stars' :
        apfieldname=apfieldname.replace('/stars/','/'+apstar_vers+'/')
        apfield = fits.open(apfieldname)[1].data
    else :
        apfield=load.apField(field)[1].data
    print('apField file: ', apfieldname)

    # add GAIA data
    apfield=gaia.add_gaia(apfield)

    aspcapfield=Table(apfield)
    if nobj is not None : aspcapfield=aspcapfield[0:nobj]
    try : test = aspcapfield['MEANFIB']
    except : aspcapfield['MEANFIB'] = 150
    # add new columns
    nparam = len(params()[0])

    # read ASPCAP configuration
    config = yaml.safe_load(open(os.environ['APOGEE_DIR']+'/config/aspcap/'+aspcap_config+'/'+instrument+'.yml','r'))
 
    # add tags to structure
    ngrids=len(config['grids'])
    aspcapfield.add_column(Column(name='CLASS',dtype='S8',length=len(aspcapfield)))
    aspcapfield.add_column(Column(name='FPARAM_CLASS',dtype=float,shape=(ngrids,nparam),length=len(aspcapfield)))
    aspcapfield.add_column(Column(name='FPARAM_COV_CLASS',dtype=float,shape=(ngrids,nparam,nparam),length=len(aspcapfield)))
    aspcapfield.add_column(Column(name='CHI2_CLASS',dtype=float,shape=(ngrids),length=len(aspcapfield)))
    aspcapfield.add_column(Column(name='FPARAM',dtype=float,shape=(nparam),length=len(aspcapfield)))
    aspcapfield.add_column(Column(name='FPARAM_COV',dtype=float,shape=(nparam,nparam),length=len(aspcapfield)))
    aspcapfield.add_column(Column(name='PARAM_CHI2',dtype=float,length=len(aspcapfield)))
    aspcapfield.add_column(Column(name='PARAM',dtype=float,shape=(nparam),length=len(aspcapfield)))
    aspcapfield.add_column(Column(name='PARAM_COV',dtype=float,shape=(nparam,nparam),length=len(aspcapfield)))
    aspcapfield.add_column(Column(name='PARAMFLAG',dtype=np.uint64,shape=(nparam),length=len(aspcapfield)))
    aspcapfield.add_column(Column(name='ASPCAPFLAG',dtype=np.uint64,length=len(aspcapfield)))
    aspcapfield.add_column(Column(name='ASPCAPFLAGS',dtype='S132',length=len(aspcapfield)))
    nelem = len(elems()[0])
    aspcapfield.add_column(Column(name='FELEM',dtype=float,shape=(nelem),length=len(aspcapfield)))
    aspcapfield.add_column(Column(name='FELEM_ERR',dtype=float,shape=(nelem),length=len(aspcapfield)))
    aspcapfield.add_column(Column(name='X_H',dtype=float,shape=(nelem),length=len(aspcapfield)))
    aspcapfield.add_column(Column(name='X_H_ERR',dtype=float,shape=(nelem),length=len(aspcapfield)))
    aspcapfield.add_column(Column(name='X_M',dtype=float,shape=(nelem),length=len(aspcapfield)))
    aspcapfield.add_column(Column(name='X_M_ERR',dtype=float,shape=(nelem),length=len(aspcapfield)))
    aspcapfield.add_column(Column(name='ELEM_CHI2',dtype=float,shape=(nelem),length=len(aspcapfield)))
    aspcapfield.add_column(Column(name='ELEMFLAG',dtype=float,shape=(nelem),length=len(aspcapfield)))

    # load spectra
    # create table for output spectral data
    aspcapspec = Table()
    nwave = nw_chip.sum()
    aspcapspec.add_column(Column(name='SPEC',dtype=float,shape=(nwave),length=len(aspcapfield)))
    aspcapspec.add_column(Column(name='SPEC_ERR',dtype=float,shape=(nwave),length=len(aspcapfield)))
    aspcapspec.add_column(Column(name='MASK',dtype=float,shape=(nwave),length=len(aspcapfield)))
    aspcapspec.add_column(Column(name='SPEC_BESTFIT',dtype=float,shape=(nwave),length=len(aspcapfield)))
    aspcapspec.add_column(Column(name='OBS',dtype=float,shape=(nwave),length=len(aspcapfield)))
    aspcapspec.add_column(Column(name='ERR',dtype=float,shape=(nwave),length=len(aspcapfield)))
    aspcapspec.add_column(Column(name='NORM',dtype=float,shape=(nwave),length=len(aspcapfield)))

    pixelmask=bitmask.PixelBitMask()
    badval=pixelmask.badval()|pixelmask.getval('SIG_SKYLINE')
    aspcapmask=bitmask.AspcapBitMask()
    for istar,star in enumerate(aspcapfield) :

        apstarfile=load.filename('Star',field=field,obj=star['APOGEE_ID'])
        if apstar_vers != 'stars' :
            apstarfile=apstarfile.replace('/stars/','/'+apstar_vers+'/')
            try: 
                hdulist=fits.open(apstarfile)
                wave=spectra.fits2vector(hdulist[1].header,1)
                apstar=apload.ApSpec(hdulist[1].data,header=hdulist[0].header,
                              err=hdulist[2].data,bitmask=hdulist[3].data,wave=wave,
                              sky=hdulist[4].data,skyerr=hdulist[5].data,
                              telluric=hdulist[6].data,telerr=hdulist[7].data)
            except: 
                print('No apStar file found: ',apstarfile)
                try: aspcapfield['ASPCAPFLAG'][istar] |= aspcapmask.getval('MISSING_APSTAR')
                except: pdb.set_trace()
                apstar=None
        else :
            apstar=load.apStar(field,star['APOGEE_ID'],load=True)
        print('apstarfile: ',apstarfile)
        if apstar is None: continue

        norm=np.nanmedian(apStar2aspcap(apstar.flux[0,:]))
        aspcapspec['OBS'][istar] = apStar2aspcap(apstar.flux[0,:])/norm
        mask= np.where((apStar2aspcap(apstar.bitmask[0,:]) & badval) > 0)[0]
        tmp = apStar2aspcap(apstar.err[0,:])
        tmp[mask] *= 100.
        # set uncertainty floor to minerr (fractional)
        ind = np.where(tmp/apStar2aspcap(apstar.flux[0,:]) < minerr)[0]
        tmp[ind] = minerr*apStar2aspcap(apstar.flux[0,:])[ind]
        aspcapspec['ERR'][istar] = tmp/norm
        bd=np.where(~np.isfinite(aspcapspec['ERR'][istar]))[0]
        if len(bd) > 0 : 
            aspcapspec['ERR'][istar,bd] = 1.e10
        bd=np.where(~np.isfinite(aspcapspec['OBS'][istar]))[0]
        if len(bd) > 0 : 
            aspcapspec['OBS'][istar,bd] = 0.
            aspcapspec['ERR'][istar,bd] = 1.e10
        aspcapspec['MASK'][istar] = apStar2aspcap(apstar.bitmask[0,:])
        aspcapspec['NORM'][istar] = 1.

    return aspcapfield, aspcapspec

def fit_params(planfile,aspcapdata=None,clobber=False,nobj=None,write=True,minerr=0.005,apstar_vers=None,plot=False,
               init='RV',fix=None,renorm=False,suffix='',html=True) :
    """ run ASPCAP on a field to get parameters
    """
    # read configuration file
    plan=yaml.safe_load(open(planfile,'r'))
    apred=plan['apred_vers']
    if apstar_vers is None : apstar_vers=plan['apstar_vers'] if plan.get('apstar_vers') else 'stars'
    aspcap_vers=plan['aspcap_vers']
    aspcap_config=plan['aspcap_config']
    instrument=plan['instrument']
    telescope=plan['telescope']
    field=plan['field']

    # get initial aspcapField if not provided
    if aspcapdata is None : 
        aspcapfield,aspcapspec=apField2aspcapField(planfile,nobj=nobj,minerr=minerr,apstar_vers=apstar_vers)
    else :
        aspcapfield = copy.deepcopy(aspcapdata[0])
        aspcapspec = copy.deepcopy(aspcapdata[1])

    # output directory
    load = apload.ApLoad(apred=apred,aspcap=aspcap_vers,telescope=telescope)
    outfield=load.filename('aspcapField',field=field)
    outdir=os.path.dirname(outfield)

    # read ASPCAP configuration
    config = yaml.safe_load(open(os.environ['APOGEE_DIR']+'/config/aspcap/'+aspcap_config+'/'+instrument+'.yml','r'))

    # pixel masking
    pixelmask=bitmask.PixelBitMask()
    badval=pixelmask.badval()|pixelmask.getval('SIG_SKYLINE')

    # set NOGRID bit until we have a grid
    aspcapmask=bitmask.AspcapBitMask()
    aspcapfield['ASPCAPFLAG'] |= aspcapmask.getval('NO_GRID')

    # reset CHI2 if we are iterating on a previous solution
    aspcapfield['PARAM_CHI2'] = 0.
    pars=params()[0]

    # loop over all grids
    param_class=[]
    spec_class=[]
    chi2_class=[]
    for igrid,grid in enumerate(config['grids']) :

        # set up output FERRE directory for this grid
        if  init == 'RV' :
            out=outdir+'/ferre/class_'+grid['name']+'/'+grid['name']+'-'+field
        else :
            out=outdir+'/ferre/class_'+grid['name']+'_'+init+'/'+grid['name']+'-'+field

        try: os.makedirs(os.path.dirname(out))
        except: pass
        link(os.environ['APOGEE_SPECLIB']+'/synth/',os.path.dirname(out)+'/lib')

        # get FERRE library information
        libfile = 'lib/'+grid['lib']+'.hdr'
        libhead0,libhead = ferre.rdlibhead(os.path.dirname(out)+'/'+libfile)
        libhead0['FILE'] = libfile
        nparams=libhead0['N_OF_DIM']
        # get the index numbers in input parameter array for correct library parameter order
        index=np.zeros(nparams,dtype=int)
        for i in range(nparams) : index[i] = np.where(params()[0] == libhead0['LABEL'][i].decode())[0]

        # select stars for this grid
        if init == 'RV' :
            gd = np.where((aspcapfield['RV_TEFF'] >= grid['teff_range'][0]) &
                          (aspcapfield['RV_TEFF'] <= grid['teff_range'][1]) &
                          (aspcapfield['RV_LOGG'] >= grid['logg_range'][0]) &
                          (aspcapfield['RV_LOGG'] <= grid['logg_range'][1]) &
                          (aspcapfield['MEANFIB'] >= grid['fibermin']) &
                          (aspcapfield['MEANFIB'] < grid['fibermax']) ) [0]
        else :
            gd = np.where(aspcapfield['CLASS'] == grid['name'])[0]

        print(grid['name'],len(gd))
        if len(gd) == 0 : continue

        # turn off NO_GRID
        aspcapfield['ASPCAPFLAG'][gd] &= ~aspcapmask.getval('NO_GRID')

        # loop over stars and accumulate input for FERRE
        inpars=[]
        flux=[]
        err=[]
        stars=[]
        for star,spec in zip(aspcapfield[gd],aspcapspec[gd]) :
            if star['ASPCAPFLAG'] & aspcapmask.getval('MISSING_APSTAR') : continue

            stars.append(star['APOGEE_ID'])
            if init == 'RV' :
                vmicro=np.array([0.372160,-0.090531,-0.000802,0.001263,-0.027321])
                logg=star['RV_LOGG']
                vm = vmicro[0]+vmicro[1]*logg+vmicro[2]*logg**2+vmicro[3]*logg**3
                inpars.append([star['RV_TEFF'],star['RV_LOGG'],vm,star['RV_FEH'],0.,0.,0.,1.,0.])
            elif init == 'FPARAM' :
                inpars.append(list(star['FPARAM']))
            elif init == 'PARAM' :
                inpars.append(list(star['PARAM']))
            if renorm :
                flux.append(spec['OBS']/spec['NORM'])
                err.append(spec['ERR']/spec['NORM'])
            else :
                flux.append(spec['OBS'])
                err.append(spec['ERR'])
      

        # write FERRE files and run FERRE
        if clobber or not os.path.exists(out+'.spm') :
            if fix is not None :
                indv=[]
                for i,par in enumerate(libhead0['LABEL']) :
                    if par.decode() not in fix : indv.append(i+1)
            else : indv=None
            ferre.writeipf(out,libfile,stars,param=np.array(inpars))
            ferre.writespec(out+'.obs',flux)
            ferre.writespec(out+'.err',err)
            ferre.writenml(out+'.nml',os.path.basename(out),libhead0,init=0,indv=indv,
                       algor=grid['algor'],ncpus=plan['ncpus'],
                       obscont=grid['obscont'],rejectcont=grid['rejectcont'],
                       renorm=abs(grid['renorm']),
                       filterfile=os.environ['APOGEE_DIR']+'/data/windows/'+grid['mask'])
            fout=open(out+'.stdout','w')
            ferr=open(out+'.stderr','w')
            print('running ferre.x: ', os.path.basename(out)+'.nml')
            start = time.time()
            subprocess.call(['ferre.x',os.path.basename(out)+'.nml'],shell=False,
                            cwd=os.path.dirname(out),stdout=fout,stderr=ferr)
            print('elapsed: ',time.time()-start)
            fout.close()
            ferr.close()

        # read FERRE output
        param,spec,wave=ferre.read(out,libfile)
        # fill in locked parameters
        fill_plock(param,grid['PLOCK'])
        param_class.append(param) 
        spec_class.append(spec) 

        # load into apcapField
        for istar,star in enumerate(aspcapfield[gd]) :
            i = np.where(param['APOGEE_ID'] == star['APOGEE_ID'].encode())[0]
            if len(i) == 0 : continue
            aspcapfield['FPARAM_CLASS'][gd[istar],igrid,:] = param['FPARAM'][i]
            aspcapfield['FPARAM_COV_CLASS'][gd[istar],igrid,:] = param['FPARAM_COV'][i]
            aspcapfield['CHI2_CLASS'][gd[istar],igrid] = param['PARAM_CHI2'][i]
            # if this is the lowest CHI2, store in FPARAM and PARAM_CHI2
            # penalize GK grid in favor of finer M grid
            chi2 = param['PARAM_CHI2'][i]
            if ('GK' in grid['name']) and (param['FPARAM'][i,0] < 3985) : 
                print('penalizing GK')
                chi2 *= 10
            if chi2 < aspcapfield['PARAM_CHI2'][gd[istar]] or aspcapfield['PARAM_CHI2'][gd[istar]]<=0 :
                aspcapfield['CLASS'][gd[istar]] = grid['name']
                aspcapfield['PARAM_CHI2'][gd[istar]] = param['PARAM_CHI2'][i]
                aspcapfield['FPARAM'][gd[istar]] = param['FPARAM'][i]
                aspcapfield['FPARAM_COV'][gd[istar]] = param['FPARAM_COV'][i]
                aspcapfield['ASPCAPFLAG'][gd[istar]] = param['ASPCAPFLAG'][i]
                aspcapfield['ASPCAPFLAGS'][gd[istar]] = aspcapmask.getname(aspcapfield['ASPCAPFLAG'][gd[istar]])
                aspcapspec['SPEC'][gd[istar]] = spec['frd'][i]
                aspcapspec['SPEC_BESTFIT'][gd[istar]] =  spec['mdl'][i]
                # continuum correct
                width=151
                p1=0
                for ichip in range(len(libhead)) :
                    npix=libhead[ichip]['NPIX']
                    obs=spec['frd'][i,p1:p1+npix].flatten()
                    # replace bad pixels with median filtered value
                    med=scipy.ndimage.filters.median_filter(obs,[5*width],mode='nearest')
                    bd=np.where(obs < 0.01)[0]
                    obs[bd] = med[bd]
                    # populate error with FERRE-normalized error
                    err= spec['err'][i,p1:p1+npix]*spec['frd'][i,p1:p1+npix]/obs
                    bd=np.where(~np.isfinite(err))[0]
                    err[bd] = 1.e10
                    aspcapspec['SPEC_ERR'][gd[istar]][p1:p1+npix] = err
                  
                    # get ratio of observed / model, make sure edge is reasonable number
                    mdl=spec['mdl'][i,p1:p1+npix].flatten()
                    ratio=obs/mdl
                    corr=scipy.ndimage.filters.median_filter(ratio,[width],mode='nearest')
                    bd=np.where(~np.isfinite(corr) | np.isclose(corr,0.))[0]
                    corr[bd] = 1.
                    if not renorm: 
                        aspcapspec['NORM'][gd[istar]][p1:p1+npix] = corr
                    p1+=npix

                if plot :
                    plt.figure()
                    plt.plot(aspcapspec['OBS'][gd[istar]])
                    plt.plot(aspcapspec['SPEC'][gd[istar]])
                    plt.plot(aspcapspec['SPEC_BESTFIT'][gd[istar]])
                    plt.plot(aspcapspec['NORM'][gd[istar]])
                    plt.plot(aspcapspec['ERR'][gd[istar]])
                    plt.ylim(0,1.5)
                    plt.draw()
                    pdb.set_trace()
                    plt.close()

                # get mask from apStar, and supplement with FERRE mask
                try:
                    filterfile=os.environ['APOGEE_DIR']+'/data/windows/'+grid['mask']
                    fmask=np.loadtxt(filterfile)
                    bd = np.where(fmask < 0.001)[0]
                    aspcapspec['MASK'][gd[istar]][bd] !=  pixelmask.getval('FERRE_MASK')
                except: pdb.set_trace()


    # Results in astropy tables
    aspcapfield=Table(aspcapfield)
    aspcapspec=Table(aspcapspec)
    wave=np.hstack(gridWave())
    aspcapkey=np.empty(1,dtype=[('WAVE','f4',len(wave)),
                          ('PARAM_SYMBOL','S16',len(params()[0])),
                          ('ELEM_SYMBOL','S5',len(elems()[0])),
                          ('ELEMTOH','i4',len(elems()[0])),
                          ('ELEM_VALUE','S12',len(elems()[2])),
                          ('CLASSES','S5',len(config['grids']))])
    aspcapkey['PARAM_SYMBOL'] = params()[1]
    aspcapkey['ELEM_SYMBOL'] = elems()[0]
    aspcapkey['ELEMTOH'] = elems()[1]
    aspcapkey['ELEM_VALUE'] = elems()[2]
    grids=[]
    for grid in config['grids'] : grids.append(grid['name'])
    aspcapkey['CLASSES'] = grids
    aspcapkey=Table(aspcapkey)

    if write :
        hdulist=fits.HDUList()
        hdulist.append(fits.table_to_hdu(aspcapfield))
        hdulist.append(fits.table_to_hdu(aspcapspec))
        hdulist.append(fits.table_to_hdu(aspcapkey))

        #output aspcapField
        outfield=load.filename('aspcapField',field=field)
        try: os.makedirs(os.path.dirname(outfield))
        except: pass
        outfile=os.path.dirname(outfield)+'/'+os.path.splitext(os.path.basename(outfield))[0]+suffix+'.fits'
        hdulist.writeto(outfile,overwrite=True)

    if html :
        # create output HTML page
        mkhtml(field,suffix='',apred=apred,aspcap_vers=aspcap_vers,telescope=telescope)

    return aspcapfield,aspcapspec,aspcapkey

def fit_elems(planfile,aspcapdata=None,clobber=False,nobj=None,write=True,calib=False,renorm=True,suffix='',html=True) :
    """ run ASPCAP on a field for elemental abundances
    """

    # read configuration file
    plan=yaml.safe_load(open(planfile,'r'))
    apred=plan['apred_vers']
    aspcap_vers=plan['aspcap_vers']
    aspcap_config=plan['aspcap_config']
    instrument=plan['instrument']
    telescope=plan['telescope']
    field=plan['field']

    # get aspcapField if not provided
    load = apload.ApLoad(apred=apred,aspcap=aspcap_vers,telescope=telescope)
    outfield=load.filename('aspcapField',field=field)
    outdir=os.path.dirname(outfield)
    if aspcapdata is None : 
        aspcapfield=Table(fits.open(outfield)[1].data)
        aspcapspec=Table(fits.open(outfield)[2].data)
        aspcapkey=Table(fits.open(outfield)[3].data)
    else :
        aspcapfield = copy.deepcopy(aspcapdata[0])
        aspcapspec = copy.deepcopy(aspcapdata[1])
        aspcapkey = copy.deepcopy(aspcapdata[2])

    # output directory

    # read ASPCAP configuration
    config = yaml.safe_load(open(os.environ['APOGEE_DIR']+'/config/aspcap/'+aspcap_config+'/'+instrument+'.yml','r'))

    if calib : useparam='PARAM'
    else : useparam ='FPARAM'

    # loop over grids
    for igrid,grid in enumerate(config['grids']) :

        # output directory for spectra
        outspec=outdir+'/ferre/spectra/'+grid['name']+'-'+field
        try: os.makedirs(os.path.dirname(outspec))
        except: pass

        libfile = 'lib_'+grid['name']+'/'+grid['lib']+'.hdr'
        # get stars for which best fit was in this grid
        gd = np.where(aspcapfield['CLASS'] == grid['name'])[0]
        print('Grid: ',grid['name'],len(gd))
        if len(gd) == 0 : continue

        if clobber or not os.path.exists(outspec+'.obs') :
            link(os.environ['APOGEE_SPECLIB']+'/synth/',os.path.dirname(outspec)+'/../lib_'+grid['name'])

            # get FERRE library information for this grid
            libhead0,libhead = ferre.rdlibhead(outdir+'/ferre/'+libfile)
            libhead0['FILE'] = libfile

            # loop over stars and accumulate input for FERRE
            inpars=[]
            flux=[]
            err=[]
            stars=[]
            for star,spec in zip(aspcapfield[gd],aspcapspec[gd]) :
                stars.append(star['APOGEE_ID'])
                inpars.append(star[useparam])
                print(star['APOGEE_ID'])
                if renorm :
                    flux.append(spec['OBS']/spec['NORM'])
                    err.append(spec['ERR']/spec['NORM'])
                else :
                    flux.append(spec['OBS'])
                    err.append(spec['ERR'])

            ferre.writeipf(outspec,outdir+'/ferre/'+libfile,stars,param=np.array(inpars))
            ferre.writespec(outspec+'.obs',flux)
            ferre.writespec(outspec+'.err',err)

        # loop over elements
        fp=open(outdir+'/ferre/'+grid['name']+'.nmlfiles','w')
        fit = False
        for ielem,elem in enumerate(config['elems']) :
        
            # set up output FERRE directory for this grid
            if calib:
                out=outdir+'/ferre/elem_'+elem['name']+'_PARAM/'+elem['name']+'-'+grid['name']+'-'+field
            else :
                out=outdir+'/ferre/elem_'+elem['name']+'/'+elem['name']+'-'+grid['name']+'-'+field
            try: os.makedirs(os.path.dirname(out))
            except: pass

            # write FERRE files
            if clobber or not os.path.exists(out+'.spm') :
                link(os.environ['APOGEE_SPECLIB']+'/synth/',os.path.dirname(out)+'/lib')
                filterfile = elem['name']+'.mask'
                shutil.copyfile(os.environ['APOGEE_DIR']+'/data/windows/'+grid['windows']+'/'+filterfile,
                                     os.path.dirname(out)+'/'+elem['name']+'.mask')

                # links for ipf, obs, and err files
                for ext in ['.ipf','.obs','.err'] :
                    try: os.remove(out+ext)
                    except: pass
                    link('../spectra/'+os.path.basename(outspec)+ext,out+ext)
                index = np.where(libhead0['LABEL'] == elem['griddim'].encode())[0][0]+1
                if elem['griddim'] == 'METALS' : 
                    ttie=[]
                    for dim in ['C','N','O Mg Si S Ca Ti'] :
                        j=np.where(libhead0['LABEL'] == dim.encode())[0]
                        if len(j) > 0 : ttie.append(j[0]+1)
                        else : ttie.append(-1)
                else : ttie=[-1,-1,-1]
                ferre.writenml(out+'.nml','elem_'+elem['name']+'/'+os.path.basename(out),libhead0,init=0,
                           nov=1,indv=[index],ttie=ttie,
                           algor=grid['algor'],ncpus=plan['ncpus'],
                           obscont=grid['obscont'],rejectcont=grid['rejectcont'],
                           renorm=abs(grid['renorm']),
                           filterfile='elem_'+elem['name']+'/'+filterfile)
                fp.write('elem_'+elem['name']+'/'+os.path.basename(out)+'.nml\n')
                fit = True

        fp.close()
        # run FERRE for all elements for this grid
        if fit :
            fout=open(outdir+'/ferre/'+grid['name']+'.stdout','w')
            ferr=open(outdir+'/ferre/'+grid['name']+'.stderr','w')
            subprocess.call(['ferre.x','-l',grid['name']+'.nmlfiles'],shell=False,
                            cwd=outdir+'/ferre',stdout=fout,stderr=ferr)
            fout.close()
            ferr.close()

        for ielem,elem in enumerate(config['elems']) :
            out=outdir+'/ferre/elem_'+elem['name']+'/'+elem['name']+'-'+grid['name']+'-'+field
            # read FERRE output
            param,spec,wave=ferre.read(out,outdir+'/ferre/'+libfile)
            # fill in locked parameters
            fill_plock(param,grid['PLOCK'])
            # location for this element in FELEM array
            jelem=np.where(elems()[0] == elem['name'])[0]

            # load into aspcapField
            for istar,star in enumerate(aspcapfield[gd]) :
                i = np.where(param['APOGEE_ID'] == star['APOGEE_ID'].encode())[0]
                if len(i) == 0 : continue
                index = np.where(params()[0] == elem['griddim'])[0]
                try:
                    aspcapfield['FELEM'][gd[istar],jelem] = param['FPARAM'][i,index]
                    if param['FPARAM_COV'][i,index,index] > 0 :
                        aspcapfield['FELEM_ERR'][gd[istar],jelem] = np.sqrt(param['FPARAM_COV'][i,index,index])
                    aspcapfield['ELEM_CHI2'][gd[istar],jelem] = param['PARAM_CHI2'][i]
                except: pdb.set_trace()

    # Results into an HDUList 

    if write :
        hdulist=fits.HDUList()
        hdulist.append(fits.table_to_hdu(aspcapfield))
        hdulist.append(fits.table_to_hdu(aspcapspec))
        hdulist.append(fits.table_to_hdu(aspcapkey))

        #output aspcapField
        outfield=load.filename('aspcapField',field=field)
        try: os.makedirs(os.path.dirname(outfield))
        except: pass
        outfile=os.path.dirname(outfield)+'/'+os.path.splitext(os.path.basename(outfield))[0]+suffix+'.fits'
        hdulist.writeto(outfile,overwrite=True)

    if html :
        # create output HTML page
        mkhtml(field,suffix='',apred=apred,aspcap_vers=aspcap_vers,telescope=telescope)

    return aspcapfield,aspcapspec,aspcapkey


def mkhtml(field,suffix='',apred='r13',aspcap_vers='l33',telescope='apo25m') :
    """ Create ASPCAP field web page and plots
    """

    matplotlib.use('Agg')

    load = apload.ApLoad(apred=apred,aspcap=aspcap_vers,telescope=telescope)
    infile=load.filename('aspcapField',field=field)
    a=fits.open(infile)
    aspcapfield=a[1].data
    aspcapspec=a[2].data

    outdir=os.path.dirname(infile)
    try: os.makedirs(outdir+'/plots/')
    except: pass
    fp=open(outdir+'/'+field+suffix+'.html','w')
    fp.write('<HTML>\n')
    fp.write('<HEAD><script type=text/javascript src=../../../html/sorttable.js></script></head>')
    fp.write('<BODY>\n')
    fp.write('<H2> Field: {:s}</H2><p>\n'.format(field))

    # parameter plots: Kiel and [alpha/M] vs [M/H]
    fig,ax=plots.multi(1,1)
    plots.plotc(ax,aspcapfield['FPARAM'][:,0],aspcapfield['FPARAM'][:,1],aspcapfield['FPARAM'][:,3],
                xr=[8000,3000],yr=[6,-1],zr=[-2,0.5],size=10,colorbar=True,xt='Teff',yt='logg',zt='[M/H]')
    fig.savefig(outdir+'/plots/'+field+'_hr.png')
    plt.close()
    fp.write('<TD><A HREF=plots/{:s}_hr.png><IMG SRC=plots/{:s}_hr.png></A>\n'.format(field,field))

    yr=[-0.3,0.75]
    fig,ax=plots.multi(1,1)
    plots.plotc(ax,aspcapfield['FPARAM'][:,3],aspcapfield['FPARAM'][:,6],aspcapfield['FPARAM'][:,0],
                zr=[3000,8000],yr=yr,xr=[-2.5,0.5],size=10,colorbar=True,xt='[M/H]',yt='[alpha/M]',zt='Teff')
    fig.savefig(outdir+'/plots/'+field+'_alpha.png')
    plt.close()
    fp.write('<TD><A HREF=plots/{:s}_alpha.png><IMG SRC=plots/{:s}_alpha.png></A>\n'.format(field,field))

    # individual element abundance plots
    fig,ax=plots.multi(1,4,hspace=0.001)
    yr=[-0.3,1.5]
    for iel,el in enumerate(['C','CI','N']) :
        jel = np.where(elems()[0] == el)[0][0]
        yt='['+el+'/M]'
        plots.plotc(ax[iel],aspcapfield['FPARAM'][:,3],aspcapfield['FELEM'][:,jel],aspcapfield['FPARAM'][:,0],
                zr=[3000,8000],yr=yr,xr=[-2.5,0.5],size=10,colorbar=True,xt='[M/H]',yt=yt,zt='Teff',label=[0.9,0.9,el])
    plots.plotc(ax[3],aspcapfield['FPARAM'][:,3],aspcapfield['FELEM'][:,0]-aspcapfield['FELEM'][:,2],aspcapfield['FPARAM'][:,0],
            zr=[3000,8000],yr=yr,xr=[-2.5,0.5],size=10,colorbar=True,xt='[M/H]',yt='[C/N]',zt='Teff',label=[0.9,0.9,'[C/N]'])
    fig.savefig(outdir+'/plots/'+field+'_cnelem.png')
    plt.close()
    fp.write('<TD><A HREF=plots/{:s}_cnelem.png><IMG SRC=plots/{:s}_cnelem.png></A>\n'.format(field,field))

    fig,ax=plots.multi(1,6,hspace=0.001)
    yr=[-0.3,0.75]
    for iel,el in enumerate(['O','Mg','Si','S','Ca','Ti']) :
        jel = np.where(elems()[0] == el)[0][0]
        yt='['+el+'/M]'
        plots.plotc(ax[iel],aspcapfield['FPARAM'][:,3],aspcapfield['FELEM'][:,jel],aspcapfield['FPARAM'][:,0],
                zr=[3000,8000],yr=yr,xr=[-2.5,0.5],size=10,colorbar=True,xt='[M/H]',yt=yt,zt='Teff',label=[0.9,0.9,el])
    fig.savefig(outdir+'/plots/'+field+'_alphaelem.png')
    plt.close()
    fp.write('<TD><A HREF=plots/{:s}_alphaelem.png><IMG SRC=plots/{:s}_alphaelem.png></A>\n'.format(field,field))

    fig,ax=plots.multi(1,4,hspace=0.001)
    for iel,el in enumerate(['Na','Al','P','K']) :
        jel = np.where(elems()[0] == el)[0][0]
        yt='['+el+'/M]'
        try :
          plots.plotc(ax[iel],aspcapfield['FPARAM'][:,3],aspcapfield['FELEM'][:,jel]-aspcapfield['FPARAM'][:,3],aspcapfield['FPARAM'][:,0],
                zr=[3000,8000],yr=yr,xr=[-2.5,0.5],size=10,colorbar=True,xt='[M/H]',yt=yt,zt='Teff',label=[0.9,0.9,el])
        except: pdb.set_trace()
    fig.savefig(outdir+'/plots/'+field+'_oddzelem.png')
    plt.close()
    fp.write('<TD><A HREF=plots/{:s}_oddzelem.png><IMG SRC=plots/{:s}_oddzelem.png></A>\n'.format(field,field))

    fig,ax=plots.multi(1,7,hspace=0.001)
    for iel,el in enumerate(['V','Cr','Mn','Fe','Co','Ni','Cu'] ) :
        jel = np.where(elems()[0] == el)[0][0]
        yt='['+el+'/M]'
        plots.plotc(ax[iel],aspcapfield['FPARAM'][:,3],aspcapfield['FELEM'][:,jel]-aspcapfield['FPARAM'][:,3],aspcapfield['FPARAM'][:,0],
                zr=[3000,8000],yr=yr,xr=[-2.5,0.5],size=10,colorbar=True,xt='[M/H]',yt=yt,zt='Teff',label=[0.9,0.9,el])
    fig.savefig(outdir+'/plots/'+field+'_feelem.png')
    plt.close()
    fp.write('<TD><A HREF=plots/{:s}_feelem.png><IMG SRC=plots/{:s}_feelem.png></A>\n'.format(field,field))

    fp.write('<BR>Click on column headers to sort by column value<BR>\n')
    fp.write('<TABLE BORDER=2 CLASS=sortable>\n')
    fp.write('<TR><TD>Obj<TD>Grid<TD>CHI2\n')
    parnames=params()[0]
    for ipar in range(8) :
        parname=parnames[ipar]
        if ipar == 2 : parname = 'VMICRO'
        if ipar == 7 : parname = 'VSINI/ VMACRO'
        fp.write('<TD>{:s}\n'.format(parname))

    w = np.hstack(gridWave())
    pix = gridPix(apStar=False)    
    for istar,star in enumerate(aspcapfield['APOGEE_ID']) :
        #fig,ax = plots.multi(1,3,hspace=0.5,figsize=(12,6))
        fig,ax = plots.multi(1,1,figsize=(18,3))
        gd = np.where(aspcapspec[istar]['err'] < 0.5)[0]
        plots.plotl(ax,w.flatten(),aspcapspec[istar]['spec'],yr=[0.,1.2],color='k')
        plots.plotl(ax,w.flatten()[gd],aspcapspec[istar]['spec'][gd],yr=[0.,1.2],color='g')
        plots.plotl(ax,w.flatten(),aspcapspec[istar]['spec_bestfit'],yr=[0.,1.2],color='r')
        plots.plotl(ax,w.flatten(),aspcapspec[istar]['err'],yr=[0.,1.2],color='m')
        #for ichip in range(3) :
        #    plots.plotl(ax,w[ichip],aspcapspec[istar]['spec'][pix[ichip][0]:pix[ichip][1]],yr=[0.,1.2],color='k')
        #    plots.plotl(ax,w[ichip],aspcapspec[istar]['spec_bestfit'][pix[ichip][0]:pix[ichip][1]])
        #    plots.plotl(ax,w[ichip],aspcapspec[istar]['err'][pix[ichip][0]:pix[ichip][1]],color='r')
        pars=r'T$_e$: {:8.1f} logg: {:5.2f} log(v$_{{micro}}$): {:5.2f}  [M/H]: {:5.2f}  [C/M]: {:5.2f}  [N/M]: {:5.2f} [$\alpha$/M]: {:5.2f}  log(vsini/v$_{{macro}}$): {:5.2f}'.format(*aspcapfield['FPARAM'][istar,0:8])
        fig.suptitle(star+' : '+pars)
        fig.savefig(outdir+'/plots/'+star+'.png')
        plt.close()

        fp.write('<TR><TD>{:s}<BR>\n'.format(star))
        fp.write('H  = {:7.2f}<br>'.format(aspcapfield['H'][istar]))
        fp.write('SNR  = {:7.2f}<br>'.format(aspcapfield['SNR'][istar]))
        fp.write('<FONT COLOR=blue> TARGFLAGS </FONT>: {:s}<BR>\n'.format(aspcapfield['TARGFLAGS'][istar]))
        fp.write('<FONT COLOR=blue> STARFLAGS </FONT>: {:s}<BR>\n'.format(aspcapfield['STARFLAGS'][istar]))
        fp.write('<FONT COLOR=blue> ASPCAPFLAGS </FONT>: {:s}<BR>\n'.format(aspcapfield['ASPCAPFLAGS'][istar]))
        fp.write('<TD>{:s}\n'.format(aspcapfield['CLASS'][istar]))
        fp.write('<TD>{:6.1f}\n'.format(aspcapfield['PARAM_CHI2'][istar]))
        for ipar in range(8) :
            p = aspcapfield['FPARAM'][istar,ipar]
            if aspcapfield['FPARAM_COV'][istar,ipar,ipar] > 0 :
                perr = np.sqrt(aspcapfield['FPARAM_COV'][istar,ipar,ipar])
            else : perr = -999.
            if ipar==2 or ipar==7 : 
                p = 10.**p
                if perr > 0 : perr = perr * p * np.log(10)
            pm = '&plusmn;'
            fp.write(r'<TD>{:8.2f}{:s}{:8.2f}'.format(p,pm,perr))
        fp.write('<TD><A HREF=plots/{:s}.png><IMG SRC=plots/{:s}.png></A>\n'.format(star,star))

    fp.write('</TABLE></BODY></HTML>')
    fp.close()

def fill_plock(a,plocks) :
    """ Fill FPARAM array in input structure with values of locked parameters
    """
    te_index=np.where(params()[0] == 'TEFF')[0]
    logg_index=np.where(params()[0] == 'LOGG')[0]
    mh_index=np.where(params()[0] == 'METALS')[0]
    for plock in plocks :
        index = np.where(params()[0] == plock['name'])[0]
        a['FPARAM'][:,index] = (plock['const']+
                               plock['te_coef']*a['FPARAM'][:,te_index]+
                               plock['logg_coef'][0]*a['FPARAM'][:,logg_index]+
                               plock['logg_coef'][1]*a['FPARAM'][:,logg_index]**2+
                               plock['logg_coef'][2]*a['FPARAM'][:,logg_index]**3+
                               plock['mh_coef']*a['FPARAM'][:,mh_index])
       
def chi2(data) :
    """ calculate chi2 and cumulative chi2
    """ 
    chi2=(data['SPEC']-data['SPEC_BESTFIT'])**2/data['ERR']**2
    cumchi2=np.cumsum(chi2,axis=1)
    return chi2, cumchi2

def comp(a,b,terange=[4500,5000],loggrange=[2.5,3.5]) :
    """ Compare two versions in region on Kiel diagram
    """
    matplotlib.use('TkAgg')
    gd=np.where( (a[1].data['FPARAM'][:,0]>=terange[0]) &
                 (a[1].data['FPARAM'][:,0]<terange[1]) &
                 (a[1].data['FPARAM'][:,1]>=loggrange[0]) &
                 (a[1].data['FPARAM'][:,1]<loggrange[1]) )[0]
    wave=np.hstack(gridWave())
    achi,acum=chi2(a[2].data)
    bchi,bcum=chi2(b[2].data)
    fig,ax=plots.multi(1,4,hspace=0.001,sharex=True) 
    plots.plotl(ax[0],wave,np.nanmedian(bcum-acum,axis=0))
    plots.plotl(ax[1],wave,np.nanmedian(b[2].data['SPEC_BESTFIT']-a[2].data['SPEC_BESTFIT'],axis=0))
    pdb.set_trace()
    for i in gd :
        print(i)
        plots.plotl(ax[2],wave,(bcum[i]-acum[i])/(bcum[i]-acum[i]).sum())
        plots.plotl(ax[3],wave,a[2].data['SPEC'][i])
        plots.plotl(ax[3],wave,a[2].data['SPEC_BESTFIT'][i])
        plots.plotl(ax[3],wave,b[2].data['SPEC_BESTFIT'][i])

