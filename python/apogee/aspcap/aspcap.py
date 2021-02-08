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
from scipy.ndimage.filters import median_filter, gaussian_filter
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table, Column, vstack
#from holtz.tools import struct
from tools import plots
from tools import match
from tools import html
from apogee.utils import apload
from apogee.utils import bitmask
try: from apogee.utils import gaia
except: print('gaia not available')
from apogee.utils import spectra
try: from apogee.aspcap import ferre
except: print('ferre not available')

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

    elems=['C','CI','N','O','Na','Mg','Al','Si','P','S','K','Ca','Ti','TiII','V','Cr','Mn','Fe','Co','Ni','Cu','Ge','Rb','Ce','Nd','Yb','C13']
    #return,['C','N','O','Na','Mg','Al','Si','S','K','Ca','Ti','V','Mn','Fe','Ni','Nd',]
    #return,['C','N','O','Na','Mg','Al','Si','S','K','Ca','Ti','V','Mn','Fe','Ni']
    elemtoh=[0,0,0,0,1,0,1,0,1,0,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0] #0,0,0,0,0,0,0,0,0]
    tagnames=[]
    elemfitnames=[]
    for i,el in enumerate(elems) :
        if el == 'Fe' : 
            tagnames.append('{:s}_H'.format(el).upper())
            elemfitnames.append('[{:s}/H]'.format(el))
        else :
            tagnames.append('{:s}_Fe'.format(el).upper())
            if elemtoh[i] :
                elemfitnames.append('[{:s}/H]'.format(el))
            else :
                elemfitnames.append('[{:s}/M]'.format(el))

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
                    diff=gaussian_filter(diff,smooth)
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

def apField2aspcapField(planfile,minerr=0.005,apstar_vers='stars',addgaia=False) :
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
    if type(field) is str : field=[field]
    visits = plan['visits'] if plan.get('visits') else 0
    nobj = plan['nobj'] if plan.get('nobj') else None
    obj = plan['obj'] if plan.get('obj') else None

    # setup reader and load apField file
    load = apload.ApLoad(apred=apred,aspcap=aspcap_vers,telescope=telescope)

    # if we have multiple fields, read them all and stack them
    apfield=[]
    for f in field :
        apfieldname=load.filename('Field',field=f)
        if apstar_vers != 'stars' :
            apfieldname=apfieldname.replace('/stars/','/'+apstar_vers+'/')
            apfield.append(Table.read(apfieldname))
        else :
            apfield.append(Table(load.apField(f)[1].data))
        print('apField file: ', apfieldname)
    apfield=vstack(apfield)

    # add GAIA data
    if addgaia: apfield=gaia.add_gaia(apfield)

    # create out output table
    aspcapfield=Table(apfield)

    # limited or single object specified?
    if nobj is not None : aspcapfield=aspcapfield[0:nobj]
    if obj is not None :
        if isinstance(obj,str) : obj=[obj]
        j=[]
        for o in obj :
            j.append(np.where(aspcapfield['APOGEE_ID'] == o)[0][0])
        aspcapfield=aspcapfield[j]

    #if we want individual visits, need to add rows to table
    if visits > 0 :
        aspcapfield.add_column(Column(name='VISIT',dtype=int,length=len(aspcapfield)))
        newfield = aspcapfield[:0].copy() 
        for row in aspcapfield :
            newfield.add_row(row)
            #if visits=1 take ALL visits, else take min(visits,row['NVISITS'])
            if visits == 1 : nvisits = row['NVISITS']
            else : nvisits = np.min([visits,row['NVISITS']])
            if nvisits > 1 :
                for ivisit in range(nvisits) :
                    row['VISIT'] = ivisit+1
                    newfield.add_row(row)
        aspcapfield = newfield

    # read ASPCAP configuration
    config = yaml.safe_load(open(os.environ['APOGEE_DIR']+'/config/aspcap/'+aspcap_config+'/'+instrument+'.yml','r'))
 
    # add new columns
    for col in ['ASPCAP_GRID','FPARAM_GRID','CHI2_GRID','FPARAM','FPARAM_COV','ASPCAP_CHI2',
                'PARAM','PARAM_COV','PARAMFLAG','ASPCAPFLAG','ASPCAPFLAGS','FRAC_BADPIX','FRAC_LOWSNR','FRAC_SIGSKY',
                'FELEM','FELEM_ERR','X_H','X_H_ERR','X_M','X_M_ERR','ELEM_CHI2','ELEMFLAG'] :
        try : aspcapfield.remove_column(col)
        except KeyError: pass

    nparam = len(params()[0])
    # add tags to structure
    ngrids=len(config['grids'])
    n=len(aspcapfield)
    aspcapfield.add_column(Column(name='ASPCAP_GRID',dtype='S8',data=np.full([n],'None')))
    aspcapfield.add_column(Column(name='FPARAM_GRID',dtype=np.float32,data=np.full([n,ngrids,nparam],np.nan)))
    aspcapfield.add_column(Column(name='CHI2_GRID',dtype=np.float32,data=np.full([n,ngrids],np.nan)))
    aspcapfield.add_column(Column(name='FPARAM',dtype=np.float32,data=np.full([n,nparam],np.nan)))
    aspcapfield.add_column(Column(name='FPARAM_COV',dtype=np.float32,data=np.full([n,nparam,nparam],np.nan)))
    aspcapfield.add_column(Column(name='ASPCAP_CHI2',dtype=np.float32,data=np.full([n],np.nan)))
    aspcapfield.add_column(Column(name='PARAM',dtype=np.float32,data=np.full([n,nparam],np.nan)))
    aspcapfield.add_column(Column(name='PARAM_COV',dtype=np.float32,data=np.full([n,nparam,nparam],np.nan)))
    aspcapfield.add_column(Column(name='PARAMFLAG',dtype=np.int,data=np.full([n,nparam],0)))
    aspcapfield.add_column(Column(name='ASPCAPFLAG',dtype=np.int64,data=np.full([n],0)))
    aspcapfield.add_column(Column(name='ASPCAPFLAGS',dtype='S256',data=np.full([n],'')))
    aspcapfield.add_column(Column(name='FRAC_BADPIX',dtype=np.float32,data=np.full([n],np.nan)))
    aspcapfield.add_column(Column(name='FRAC_LOWSNR',dtype=np.float32,data=np.full([n],np.nan)))
    aspcapfield.add_column(Column(name='FRAC_SIGSKY',dtype=np.float32,data=np.full([n],np.nan)))
    nelem = len(elems()[0])
    aspcapfield.add_column(Column(name='FELEM',dtype=np.float32,data=np.full([n,nelem],np.nan)))
    aspcapfield.add_column(Column(name='FELEM_ERR',dtype=float,data=np.full([n,nelem],np.nan)))
    aspcapfield.add_column(Column(name='X_H',dtype=np.float32,data=np.full([n,nelem],np.nan)))
    aspcapfield.add_column(Column(name='X_H_ERR',dtype=np.float32,data=np.full([n,nelem],np.nan)))
    aspcapfield.add_column(Column(name='X_M',dtype=np.float32,data=np.full([n,nelem],np.nan)))
    aspcapfield.add_column(Column(name='X_M_ERR',dtype=np.float32,data=np.full([n,nelem],np.nan)))
    aspcapfield.add_column(Column(name='ELEM_CHI2',dtype=np.float32,data=np.full([n,nelem],np.nan)))
    aspcapfield.add_column(Column(name='ELEMFRAC',dtype=np.float32,data=np.full([n,nelem],np.nan)))
    aspcapfield.add_column(Column(name='ELEMFLAG',dtype=np.int,data=np.full([n,nelem],np.nan)))

    # load spectra
    # create table for output spectral data
    aspcapspec = Table()
    nwave = nw_chip.sum()
    aspcapspec.add_column(Column(name='SPEC',dtype=np.float32,shape=(nwave,),length=len(aspcapfield)))
    aspcapspec.add_column(Column(name='SPEC_ERR',dtype=np.float32,shape=(nwave,),length=len(aspcapfield)))
    aspcapspec.add_column(Column(name='MASK',dtype=np.int,shape=(nwave,),length=len(aspcapfield)))
    aspcapspec.add_column(Column(name='SPEC_BESTFIT',dtype=np.float32,shape=(nwave,),length=len(aspcapfield)))
    aspcapspec.add_column(Column(name='OBS',dtype=np.float32,shape=(nwave,),length=len(aspcapfield)))
    aspcapspec.add_column(Column(name='ERR',dtype=np.float32,shape=(nwave,),length=len(aspcapfield)))
    aspcapspec.add_column(Column(name='NORM',dtype=np.float32,shape=(nwave,),length=len(aspcapfield)))

    pixelmask=bitmask.PixelBitMask()
    badval=pixelmask.badval()|pixelmask.getval('SIG_SKYLINE')
    aspcapmask=bitmask.AspcapBitMask()
    #skip=open('skip.txt','a+')
    for istar,star in enumerate(aspcapfield) :

        apstarfile=load.filename('Star',field=star['FIELD'],obj=star['APOGEE_ID'])
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
                apstar=None
        else :
            apstar=load.apStar(star['FIELD'],star['APOGEE_ID'],load=True)
        print('apstarfile: ',apstarfile)
        if apstar is None: 
            print('No apStar file found: ',apstarfile)
            aspcapfield['ASPCAPFLAG'][istar] |= aspcapmask.getval('MISSING_APSTAR')
            aspcapfield['ASPCAPFLAG'][istar] |= aspcapmask.getval('NO_ASPCAP_RESULT')
            continue

        # if we have an individual visit, use it, and load appropriate SNR
        row = 0
        if visits > 0 :
            if star['VISIT'] == 0 : 
                row=0
            else : 
                row=star['VISIT']+1
                star['SNR'] = apstar.header['SNRVIS{:d}'.format(star['VISIT'])]
                star['APOGEE_ID'] = '{:s}.{:d}'.format(star['APOGEE_ID'],star['VISIT'])
        print(star['APOGEE_ID'])

        norm=np.nanmedian(apStar2aspcap(apstar.flux[row,:]))
        aspcapspec['OBS'][istar] = apStar2aspcap(apstar.flux[row,:])/norm
        aspcapspec['MASK'][istar] = apStar2aspcap(apstar.bitmask[row,:])
        aspcapspec['NORM'][istar] = 1.

        # enhance error around sky lines
        tmp = apStar2aspcap(apstar.err[row,:])
        mask= np.where((aspcapspec['MASK'][istar] & pixelmask.getval('SIG_SKYLINE')) > 0)[0]
        tmp[mask] *= 100.
        aspcapfield['FRAC_SIGSKY'][istar] = len(mask) / len(aspcapspec['OBS'][istar])

        # set uncertainty floor to minerr (fractional): check for negative fluxes!
        ind = np.where(tmp/np.abs(apStar2aspcap(apstar.flux[row,:])) < minerr)[0]
        #cont = gaussian_filter(median_filter(apstar.flux[row,:],[501],mode='reflect'),100)
        ##tmp[ind] = minerr*np.abs(apStar2aspcap(apstar.flux[row,:]))[ind]
        #tmp[ind] = minerr*np.abs(cont)[ind]
        #aspcapspec['ERR'][istar] = tmp/norm
        tmp /= norm
        ind = np.where(tmp < minerr)[0]
        tmp[ind] = minerr
        aspcapspec['ERR'][istar] = tmp
  
        # fraction of low S/N pixels 
        fracerr = aspcapspec['ERR'][istar]/aspcapspec['OBS'][istar] 
        aspcapfield['FRAC_LOWSNR'][istar] = len(np.where(fracerr > 0.2)[0]) / len(aspcapspec['OBS'][istar])
        # skip if lowsnr fraction > 0.5 
        if aspcapfield['FRAC_LOWSNR'][istar] > 0.5 :
            aspcapfield['ASPCAPFLAG'][istar] |= aspcapmask.getval('BAD_FRAC_BADPIX')
            aspcapfield['ASPCAPFLAG'][istar] |= aspcapmask.getval('NO_ASPCAP_RESULT')
            #skip.write('{:s} {:s} FRAC_LOWSNR: {:8.2f}\n'.format(aspcapfield['APOGEE_ID'][istar],aspcapfield['FIELD'][istar],aspcapfield['FRAC_LOWSNR'][istar]))

        # make sure bad pixels are not included
        bd=np.where(~np.isfinite(aspcapspec['ERR'][istar]) | 
                    ~np.isfinite(aspcapspec['OBS'][istar]) |
                    (aspcapspec['MASK'][istar] & pixelmask.badval()) > 0 |
                    (aspcapspec['OBS'][istar] < 0.) |
                    (aspcapspec['ERR'][istar] < 0.) )[0]
        if len(bd) > 0 : 
            # make sure there are no NaNs input to FERRE
            aspcapspec['OBS'][istar,bd] = 0.0001
            aspcapspec['ERR'][istar,bd] = 1.e10
        # fraction of bad pixels
        aspcapfield['FRAC_BADPIX'][istar] = len(bd) / len(aspcapspec['OBS'][istar])
        # skip if bad pixel fraction > 0.5 or > 0.33 in any individual chip
        if aspcapfield['FRAC_BADPIX'][istar] > 0.5 :
            aspcapfield['ASPCAPFLAG'][istar] |= aspcapmask.getval('BAD_FRAC_BADPIX')
            aspcapfield['ASPCAPFLAG'][istar] |= aspcapmask.getval('NO_ASPCAP_RESULT')
        p0=0
        for pix in nw_chip :
            badfrac = len(np.where((bd>=p0)&(bd<p0+pix))[0]) / pix
            if badfrac > 0.33 :
                #skip.write('{:s} {:s} FRAC_BADPIX: {:d} {:8.2f}\n'.format(aspcapfield['APOGEE_ID'][istar],aspcapfield['FIELD'][istar],p0,badfrac))
                aspcapfield['ASPCAPFLAG'][istar] |= aspcapmask.getval('BAD_FRAC_LOWSNR')
                aspcapfield['ASPCAPFLAG'][istar] |= aspcapmask.getval('NO_ASPCAP_RESULT')
                p0+=pix

    #skip.close()

    # add element window pixel fractions
    for iel,el in enumerate(elems()[0]) :
        aspcapfield['ELEMFRAC'][:,iel] = elemfrac(aspcapspec,el)

    return aspcapfield, aspcapspec

def fit_params(planfile,aspcapdata=None,clobber=False,write=True,minerr=0.005,apstar_vers=None,plot=False,
               init='RV',coarse=False,fix=None,renorm=False,suffix='',html=True,mult=False,dostar=None) :
    """ run ASPCAP on a field to get parameters
    """
    # read configuration file
    plan=yaml.safe_load(open(planfile,'r'))
    if plan['apogee_ver'] != os.environ['APOGEE_VER'] :
        print('apogee_ver {:s} does not match running version {:s}'.format(plan['apogee_ver'],os.environ['APOGEE_VER']))
        pdb.set_trace()

    apred=plan['apred_vers']
    if apstar_vers is None : apstar_vers=plan['apstar_vers'] if plan.get('apstar_vers') else 'stars'
    aspcap_vers=plan['aspcap_vers']
    aspcap_config=plan['aspcap_config']
    instrument=plan['instrument']
    telescope=plan['telescope']
    field=plan['field']
    if 'outfield' not in plan.keys() : outfield = field
    else : outfield = plan['outfield']

    # output directory
    load = apload.ApLoad(apred=apred,aspcap=aspcap_vers,telescope=telescope)
    if isinstance(field,list) : bundleload = apload.ApLoad(apred=apred,aspcap=aspcap_vers,telescope='bundle_'+telescope)
    else : bundleload = load
    outfile=bundleload.filename('aspcapField',field=outfield)
    outdir=os.path.dirname(outfile)

    # get initial aspcapField if not provided
    if aspcapdata is None : 
        aspcapfield,aspcapspec=apField2aspcapField(planfile,minerr=minerr,apstar_vers=apstar_vers)
    elif aspcapdata=='read' :
        aspcapfield = fits.open(outfile)[1].data
        aspcapspec = fits.open(outfile)[2].data
    else :
        aspcapfield = copy.deepcopy(aspcapdata[0])
        aspcapspec = copy.deepcopy(aspcapdata[1])

    # read ASPCAP configuration
    config = yaml.safe_load(open(os.environ['APOGEE_DIR']+'/config/aspcap/'+aspcap_config+'/'+instrument+'.yml','r'))

    # pixel masking
    pixelmask=bitmask.PixelBitMask()
    badval=pixelmask.badval()|pixelmask.getval('SIG_SKYLINE')

    # set NOGRID bit until we have a grid
    aspcapmask=bitmask.AspcapBitMask()
    parammask=bitmask.ParamBitMask()
    aspcapfield['ASPCAPFLAG'] |= aspcapmask.getval('NO_GRID')

    # reset CHI2 if we are iterating on a previous solution
    aspcapfield['ASPCAP_CHI2'] = 0.
    pars=params()[0]

    if mult :
        aspcapfield=Table(aspcapfield)
        aspcapspec=Table(aspcapspec)
        nparam = len(params()[0])
        aspcapfield.add_column(Column(name='FPARAM_MULT',dtype=np.float32,shape=(27,nparam),length=len(aspcapfield)))
        nwave = nw_chip.sum()
        aspcapspec.add_column(Column(name='SPEC_BESTFIT_MULT',dtype=np.float32,shape=(27,nwave),length=len(aspcapfield)))

    # loop over all grids
    for igrid,grid in enumerate(config['grids']) :

        # set up output FERRE directory for this grid
        if  init == 'RV' :
            out=outdir+'/ferre/class_'+grid['name']+'/'+grid['name']+'-'+outfield
        else :
            out=outdir+'/ferre/class_'+grid['name']+'_'+init+'/'+grid['name']+'-'+outfield

        if mult: out=out+'_mult'

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
            gd = np.where(aspcapfield['ASPCAP_GRID'] == grid['name'])[0]

        if dostar is not None :
            gdstar=np.where(aspcapfield['APOGEE_ID'][gd] == dostar)[0]
            gd=gd[gdstar]

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
            if star['ASPCAPFLAG'] & aspcapmask.getval('NO_ASPCAP_RESULT') : continue

            # append field in case we have duplicates
            stars.append(star['APOGEE_ID']+'__'+star['FIELD'])

      
        # write FERRE files and run FERRE
        if clobber or not os.path.exists(out+'.spm') or \
               (os.path.exists(out+'.spm') and len(open(out+'.spm').readlines()) < len(stars) ) :

            # load up start parameters, fluxes, and errors
            for star,spec in zip(aspcapfield[gd],aspcapspec[gd]) :
                if star['ASPCAPFLAG'] & aspcapmask.getval('NO_ASPCAP_RESULT') : continue
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
                    normflux=spec['OBS']/spec['NORM']
                    normerr=spec['ERR']/spec['NORM']
                    bd=np.where(spec['ERR']>0.99e10)[0]
                    if len(bd) > 0 : 
                        normflux[bd] = 0.0001
                        normerr[bd] = 1.e10
                else :
                    normflux=spec['OBS']
                    normerr=spec['ERR']
                flux.append(normflux)
                err.append(normerr)
                if mult :
                    imult=0
                    for dlogg in np.linspace(-0.2,0.2,3) :
                        for dam in np.linspace(-0.2,0.2,3) :
                            for dcm in np.linspace(-0.2,0.2,3) :
                                stars.append(star['APOGEE_ID']+'.{:d}'.format(imult))
                                cent=copy.copy(star['FPARAM'])
                                cent[1]+=dlogg
                                cent[6]+=dam
                                cent[4]+=dcm
                                inpars.append(list(cent))
                                flux.append(normflux)
                                err.append(normerr)
                                imult+=1
            #exclude dimensions if coarse run or fixed parameters specified
            exclude=[]
            if coarse : exclude.extend(grid['fix_coarse'])
            if fix is not None : exclude.extend(fix)
            if len(exclude) > 0 :    
                indv=[]
                for i,par in enumerate(libhead0['LABEL']) :
                    if par.decode() not in exclude : indv.append(i+1)
            else : indv=None

            # write ferre files and run ferre
            ferre.writeipf(out,os.path.dirname(out)+'/'+libfile,stars,param=np.array(inpars))
            ferre.writespec(out+'.obs',flux)
            ferre.writespec(out+'.err',err)
            ferre.writenml(out+'.nml',os.path.basename(out),libhead0,init=0,indv=indv,f_sort=0,
                       algor=grid['algor'],ncpus=plan['ncpus'],
                       obscont=grid['obscont'],rejectcont=grid['rejectcont'],
                       renorm=abs(grid['renorm']),
                       filterfile=os.environ['APOGEE_DIR']+'/data/windows/'+grid['mask'])
            fout=open(out+'.stdout','w')
            ferr=open(out+'.stderr','w')
            print('running ferre.x: ', os.path.basename(out)+'.nml')
            start = time.time()
            ret = subprocess.call(['ferre.x',os.path.basename(out)+'.nml'],shell=False,
                            cwd=os.path.dirname(out),stdout=fout,stderr=ferr)
            print('elapsed: ',ret, time.time()-start)
            fout.close()
            ferr.close()

        # read FERRE output
        print('reading FERRE output...')
        param,spec,wave=ferre.read(out,os.path.dirname(out)+'/'+libfile)
        # fill in locked parameters
        fill_plock(param,grid['PLOCK'])

        # load into apcapField
        print('loading FERRE output...')
        for istar,star in enumerate(aspcapfield[gd]) :
            i = np.where(param['APOGEE_ID'] == (star['APOGEE_ID']+'__'+star['FIELD']).encode())[0]
            if len(i) == 0 : continue
            aspcapfield['FPARAM_GRID'][gd[istar],igrid,:] = param['FPARAM'][i]
            aspcapfield['CHI2_GRID'][gd[istar],igrid] = param['ASPCAP_CHI2'][i]
            # if this is the lowest CHI2, store in FPARAM and ASPCAP_CHI2
            # penalize GK grid in favor of finer M grid
            chi2 = param['ASPCAP_CHI2'][i]
            if ('GK' in grid['name']) and (param['FPARAM'][i,0] < 3900) : 
                print('penalizing GK')
                chi2 *= 10
            #penalize stars within 1 step of grid edge in TEFF or LOGG
            if param['PARAMFLAG'][i,0] & (parammask.getval('GRIDEDGE_WARN')|parammask.getval('GRIDEDGE_BAD')) > 0 : chi2*=5
            if param['PARAMFLAG'][i,1] & (parammask.getval('GRIDEDGE_WARN')|parammask.getval('GRIDEDGE_BAD')) > 0 : chi2*=5

            # load star for this grid if its the best fit
            if chi2 < aspcapfield['ASPCAP_CHI2'][gd[istar]] or aspcapfield['ASPCAP_CHI2'][gd[istar]]<=0 :
                aspcapfield['ASPCAP_GRID'][gd[istar]] = grid['name']
                aspcapfield['ASPCAP_CHI2'][gd[istar]] = chi2
                aspcapfield['FPARAM'][gd[istar]] = param['FPARAM'][i]
                aspcapfield['FPARAM_COV'][gd[istar]] = param['FPARAM_COV'][i]
                aspcapfield['PARAMFLAG'][gd[istar]] = param['PARAMFLAG'][i]
                aspcapfield['ASPCAPFLAG'][gd[istar]] = param['ASPCAPFLAG'][i]
                aspcapspec['SPEC'][gd[istar]] = spec['frd'][i]
                aspcapspec['SPEC_BESTFIT'][gd[istar]] =  spec['mdl'][i]
                # continuum correct
                width=151
                p1=0
                for ichip in range(len(libhead)) :
                    npix=libhead[ichip]['NPIX']
                    obs=spec['frd'][i,p1:p1+npix].flatten()
                    # replace bad pixels with median filtered value
                    med=median_filter(obs,[5*width],mode='nearest')
                    bd=np.where(obs < 0.01)[0]
                    obs[bd] = med[bd]
                    # populate SPEC_ERR with FERRE-normalized error
                    err= spec['err'][i,p1:p1+npix]*spec['frd'][i,p1:p1+npix]/obs
                    bd=np.where(~np.isfinite(err))[0]
                    err[bd] = 1.e10
                    aspcapspec['SPEC_ERR'][gd[istar]][p1:p1+npix] = err
                  
                    # get ratio of observed / model, make sure edge is reasonable number
                    mdl=spec['mdl'][i,p1:p1+npix].flatten()
                    ratio=obs/mdl
                    corr=median_filter(ratio,[width],mode='nearest')
                    bd=np.where(~np.isfinite(corr) | (corr < 0.1) | (corr > 10.) )[0]
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
                    aspcapspec['MASK'][gd[istar]][bd] |=  pixelmask.getval('FERRE_MASK')
                except: pdb.set_trace()
            if mult :
                for imult in range(27) :
                    i = np.where(param['APOGEE_ID'] == '{:s}.{:d}'.format(star['APOGEE_ID'],imult).encode())[0]
                    if len(i) == 0 : continue
                    aspcapfield['FPARAM_MULT'][gd[istar],imult,:] = param['FPARAM'][i]
                    aspcapspec['SPEC_BESTFIT_MULT'][gd[istar],imult,:] =  spec['mdl'][i]

    # populate ASPCAPFLAGS
    for i,star in enumerate(aspcapfield) :
        aspcapfield['ASPCAPFLAGS'][i] = aspcapmask.getname(star['ASPCAPFLAG'])

    # Results in astropy tables
    aspcapfield=Table(aspcapfield)
    aspcapspec=Table(aspcapspec)
    wave=np.hstack(gridWave())
    aspcapkey=np.empty(1,dtype=[('WAVE','f4',len(wave)),
                          ('PARAM_SYMBOL','S16',len(params()[0])),
                          ('ELEM_SYMBOL','S5',len(elems()[0])),
                          ('ELEMTOH','i4',len(elems()[0])),
                          ('ELEM_VALUE','S12',len(elems()[2])),
                          ('GRIDS','S5',len(config['grids']))])
    aspcapkey['WAVE'] = wave
    aspcapkey['PARAM_SYMBOL'] = params()[1]
    aspcapkey['ELEM_SYMBOL'] = elems()[0]
    aspcapkey['ELEMTOH'] = elems()[1]
    aspcapkey['ELEM_VALUE'] = elems()[2]
    grids=[]
    for grid in config['grids'] : grids.append(grid['name'])
    aspcapkey['GRIDS'] = grids
    aspcapkey=Table(aspcapkey)

    if write :
        if mult : suffix='_mult'
        print('writing master aspcapField...')
        writefiles(bundleload,outfield,aspcapfield,aspcapspec,aspcapkey,suffix=suffix,aspcapstar=False)
        if isinstance(field,list) :
            for f in field :
                print('writing aspcapField...',f)
                gd = np.where(aspcapfield['FIELD'] == f)[0]
                if len(gd) > 0 : writefiles(load,f,aspcapfield[gd],aspcapspec[gd],aspcapkey,suffix=suffix)

    if html :
        # create output HTML page
        if isinstance(field,list) :
            for f in field :
                gd = np.where(aspcapfield['FIELD'] == f)[0]
                if len(gd) > 0 : mkhtml(load,f,suffix='')
        else :
            mkhtml(load,outfield,suffix='')

    return aspcapfield,aspcapspec,aspcapkey

def fit_elems(planfile,aspcapdata=None,clobber=False,nobj=None,write=True,calib=False,renorm=True,suffix='',html=True) :
    """ run ASPCAP on a field for elemental abundances
    """

    # read configuration file
    plan=yaml.safe_load(open(planfile,'r'))
    if plan['apogee_ver'] != os.environ['APOGEE_VER'] :
        print('apogee_ver {:s} does not match running version {:s}'.format(plan['apogee_ver'],os.environ['APOGEE_VER']))
        pdb.set_trace()

    apred=plan['apred_vers']
    aspcap_vers=plan['aspcap_vers']
    aspcap_config=plan['aspcap_config']
    instrument=plan['instrument']
    telescope=plan['telescope']
    field=plan['field']
    if 'outfield' not in plan.keys() : outfield = field
    else : outfield = plan['outfield']

    # get aspcapField if not provided
    load = apload.ApLoad(apred=apred,aspcap=aspcap_vers,telescope=telescope)
    if isinstance(field,list) : bundleload = apload.ApLoad(apred=apred,aspcap=aspcap_vers,telescope='bundle_'+telescope)
    else : bundleload = load
    outfile=bundleload.filename('aspcapField',field=outfield)
    outdir=os.path.dirname(outfile)

    if aspcapdata is None : 
        aspcapfield=Table(fits.open(outfile)[1].data)
        aspcapspec=Table(fits.open(outfile)[2].data)
        aspcapkey=Table(fits.open(outfile)[3].data)
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
        outspec=outdir+'/ferre/spectra/'+grid['name']+'-'+outfield
        try: os.makedirs(os.path.dirname(outspec))
        except: pass

        # get FERRE library information for this grid
        link(os.environ['APOGEE_SPECLIB']+'/synth/',os.path.dirname(outspec)+'/../lib_'+grid['name'])
        libfile = 'lib_'+grid['name']+'/'+grid['lib']+'.hdr'
        libhead0,libhead = ferre.rdlibhead(outdir+'/ferre/'+libfile)
        libhead0['FILE'] = libfile

        # get stars for which best fit was in this grid
        gd = np.where(aspcapfield['ASPCAP_GRID'] == grid['name'])[0]
        print('Grid: ',grid['name'],len(gd))
        if len(gd) == 0 : continue

        # loop over stars and accumulate input for FERRE
        inpars=[]
        flux=[]
        err=[]
        stars=[]
        for star,spec in zip(aspcapfield[gd],aspcapspec[gd]) :
            stars.append(star['APOGEE_ID']+'__'+star['FIELD'])
            inpars.append(star[useparam])
            print(star['APOGEE_ID'])
            if renorm :
                normflux=spec['OBS']/spec['NORM']
                normerr=spec['ERR']/spec['NORM']
                bd=np.where(spec['ERR']>0.99e10)[0]
                if len(bd) > 0 : 
                    normflux[bd] = 0.0001
                    normerr[bd] = 1.e10
            else :
                normflux=spec['OBS']
                normerr=spec['ERR']
            flux.append(normflux)
            err.append(normerr)

        if clobber or not os.path.exists(outspec+'.obs') :

            ferre.writespec(outspec+'.obs',flux)
            ferre.writespec(outspec+'.err',err)

        # loop over elements
        fp=open(outdir+'/ferre/'+grid['name']+'.nmlfiles','w')
        fit = False
        for ielem,elem in enumerate(config['elems']) :
            print('Element: ', elem) 
            # set up output FERRE directory for this grid
            dirname='elem_'+elem['name']
            if calib: dirname=dirname+'_PARAM'
            out=outdir+'/ferre/'+dirname+'/'+elem['name']+'-'+grid['name']+'-'+outfield
            try: os.makedirs(os.path.dirname(out))
            except: pass

            # write FERRE files
            if clobber or not os.path.exists(out+'.spm') or \
               (os.path.exists(out+'.spm') and len(open(out+'.spm').readlines()) < len(stars) ) :
                link(os.environ['APOGEE_SPECLIB']+'/synth/',os.path.dirname(out)+'/lib')
                filterfile = elem['name']+'.mask'
                shutil.copyfile(os.environ['APOGEE_DIR']+'/data/windows/'+grid['windows']+'/'+filterfile,
                                     os.path.dirname(out)+'/'+elem['name']+'.mask')

                # links for obs, and err files
                for ext in ['.obs','.err'] :
                    try: os.remove(out+ext)
                    except: pass
                    link('../spectra/'+os.path.basename(outspec)+ext,out+ext)
                ferre.writeipf(out,outdir+'/ferre/'+libfile,stars,param=np.array(inpars))
                index = np.where(libhead0['LABEL'] == elem['griddim'].encode())[0][0]+1
                if elem['griddim'] == 'METALS' : 
                    ttie=[]
                    for dim in ['C','N','O Mg Si S Ca Ti'] :
                        j=np.where(libhead0['LABEL'] == dim.encode())[0]
                        if len(j) > 0 : ttie.append(j[0]+1)
                        else : ttie.append(-1)
                else : ttie=[-1,-1,-1]
                ferre.writenml(out+'.nml',dirname+'/'+os.path.basename(out),libhead0,init=0,f_sort=0,
                           nov=1,indv=[index],ttie=ttie,
                           algor=grid['algor'],ncpus=plan['ncpus'],
                           obscont=grid['obscont'],rejectcont=grid['rejectcont'],
                           renorm=abs(grid['renorm']),
                           filterfile='elem_'+elem['name']+'/'+filterfile)
                fp.write(dirname+'/'+os.path.basename(out)+'.nml\n')
                fit = True

        fp.close()
        # run FERRE for all elements for this grid
        if fit :
            fout=open(outdir+'/ferre/'+grid['name']+'.stdout','w')
            ferr=open(outdir+'/ferre/'+grid['name']+'.stderr','w')
            print('running ferre.x: -l', grid['name']+'.nmlfiles')
            start = time.time()
            ret = subprocess.call(['ferre.x','-l',grid['name']+'.nmlfiles'],shell=False,
                            cwd=outdir+'/ferre',stdout=fout,stderr=ferr)
            print('return, elapsed: ',ret, time.time()-start)
            fout.close()
            ferr.close()

        for ielem,elem in enumerate(config['elems']) :
            dirname='elem_'+elem['name']
            if calib: dirname=dirname+'_PARAM'
            out=outdir+'/ferre/'+dirname+'/'+elem['name']+'-'+grid['name']+'-'+outfield
            # read FERRE output
            param,spec,wave=ferre.read(out,outdir+'/ferre/'+libfile)
            # fill in locked parameters
            fill_plock(param,grid['PLOCK'])
            # location for this element in FELEM array
            jelem=np.where(elems()[0] == elem['name'])[0]

            # load into aspcapField
            for istar,star in enumerate(aspcapfield[gd]) :
                i = np.where(param['APOGEE_ID'] == (star['APOGEE_ID']+'__'+star['FIELD']).encode())[0]
                if len(i) == 0 : continue
                index = np.where(params()[0] == elem['griddim'])[0]
                try:
                    aspcapfield['FELEM'][gd[istar],jelem] = param['FPARAM'][i,index]
                    if param['FPARAM_COV'][i,index,index] > 0 :
                        aspcapfield['FELEM_ERR'][gd[istar],jelem] = np.sqrt(param['FPARAM_COV'][i,index,index])
                    aspcapfield['ELEM_CHI2'][gd[istar],jelem] = param['ASPCAP_CHI2'][i]
                    aspcapfield['ELEMFLAG'][gd[istar],jelem] = param['PARAMFLAG'][i,index]
                except: pdb.set_trace()

    # Results into an HDUList 

    if write :
        writefiles(bundleload,outfield,aspcapfield,aspcapspec,aspcapkey,suffix=suffix,aspcapstar=False)
        if isinstance(field,list) :
            for f in field :
                gd = np.where(aspcapfield['FIELD'] == f)[0]
                if len(gd) > 0 : writefiles(load,f,aspcapfield[gd],aspcapspec[gd],aspcapkey,suffix=suffix)

    if html :
        # create output HTML page, but if we are bundling, only for the individual fields
        if isinstance(field,list) :
            for f in field :
                gd = np.where(aspcapfield['FIELD'] == f)[0]
                if len(gd) > 0 : mkhtml(load,f,suffix='')
        else :
            mkhtml(load,outfield,suffix='')

    return aspcapfield,aspcapspec,aspcapkey

def writefiles(load,field,aspcapfield,aspcapspec,aspcapkey,suffix='',aspcapstar=True) :

    """ Write aspcapField and aspcapStar files
    """

    hdulist=fits.HDUList()
    hdulist.append(fits.table_to_hdu(aspcapfield))
    hdulist.append(fits.table_to_hdu(aspcapspec))
    hdulist.append(fits.table_to_hdu(aspcapkey))

    #output aspcapField
    outfield=load.filename('aspcapField',field=field)
    try: os.makedirs(os.path.dirname(outfield))
    except: pass
    outfile=os.path.dirname(outfield)+'/'+os.path.splitext(os.path.basename(outfield))[0]+suffix+'.fits'
    hdulist[0].header['VERSION'] = (os.environ['APOGEE_VER'],'APOGEE software version APOGEE_VER')
    hdulist.writeto(outfile,overwrite=True)

    #output aspcapStar
    aspcapmask=bitmask.AspcapBitMask()
    if aspcapstar :
        for star,spec in zip(aspcapfield,aspcapspec) :
          # don't write aspcapStar if no ASPCAP result
          if not (star['ASPCAPFLAG'] & aspcapmask.getval('NO_ASPCAP_RESULT')) :
            outfile=load.filename('aspcapStar',field=field,obj=star['APOGEE_ID'])
            hdulist=fits.HDUList()
            hdulist.append(fits.PrimaryHDU())
            hdu=fits.ImageHDU(aspcap2apStar(spec['SPEC']))
            add_header(hdu)
            hdulist.append(hdu)
            hdu=fits.ImageHDU(aspcap2apStar(spec['SPEC_ERR']))
            add_header(hdu)
            hdulist.append(hdu)
            hdu=fits.ImageHDU(aspcap2apStar(spec['SPEC_BESTFIT']))
            add_header(hdu)
            hdulist.append(hdu)
            hdulist.append(fits.table_to_hdu(Table(star)))
            hdulist[0].header['VERSION'] = (os.environ['APOGEE_VER'],'APOGEE software version APOGEE_VER')
            hdulist.writeto(outfile,overwrite=True)

    return

def add_header(hdu) :
    hdu.header['CRVAL1']=logw0
    hdu.header['CDELT1']=dlogw
    hdu.header['CRPIX1']=1
    hdu.header['CTYPE1'] = 'LOG-LINEAR'
    hdu.header['DC-FLAG'] = 1

def mkhtml(load,field,suffix='') :
    """ Create ASPCAP field web page and plots
    """

    matplotlib.use('Agg')

    infile=load.filename('aspcapField',field=field)
    infile=os.path.dirname(infile)+'/'+os.path.splitext(os.path.basename(infile))[0]+suffix+'.fits'
    a=fits.open(infile)
    aspcapfield=a[1].data
    aspcapspec=a[2].data

    outdir=os.path.dirname(infile)
    plotdir=outdir+'/plots'+suffix+'/'
    reldir='plots'+suffix+'/'
    try: os.makedirs(plotdir)
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
    fig.savefig(plotdir+field+'_hr.png')
    plt.close()
    fp.write('<TD><A HREF={:s}/{:s}_hr.png><IMG SRC={:s}/{:s}_hr.png></A>\n'.format(reldir,field,reldir,field))

    yr=[-0.3,0.75]
    fig,ax=plots.multi(1,1)
    plots.plotc(ax,aspcapfield['FPARAM'][:,3],aspcapfield['FPARAM'][:,6],aspcapfield['FPARAM'][:,0],
                zr=[3000,8000],yr=yr,xr=[-2.5,0.5],size=10,colorbar=True,xt='[M/H]',yt='[alpha/M]',zt='Teff')
    fig.savefig(plotdir+field+'_alpha.png')
    plt.close()
    fp.write('<TD><A HREF={:s}/{:s}_alpha.png><IMG SRC={:s}/{:s}_alpha.png></A>\n'.format(reldir,field,reldir,field))

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
    fig.savefig(plotdir+field+'_cnelem.png')
    plt.close()
    fp.write('<TD><A HREF={:s}/{:s}_cnelem.png><IMG SRC={:s}/{:s}_cnelem.png></A>\n'.format(reldir,field,reldir,field))

    fig,ax=plots.multi(1,6,hspace=0.001)
    yr=[-0.3,0.75]
    for iel,el in enumerate(['O','Mg','Si','S','Ca','Ti']) :
        jel = np.where(elems()[0] == el)[0][0]
        yt='['+el+'/M]'
        plots.plotc(ax[iel],aspcapfield['FPARAM'][:,3],aspcapfield['FELEM'][:,jel],aspcapfield['FPARAM'][:,0],
                zr=[3000,8000],yr=yr,xr=[-2.5,0.5],size=10,colorbar=True,xt='[M/H]',yt=yt,zt='Teff',label=[0.9,0.9,el])
    fig.savefig(plotdir+field+'_alphaelem.png')
    plt.close()
    fp.write('<TD><A HREF={:s}/{:s}_alphaelem.png><IMG SRC={:s}/{:s}_alphaelem.png></A>\n'.format(reldir,field,reldir,field))

    fig,ax=plots.multi(1,4,hspace=0.001)
    for iel,el in enumerate(['Na','Al','P','K']) :
        jel = np.where(elems()[0] == el)[0][0]
        yt='['+el+'/M]'
        try :
          plots.plotc(ax[iel],aspcapfield['FPARAM'][:,3],aspcapfield['FELEM'][:,jel]-aspcapfield['FPARAM'][:,3],aspcapfield['FPARAM'][:,0],
                zr=[3000,8000],yr=yr,xr=[-2.5,0.5],size=10,colorbar=True,xt='[M/H]',yt=yt,zt='Teff',label=[0.9,0.9,el])
        except: pdb.set_trace()
    fig.savefig(plotdir+field+'_oddzelem.png')
    plt.close()
    fp.write('<TD><A HREF={:s}/{:s}_oddzelem.png><IMG SRC={:s}/{:s}_oddzelem.png></A>\n'.format(reldir,field,reldir,field))

    fig,ax=plots.multi(1,7,hspace=0.001)
    for iel,el in enumerate(['V','Cr','Mn','Fe','Co','Ni','Cu'] ) :
        jel = np.where(elems()[0] == el)[0][0]
        yt='['+el+'/M]'
        plots.plotc(ax[iel],aspcapfield['FPARAM'][:,3],aspcapfield['FELEM'][:,jel]-aspcapfield['FPARAM'][:,3],aspcapfield['FPARAM'][:,0],
                zr=[3000,8000],yr=yr,xr=[-2.5,0.5],size=10,colorbar=True,xt='[M/H]',yt=yt,zt='Teff',label=[0.9,0.9,el])
    fig.savefig(plotdir+field+'_feelem.png')
    plt.close()
    fp.write('<TD><A HREF={:s}/{:s}_feelem.png><IMG SRC={:s}/{:s}_feelem.png></A>\n'.format(reldir,field,reldir,field))

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
        fig.savefig(plotdir+star+'.png')
        plt.close()

        fp.write('<TR><TD>{:s}<BR>\n'.format(star))
        fp.write('H  = {:7.2f}<br>'.format(aspcapfield['H'][istar]))
        fp.write('SNR  = {:7.2f}<br>'.format(aspcapfield['SNR'][istar]))
        fp.write('<FONT COLOR=blue> TARGFLAGS </FONT>: {:s}<BR>\n'.format(aspcapfield['TARGFLAGS'][istar]))
        fp.write('<FONT COLOR=blue> STARFLAGS </FONT>: {:s}<BR>\n'.format(aspcapfield['STARFLAGS'][istar]))
        fp.write('<FONT COLOR=blue> ASPCAPFLAGS </FONT>: {:s}<BR>\n'.format(aspcapfield['ASPCAPFLAGS'][istar]))
        fp.write('<TD>{:s}\n'.format(aspcapfield['ASPCAP_GRID'][istar]))
        fp.write('<TD>{:6.1f}\n'.format(aspcapfield['ASPCAP_CHI2'][istar]))
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
        fp.write('<TD><A HREF={:s}/{:s}.png><IMG SRC={:s}/{:s}.png></A>\n'.format(reldir,star,reldir,star))

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

def comp(a,b,terange=[4500,5000],loggrange=[2.5,3.5],mhrange=[-2.5,1.0],amrange=[-1,1],spec=True) :
    """ Compare two versions in region on Kiel diagram
    """
    matplotlib.use('Agg')
    gd=np.where( (a[1].data['FPARAM'][:,0]>=terange[0]) &
                 (a[1].data['FPARAM'][:,0]<terange[1]) &
                 (a[1].data['FPARAM'][:,1]>=loggrange[0]) &
                 (a[1].data['FPARAM'][:,1]<loggrange[1]) &
                 (a[1].data['FPARAM'][:,6]>=amrange[0]) &
                 (a[1].data['FPARAM'][:,6]<amrange[1]) &
                 (a[1].data['FPARAM'][:,3]>=mhrange[0]) &
                 (a[1].data['FPARAM'][:,3]<mhrange[1]) )[0]
    print('{:d} stars'.format(len(gd)))

    fp=html.head('comp.html')

    grid=[]
    for iplot in range(3) :
        if iplot == 0 :
            a_z = a[1].data['FPARAM'][:,3]
            b_z = b[1].data['FPARAM'][:,3]
            zr=[-2.5,0.5]
            zt='[M/H]'
        elif iplot == 1 :
            a_z = a[1].data['FPARAM'][:,6]
            b_z = b[1].data['FPARAM'][:,6]
            zr=[-0.25,0.5]
            zt='[alpha/M]'
        elif iplot == 2 :
            a_z = a[1].data['FPARAM'][:,4]-a[1].data['FPARAM'][:,5]
            b_z = b[1].data['FPARAM'][:,6]-b[1].data['FPARAM'][:,5]
            zr=[-0.5,1.0]
            zt='[C/N]'
        fig,ax=plots.multi(2,1,figsize=(4,3),wspace=0.001)
        plots.plotc(ax[0],a[1].data['FPARAM'][:,0],a[1].data['FPARAM'][:,1],a_z,
                    xr=[8000,3000],yr=[6,-1],zr=zr,size=1,colorbar=False,zt=zt,yt='log g',xt='Teff')
        plots.plotc(ax[0],a[1].data['FPARAM'][gd,0],a[1].data['FPARAM'][gd,1],a_z[gd],
                    xr=[8000,3000],yr=[6,-1],zr=zr,size=10,colorbar=False,zt=zt,yt='log g',xt='Teff')
        plots.plotc(ax[1],b[1].data['FPARAM'][:,0],b[1].data['FPARAM'][:,1],b_z,
                    xr=[8000,3000],yr=[6,-1],zr=zr,size=1,colorbar=False,zt=zt,xt='Teff')
        im=plots.plotc(ax[1],b[1].data['FPARAM'][gd,0],b[1].data['FPARAM'][gd,1],b_z[gd],
                    xr=[8000,3000],yr=[6,-1],zr=zr,size=10,colorbar=False,zt=zt,xt='Teff')
        cbar = fig.colorbar(im, ax=ax.ravel().tolist(), shrink=0.95)
        cbar.set_label(zt)
        fig.savefig('kiel_{:d}.png'.format(iplot))
        plt.close()
    fp.write(html.table([['kiel_0.png','kiel_1.png','kiel_2.png']]))

    fig,ax=plots.multi(1,2,figsize=(4,3),hspace=0.001)
    zr=[3000,8000]
    zt='Teff'
    plots.plotc(ax[0],a[1].data['FPARAM'][:,3],a[1].data['FPARAM'][:,6],a[1].data['FPARAM'][:,0],
                xr=[-2.5,1],yr=[-0.5,0.75],zr=zr,size=1,colorbar=False,zt=zt,yt='[alpha/M]',xt='[M/H]')
    plots.plotc(ax[0],a[1].data['FPARAM'][gd,3],a[1].data['FPARAM'][gd,6],a[1].data['FPARAM'][gd,0],
                xr=[-2.5,1],yr=[-0.5,0.75],zr=zr,size=10,colorbar=False,zt=zt,yt='[alpha/M]',xt='[M/H]')
    plots.plotc(ax[1],b[1].data['FPARAM'][:,3],b[1].data['FPARAM'][:,6],b[1].data['FPARAM'][:,0],
                xr=[-2.5,1],yr=[-0.5,0.75],zr=zr,size=1,colorbar=False,zt=zt,yt='[alpha/M]',xt='[M/H]')
    im=plots.plotc(ax[1],b[1].data['FPARAM'][gd,3],b[1].data['FPARAM'][gd,6],b[1].data['FPARAM'][gd,0],
                xr=[-2.5,1],yr=[-0.5,0.75],zr=zr,size=10,colorbar=False,zt=zt,yt='[alpha/M]',xt='[M/H]')
    cbar = fig.colorbar(im, ax=ax.ravel().tolist(), shrink=0.95)
    cbar.set_label(zt)
    fig.savefig('alpha.png')
    plt.close()
    fig,ax=plots.multi(1,2,figsize=(4,3),hspace=0.001)
    zr=[-2.5,0.5]
    zt='[M/H]'
    plots.plotc(ax[0],a[1].data['FPARAM'][:,1],a[1].data['FPARAM'][:,4]-a[1].data['FPARAM'][:,5],a[1].data['FPARAM'][:,3],
                xr=[5,0],yr=[-0.5,1.],zr=zr,size=1,colorbar=False,zt=zt,yt='[C/N]',xt='log g')
    plots.plotc(ax[0],a[1].data['FPARAM'][gd,1],a[1].data['FPARAM'][gd,4]-a[1].data['FPARAM'][gd,5],a[1].data['FPARAM'][gd,3],
                xr=[5,0],yr=[-0.5,1.],zr=zr,size=10,colorbar=False,zt=zt,yt='[C/N]',xt='log g')
    plots.plotc(ax[1],b[1].data['FPARAM'][:,1],b[1].data['FPARAM'][:,4]-b[1].data['FPARAM'][:,5],b[1].data['FPARAM'][:,3],
                xr=[5,0],yr=[-0.5,1.],zr=zr,size=1,colorbar=False,zt=zt,yt='[C/N]',xt='log g')
    im=plots.plotc(ax[1],b[1].data['FPARAM'][gd,1],b[1].data['FPARAM'][gd,4]-b[1].data['FPARAM'][gd,5],b[1].data['FPARAM'][gd,3],
                xr=[5,0],yr=[-0.5,1.],zr=zr,size=10,colorbar=False,zt=zt,yt='[C/N]',xt='log g')
    cbar = fig.colorbar(im, ax=ax.ravel().tolist(), shrink=0.95)
    cbar.set_label(zt)
    fig.savefig('cn.png')
    plt.close()
    fp.write(html.table([['alpha.png','cn.png']]))

    fig,ax=plots.multi(1,7,hspace=0.001)
    for i in range(7) :
        if i == 0 : yr=[-100,100]
        else : yr=[-0.25,0.25]
        xt = 'Delta log g'
        plots.plotc(ax[i],b[1].data['FPARAM'][:,1]-a[1].data['FPARAM'][:,1],b[1].data['FPARAM'][:,i]-a[1].data['FPARAM'][:,i],
                    a[1].data['FPARAM'][:,3],zr=[-2.5,0.5],yt=a[3].data['PARAM_SYMBOL'][0][i],size=1,xr=[-0.5,0.5],yr=yr,xt=xt)
        plots.plotc(ax[i],b[1].data['FPARAM'][gd,1]-a[1].data['FPARAM'][gd,1],b[1].data['FPARAM'][gd,i]-a[1].data['FPARAM'][gd,i],
                    a[1].data['FPARAM'][gd,3],zr=[-2.5,0.5],yt=a[3].data['PARAM_SYMBOL'][0][i],xr=[-0.5,0.5],yr=yr)
    fig.savefig('dparam_dlogg.png')
    plt.close()

    fig,ax=plots.multi(1,7,hspace=0.001)
    for i in range(7) :
        if i == 0 : yr=[-100,100]
        else : yr=[-0.25,0.25]
        xt = 'Delta Teff'
        plots.plotc(ax[i],b[1].data['FPARAM'][:,0]-a[1].data['FPARAM'][:,0],b[1].data['FPARAM'][:,i]-a[1].data['FPARAM'][:,i],
                    a[1].data['FPARAM'][:,3],zr=[-2.5,0.5],yt=a[3].data['PARAM_SYMBOL'][0][i],size=1,xr=[-200,200],yr=yr,xt=xt)
        plots.plotc(ax[i],b[1].data['FPARAM'][gd,0]-a[1].data['FPARAM'][gd,0],b[1].data['FPARAM'][gd,i]-a[1].data['FPARAM'][gd,i],
                    a[1].data['FPARAM'][gd,3],zr=[-2.5,0.5],yt=a[3].data['PARAM_SYMBOL'][0][i],xr=[-200,200],yr=yr)
    fig.savefig('dparam_dteff.png')
    plt.close()
    fp.write(html.table([['dparam_dlogg.png','dparam_dteff.png']]))

    for iplot in range(2) :
        fig,ax=plots.multi(3,6,hspace=0.001,wspace=0.001,figsize=(12,8))
        yr=[-0.14,0.14]
        x=a[1].data['FPARAM'][gd,3]
        xr=[8000,3000]
        if iplot == 0 :
            x=b[1].data['FPARAM'][:,1]-a[1].data['FPARAM'][:,1]
            xr=[-0.5,0.5]
            xt = 'Delta log g'
        else :
            x=b[1].data['FPARAM'][:,0]-a[1].data['FPARAM'][:,0]
            xr=[-200,200]
            xt = 'Delta Teff'
        for iel,el in enumerate(['O','Mg','Si','S','Ca','Ti']) :
            jel = np.where(elems()[0] == el)[0][0]
            print(jel)
            yt='['+el+'/M]'
            x_m_a=a[1].data['FELEM'][:,jel]
            x_m_b=b[1].data['FELEM'][:,jel]
            plots.plotc(ax[iel,0],x,x_m_b-x_m_a,a[1].data['FPARAM'][:,1],
                    xr=xr,yr=yr,zr=[0,5],size=1,colorbar=False,zt='[M/H]',xt=xt,label=[0.9,0.9,el])
            plots.plotc(ax[iel,0],x[gd],x_m_b[gd]-x_m_a[gd],a[1].data['FPARAM'][gd,1],
                    xr=xr,yr=yr,zr=[0,5],size=10,colorbar=False,zt='[M/H]',xt=xt,label=[0.9,0.9,el])
        for iel,el in enumerate(['Na','Al','P','K']) :
            jel = np.where(elems()[0] == el)[0][0]
            print(jel)
            yt='['+el+'/M]'
            x_m_a=a[1].data['FELEM'][:,jel]-a[1].data['FPARAM'][:,3]
            x_m_b=b[1].data['FELEM'][:,jel]-b[1].data['FPARAM'][:,3]
            plots.plotc(ax[iel,1],x,x_m_b-x_m_a,a[1].data['FPARAM'][:,1],
                    xr=xr,yr=yr,zr=[0,5],size=1,colorbar=False,zt='[M/H]',xt=xt,label=[0.9,0.9,el])
            plots.plotc(ax[iel,1],x[gd],x_m_b[gd]-x_m_a[gd],a[1].data['FPARAM'][gd,1],
                    xr=xr,yr=yr,zr=[0,5],size=10,colorbar=False,zt='[M/H]',xt=xt,label=[0.9,0.9,el])
        for iel,el in enumerate(['V','Cr','Mn','Fe','Co','Ni'] ) :
            jel = np.where(elems()[0] == el)[0][0]
            print(jel)
            yt='['+el+'/M]'
            x_m_a=a[1].data['FELEM'][:,jel]-a[1].data['FPARAM'][:,3]
            x_m_b=b[1].data['FELEM'][:,jel]-b[1].data['FPARAM'][:,3]
            plots.plotc(ax[iel,2],x,x_m_b-x_m_a,a[1].data['FPARAM'][:,1],
                    xr=xr,yr=yr,zr=[0,5],size=1,colorbar=False,zt='[M/H]',xt=xt,label=[0.9,0.9,el])
            plots.plotc(ax[iel,2],x[gd],x_m_b[gd]-x_m_a[gd],a[1].data['FPARAM'][gd,1],
                    xr=xr,yr=yr,zr=[0,5],size=10,colorbar=False,zt='[M/H]',xt=xt,label=[0.9,0.9,el])
        fig.savefig('elem_{:d}.png'.format(iplot))
        plt.close()
    fp.write(html.table([['elem_0.png','elem_1.png']]))

    if spec :
        wave=np.hstack(gridWave())
        achi,acum=chi2(a[2].data)
        bchi,bcum=chi2(b[2].data)

        fig,ax=plots.multi(1,4,hspace=0.001,sharex=True) 
        plots.plotl(ax[0],wave,np.nanmedian(bcum-acum,axis=0))
        plots.plotl(ax[1],wave,np.nanmedian(b[2].data['SPEC_BESTFIT'][gd]-a[2].data['SPEC_BESTFIT'][gd],axis=0))
        for i in gd :
            plots.plotl(ax[2],wave,(bcum[i]-acum[i])/(bcum[i]-acum[i]).sum())
            #plots.plotl(ax[3],wave,a[2].data['SPEC'][i])
            #plots.plotl(ax[3],wave,a[2].data['SPEC_BESTFIT'][i])
            plots.plotl(ax[3],wave,b[2].data['SPEC_BESTFIT'][i]-a[2].data['SPEC_BESTFIT'][i])
        fig.savefig('spec.png')
        plt.close()
        fp.write(html.table([['spec.png']]))

    html.tail(fp)

def elemfrac(spec,el) :
    """ calculate fractional contribution of good pixels
    """
    print('elemfrac: ',el)
    mask=np.loadtxt(os.environ['APOGEE_DIR']+'/data/windows/dr17/'+el+'.filt')
    rat=[]
    pixmask=bitmask.PixelBitMask()
    for i in range(len(spec)) :
        gd=np.where((spec['MASK'][i] & (pixmask.badval()|pixmask.getval('SIG_SKYLINE'))) == 0)[0]
        rat.append(mask[gd].sum() / mask.sum())
        #err=np.median(spec['ERR'][i,:])
        #rat.append((mask/spec['ERR'][i,:]**2).sum()/(mask/spec['RAWERR']**2).sum())
    rat=np.array(rat)

    return rat
