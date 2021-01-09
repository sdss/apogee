import numpy as np
from esutil import htm
import astropy
from astropy.table import Table, Column, vstack
from astropy.io import fits
from apogee.utils import bitmask
from apogee.apred import apstar
from apogee.aspcap import aspcap, teff, logg, cal
from tools import match
import os
import pdb
import glob
import yaml

def all(planfile,dofix=False,suffix=None,allvisit=True,allplate=True, calsample=False) :
   """ Summary files after everything is run
   """
   plan=yaml.safe_load(open(planfile,'r'))
   apred_vers=plan['apred_vers']
   apstar_vers=plan['apstar_vers']
   aspcap_vers=plan['aspcap_vers']
   if suffix==None : suffix=plan['suffix']

   # allStar file
   apstar_dir=os.environ['APOGEE_REDUX']+'/'+apred_vers+'/'+apstar_vers+'/'
   aspcap_dir=os.environ['APOGEE_ASPCAP']+'/'+apred_vers+'/'+aspcap_vers+'/'
   print('Create allStar file')
   if calsample : skip=['Field-apo25m_','Field-lco25m_','Field-apo1m_','apo25m.','lco25m.']
   else : skip = None
   tab,dat=allStar(search=[aspcap_dir+'apo*/*/aspcapField-*.fits',aspcap_dir+'lco*/*/aspcapField-*.fits'],skip=skip,out=None,dofix=dofix)

   # allVisit file
   if allvisit :
       print('Create allVisit file')
       allvisit = apstar.allFieldVisit(search=[apstar_dir+'apo*/*/a?FieldVisits-*.fits',apstar_dir+'lco*/*/a?FieldVisits-*.fits'],
                  out=aspcap_dir+'allVisit-'+apred_vers+'-'+aspcap_vers+suffix+'.fits')

   # fix up issues in allStar file that require allVisit (SB2s)
   if dofix :
       print('fix allStar file')
       fix(tab,allvisit)
       # write it out
       hdulist=fits.HDUList()
       hdulist.append(fits.BinTableHDU(tab))
       hdulist.append(fits.BinTableHDU(dat))
       hdulist.append(fits.BinTableHDU(dat))
       hdulist.writeto(aspcap_dir+'allStar-'+apred_vers+'-'+aspcap_vers+suffix+'.fits',overwrite=True)
 
   # allPlate file
   if allplate :
       print('Create allPlate file')
       allplate = apstar.allPlate(allvisit,
                  out=aspcap_dir+'allPlate-'+apred_vers+'-'+aspcap_vers+suffix+'.fits')

def allStar(search=['apo*/*/aspcapField-*.fits','lco*/*/aspcapField-*.fits'],out='allStar.fits',
            skip=['Field-cal_','Field-apo25m_','Field-lco25m_','Field-apo1m_','apo25m.','lco25m.'], dofix=False) :
    '''
    Concatenate set of aspcapField files, and add named_tags, extratarg
    '''

    # search for input files
    if type(search) == str:
        search=[search]
    allfiles=[]
    for path in search :
        allfiles.extend(glob.glob(path))

    # read and append all of the individual field tables
    a=[]
    for file in allfiles :
        if skip is not None and doskip(file,skip) : continue
        print(file)
        dat=Table.read(file,hdu=1)
        try : dat.remove_column('FPARAM_COV_CLASS')
        except : pass
        a.append(dat)
    # stack them
    tab =vstack(a)
    del(a)
    tab.sort(['RA'])

    if dofix : fix(tab)

    # add named tags
    add_named_tags(tab)

    # add EXTRATARG, H_MIN, H_MAX, JKMIN, JKMAX
    add_extratarg(tab)

    # set ASPCAPFLAG bits and ASPCAPFLAGS
    allcal(tab)
    aspcapflag(tab)

    # construct HDUList
    hdulist=fits.HDUList()
    hdulist.append(fits.BinTableHDU(tab))
    dat=Table.read(file,hdu=3)
    hdulist.append(fits.BinTableHDU(dat))
    hdulist.append(fits.BinTableHDU(dat))
    # write out the file
    if out is not None:
        print('writing',out)
        hdulist[0].header['APOGEE_VER'] = os.environ['APOGEE_VER']
        hdulist.writeto(out,overwrite=True)

    return tab, dat

def doskip(file,skip) :
    for sk in skip : 
        if sk in file : return True
    return False

def add_named_tags(tab) :
    """ Add abundance named tags
    """
    if not isinstance(tab,astropy.table.table.Table) : tab=Table(tab)
    # named param flags
    for name in ['TEFF', 'LOGG', 'M_H', 'ALPHA_M'] :
        col = Column(np.full([len(tab)],np.nan),name=name,dtype=np.float32)
        try : tab.remove_column(name)
        except: pass
        tab.add_column(col)
        col = Column(np.full([len(tab)],np.nan),name=name+'_ERR',dtype=np.float32)
        try : tab.remove_column(name+'_ERR')
        except: pass
        tab.add_column(col)
    for name in ['VMICRO', 'VMACRO', 'VSINI', 'TEFF_SPEC', 'LOGG_SPEC'] :
        col = Column(np.full([len(tab)],np.nan),name=name,dtype=np.float32)
        try : tab.remove_column(name)
        except: pass
        tab.add_column(col)
    tab['TEFF'] = tab['PARAM'][:,0].astype(np.float32)
    tab['TEFF_ERR'] = np.sqrt(tab['PARAM_COV'][:,0,0]).astype(np.float32)
    tab['LOGG'] = tab['PARAM'][:,1].astype(np.float32)
    tab['LOGG_ERR'] = np.sqrt(tab['PARAM_COV'][:,1,1]).astype(np.float32)
    tab['M_H'] = tab['PARAM'][:,3].astype(np.float32)
    tab['M_H_ERR'] = np.sqrt(tab['PARAM_COV'][:,3,3]).astype(np.float32)
    tab['ALPHA_M'] = tab['PARAM'][:,6].astype(np.float32)
    tab['ALPHA_M_ERR'] = np.sqrt(tab['PARAM_COV'][:,6,6]).astype(np.float32)
    tab['TEFF_SPEC'] = tab['FPARAM'][:,0].astype(np.float32)
    tab['LOGG_SPEC'] = tab['FPARAM'][:,1].astype(np.float32)
    tab['VMICRO'] = 10.**tab['FPARAM'][:,2].astype(np.float32)
    dw=np.where(np.core.defchararray.find(tab['ASPCAP_GRID'],b'BA') |
                np.core.defchararray.find(tab['ASPCAP_GRID'],b'GKd') |
                np.core.defchararray.find(tab['ASPCAP_GRID'],b'Fd') |
                np.core.defchararray.find(tab['ASPCAP_GRID'],b'Md') ) [0]
    tab['VMACRO'][dw] = 0.
    tab['VSINI'][dw] = 10.**tab['FPARAM'][:,7].astype(np.float32)
    giant=np.where(np.core.defchararray.find(tab['ASPCAP_GRID'],b'GKg') |
                np.core.defchararray.find(tab['ASPCAP_GRID'],b'Mg') ) [0]
    tab['VMACRO'][giant] = 10.**tab['FPARAM'][:,7].astype(np.float32)

    # named element flags
    elems, elemtoh, tagnames, elemfitnames = aspcap.elems()
    newtags=[]
    for name in tagnames :
        col = Column(np.full([len(tab)],np.nan),name=name,dtype=np.float32)
        try : tab.remove_column(name)
        except: pass
        tab.add_column(col)
        col = Column(np.full([len(tab)],np.nan),name=name+'_ERR',dtype=np.float32)
        try : tab.remove_column(name+'_ERR')
        except: pass
        tab.add_column(col)
        col = Column(np.full([len(tab)],np.uint32(0)),name=name+'_FLAG',dtype=np.uint32)
        try : tab.remove_column(name+'_FLAG')
        except: pass
        tab.add_column(col)
    ife =np.where(elems == 'Fe')[0][0]
    for i,(el,tag) in enumerate(zip(elems,tagnames)) :
        if el == 'Fe' :
            tab[tag] = tab['X_H'][:,i].astype(np.float32)
        else :
            tab[tag] = tab['X_H'][:,i].astype(np.float32) - tab['X_H'][:,ife].astype(np.float32)
        tab[tag+'_ERR'] = tab['X_H_ERR'][:,i].astype(np.float32)
        tab[tag+'_FLAG'] = tab['ELEMFLAG'][:,i].astype(np.uint32)

def add_extratarg(tab) :
    """ add EXTRATARG, MIN_H, MAX_H, MIN_JK, MAX_JK
    """
    if not isinstance(tab,astropy.table.table.Table) : tab=Table(tab)
    col = Column(np.full([len(tab)],np.uint32(0)),name='EXTRATARG',dtype=np.uint32)
    try : tab.remove_column('EXTRATARG')
    except: pass
    # start with OTHER set
    tab.add_column(col)
    apogee_targ1 = bitmask.ApogeeTarget1()
    apogee_targ2 = bitmask.ApogeeTarget2()
    apogee2_targ1 = bitmask.Apogee2Target1()
    apogee2_targ2 = bitmask.Apogee2Target2()
    
    # main : turn off OTHER
    main = np.where( tab['APOGEE_TARGET1'] & (apogee_targ1.getval('APOGEE_SHORT')|apogee_targ1.getval('APOGEE_MEDIUM')|apogee_targ1.getval('APOGEE_LONG')) ) [0]
    tab['EXTRATARG'][main] = 0
    print('APOGEE main',len(main))
    main = np.where( tab['APOGEE2_TARGET1'] & (apogee2_targ1.getval('APOGEE2_SHORT')|apogee2_targ1.getval('APOGEE2_MEDIUM')|apogee2_targ1.getval('APOGEE2_LONG')) ) [0]
    tab['EXTRATARG'][main] = 0
    print('APOGEE2 main',len(main))

    #telluric
    tell = np.where( tab['APOGEE_TARGET2'] & (apogee_targ2.getval('APOGEE_TELLURIC')|apogee_targ2.getval('APOGEE_TELLURIC_BAD')) )[0]
    tab['EXTRATARG'][tell] |= 4
    print('APOGEE telluric',len(tell))
    tell = np.where( tab['APOGEE2_TARGET2'] & (apogee2_targ2.getval('APOGEE2_TELLURIC')|apogee2_targ2.getval('APOGEE2_TELLURIC_BAD')) )[0]
    tab['EXTRATARG'][tell] |= 4
    print('APOGEE2 telluric',len(tell))

    # apo1m
    j= np.where(tab['TELESCOPE'] == 'apo1m')[0]
    tab['EXTRATARG'][j] |= 8
    print('apo1m',len(j))

def allcal(aspcapfield) :
    """ Apply the calibrations
    """
    teff.cal(aspcapfield)
    logg.cal(aspcapfield)
    cal.elem(aspcapfield)

def aspcapflag(aspcapfield) :
    """ Set bits in ASPCAPFLAG
    """

    parambitmask=bitmask.ParamBitMask()
    aspcapbitmask=bitmask.AspcapBitMask()

    gd=np.where((aspcapfield['ASPCAPFLAG'] & aspcapbitmask.getval('NO_ASPCAP_RESULT') ==0) &
                (aspcapfield['ASPCAPFLAG'] & aspcapbitmask.getval('NO_GRID') == 0))  [0]

    # set ASPCAPFLAG bits for grid edge with final adopted grid
    for iparam,flagname in enumerate(aspcap.params()[2]) :
        j=np.where(aspcapfield['PARAMFLAG'][gd,iparam] & parambitmask.getval('GRIDEDGE_BAD') )[0]
        aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval(flagname+'_BAD')
        print(flagname+'_BAD',len(j))
        j=np.where(aspcapfield['PARAMFLAG'][gd,iparam] & parambitmask.getval('GRIDEDGE_WARN') )[0]
        aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval(flagname+'_WARN') 
        print(flagname+'_WARN',len(j))
    
    # chi**2 
    tmp= aspcapfield['ASPCAP_CHI2']/(aspcapfield['SNR']/100.)**2
    j=np.where(tmp[gd] > 50)[0]
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('CHI2_BAD')
    print('CHI2_BAD',len(j))
    j=np.where((tmp[gd] > 30) & (aspcapfield['ASPCAPFLAG'][gd]&aspcapbitmask.getval('CHI2_BAD')==0) )[0]
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('CHI2_WARN')
    print('CHI2_WARN',len(j))

    # rotation in giant grids
    j=np.where( ( (np.core.defchararray.find(aspcapfield['ASPCAP_GRID'][gd].astype(str),'GKg') >=0)  |
                  (np.core.defchararray.find(aspcapfield['ASPCAP_GRID'][gd].astype(str),'Mg') >=0) ) &
                (aspcapfield['RV_CCFWHM'][gd]/aspcapfield['RV_AUTOFWHM'][gd] > 2.0 ) ) [0]
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('ROTATION_BAD')
    print('ROTATION_BAD',len(j))
    j=np.where( ( (np.core.defchararray.find(aspcapfield['ASPCAP_GRID'][gd].astype(str),'GKg') >=0)  |
                  (np.core.defchararray.find(aspcapfield['ASPCAP_GRID'][gd].astype(str),'Mg') >=0) ) &
                (aspcapfield['RV_CCFWHM'][gd]/aspcapfield['RV_AUTOFWHM'][gd] > 1.5 ) &
                (aspcapfield['ASPCAPFLAG'][gd]&aspcapbitmask.getval('ROTATION_BAD')==0) ) [0]
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('ROTATION_WARN')
    print('ROTATION_WARN',len(j))

    # S/N
    j=np.where(aspcapfield['SNR'][gd] < 30)[0]  
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('SN_BAD')
    print('SNR_BAD',len(j))
    j=np.where((aspcapfield['SNR'][gd] < 70) & ( aspcapfield['ASPCAPFLAG'][gd]&aspcapbitmask.getval('SN_BAD')) )[0]  
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('SN_WARN')
    print('SNR_WARN',len(j))

    #color-Teff
    jk0=(aspcapfield['J']-aspcapfield['K'])-1.5*np.clip(aspcapfield['AK_TARG'],0.,None)
    j = np.where( (aspcapfield['H'][gd] < 90) & (aspcapfield['AK_TARG'][gd] != 0.) &
                  (np.abs(aspcapfield['FPARAM'][gd,0]-teff.cte_ghb(jk0[gd],aspcapfield['FPARAM'][gd,3])[0]) > 1000) )[0]
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('COLORTE_BAD')
    print('COLORTE_BAD',len(j))
    j = np.where( (aspcapfield['H'][gd] < 90) & (aspcapfield['AK_TARG'][gd] != 0.) &
                  (np.abs(aspcapfield['FPARAM'][gd,0]-teff.cte_ghb(jk0[gd],aspcapfield['FPARAM'][gd,3])[0]) > 500) &
                  (aspcapfield['ASPCAPFLAG'][gd] & aspcapbitmask.getval('COLORTE_BAD')==0) )[0]
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('COLORTE_WARN')
    print('COLORTE_WARN',len(j))

    # bad targets
    j = np.where( (np.core.defchararray.find(aspcapfield['TARGFLAGS'][gd].astype(str),'EXTENDED') >=0) |
                  (np.core.defchararray.find(aspcapfield['TARGFLAGS'][gd].astype(str),'EMBEDDED') >=1) |
                  (np.core.defchararray.find(aspcapfield['TARGFLAGS'][gd].astype(str),'APOGEE2_M31') >=1) |
                  (np.core.defchararray.find(aspcapfield['TARGFLAGS'][gd].astype(str),'APOGEE2_M33') >=1) |
                  (np.core.defchararray.find(aspcapfield['TARGFLAGS'][gd].astype(str),'APOGEE2_QSO') >=1) ) [0]
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('PROBLEM_TARGET')
    print('PROBLEM_TARGET',len(j))

    # star level flag
    j = np.where( aspcapfield['ASPCAPFLAG']&aspcapbitmask.badval() > 0)[0]
    aspcapfield['ASPCAPFLAG'][j] |= aspcapbitmask.getval('STAR_BAD')
    print('STAR_BAD',len(j))
    j = np.where( (aspcapfield['ASPCAPFLAG']&aspcapbitmask.warnval() > 0) &
                  (aspcapfield['ASPCAPFLAG']&aspcapbitmask.getval('STAR_BAD') == 0) )[0]
    aspcapfield['ASPCAPFLAG'][j] |= aspcapbitmask.getval('STAR_WARN')
    print('STAR_WARN',len(j))

    #populate character ASPCAPFLAGS
    maxlen=0
    for istar in range(len(aspcapfield)) :
        # set ASPCAPFLAGS character string
        aspcapflags=aspcapbitmask.getname(aspcapfield['ASPCAPFLAG'][istar])
        if len(aspcapflags) > maxlen: 
            maxlen=len(aspcapflags)
            print(aspcapflags)
        aspcapfield['ASPCAPFLAGS'][istar] = aspcapflags
    print('maxlen: ', maxlen)


def fix(tab,visit=None) :
    """ Fix up broken things
    """
    aspcapmask=bitmask.AspcapBitMask()
    parammask=bitmask.ParamBitMask()
    starmask=bitmask.StarBitMask()
    rvmask=bitmask.RVMask()

    # add in SNREV
    try: tab.add_column(Column(tab['SNR'],name='SNREV'))
    except ValueError: print('SNREV already exists...')

    #rename column CLASS to ASPCAP_CLASS
    try: tab.rename_column('CLASS','ASPCAP_GRID')
    except KeyError: print('ASPCAP_GRID already exists')
    try: tab.rename_column('ASPCAP_CLASS','ASPCAP_GRID')
    except KeyError: print('ASPCAP_GRID already exists')
    
    # flags for no good visits for RV
    j=np.where(tab['NVISITS'] == 0)[0]
    tab['STARFLAG'][j] |= starmask.getval('RV_FAIL')
    for jj in j : tab['STARFLAGS'][jj] = starmask.getname(tab['STARFLAG'][jj])
    tab['RV_FLAG'][j] |= rvmask.getval('NO_GOOD_VISITS')

    # FERRE_FAIL (mostly fixed with edge issues?
    j=np.where(tab['FPARAM'][:,0]<-999)[0]
    tab['ASPCAPFLAG'][j] |= aspcapmask.getval('FERRE_FAIL')
    tab['FPARAM'][j,:] = np.nan

    # grid edge flag repair (could have been set for non-adopted grid or first pass)
    for par in ['TEFF','LOGG','M_H','ALPHA_M','C_M','N_M'] :
        tab['ASPCAPFLAG'] &= ~aspcapmask.getval(par+'_BAD')
        tab['ASPCAPFLAG'] &= ~aspcapmask.getval(par+'_WARN')
    
    grids=['GKg','GKd','Mg','Md','Fd']
    teff_pars=[[3500,6000,250],[3500,6000,250],[3000,4000,100],[3000,4000,100],[5500,8000,250]]
    logg_pars=[[0.,4.5,0.5],[2.5,5.5,0.5],[-0.5,3,0.5],[2.5,5.5,0.5],[2.5,5.5,0.5]]
    mh_pars=[[-2.5,1.0,0.25],[-2.5,1.0,0.25],[-2.5,1.0,0.25],[-2.5,1.0,0.25],[-2.5,1.0,0.25]]
    cm_pars=[[-1.5,1.0,0.25],[-1.5,1.0,0.25],[-0.5,0.5,0.25],[-0.5,0.5,0.25],[-0.5,0.0,0.25]]
    nm_pars=[[-0.5,2.0,0.50],[-0.5,2.0,0.50],[-0.5,1.5,0.5],[-0.5,1.5,0.5],[-0.5,1.5,0.25]]
    am_pars=[[-0.75,1.0,0.25],[-0.75,1.0,0.25],[-0.75,1.0,0.25],[-0.75,1.0,0.25],[-0.75,1.0,0.25]]
    
    for grid,teff_par,logg_par,mh_par,cm_par,nm_par,am_par in zip(grids,teff_pars,logg_pars,mh_pars,cm_pars,nm_pars,am_pars) :
        j=np.where((np.core.defchararray.find(tab['ASPCAP_GRID'].astype(str),grid) >= 0))[0]
        print(grid,len(j))
        j=np.where((np.core.defchararray.find(tab['ASPCAP_GRID'].astype(str),grid) >= 0) &
            ((tab['FPARAM'][:,0] < teff_par[0]+teff_par[2]/8.) | 
             (tab['FPARAM'][:,0] > teff_par[1]-teff_par[2]/8.) ) )[0]
        tab['ASPCAPFLAG'][j] |= aspcapmask.getval('TEFF_BAD')
        tab['PARAMFLAG'][j,0] |= parammask.getval('GRIDEDGE_BAD')
        print('TEFF: ',len(j))
        j=np.where((np.core.defchararray.find(tab['ASPCAP_GRID'].astype(str),grid) >= 0) &
            ((tab['FPARAM'][:,0] < teff_par[0]+teff_par[2]/2.) | 
             (tab['FPARAM'][:,0] > teff_par[1]-teff_par[2]/2.) ) )[0]
        tab['ASPCAPFLAG'][j] |= aspcapmask.getval('TEFF_WARN')
        tab['PARAMFLAG'][j,0] |= parammask.getval('GRIDEDGE_WARN')
        print('TEFF: ',len(j))

        j=np.where((np.core.defchararray.find(tab['ASPCAP_GRID'].astype(str),grid) >= 0) &
            ((tab['FPARAM'][:,1] < logg_par[0]+logg_par[2]/8.) | 
             (tab['FPARAM'][:,1] > logg_par[1]+logg_par[2]/8.) ) )[0]
        tab['ASPCAPFLAG'][j] |= aspcapmask.getval('LOGG_BAD')
        tab['PARAMFLAG'][j,1] |= parammask.getval('GRIDEDGE_BAD')
        print('LOGG: ',len(j))
        j=np.where((np.core.defchararray.find(tab['ASPCAP_GRID'].astype(str),grid) >= 0) &
            ((tab['FPARAM'][:,1] < logg_par[0]+logg_par[2]/2.) | 
             (tab['FPARAM'][:,1] > logg_par[1]+logg_par[2]/2.) ) )[0]
        tab['ASPCAPFLAG'][j] |= aspcapmask.getval('LOGG_WARN')
        tab['PARAMFLAG'][j,1] |= parammask.getval('GRIDEDGE_WARN')
        print('LOGG: ',len(j))

        j=np.where((np.core.defchararray.find(tab['ASPCAP_GRID'].astype(str),grid) >= 0) &
            ((tab['FPARAM'][:,3] < mh_par[0]+mh_par[2]/8.) | 
             (tab['FPARAM'][:,3] > mh_par[1]+mh_par[2]/8.) ) )[0]
        tab['ASPCAPFLAG'][j] |= aspcapmask.getval('M_H_BAD')
        tab['PARAMFLAG'][j,3] |= parammask.getval('GRIDEDGE_BAD')
        print('MH: ',len(j))
        j=np.where((np.core.defchararray.find(tab['ASPCAP_GRID'].astype(str),grid) >= 0) &
            ((tab['FPARAM'][:,3] < mh_par[0]+mh_par[2]/2.) | 
             (tab['FPARAM'][:,3] > mh_par[1]+mh_par[2]/2.) ) )[0]
        tab['ASPCAPFLAG'][j] |= aspcapmask.getval('M_H_WARN')
        tab['PARAMFLAG'][j,3] |= parammask.getval('GRIDEDGE_WARN')
        print('MH: ',len(j))

        j=np.where((np.core.defchararray.find(tab['ASPCAP_GRID'].astype(str),grid) >= 0) &
            ((tab['FPARAM'][:,4] < cm_par[0]+cm_par[2]/8.) | 
             (tab['FPARAM'][:,4] > cm_par[1]+cm_par[2]/8.) ) )[0]
        tab['ASPCAPFLAG'][j] |= aspcapmask.getval('C_M_BAD')
        tab['PARAMFLAG'][j,4] |= parammask.getval('GRIDEDGE_BAD')
        print('CM: ',len(j))
        j=np.where((np.core.defchararray.find(tab['ASPCAP_GRID'].astype(str),grid) >= 0) &
            ((tab['FPARAM'][:,4] < cm_par[0]+cm_par[2]/2.) | 
             (tab['FPARAM'][:,4] > cm_par[1]+cm_par[2]/2.) ) )[0]
        tab['ASPCAPFLAG'][j] |= aspcapmask.getval('C_M_WARN')
        tab['PARAMFLAG'][j,4] |= parammask.getval('GRIDEDGE_WARN')
        print('CM: ',len(j))

        j=np.where((np.core.defchararray.find(tab['ASPCAP_GRID'].astype(str),grid) >= 0) &
            ((tab['FPARAM'][:,5] < nm_par[0]+nm_par[2]/8.) | 
             (tab['FPARAM'][:,5] > nm_par[1]+nm_par[2]/8.) ) )[0]
        tab['ASPCAPFLAG'][j] |= aspcapmask.getval('N_M_BAD')
        tab['PARAMFLAG'][j,5] |= parammask.getval('GRIDEDGE_BAD')
        print('NM: ',len(j))
        j=np.where((np.core.defchararray.find(tab['ASPCAP_GRID'].astype(str),grid) >= 0) &
            ((tab['FPARAM'][:,5] < nm_par[0]+nm_par[2]/2.) | 
             (tab['FPARAM'][:,5] > nm_par[1]+nm_par[2]/2.) ) )[0]
        tab['ASPCAPFLAG'][j] |= aspcapmask.getval('N_M_WARN')
        tab['PARAMFLAG'][j,5] |= parammask.getval('GRIDEDGE_WARN')
        print('NM: ',len(j))

        j=np.where((np.core.defchararray.find(tab['ASPCAP_GRID'].astype(str),grid) >= 0) &
            ((tab['FPARAM'][:,6] < am_par[0]+am_par[2]/8.) | 
             (tab['FPARAM'][:,6] > am_par[1]+am_par[2]/8.) ) )[0]
        tab['ASPCAPFLAG'][j] |= aspcapmask.getval('ALPHA_M_BAD')
        tab['PARAMFLAG'][j,6] |= parammask.getval('GRIDEDGE_BAD')
        print('AM: ',len(j))
        j=np.where((np.core.defchararray.find(tab['ASPCAP_GRID'].astype(str),grid) >= 0) &
            ((tab['FPARAM'][:,6] < am_par[0]+am_par[2]/2.) | 
             (tab['FPARAM'][:,6] > am_par[1]+am_par[2]/2.) ) )[0]
        tab['ASPCAPFLAG'][j] |= aspcapmask.getval('ALPHA_M_WARN')
        tab['PARAMFLAG'][j,6] |= parammask.getval('GRIDEDGE_WARN')
        print('AM: ',len(j))

    # SB2 flagging revision
    if visit is not None :
        sb2 = np.where(np.abs(tab['N_COMPONENTS']) > 1)[0]
        print('SB2s: ', len(sb2))
        nbd=0
        for i,s in enumerate(sb2) :
            j=np.where(visit['APOGEE_ID'] == tab['APOGEE_ID'][s])[0]
            n=np.where((visit['N_COMPONENTS'][j]>1) & (visit['SNR'][j]>10) )[0]
            if len(n)/len(j) < 0.1 : 
                print(tab['APOGEE_ID'][s],i,len(n),len(j))
                tab['N_COMPONENTS'][s] = -1*np.abs(tab['N_COMPONENTS'][s])
                tab['STARFLAG'][s] &= ~starmask.getval('MULTIPLE_SUSPECT')
                tab['STARFLAGS'][s] = starmask.getname(tab['STARFLAG'][s])
                nbd+=1
        print('SB2s turned off: ', nbd)


# routines here used for DR16 only

def mkcoord(file='allStar-r12-l33-58358.fits') :
    """ Create coordinate CSV from allStar file, to use for GAIA cross-match
    """
    hdulist=fits.open(file)
    data=hdulist[1].data
    coords=np.vstack((data['RA'],data['DEC'])).T
    np.savetxt('coords.csv',coords,delimiter=',',header='ra,dec')
    # now get the 2MASS object IDs for GAIA crossmatch
    j=np.where(np.core.defchararray.find(data['APOGEE_ID'],'2M') == 0)[0]
    out=np.unique(np.core.defchararray.replace(data['APOGEE_ID'][j],'2M',''))
    np.savetxt('twomass.txt',out,fmt='%s')


def add_gaia(data,gaia_1='gaia_2mass_xmatch.fits.gz', gaia_2='gaia_posn_xmatch.fits.gz') :
    """ Add GAIA data to allStar file, with coordinate match to (cross-matched) GAIA reference file
    """
    tab=Table(data)
    in_names=('source_id','parallax','parallax_error','pmra','pmra_error','pmdec','pmdec_error',
              'phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag','a_g_val', 'e_bp_min_rp_val',
              'radial_velocity','radial_velocity_error', 'r_est','r_lo','r_hi')
    dtypes=('i8','f8','f8','f8','f8','f8','f8','f4','f4','f4','f4','f4','f8','f8','f8','f8','f8')
    out_names=[]
    for name in in_names: out_names.append(('gaia_'+name).upper())
    newcols=Table(np.zeros([len(tab),len(out_names)])-9999.,names=out_names,dtype=dtypes)
    # for source_id, default to 0, not -9999.
    newcols['GAIA_SOURCE_ID'] = 0
    # get rid of targetting proper motions to avoid confusion!
    #tab.remove_columns(['PMRA','PMDEC','PM_SRC'])
    # change to rename
    tab.rename_column('PMRA','TARG_PMRA')
    tab.rename_column('PMDEC','TARG_PMDEC')
    tab.rename_column('PM_SRC','TARG_PM_SRC')
    # add unpopulated columns
    tab.add_columns(newcols.columns.values())

    # read gaia 2MASS matched file, match by 2MASS ID, and populate
    gaia=fits.open(gaia_1)[1].data
    print('number in GAIA-2MASS xmatch catalog: ',len(gaia),len(set(gaia['original_ext_source_id'])))
    while True :
        # loop for matches since we have repeats and want them all matched
        j=np.where(tab['GAIA_SOURCE_ID'] == 0)[0]
        print('Number missing gaia_source_id: ', len(j))
        m1,m2=match.match(np.core.defchararray.replace(tab['APOGEE_ID'][j],'2M',''),gaia['original_ext_source_id'])
        print('Number matched by 2MASS: ', len(m1))
        if len(m1) == 0 : break
        for inname,outname in zip(in_names,out_names) :
            tab[outname][j[m1]] = gaia[inname][m2]
        #if len(m1) < 100 : 
        #    for i in m1 : print(tab['APOGEE_ID'][j[m1]])
    j=np.where(tab['GAIA_SOURCE_ID'] > 0)[0]
    print('number of unique APOGEE_ID matches: ',len(set(tab['APOGEE_ID'][j])))

    j=np.where(tab['GAIA_SOURCE_ID'] == 0)[0]
    print('missing sources after 2MASS matches: ',len(j))
    gaia=fits.open(gaia_2)[1].data
    h=htm.HTM()
    # now do a positional match, take the brightest object within 3 arcsec (which is the max from the GAIA crossmatch)
    maxrad=3./3600.
    m1,m2,rad=h.match(tab['RA'][j],tab['DEC'][j],gaia['RA'],gaia['DEC'],maxrad,maxmatch=10)
    for m in set(m1) :
        jj=np.where(m1 == m)[0]
        ii=np.argsort(gaia['phot_rp_mean_mag'][m2[jj]])
        #print(tab['RA'][j[m]],tab['DEC'][j[m]],tab['TARGET_ID'][j[m]])
        #print(gaia['RA'][m2[jj]],gaia['DEC'][m2[jj]],gaia['PHOT_RP_MEAN_MAG'][m2[jj]])
        #print(ii)
        for inname,outname in zip(in_names,out_names) :
            tab[outname][j[m]] = gaia[inname][m2[jj[ii[0]]]]
    j=np.where(tab['GAIA_SOURCE_ID'] == 0)[0]
    print('missing sources after second match: ',len(j))

    # replace NaNs
    for name in out_names :
        bd = np.where(np.isnan(tab[name]))[0]
        tab[name][bd] = -9999.

    return tab

        
def trimfile(data) :
    """ Write a 'lite' allStar file, removing some of the big space users
    """
    tab=Table(data)
    remove=['ALL_VISITS','VISITS','ALL_VISIT_PK','VISIT_PK','FPARAM_CLASS','CHI2_CLASS',
            'FPARAM','FELEM','FPARAM_COV','PARAM','PARAM_COV','ELEMFLAG','FELEM_ERR',
            'APSTAR_ID','TARGET_ID','ASPCAP_ID','FILE','LOCATION_ID']
    out=fits.HDUList()
    tab.remove_columns(remove)

    return tab

def new(infile='allStar-r12-l33-58358.fits',new='allStar-r12-l33.fits',trim='allStarLite-r12-l33.fits',
        gaia_1='gaia/gaia_2mass_xmatch.fits.gz', gaia_2='gaia/gaia_posn_xmatch.fits.gz',dr16new=False) :
    """ take allStar file, add GAIA info and new _spec columns, outpu
        also output allStarLite version
    """
    if dr16new :
        hdulist=fits.open(infile)
        tab=dr16fix(hdulist[1].data,dr16file='sav/allStar-r12-l33-58358.fits')
    else :
        hdulist=fits.open(infile)
        tab=hdulist[1].data
    gd=np.where(np.core.defchararray.strip(tab['APOGEE_ID']) != '')[0]
    print('adding gaia....')
    tab=add_gaia(tab[gd],gaia_1=gaia_1,gaia_2=gaia_2)
    print('adding _SPEC columns....')
    tab=add_spec(tab)
    # populate TARGFLAG for 1m observations
    t2a=bitmask.ApogeeTarget2()
    t2b=bitmask.Apogee2Target2()
    j=np.where(tab['TELESCOPE'] == 'apo1m')[0]
    tab['APOGEE2_TARGET2'][j] |= t2a.getval('APOGEE_1MTARGET')
    tab['APOGEE_TARGET2'][j] |= t2b.getval('APOGEE2_1MTARGET')
    # repopulate TARGFLAGS with latest bitmask info
    for i in range(len(tab)) :
        if i%10000 == 0 : print('repopulating targflags: ',i)
        if 'apogee2' in tab['SURVEY'][i] :
            tab['TARGFLAGS'][i] = bitmask.targflags(
                                      tab['APOGEE2_TARGET1'][i],
                                      tab['APOGEE2_TARGET2'][i],
                                      tab['APOGEE2_TARGET3'][i],survey='apogee2')
        else :
            tab['TARGFLAGS'][i] = bitmask.targflags(
                                      tab['APOGEE_TARGET1'][i],
                                      tab['APOGEE_TARGET2'][i],0,survey='apogee')

    # write out the modified file
    print('writing file ',new)
    out=fits.HDUList()
    hdu=fits.PrimaryHDU()
    hdu.header.add_comment('APOGEE_VER:'+os.environ['APOGEE_VER'])
    out.append(hdu)
    out.append(fits.BinTableHDU(tab))
    out.append(hdulist[2])
    out.append(hdulist[3])
    out.writeto(new,overwrite=True)

    print('writing file ',trim)
    out=fits.HDUList()
    out.append(hdu)
    out.append(fits.BinTableHDU(trimfile(tab)))
    out.writeto(trim,overwrite=True)


def dr16fix(a,dr16file='allStar-r12-l33-58358.fits') :
    """  Routine to take a modified allStar file plus the original DR16 allStarfile
         and return the new version, but with the stars in the exact same order/location
         as the original DR16 file. Written after fixes for better apogeeObject target
         matching (which led to incorrect matches in original DR16 file), and changes
         for TARGET_ID and some of the selection function tags
    """

    # get old and new files
    old=fits.open(dr16file)[1].data

    # loop over the old file entries, and find the matching entry in new file
    #   based on a number of tags (but not necessarily APOGEE_ID
    ind=[]
    for i in range(len(old)) :
        if ((a['APOGEE_ID'][i] != old['APOGEE_ID'][i]) or
            (a['FIELD'][i] != old['FIELD'][i]) or
            (a['TEFF'][i] != old['TEFF'][i]) or
            (a['FPARAM'][i,0] != old['FPARAM'][i,0]) or
            (a['FPARAM'][i,1] != old['FPARAM'][i,1]) or
            (a['STARFLAG'][i] != old['STARFLAG'][i]) or
            (a['SNR'][i] != old['SNR'][i]) or
            (a['TELESCOPE'][i] != old['TELESCOPE'][i]) or
            (a['MEANFIB'][i] != old['MEANFIB'][i]) 
            ) :
            match=-1
            for j in range(-10,10) :
                if ((a['FIELD'][i+j] == old['FIELD'][i]) and
                    (a['TEFF'][i+j] == old['TEFF'][i]) and
                    (a['FPARAM'][i+j,0] == old['FPARAM'][i,0]) and
                    (a['FPARAM'][i+j,1] == old['FPARAM'][i,1]) and
                    (a['STARFLAG'][i+j] == old['STARFLAG'][i]) and
                    (a['SNR'][i+j] == old['SNR'][i]) and
                    (a['TELESCOPE'][i+j] == old['TELESCOPE'][i]) and
                    (a['MEANFIB'][i+j] == old['MEANFIB'][i])
                    ) :
                    #print(i,i+j,a['APOGEE_ID'][i],old['APOGEE_ID'][i],a['TEFF'][i],old['TEFF'][i])
                    match=i+j
                    break
            if match < 0 : print('no match: ',old['APOGEE_ID'][i])
        else : match = i
        ind.append(match)
    print(len(ind),len(set(ind)))
    print(set(range(len(a)))-set(ind))
    return a[ind]

def dr16check(file='newStar-r12-l33.fits',dr16file='allStar-r12-l33.fits') :

    new=fits.open(file)[1].data
    old=fits.open(dr16file)[1].data
    print(len(old),len(new))
    nbad=0
    fp=open('dr16check.txt','w')
    for i in range(len(old)):
        if i%1000 == 0 : print(i)
        for col in old.columns.names :
            if col == 'PARAM_COV' : continue
            elif col == 'FPARAM_COV' : continue
            elif col == 'FPARAM_CLASS' : continue
            elif col == 'REDUCTION_ID' : continue
            elif col == 'TARGET_ID' : continue
            #elif col == 'MIN_H' : continue
            #elif col == 'MAX_H' : continue
            #elif col == 'MIN_JK' : continue
            #elif col == 'MAX_JK' : continue
            elif col == 'GAIA_PARALLAX' and np.isnan(old[col][i]) : continue
            elif col == 'GAIA_PARALLAX_ERROR' and np.isnan(old[col][i]) : continue
            elif col == 'GAIA_PMRA' and np.isnan(old[col][i]) : continue
            elif col == 'GAIA_PMRA_ERROR' and np.isnan(old[col][i]) : continue
            elif col == 'GAIA_PMDEC' and np.isnan(old[col][i]) : continue
            elif col == 'GAIA_PMDEC_ERROR' and np.isnan(old[col][i]) : continue
            elif col == 'GAIA_RADIAL_VELOCITY' and np.isnan(old[col][i]) : continue
            elif col == 'GAIA_RADIAL_VELOCITY_ERROR' and np.isnan(old[col][i]) : continue
            elif col == 'GAIA_R_EST' and np.isnan(old[col][i]) : continue
            elif col == 'GAIA_R_LO' and np.isnan(old[col][i]) : continue
            elif col == 'GAIA_R_HI' and np.isnan(old[col][i]) : continue
            elif col == 'GAIA_PHOT_BP_MEAN_MAG' and np.isnan(old[col][i]) : continue
            elif col == 'GAIA_PHOT_RP_MEAN_MAG' and np.isnan(old[col][i]) : continue
            try: 
                for j in range(len(new[col][i])) :
                    if new[col][i,j] != old[col][i,j] :
                        print(i,col,new[col][i,j],old[col][i,j],new['FIELD'][i],old['FIELD'][i])
                        fp.write('{:d} {:s} {} {} {:s} {:s}\n'.format(i,col,new[col][i,j],old[col][i,j],new['FIELD'][i],old['FIELD'][i]))
            except:
                if new[col][i] != old[col][i] :
                    print(i,col,new[col][i],old[col][i],new['FIELD'][i],old['FIELD'][i])
                    fp.write('{:d} {:s} {} {} {:s} {:s}\n'.format(i,col,new[col][i],old[col][i],new['FIELD'][i],old['FIELD'][i]))

        #if ((new['FIELD'][i] != old['FIELD'][i]) or
        #    (new['TEFF'][i] != old['TEFF'][i]) or
        #    (new['FPARAM'][i,0] != old['FPARAM'][i,0]) or
        #    (new['FPARAM'][i,1] != old['FPARAM'][i,1]) or
        #    (new['STARFLAG'][i] != old['STARFLAG'][i]) or
        #    (new['TARGFLAGS'][i] != old['TARGFLAGS'][i]) or
        #    (new['FE_H'][i] != old['FE_H'][i]) or
        #    (new['MG_FE'][i] != old['MG_FE'][i]) or
        #    (new['SNR'][i] != old['SNR'][i]) or
        #    (new['TELESCOPE'][i] != old['TELESCOPE'][i]) or
        #    (new['H'][i] != old['H'][i]) or
        #    (new['AK_TARG'][i] != old['AK_TARG'][i]) or
        #    (new['AK_WISE'][i] != old['AK_WISE'][i]) or
        #    (new['APOGEE_ID'][i] != old['APOGEE_ID'][i]) or
        #    (new['MEANFIB'][i] != old['MEANFIB'][i]) 
        #   ) :
        #     print(i,new['APOGEE_ID'][i],old['APOGEE_ID'][i],new['H'][i],old['H'][i])
        #     print('   ',new['AK_TARG'][i],old['AK_TARG'][i],new['FIELD'][i],old['FIELD'][i])
        #     print('   ',new['AK_TARG_METHOD'][i],old['AK_TARG_METHOD'][i],new['SURVEY'][i],old['SURVEY'][i])
        #     nbad+=1

    print('nbad:',nbad)
    fp.close()
 
