import matplotlib.pyplot as plt
import numpy as np
from esutil import htm
import astropy
from astropy.table import Table, Column, vstack
from astropy.io import fits, ascii
from apogee.utils import bitmask, apselect,gaia, apload
from apogee.apred import apstar
from apogee.aspcap import aspcap, teff, logg, cal, err, elem, qa
from tools import match, struct, html, plots
import os
import pdb
import glob
import yaml
import shutil
import time
from multiprocessing import Process

def all(planfile,dofix=False,suffix=None,allplate=True, calsample=False) :
   """ Summary files after everything is run
   """
   plan=yaml.safe_load(open(planfile,'r'))
   apred_vers=plan['apred_vers']
   apstar_vers=plan['apstar_vers']
   aspcap_vers=plan['aspcap_vers']
   if suffix==None : suffix=plan['suffix']
   outdir='allStar-'+apred_vers+'-'+aspcap_vers+suffix
   try: os.mkdir(outdir)
   except: pass

   # allStar file
   apstar_dir=os.environ['APOGEE_REDUX']+'/'+apred_vers+'/'+apstar_vers+'/'
   aspcap_dir=os.environ['APOGEE_ASPCAP']+'/'+apred_vers+'/'+aspcap_vers+'/'
   print('Create allStar file')
   if calsample : 
       tab,dat=allStar(search=[aspcap_dir+'apo*/*/aspcapField-*.fits',aspcap_dir+'lco*/*/aspcapField-*.fits'],
               skip=['Field-apo25m_','Field-lco25m_','Field-apo1m_','apo25m.','lco25m.'],out=None,dofix=dofix,outdir=outdir)
   else : 
       tab,dat=allStar(search=[aspcap_dir+'apo*/*/aspcapField-*.fits',aspcap_dir+'lco*/*/aspcapField-*.fits'],out=None,dofix=dofix,outdir=outdir)

   # allVisit file
   print('Create allVisit file')
   allvisit = apstar.allFieldVisit(search=[apstar_dir+'apo*/*/a?FieldVisits-*.fits',apstar_dir+'lco*/*/a?FieldVisits-*.fits'],out=None,apred_vers=apred_vers)
   allvisit.write(aspcap_dir+'allVisit-'+apred_vers+'-'+aspcap_vers+suffix+'.fits',overwrite=True)

   # add VISIT_PK
   nvisits=add_visitpk(tab,allvisit)

   # fix up issues in allStar file that require allVisit (SB2s)
   if dofix :
       print('fix allStar file')
       fix(tab,allvisit)

   # write allStar out
   write(tab, dat, dat, aspcap_dir+'allStar-'+apred_vers+'-'+aspcap_vers+suffix+'.fits') 
 
   # allPlate file
   if allplate :
       print('Create allPlate file')
       allplate = apstar.allPlate(allvisit,
                  out=aspcap_dir+'allPlate-'+apred_vers+'-'+aspcap_vers+suffix+'.fits')

   return nvisits

def allStar(search=['apo*/*/aspcapField-*.fits','lco*/*/aspcapField-*.fits'],out='allStar.fits',
            skip=['Field-redo','Field-cal_','Field-apo25m_','Field-lco25m_','Field-apo1m_','apo25m.','lco25m.','apo1m.'], 
            doerr=True,docal=True,doextratarg=True,dofix=False,addgaia=False,outdir=None) :
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
        if addgaia: 
            try: dat=gaia.add_gaia(dat)
            except: print('failed to add_gaia: ',file)

        a.append(dat)
    # stack them
    tab =vstack(a)
    del(a)
    tab.sort(['RA','DEC','FIELD'])
    tab3=Table.read(file,hdu=3)

    # fixes
    #fix(tab)
    washphotfix(tab)
    targflagfix(tab)
    elemflagfix(tab)
    nanfix(tab)

    # empirical uncertainties
    if doerr : repeat_err(tab,tab3,outdir=outdir)

    # add EXTRATARG, H_MIN, H_MAX, JKMIN, JKMAX
    if doextratarg: add_extratarg(tab)

    # add MEMBERS
    add_members(tab)

    # set ASPCAPFLAG bits and ASPCAPFLAGS
    aspcapflag(tab)

    # calibration depends on having STAR_BAD set
    if docal : allcal(tab,tab3,outdir=outdir)

    # add named tags
    add_named_tags(tab)

    # add unique identifiers for database
    col = Column(aspcap_id(tab),name='ASPCAP_ID')
    tab.add_column(col,index=2)
    col = Column(apstar_id(tab),name='APSTAR_ID')
    tab.add_column(col,index=2)
    col = Column(target_id(tab),name='TARGET_ID')
    tab.add_column(col,index=2)

    # write out the file
    if out is not None: hdulist = write(tab,tab3,tab3,out)
   
    outdir=outdir+'/' 
    qa_plots(tab,tab3,prefix=outdir)
    qa_html(tab,tab3,prefix=outdir)

    return tab, tab3

def rewrite(planfile) :
    """ rewrite aspcapField and aspcapStar files with allStar record
    """      
    plan=yaml.safe_load(open(planfile,'r'))
    apred_vers=plan['apred_vers']
    apstar_vers=plan['apstar_vers']
    aspcap_vers=plan['aspcap_vers']
    telescope=plan['telescope']
    aspcap_dir=os.environ['APOGEE_ASPCAP']+'/'+apred_vers+'/'+aspcap_vers+'/'
 
    tab=Table.read(aspcap_dir+'/allStar-'+apred_vers+'-'+aspcap_vers+'.fits')

    load=apload.ApLoad(apred=apred_vers,aspcap=aspcap_vers,telescope=telescope)
    for field in plan['field'] : 
        print(field)
        aspcapfield=load.aspcapField(field)
        dat=Table(aspcapfield[1].data)
        spec=Table(aspcapfield[2].data)
        key=Table(aspcapfield[3].data)
        j=np.where(( tab['FIELD'].astype(str) == field) & (tab['TELESCOPE'].astype(str) == telescope) )[0]
        i1,i2=match.match(tab['APOGEE_ID'][j],dat['APOGEE_ID'])
        if len(i1) != len(j) or  len(i1) != len(dat) : pdb.set_trace()
        name=os.path.dirname(load.filename('aspcapField',field=field))
        shutil.move(name,name+'.orig')
        aspcap.writefiles(load,field,tab[j[i1]],spec[i2],key,version=aspcapfield[0].header['VERSION'])


def write(tab, tab2, tab3, out) :
    """ Write out allStar file from input HDU tables
    """
    # construct HDUList
    hdulist=fits.HDUList()
    hdulist.append(fits.BinTableHDU(tab))
    hdulist.append(fits.BinTableHDU(tab2))
    hdulist.append(fits.BinTableHDU(tab3))
    # write out the file
    if out is not None:
        print('writing',out)
        hdulist[0].header['VERSION'] = (os.environ['APOGEE_VER'],'APOGEE software version APOGEE_VER')
        hdulist.writeto(out,overwrite=True)

    return hdulist

def doskip(file,skip) :
    for sk in skip : 
        if sk in file : return True
    return False

def repeat_err(tab,tab3,outdir=None) :
    """ get scatter from repeat observations
    """
    start = time.time()
    print('running repeat....')
    try: os.makedirs(outdir+'/repeat')
    except FileExistsError : pass

    procs=[]
    kw={'data' : tab, 'out' : outdir+'/repeat/giant_', 'params' : tab3['PARAM_SYMBOL'][0].astype(str), 
        'elems' : tab3['ELEM_SYMBOL'][0].astype(str), 'logg' : [-1,3.8], 'teff' : [3000,6000]}
    procs.append(Process(target=err.repeat,kwargs=kw))
    kw={'data' : tab, 'out' : outdir+'/repeat/dwarf_', 'params' : tab3['PARAM_SYMBOL'][0].astype(str), 
        'elems' : tab3['ELEM_SYMBOL'][0].astype(str), 'logg' : [3.8,5.5]}
    procs.append(Process(target=err.repeat,kwargs=kw))
    for proc in procs : proc.start()
    for proc in procs : proc.join()
    print('repeat elapsed: ',time.time()-start)

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

    aspcapmask=bitmask.AspcapBitMask()
    parammask=bitmask.ParamBitMask()
    gd=np.where(((tab['ASPCAPFLAG'] & aspcapmask.badval()) == 0) &
                ((tab['ASPCAPFLAG'] & aspcapmask.getval('CHI2_BAD')) == 0) )[0]

    tab['TEFF'][gd] = tab['PARAM'][gd,0].astype(np.float32)
    tab['TEFF_ERR'][gd] = np.sqrt(tab['PARAM_COV'][gd,0,0]).astype(np.float32)
    tab['TEFF_SPEC'][gd] = tab['FPARAM'][gd,0].astype(np.float32)

    tab['LOGG'][gd] = tab['PARAM'][gd,1].astype(np.float32)
    tab['LOGG_ERR'][gd] = np.sqrt(tab['PARAM_COV'][gd,1,1]).astype(np.float32)
    tab['LOGG_SPEC'][gd] = tab['FPARAM'][gd,1].astype(np.float32)

    # populated named abundance tags only off grid edge
    gdel = np.where((tab['PARAMFLAG'][gd,3]&parammask.badval()) == 0)[0]
    tab['M_H'][gd[gdel]] = tab['PARAM'][gd[gdel],3].astype(np.float32)
    tab['M_H_ERR'][gd[gdel]] = np.sqrt(tab['PARAM_COV'][gd[gdel],3,3]).astype(np.float32)

    gdel = np.where((tab['PARAMFLAG'][gd,6]&parammask.badval()) == 0)[0]
    tab['ALPHA_M'][gd[gdel]] = tab['PARAM'][gd[gdel],6].astype(np.float32)
    tab['ALPHA_M_ERR'][gd[gdel]] = np.sqrt(tab['PARAM_COV'][gd[gdel],6,6]).astype(np.float32)

    tab['VMICRO'][gd] = 10.**tab['FPARAM'][gd,2].astype(np.float32)
    # populate VSINI and VMACRO for dwarfs (latter is fixed)
    dw=np.where((np.core.defchararray.find(tab['ASPCAP_GRID'][gd].astype(str),'BA') >= 0) |
                (np.core.defchararray.find(tab['ASPCAP_GRID'][gd].astype(str),'GKd')>= 0)  |
                (np.core.defchararray.find(tab['ASPCAP_GRID'][gd].astype(str),'Fd')>= 0)  |
                (np.core.defchararray.find(tab['ASPCAP_GRID'][gd].astype(str),'Md')>= 0)  ) [0]
    tab['VMACRO'][gd[dw]] = 0.
    tab['VSINI'][gd[dw]] = 10.**tab['FPARAM'][gd[dw],7].astype(np.float32)

    # don't populate VSINI for giants (leave as NaN) since these are not derived
    giant=np.where((np.core.defchararray.find(tab['ASPCAP_GRID'][gd].astype(str),'GKg') >=0) |
                   (np.core.defchararray.find(tab['ASPCAP_GRID'][gd].astype(str),'Mg')  >=0)) [0]
    # VMACRO for giants is also fixed, at the adopted relation
    tab['VMACRO'][gd[giant]] = 10.**tab['FPARAM'][gd[giant],7].astype(np.float32)

    # named element flags
    elems, elemtoh, tagnames, elemfitnames = aspcap.elems()
    newtags=[]
    for name in tagnames :
        col = Column(np.full([len(tab)],np.nan),name=name,dtype=np.float32)
        try : tab.remove_column(name)
        except: pass
        tab.add_column(col)
        col = Column(np.full([len(tab)],np.nan),name=name+'_SPEC',dtype=np.float32)
        try : tab.remove_column(name+'_SPEC')
        except: pass
        tab.add_column(col)
        col = Column(np.full([len(tab)],np.nan),name=name+'_ERR',dtype=np.float32)
        try : tab.remove_column(name+'_ERR')
        except: pass
        tab.add_column(col)
        col = Column(np.full([len(tab)],np.int32(0)),name=name+'_FLAG',dtype=np.int32)
        try : tab.remove_column(name+'_FLAG')
        except: pass
        tab.add_column(col)
    ife =np.where(elems == 'Fe')[0][0]
    for i,(el,tag) in enumerate(zip(elems,tagnames)) :
        gdel = np.where((tab['ELEMFLAG'][gd,i]&parammask.badval()) == 0)[0]
        if el == 'Fe' :
            tab[tag][gd[gdel]] = tab['X_H'][gd[gdel],i].astype(np.float32)
            tab[tag+'_SPEC'][gd[gdel]] = tab['X_H_SPEC'][gd[gdel],i].astype(np.float32)
        else :
            tab[tag][gd[gdel]] = tab['X_H'][gd[gdel],i].astype(np.float32) - tab['X_H'][gd[gdel],ife].astype(np.float32)
            tab[tag+'_SPEC'][gd[gdel]] = tab['X_H_SPEC'][gd[gdel],i].astype(np.float32) - tab['X_H_SPEC'][gd[gdel],ife].astype(np.float32)
        tab[tag+'_ERR'][gd[gdel]] = tab['X_H_ERR'][gd[gdel],i].astype(np.float32)
        tab[tag+'_FLAG'] = tab['ELEMFLAG'][:,i].astype(np.int32)

def add_extratarg(tab) :
    """ add EXTRATARG, MIN_H, MAX_H, MIN_JK, MAX_JK
    """
    if not isinstance(tab,astropy.table.table.Table) : tab=Table(tab)
    # start with OTHER set
    col = Column(np.full([len(tab)],np.int32(1)),name='EXTRATARG',dtype=np.int32)
    try : tab.remove_column('EXTRATARG')
    except: pass
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

    # duplicates
    # get indices in file for each degree of RA, to shorten search
    ind = np.zeros(360,dtype=int)
    for i in range(1,360) : 
        ind[i] = np.where(tab['RA']>i)[0][0] 
        print(i,ind[i])

    print('duplicates: ')
    jall=[]
    for i,star in enumerate(tab) :
        ira=int(star['RA'])
        i1 = ind[ira]
        if ira < 359 : i2 = ind[ira+1 ]
        else : i2 = len(tab)
        j = i1 + np.where(tab['APOGEE_ID'][i1:i2] == star['APOGEE_ID'])[0]
        print(i,len(j))
        if len(j) > 1 :
            jall.extend(j)
            jsort = np.argsort(tab['SNR'][j])
            tab['EXTRATARG'][j[jsort[0:-1]]] |= 16
    dup = np.where(tab['EXTRATARG']&16)[0]
    print('duplicates: ',len(dup))
    return jall

def allcal(tab,tab3,outdir='./',doteff=True,dologg=True,doelem=True, doerr=True, calvers='dr17', calib=False, caldir='calib/') :
    """ Apply the calibrations
    """

    caldir=outdir+'/'+caldir
    try: os.makedirs(caldir)
    except FileExistsError : pass

    #uncertainties
    if doerr :
        err.apply(tab,caldir=outdir+'/repeat/')

    # Teff
    if doteff: 
        ebvmax=0.03
        teffcal = teff.ghb(tab,ebvmax=ebvmax,glatmin=10,out=caldir+'/tecal',yr=[-750,750],trange=[4500,7000],loggrange=[-1,6],calib=calib,doerr=False)
        Table(struct.dict2struct(teffcal)).write(caldir+'/tecal.fits',overwrite=True)
        if not calib : teff.cal(tab,caldir=caldir)
        # add calibration correction with same color scheme
        fig,ax=plots.multi(1,1)
        plots.plotc(ax,tab['FPARAM'][:,0],tab['FPARAM'][:,1],tab['FPARAM'][:,0]-tab['PARAM'][:,0],
                    xr=[8000,3000],xt='Teff',yt='log g',zt=r'$\Delta$ Teff',yr=[6,-1],zr=[-250,250],colorbar=True)
        outfile=caldir+'/calib_Teff_hr.png'
        fig.savefig(outfile)
        plt.close(fig)

        figs=[['tecal.png','tecal_berger_mh.png'],['tecal_ghb_hr.png','tecal_berger_hr.png'],['tecal_b.png','calib_Teff_hr.png']]
        ytitle=['Teff all stars','GHB&Berger HR','GHB']
        html.htmltab(figs,ytitle=ytitle,file=caldir+'/teff.html')

    # log g
    if dologg: 
        logg.nn_train(tab,out=caldir,calib=calib)
        if not calib : logg.nn_cal(tab,caldir=caldir,out=caldir)
        grid=[]
        grid.append(['nn_logg_cal.png','nn_logg_cal_hr.png'])
        grid.append(['nn_logg_hr.png','nn_logg_scatter.png'])
        grid.append(['logg_correction.png','logg_correction_hr.png'])
        html.htmltab(grid,caldir+'logg.html')
        # add RC/RGB split for information purposes only in PARAMFLAG[1]
        rcrgb = logg.rcrgb(tab,rclim=np.array([2.2,3.5]),out=caldir+'/rcrgb')
        Table(struct.dict2struct(rcrgb)).write(caldir+'/rcrgb.fits',overwrite=True)
        logg.rcrgb_class(tab,rcrgb,out=caldir)

    # abundances
    if doelem: 
        #cal.elemcal(tab,caldir=caldir)
        elems=tab3['ELEM_SYMBOL'][0].astype(str)
        elemtoh=tab3['ELEMTOH'][0]

        for col in ['GIANT_SOLARNEIGH_ZERO','DWARF_SOLARNEIGH_ZERO', 'ALLDWARF_SOLARNEIGH_ZERO', 'SOLAR_ZERO'] :
            try: tab3.remove_column(col)
            except: pass
            try: tab3.remove_column(col+'_ERR')
            except: pass

        # need to calibrate M before elems to populate X_H
        doels = np.append(np.array(['M','alpha']),elems)

        # solar neighborhood giants
        solar=apselect.solar(tab,logg=[-1,3.8])
        elemcal=elem.zerocal(tab,solar,elems,elemtoh,doels,calvers=calvers,calib=calib,extfit=4,dwarfs=False)
        tab3.add_column(Column([elemcal['extpar'][:,0]],name='GIANT_SOLARNEIGH_ZERO'))
        tab3.add_column(Column([elemcal['exterr']],name='GIANT_SOLARNEIGH_ZERO_ERR'))
        Table(elemcal).write(caldir+'/giant_abuncal.fits',overwrite=True)

        # solar neighborhood 4500-5000 dwarfs
        solar=apselect.solar(tab,logg=[4,6],teff=[4500,5000])
        elemcal=elem.zerocal(tab,solar,elems,elemtoh,doels,calvers=calvers,calib=calib,extfit=4,dwarfs=True)
        tab3.add_column(Column([elemcal['extpar'][:,0]],name='DWARF_SOLARNEIGH_ZERO'))
        tab3.add_column(Column([elemcal['exterr']],name='DWARF_SOLARNEIGH_ZERO_ERR'))
        Table(elemcal).write(caldir+'/dwarf_abuncal.fits',overwrite=True)

        # solar neighborhood all dwarfs
        solar=apselect.solar(tab,logg=[4,6])
        elemcal=elem.zerocal(tab,solar,elems,elemtoh,doels,calvers=calvers,calib=calib,extfit=4,dwarfs=True)
        tab3.add_column(Column([elemcal['extpar'][:,0]],name='ALLDWARF_SOLARNEIGH_ZERO'))
        tab3.add_column(Column([elemcal['exterr']],name='ALLDWARF_SOLARNEIGH_ZERO_ERR'))

        # VESTA
        j=np.where(tab['APOGEE_ID'] == 'VESTA')[0]
        elemcal=elem.zerocal(tab,j,elems,elemtoh,doels,calvers=calvers,calib=calib,extfit=2,dwarfs=True)
        tab3.add_column(Column([elemcal['extpar'][:,0]],name='SOLAR_ZERO'))
        tab3.add_column(Column([elemcal['exterr']],name='SOLAR_ZERO_ERR'))

        # populate calibrated quantities 
        elem.cal(tab,tab3,caldir=caldir)

    # VMICRO, VSINI, O copy from FPARAM
    for i in [2,7,8] :
        tab['PARAM'][:,i] = tab['FPARAM'][:,i]

def aspcapflag(aspcapfield) :
    """ Set bits in ASPCAPFLAG
    """

    starbitmask=bitmask.StarBitMask()
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
    aspcapfield['ASPCAPFLAG'] &= ~aspcapbitmask.getval('CHI2_BAD')
    j=np.where(aspcapfield['ASPCAP_CHI2'][gd] > 50*(aspcapfield['SNR'][gd]/100.)**2 + 2) [0]
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('CHI2_BAD')
    print('CHI2_BAD',len(j))
    aspcapfield['ASPCAPFLAG'] &= ~aspcapbitmask.getval('CHI2_WARN')
    j=np.where((aspcapfield['ASPCAP_CHI2'][gd] > 30*(aspcapfield['SNR'][gd]/100.)**2 + 2) &
               (aspcapfield['ASPCAPFLAG'][gd]&aspcapbitmask.getval('CHI2_BAD')==0) )[0]
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('CHI2_WARN')
    print('CHI2_WARN',len(j))

    # rotation in giant grids, warn only
    #j=np.where( ( (np.core.defchararray.find(aspcapfield['ASPCAP_GRID'][gd].astype(str),'GKg') >=0)  |
    #              (np.core.defchararray.find(aspcapfield['ASPCAP_GRID'][gd].astype(str),'Mg') >=0) ) &
    #            (aspcapfield['RV_CCFWHM'][gd]/aspcapfield['RV_AUTOFWHM'][gd] > 4.0 ) &
    #            (aspcapfield['SNR'][gd] > 10) ) [0]
    #aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('ROTATION_BAD')
    #print('ROTATION_BAD',len(j))
    j=np.where( ( (np.core.defchararray.find(aspcapfield['ASPCAP_GRID'][gd].astype(str),'GKg') >=0)  |
                  (np.core.defchararray.find(aspcapfield['ASPCAP_GRID'][gd].astype(str),'Mg') >=0) ) &
                (aspcapfield['RV_CCFWHM'][gd]/aspcapfield['RV_AUTOFWHM'][gd] > 2.0 ) &
                (aspcapfield['ASPCAPFLAG'][gd]&aspcapbitmask.getval('ROTATION_BAD')==0) ) [0]
    aspcapfield['ASPCAPFLAG'] &= ~aspcapbitmask.getval('ROTATION_WARN')
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('ROTATION_WARN')
    print('ROTATION_WARN',len(j))

    # S/N
    aspcapfield['ASPCAPFLAG'] &= ~aspcapbitmask.getval('SN_BAD')
    j=np.where(aspcapfield['SNR'][gd] < 30)[0]  
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('SN_BAD')
    print('SNR_BAD',len(j))
    aspcapfield['ASPCAPFLAG'] &= ~aspcapbitmask.getval('SN_WARN')
    j=np.where((aspcapfield['SNR'][gd] < 70) & ( aspcapfield['ASPCAPFLAG'][gd]&aspcapbitmask.getval('SN_BAD')) )[0]  
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('SN_WARN')
    print('SNR_WARN',len(j))

    #color-Teff, but only make this a warning, and only for stars with AK_TARG<1
    aspcapfield['ASPCAPFLAG'] &= ~aspcapbitmask.getval('COLORTE_BAD')
    jk0=(aspcapfield['J']-aspcapfield['K'])-1.5*np.clip(aspcapfield['AK_TARG'],0.,None)
    j = np.where( (aspcapfield['H'][gd] < 90) & (aspcapfield['AK_TARG'][gd] != 0.) &
                  (aspcapfield['AK_TARG'][gd] > -1) & (aspcapfield['AK_TARG'][gd] < 1) &
                  (np.abs(aspcapfield['FPARAM'][gd,0]-teff.cte_ghb(jk0[gd],aspcapfield['FPARAM'][gd,3])[0]) > 1000) )[0]
    #aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('COLORTE_BAD')
    #print('COLORTE_BAD',len(j))
    aspcapfield['ASPCAPFLAG'] &= ~aspcapbitmask.getval('COLORTE_WARN')
    j = np.where( (aspcapfield['H'][gd] < 90) & (aspcapfield['AK_TARG'][gd] != 0.) &
                  (aspcapfield['AK_TARG'][gd] > -1) & (aspcapfield['AK_TARG'][gd] < 1) &
                  (np.abs(aspcapfield['FPARAM'][gd,0]-teff.cte_ghb(jk0[gd],aspcapfield['FPARAM'][gd,3])[0]) > 500) &
                  (aspcapfield['ASPCAPFLAG'][gd] & aspcapbitmask.getval('COLORTE_BAD')==0) )[0]
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('COLORTE_WARN')
    print('COLORTE_WARN',len(j))

    # MULTIPLE_SUSPECT transfer to ASPCAPFLAG (where it will trigger STAR_BAD)
    aspcapfield['ASPCAPFLAG'] &= ~aspcapbitmask.getval('MULTIPLE_SUSPECT')
    j = np.where( (aspcapfield['STARFLAG'][gd] & starbitmask.getval('MULTIPLE_SUSPECT')) >0 )[0]
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('MULTIPLE_SUSPECT')

    # problematic targets
    aspcapfield['ASPCAPFLAG'] &= ~aspcapbitmask.getval('PROBLEM_TARGET')
    j = np.where( (np.core.defchararray.find(aspcapfield['TARGFLAGS'][gd].astype(str),'EXTENDED') >=0) |
                  (np.core.defchararray.find(aspcapfield['TARGFLAGS'][gd].astype(str),'EMBEDDED') >=0) |
                  (np.core.defchararray.find(aspcapfield['TARGFLAGS'][gd].astype(str),'EMISSION') >=0) |
                  (np.core.defchararray.find(aspcapfield['TARGFLAGS'][gd].astype(str),'APOGEE2_M31') >=0) |
                  (np.core.defchararray.find(aspcapfield['TARGFLAGS'][gd].astype(str),'APOGEE2_M33') >=0) |
                  (np.core.defchararray.find(aspcapfield['TARGFLAGS'][gd].astype(str),'APOGEE2_QSO') >=0) ) [0]
    aspcapfield['ASPCAPFLAG'][gd[j]] |= aspcapbitmask.getval('PROBLEM_TARGET')
    print('PROBLEM_TARGET',len(j))

    # star level flag
    aspcapfield['ASPCAPFLAG'] &= ~aspcapbitmask.getval('STAR_BAD')
    j = np.where( aspcapfield['ASPCAPFLAG']&aspcapbitmask.badval() > 0)[0]
    aspcapfield['ASPCAPFLAG'][j] |= aspcapbitmask.getval('STAR_BAD')
    print('STAR_BAD',len(j))
    aspcapfield['ASPCAPFLAG'] &= ~aspcapbitmask.getval('STAR_WARN')
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


def add_members(tab) :
    """ Add membership bitmask
    """
    col = Column(np.full([len(tab)],np.int64(0)),name='MEMBERFLAG',dtype=np.int64)
    try : tab.remove_column('MEMBERFLAG')
    except: pass
    tab.add_column(col)
    col = Column(np.full([len(tab)],''),name='MEMBER',dtype='S10')
    try : tab.remove_column('MEMBER')
    except: pass
    tab.add_column(col)

    membermask=bitmask.MembersBitMask()
    for name in membermask.name :
        j=apselect.clustmember(tab,name,usememberflag=False)
        for jj in j :
            tab['MEMBERFLAG'][jj] |= membermask.getval(name)
        j=apselect.clustmember(tab,name,dsph=True,usememberflag=False)
        for jj in j :
            tab['MEMBERFLAG'][jj] |= membermask.getval(name)

    # populate text MEMBER
    j=np.where(tab['MEMBERFLAG'] != 0)[0]
    for jj in j :
        tab['MEMBER'][jj] = membermask.getname(tab['MEMBERFLAG'][jj])
  
def add_visitpk(allstar, allvisit ) :
    """ Add reference indices from allStar into allVisit file
    """

    # add VISIT_PK column, initialize with out of range indices
    if not isinstance(allstar,astropy.table.table.Table) : allstar=Table(allstar)
    maxvisit=100
    col = Column(np.full([len(allstar),maxvisit],2**31),name='VISIT_PK',dtype=np.int32)
    try : allstar.remove_column('VISIT_PK')
    except: pass
    allstar.add_column(col)

    # get indices in allVisit file for each degree of RA, to shorten search
    ind = np.zeros(360,dtype=int)
    for i in range(1,360) : 
        ind[i] = np.where(allvisit['RA']>i)[0][0] 
        print(i,ind[i])

    # loop over all stars and get corresponding visits
    nvisits=[]
    nmax=0
    for i,star in enumerate(allstar) :
        ira=int(star['RA'])
        i1 = ind[ira]
        if ira < 359 : i2 = ind[ira+1 ]
        else : i2 = len(allvisit)
        j = np.where(allvisit['APOGEE_ID'][i1:i2] == star['APOGEE_ID'])[0]
        allstar['VISIT_PK'][i,0:min([maxvisit,len(j)])] = (i1+j)[0:min([maxvisit,len(j)])]
        print(i,len(j))
        nvisits.append(len(j))
        if len(j) > nmax : nmax = len(j)
    print('nmax: ', nmax)
    return nvisits

def washphotfix(tab):
    """ Add washington photometry that was missing
    """
    wash=ascii.read(os.environ['APOGEE_DIR']+'/data/allwash_dsph_photonly.csv')
    nmiss=0
    for obj in set(wash['APOGEE_ID']) :
        j=np.where((wash['APOGEE_ID'] == obj) & (wash['DDO51'] > 0.) )[0]
        if len(j) > 0 :
            i=np.where(tab['APOGEE_ID'] == obj)[0]
            print(obj,len(i))
            for col,wcol in zip(['WASH_M','WASH_M_ERR','WASH_T2','WASH_T2_ERR','DDO51','DDO51_ERR'],
                                ['M','MERR','T2','T2ERR','DDO51','DDO51ERR'] ) :
                tab[col][i] = wash[wcol][j[0]]
        else : 
            print('no washphot in file for ',obj)
            nmiss+=1
    print('missing: ', nmiss)


def targflagfix(tab):
    """ fix erroneous bit in K2_C12 fields and fix up formatting issues in TARGFLAGS
    """
    j = np.where(np.char.find(tab['FIELD'].astype(str),'K2_C12') >=0)[0]
    targmask=bitmask.Apogee2Target3()
    tab['APOGEE2_TARGET3'][j] &= ~targmask.getval('APOGEE2_K2_MDWARF')

    # repopulate TARGFLAGS
    for i,star in enumerate(tab) :
        newflags = (bitmask.targflags(star['APOGEE_TARGET1'],star['APOGEE_TARGET2'],0,0,survey='apogee')+
                    bitmask.targflags(star['APOGEE2_TARGET1'],star['APOGEE2_TARGET2'],star['APOGEE2_TARGET3'],star['APOGEE2_TARGET4'],survey='apogee2'))
        tab['TARGFLAGS'][i] = newflags

def elemflagfix(tab) :
    """ fix 'NaN' in ELEMFLAG
    """
    n,nel=tab['ELEMFLAG'].shape
    for i in range(nel) :
        j = np.where(tab['ELEMFLAG'][:,i] == -2**63)[0]
        tab['ELEMFLAG'][j,i] = 0

def nanfix(tab) :
    """ change -9999 to NaNs
    """
    for col in ['RA','DEC','AK_TARG','AK_WISE','SFD_EBV','J','J_ERR','H','H_ERR','K','K_ERR',
                'WASH_M','WASH_M_ERR','WASH_T2','WASH_T2_ERR','DDO51','DDO51_ERR',
                'IRAC_3_6','IRAC_3_6_ERR','IRAC_4_5','IRAC_4_5_ERR','IRAC_5_8','IRAC_5_8_ERR',
                'WISE_4_5','WISE_4_5_ERR','TARG_4_5','TARG_4_5_ERR'] :
        j=np.where((tab[col]<-9998))[0]
        print(col,len(j))
        tab[col][j] = np.nan

    # for PMRA, PMDEC, change to 0
    for col in ['TARG_PMRA','TARG_PMDEC'] :
        j=np.where((tab[col]<-9998))[0]
        print(col,len(j))
        tab[col][j] = 0

    # change VESTA RA,DEC
    j=np.where(tab['APOGEE_ID'] == 'VESTA')[0]
    tab['RA'][j] = np.nan
    tab['DEC'][j] = np.nan


def fiber300fix() :
    """ Replace records in aspcapField files that were missing for MEANFIB==300 with
       supplemental runs from redo
    """

    for telescope in ['apo25m','lco25m'] :
        a=fits.open(telescope+'/redo_000/aspcapField-redo_000.fits')
        objs=np.char.split(a[1].data['APOGEE_ID'],'.')
        for i,obj in enumerate(objs) :
            name=telescope+'/'+obj[1]+'/aspcapField-'+obj[1]+'.fits'
            old=fits.open(name)
            try: j=np.where(old[1].data['APOGEE_ID'] == obj[0])[0][0]
            except : j=np.where(np.char.find(old[1].data['APOGEE_ID'],obj[0]) >=0)[0][0]
            print(old[1].data['APOGEE_ID'][j],old[1].data['FPARAM'][j])
            print(a[1].data['APOGEE_ID'][i],a[1].data['FPARAM'][i])
            if old[1].data['APOGEE_ID'][j] != a[1].data['APOGEE_ID'][i].split('.')[0] : pdb.set_trace()
            a[1].data['APOGEE_ID'][i] = obj[0]
            a[1].data['FIELD'][i] = obj[1]
            old[1].data[j] = a[1].data[i]
            old[2].data[j] = a[2].data[i]
            try:os.makedirs('orig/'+os.path.dirname(name))
            except FileExistsError: pass
            # if we've already copied original file, don't do so again
            if not os.path.exists('orig/'+name) : os.rename(name,'orig/'+name)
            else : print('already copied orig')
            hdulist=fits.HDUList()
            hdulist.append(fits.BinTableHDU(old[1].data))
            hdulist.append(fits.BinTableHDU(old[2].data))
            hdulist.append(fits.BinTableHDU(old[3].data))
            hdulist.writeto(name,overwrite=True)


def fix(tab,visit=None) :
    """ Fix up broken things in pre-releases, e.g. dr17alpha
    """
    aspcapmask=bitmask.AspcapBitMask()
    parammask=bitmask.ParamBitMask()
    starmask=bitmask.StarBitMask()
    rvmask=bitmask.RVBitMask()

    # add in SNREV
    try: tab.add_column(Column(tab['SNR'],name='SNREV'))
    except ValueError: print('SNREV already exists...')

    #rename column CLASS to ASPCAP_CLASS
    try: tab.rename_column('CLASS','ASPCAP_GRID')
    except KeyError: print('ASPCAP_GRID already exists')
    try: tab.rename_column('ASPCAP_CLASS','ASPCAP_GRID')
    except KeyError: print('ASPCAP_GRID already exists')
    try: tab.rename_column('FPARAM_CLASS','FPARAM_GRID')
    except KeyError: print('FPARAM_GRID already exists')
    try: tab.rename_column('CHI2_CLASS','FPARAM_GRID')
    except KeyError: print('CHI2_GRID already exists')

    # replace 0. with NaN    
    j=np.where(np.isclose(tab['FPARAM'][:,0],0.))[0]
    tab['FPARAM'][j,:] = np.nan
    tab['FPARAM_COV'][j,:,:] = np.nan
    tab['ASPCAP_CHI2'][j] = np.nan

    # flags for no good visits for RV
    j=np.where(tab['NVISITS'] == 0)[0]
    tab['STARFLAG'][j] |= starmask.getval('RV_FAIL')
    for jj in j : tab['STARFLAGS'][jj] = starmask.getname(tab['STARFLAG'][jj])
    tab['RV_FLAG'][j] |= rvmask.getval('NO_GOOD_VISITS')

    # FERRE_FAIL (mostly fixed with edge issues?)
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

def apstar_id(data,apstar_vers='stars') :
    """ Unique APSTAR identifier
    """
    id = np.core.defchararray.add('apogee.',data['TELESCOPE'].astype(str))
    id = np.core.defchararray.add(id,'.')
    id = np.core.defchararray.add(id,apstar_vers)
    id = np.core.defchararray.add(id,'.')
    id = np.core.defchararray.add(id,data['FIELD'].astype(str))
    id = np.core.defchararray.add(id,'.')
    id = np.core.defchararray.add(id,data['APOGEE_ID'].astype(str))
    return id

def target_id(data) :
    """ Unique target identifier
    """
    id = np.core.defchararray.add(data['TELESCOPE'].astype(str),'.')
    id = np.core.defchararray.add(id,data['FIELD'].astype(str))
    id = np.core.defchararray.add(id,'.')
    id = np.core.defchararray.add(id,data['APOGEE_ID'].astype(str))
    return id

def aspcap_id(data,aspcap_vers='l33') :
    """ Unique ASPCAP identifier
    """
    id = np.core.defchararray.add('apogee.',data['TELESCOPE'].astype(str))
    id = np.core.defchararray.add(id,'.')
    id = np.core.defchararray.add(id,aspcap_vers)
    id = np.core.defchararray.add(id,'.')
    id = np.core.defchararray.add(id,data['FIELD'].astype(str))
    id = np.core.defchararray.add(id,'.')
    id = np.core.defchararray.add(id,data['APOGEE_ID'].astype(str))
    return id

def qa_plots(tab,tab3,prefix='allStar/',hr=True, doqa=True, doelem=True, drcomp='dr16') :
    """ Make a series of QA plots
    """
    procs=[]
    if hr :
        # HR diagrams
        try: os.makedirs(prefix+'/hr/')
        except: pass

        print('making HR diagrams....')
        kw={'a' : tab,'hard' : prefix+'hr/hr.png','xr' : [8000,3000],'grid' : True,
            'iso' : [8.0,9.0,10.0],'alpha' : 0.7,'snrbd' : 5,'target' : prefix+'hr/hr','size' : 1}
        procs.append(Process(target=qa.hr,kwargs=kw))
        kw={'a' :tab, 'hard' : prefix+'hr/hrhot.png','xr' : [20000,3000],'iso' : [8.0,9.0,10.0],'snrbd' : 5,'size' : 1,'alpha' : 0.7}
        procs.append(Process(target=qa.hr,kwargs=kw))
        kw={'a' : tab,'hard' : prefix+'hr/multihr.png','size' : 1}
        procs.append(Process(target=qa.multihr,kwargs=kw))
        kw = {'a' : tab,'hard' : prefix+'hr/hr_cal.png','xr' : [8000,3000],'grid' : True,
              'iso' : [8.0,9.0,10.0],'alpha' : 0.7,'snrbd' : 5,'param' : 'PARAM','target' : prefix+'hr/hr_cal','size' : 1}
        procs.append(Process(target=qa.hr,kwargs=kw))

    if doqa :
        try: os.makedirs(prefix+'/qa/')
        except: pass
        start = time.time()
        kw= {'tab' : tab,'out' : prefix+'qa/'}
        procs.append(Process(target=qa.chi2,kwargs=kw))
        kw = {'tab' : tab,'tab3': tab3,'out' : prefix+'qa/'}
        procs.append(Process(target=qa.apolco,kwargs=kw))
        kw = {'tab' : tab,'tab3' : tab3, 'out' : prefix+'qa/'}
        procs.append(Process(target=qa.flags,kwargs=kw))
        if drcomp is not None:
            kw={'tab' : tab,'tab3' : tab3,'dr' : drcomp,'out' : prefix+'qa/'}
            procs.append(Process(target=qa.drcomp,kwargs=kw))

    if doelem :
        try: os.makedirs(prefix+'/clust/')
        except: pass
        try: os.makedirs(prefix+'/qa/')
        except: pass
        try: os.makedirs(prefix+'/qa_calibrated/')
        except: pass

        # solar neighborhood solar metallicity and cluster star plots
        for xaxis in ['Mabs','logg','Teff'] :
            kw = {'all' : [tab],'names' : [''],'hard' : prefix+'qa/','out' : 'solarz'+'_'+xaxis,
                  'clusters' : ['solar','inner'],'xaxis' : xaxis}
            procs.append(Process(target=cal.allclust,kwargs=kw))
            kw = { 'all' : [tab],'names' : [''],'hard' : prefix+'qa/','xaxis' : xaxis}
            procs.append(Process(target=cal.allclust,kwargs=kw))
        for param in ['hr','rms','chi2','M','Cpar','Npar','alpha','Cpar_Npar','C_N'] :
            row=[]
            for xaxis in ['Mabs','logg','Teff'] :
                fig='solar'+'_'+xaxis+'_'+param+'.png'
        for elem in tab3['ELEM_SYMBOL'][0].astype(str) :
            row=[]
            for xaxis in ['Mabs','logg','Teff'] :
                fig='solar'+'_'+xaxis+'_'+elem+'.png'


        # misc QA plots
        kw= {'tab' : tab, 'tab3' : tab3, 'out' : prefix+'qa/' }
        procs.append(Process(target=qa.calib,kwargs=kw))
        kw= {'tab' : tab,'tab3' : tab3, 'out' : prefix+'qa/' }
        procs.append(Process(target=qa.m67,kwargs=kw))

        # elemental abundances
        kw= {'a' : tab,'tab3' : tab3, 'out' : prefix+'qa/'}
        procs.append(Process(target=qa.plotelems,kwargs=kw))
        kw= {'tab' : tab,'tab3' : tab3, 'out' : prefix+'qa/'}
        procs.append(Process(target=qa.plotelem_errs,kwargs=kw))
        kw= {'a' : tab,'tab3' : tab3, 'calib' : True, 'out' : prefix+'qa_calibrated/'}
        procs.append(Process(target=qa.plotelems,kwargs=kw))
        kw= {'a' : tab,'tab3' : tab3, 'out' : prefix+'qa_calibrated/all_','main' : False}
        procs.append(Process(target=qa.plotelems,kwargs=kw))
        kw= {'a' : tab,'tab3' : tab3, 'out' : prefix+'qa_calibrated/named_','named' : True}
        procs.append(Process(target=qa.plotelems,kwargs=kw))
        kw= {'a' : tab,'out' : prefix+'qa/' }
        procs.append(Process(target=qa.plotcn,kwargs=kw))

    for proc in procs : proc.start()
    for proc in procs : proc.join()

def qa_html(tab,tab3,prefix='allStar/',drcomp='dr16') :
    """ Master HTML pages for qa plots
    """
    try: os.makedirs(prefix+'/qa/')
    except: pass

    grid=[['hr/hr.png','hr/multihr.png','hr/hrhot.png'],
          ['hr/hr_main.png','hr/hr_targ.png',''],
          ['hr/hr_cal_main.png','hr/hr_cal.png','']]

    # Master summary HTML file
    f=html.head(file=prefix+'/index.html')
    f.write(html.table(grid,ytitle=['uncalibrated','uncalibrated','calibrated']))
    f.write('<br>Uncalibrated parameters:<br>')
    ids = ['VESTA','2M14153968+1910558']
    j=[]
    try:
        for id in ids: j.extend( np.where( (np.core.defchararray.strip(tab['APOGEE_ID'].astype(str)) == id) & (tab['VISIT'] == 0)) [0] )
    except:
        for id in ids: j.extend( np.where( (np.core.defchararray.strip(tab['APOGEE_ID'].astype(str)) == id) ) [0] )
    ids = ['VESTA','Arcturus']
    f.write(html.table(tab['FPARAM'][j],plots=False,ytitle=ids,xtitle=aspcap.params()[1]))
    f.write('<br>calibrated parameters:<br>')
    f.write(html.table(tab['PARAM'][j],plots=False,ytitle=ids,xtitle=aspcap.params()[1]))
    # table of abundances (relative to M)
    f.write('<br>Uncalibrated abundances:<br>')
    try: abun=tab['FELEM'][j,0,:]
    except: abun=tab['FELEM'][j,:]
    xtit=[]
    for i in range(len(tab3['ELEM_SYMBOL'][0])) :
        if tab3['ELEMTOH'][0][i] == 1 : abun[:,i]-=tab['FPARAM'][j,3]
        xtit.append('['+tab3['ELEM_SYMBOL'][0].astype(str)[i]+'/M]')
    f.write(html.table(abun,plots=False,ytitle=ids,xtitle=xtit))

    f.write('<br>calibrated abundances:<br>')
    xtit=[]
    for i in range(len(tab3['ELEM_SYMBOL'][0])) :
        xtit.append('['+tab3['ELEM_SYMBOL'][0].astype(str)[i]+'/M]')
    f.write(html.table(tab['X_M'][j],plots=False,ytitle=ids,xtitle=xtit))

    f.write('<p> Calibration relations<ul>\n')
    f.write('<li> Calibration plots\n')
    f.write('<ul>')
    f.write('<li> <a href='+'calib/logg.html> log g </a>\n')
    f.write('<li> <a href='+'calib/teff.html> Teff </a>\n')
    f.write('<li> <a href='+'qa/elem.html> Abundances </a>\n')
    f.write('</ul>')
    f.write('<li> Calibration check (calibration plots from calibrated values)\n')
    f.write('<ul>')
    f.write('<li> <a href='+'calibrated/logg.html> log g </a>\n')
    f.write('<li> <a href='+'calibrated/teff.html> Teff </a>\n')
    f.write('</ul>')
    f.write('<li> <a href='+'qa/calib.html> Calibrated-uncalibrated plots</a>\n')
    f.write('</ul>\n')
    f.write('<p> Comparisons<ul>\n')
    f.write('<li> <a href='+'optical/optical.html> Comparison with optical abundances</a>\n')
    f.write('<li> <a href='+'qa/apolco.html> APO-LCO comparison</a>\n')
    f.write('<li> <a href='+'qa/clust.html> Cluster abundances</a>\n')
    if drcomp is not None : f.write('<li> <a href='+'qa/'+drcomp+'_diffs.html> '+drcomp.upper()+' comparison plots</a>')
    f.write('</ul>\n')
    f.write('<p> Chemistry plots<ul>\n')
    f.write('<li> <a href='+'qa/elem.html> Abundances: solar neighborhood at solar metallicity and main sample</a>\n')
    f.write('<li> <a href='+'qa/elem_chem.html> Chemistry plots with uncalibrated abundances</a>\n')
    f.write('<li> <a href='+'qa_calibrated/elem_chem.html> Chemistry plots with calibrated abundances, main sample</a>\n')
    f.write('<li> <a href='+'qa_calibrated/all_elem_chem.html> Chemistry plots with calibrated abundances, all</a>\n')
    f.write('<li> <a href='+'qa_calibrated/named_elem_chem.html> Chemistry plots with calibrated abundances, named tags</a>\n')
    f.write('</ul>\n')
    f.write('<p> QA checks<ul>\n')
    f.write('<li> <a href='+'qa/chi2.png> CHI2 vs Teff</a>\n')
    f.write('<li> <a href='+'qa/cn.html> C,N parameters vs abundances</a>\n')
    f.write('<li> <a href='+'qa/flags.html> Bitmasks</a>\n')
    f.write('<li> <a href='+'qa/elem_errs.html> Abundance uncertainties</a>\n')
    f.write('</ul>\n')
    f.write('<br> Duplicates/repeats, for empirical uncertainties: <ul>\n')
    f.write('<li> <a href='+'repeat/giant_repeat_elem.html> Elemental abundances, giants</a>\n')
    f.write('<li> <a href='+'repeat/giant_repeat_param.html> Parameters,  giants</a>\n')
    f.write('<li> <a href='+'repeat/dwarf_repeat_elem.html> Elemental abundances, dwarfs</a>\n')
    f.write('<li> <a href='+'repeat/dwarf_repeat_param.html> Parameters,  dwarfs</a>\n')
    f.write('</ul>\n')
    html.tail(f)

    # optical comparison index
    grid=[]
    grid.append(['paramcomp_synspec.png',''])
    yt=['Parameters']
    for el in tab3['ELEM_SYMBOL'][0].astype(str) :
        grid.append(['{:s}_abundcomp_synspec.png'.format(el),''])
        yt.append(el)
    try: os.makedirs(prefix+'/optical/')
    except: pass
    html.htmltab(grid,file=prefix+'optical/optical.html',ytitle=yt,xtitle=['uncalibrated','calibrated'])

    # solar neighborhood
    grid=[]
    yt=[]
    for param in ['hr','rms','chi2','M','Cpar','Npar','alpha','Cpar_Npar','C_N'] :
        row=[]
        for xaxis in ['Mabs','logg','Teff'] :
            fig='solar'+'_'+xaxis+'_'+param+'.png'
            row.append(fig)
        if param == 'hr' : row.extend(['../qa/calib_Teff_hr.png','','','','','','','',''])
        elif param == 'rms': row.extend(['../qa/calib_logg_hr.png','','','','','','','',''])
        elif param == 'M': 
            row.extend(['','M_H_ERR'])
            for i in range(5) : row.append('../qa/M_H_{:d}.png'.format(i))
            row.extend(['../qa/M_H_teff.png','../qa/M_H_dwarfteff.png'])
        elif param == 'Cpar': 
            row.extend(['','C_M_ERR'])
            for i in range(5) : row.append('../qa/C_M_{:d}.png'.format(i))
            row.extend(['../qa/C_M_teff.png','../qa/C_M_dwarfteff.png'])
        elif param == 'Npar': 
            row.extend(['','N_M_ERR'])
            for i in range(5) : row.append('../qa/N_M_{:d}.png'.format(i))
            row.extend(['../qa/N_M_teff.png','../qa/N_M_dwarfteff.png'])
        elif param == 'alpha': 
            row.extend(['','ALPHA_M_ERR'])
            for i in range(5) : row.append('../qa/ALPHA_M_{:d}.png'.format(i))
            row.extend(['../qa/ALPHA_M_teff.png','../qa/ALPHA_M_dwarfteff.png'])
        else :  row.extend(['','','','','','','','',''])
        yt.append(param)
        grid.append(row)
    for elem in aspcap.elems()[0] :
        #header line for this element
        grid.append(['{:s}'.format(elem)])
        yt.append('')
        row=[]
        elemrow=[]
        xt=[]
        elemgrid=[]
        for xaxis in ['Mabs','logg','Teff'] :
            fig='solar'+'_'+xaxis+'_'+elem+'.png'
            row.append(fig)
            elemrow.append(fig)
            xt.append(xaxis)
        row.extend(['../qa/calib_'+elem+'_hr.png','../qa/'+elem+'_hr_err.png'])
        for i in range(5) : row.append('../qa/{:s}_{:d}.png'.format(elem,i))
        row.append('../qa/{:s}_teff.png'.format(elem))
        row.append('../qa/{:s}_dwarfteff.png'.format(elem))
        xt.extend(['calibration','uncertainty','giants 3000-4000','giants 4000-4500','giants 4500-8000','dwarfs 3000-4000','dwarfs 4000-8000','giants f(Teff)','dwarfs f(Teff)'])
        yt.append('<a href={:s}.html>{:s} FPARAM</a>'.format(elem, elem))
        grid.append(row)
        row=['','','','','']
        for i in range(5) : row.append('../qa_calibrated/{:s}_{:d}.png'.format(elem,i))
        row.append('../qa_calibrated/{:s}_teff.png'.format(elem))
        row.append('../qa_calibrated/{:s}_dwarfteff.png'.format(elem))
        grid.append(row)
        yt.append('<a href={:s}.html>{:s} X_M</a>'.format(elem, elem))
        row=['','','','','']
        for i in range(5) : row.append('../qa_calibrated/named_{:s}_{:d}.png'.format(elem,i))
        row.append('../qa_calibrated/named_{:s}_teff.png'.format(elem))
        row.append('../qa_calibrated/named_{:s}_dwarfteff.png'.format(elem))
        grid.append(row)
        yt.append('<a href={:s}.html>{:s} {:s}_FE</a>'.format(elem, elem,elem.upper()))

        # individual element pages
        elemrow.extend(['','','',''])
        elemgrid.append(elemrow)
        elemgrid.append(['../qa/calib_'+elem+'_hr.png','../qa/'+elem+'_hr_err.png','','','','',''])
        for j in range(3) :
            # do uncal, calibrated, and named chemistry plots for indiv elem pages
            elemrow=[]
            if j==0 : 
                dir='qa'
                pre=''
            elif j==1 : 
                dir='qa_calibrated'
                pre=''
            else : 
                dir='qa_calibrated'
                pre='named_'
            for i in range(5) : 
                elemrow.append('../{:s}/{:s}{:s}_{:d}.png'.format(dir,pre,elem,i))
            elemrow.append('../{:s}/{:s}{:s}_teff.png'.format(dir,pre,elem))
            elemrow.append('../{:s}/{:s}{:s}_dwarfteff.png'.format(dir,pre,elem))
            elemgrid.append(elemrow)
        html.htmltab(elemgrid,file=prefix+'qa/{:s}.html'.format(elem))
    html.htmltab(grid,file=prefix+'qa/elem.html',ytitle=yt,xtitle=xt)

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
 
