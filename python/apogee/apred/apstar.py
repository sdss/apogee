# routines for assessing RVs from pipeline

import os
import copy
import glob
import pdb
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import esutil
import pickle
import yaml
from astropy.io import fits
from apogee.utils import apload, applot, bitmask, gaia, spectra
from apogee.aspcap import norm, aspcap
from apogee.apred import wave, sincint, bc, target
from apogee.speclib import lsf
from tools import plots, html, match, struct
from sdss import yanny
from scipy import interpolate
from scipy.signal import correlate
from scipy.signal.windows import tukey
import scipy.linalg
from scipy.ndimage.filters import median_filter, gaussian_filter

import doppler 
import multiprocessing as mp
from astropy.table import Table, Column, vstack

colors=['r','g','b','c','m','y','k']
chips=['a','b','c']

def allField(search=['apo*/*/a?Field-*.fits','apo*/*/a?FieldC-*.fits','lco*/*/a?Field-*.fits'],out='allField.fits',verbose=False) :
    '''
    Concatenate set of apField files
    '''

    if type(search) == str:
        search=[search]
    allfiles=[]
    for path in search :
        allfiles.extend(glob.glob(path))

    a=[]
    for file in allfiles :
        if 'Field-cal_' not in file :
            print(file)
            dat=Table.read(file)
            a.append(dat)
    all =vstack(a)
    del(a)
    all.sort(['RA','DEC','FIELD'])

    # write out the file
    if out is not None:
        print('writing',out)
        all.write(out,overwrite=True)

    return all

def allFieldVisit(search=['apo*/*/a?FieldVisits-*.fits','apo*/*/a?FieldC-*.fits','lco*/*/a?FieldVisits-*.fits'],
                  apred_vers='dr17',out='allFieldVisit.fits',verbose=False) :
    """ Concatenate set of apFieldVisits files
    """

    if type(search) == str:
        search=[search]
    allfiles=[]
    for path in search :
        allfiles.extend(glob.glob(path))

    a=[]
    for file in allfiles :
        if 'Field-cal_' not in file :
            print(file)
            dat=Table.read(file)
            a.append(dat)
    allvisit =vstack(a)
    del(a)
    allvisit.sort(['RA','DEC','FIELD'])

    # add VISIT_ID
    col = Column(visit_id(allvisit,apred_vers=apred_vers),name='VISIT_ID',dtype='S64')
    allvisit.add_column(col,index=2)

    # repopulate STARFLAGS to allow for longer strings
    allvisit.remove_column('STARFLAGS')
    for ind,col in enumerate(allvisit.columns) :
        if col == 'STARFLAG' : break
    col = Column(np.full([len(allvisit)],''),name='STARFLAGS',dtype='S132')
    allvisit.add_column(col,index=ind+1)
    starmask=bitmask.StarBitMask()
    for i,visit in enumerate(allvisit) : allvisit['STARFLAGS'][i] = starmask.getname(visit['STARFLAG'])

    # put APOGEE2 target flags in APOGEE2_
    for ind,col in enumerate(allvisit.columns) :
        if col == 'APOGEE_TARGET4' : break
    for name in ['APOGEE2_TARGET4','APOGEE2_TARGET3','APOGEE2_TARGET2','APOGEE2_TARGET1'] :
        col = Column(np.full([len(allvisit)],0),name=name,dtype=np.int32)
        allvisit.add_column(col,index=ind+1)
    for i,visit in enumerate(allvisit) :
        if visit['SURVEY'].find('apogee2') >=0 :
            for name in ['TARGET1','TARGET2','TARGET3','TARGET4'] :
                allvisit['APOGEE2_'+name][i] = visit['APOGEE_'+name]
                allvisit['APOGEE_'+name][i] = 0
        elif visit['SURVEY'] == 'apo1m' :
            allvisit['APOGEE2_TARGET2'][i] = allvisit['APOGEE_TARGET2'][i]
    allvisit.remove_column('APOGEE_TARGET3')
    allvisit.remove_column('APOGEE_TARGET4')

    # write out the file
    if out is not None:
        print('writing',out)
        allvisit.write(out,overwrite=True)

    return allvisit

def visit_id(data,apred_vers='dr17') :
    """ Unique visit identifier
    """
    id = np.core.defchararray.add('apogee.',data['TELESCOPE'].astype(str))
    id = np.core.defchararray.add(id,'.')
    id = np.core.defchararray.add(id,apred_vers)
    id = np.core.defchararray.add(id,'.')
    id = np.core.defchararray.add(id,np.char.strip(data['PLATE'].astype(str)))
    id = np.core.defchararray.add(id,'.')

    tmp=[]
    for d in data :
        if d['TELESCOPE'] == 'apo1m'  :
            tmp.append(d['FILE'].split('-')[2]+'.'+d['APOGEE_ID'])
        else :
            tmp.append(str(d['MJD'])+'.'+str(d['FIBERID']))
    id = np.core.defchararray.add(id,np.array(tmp))

    return id


def allPlate(all,out=None) :
    """ Create allPlate file from allFieldVisit
    """

    platetype = np.dtype([('PLATE_VISIT_ID','S64'),('LOCATION_ID',int),('PLATE',int),('MJD',int),
                          ('APRED_VERSION','S10'),('NAME','S64'),
                          ('RACEN',float),('DECCEN',float),('RADIUS',float),('SHARED',int),
                          ('FIELD_TYPE',int),('SURVEY','S24'),('PROGRAMNAME','S24'),('PLATERUN','S24'),('CHUNK','S24'),
                          ('HA',(float,6)),('DESIGNID',int),('NSTANDARD',int),('NSCIENCE',int),('NSKY',int),
                          ('PLATEDESIGN_VERSION','S24')])

    # get unique PLATE-MJD
    gd = np.where(all['TELESCOPE'] != 'apo1m') [0]
    plate_mjd = set(
        np.core.defchararray.add(np.core.defchararray.add( np.core.defchararray.strip(all['PLATE'][gd]),'_'),
                                 all['MJD'][gd].astype(str)))

    allplate = np.zeros(len(plate_mjd),dtype=platetype)
    plans = yanny.yanny(os.environ['PLATELIST_DIR']+'/platePlans.par')['PLATEPLANS']
    for i,p_m in enumerate(plate_mjd) :
        plate,mjd = p_m.split('_')
        print(i,plate,mjd)
        allplate['PLATE'][i] = int(plate)
        allplate['MJD'][i] = int(mjd)
        allplate['APRED_VERSION'][i] =  '?'
        iplan = np.where(np.array(plans['plateid']) == int(plate))[0]
        if len(iplan) != 1 :
            print('missing or multiple plate in platePlans!')
            pdb.set_trace()
        else :
            iplan=iplan[0]
        allplate['LOCATION_ID'][i] = plans['locationid'][iplan]
        for key in ['raCen','decCen','survey','programname','platerun','chunk','ha','designid'] :
            allplate[key.upper()][i] = plans[key][iplan]

        # get some information from plateHoles
        holefile='{:s}/plates/{:04d}XX/{:06d}/plateHolesSorted-{:06d}.par'.format(
                 os.environ['PLATELIST_DIR'],int(plate)//100,int(plate),int(plate))
        #holes=yanny.yanny(holefile)
        fp = open(holefile,'r')
        fout = open('tmp.yml','w')
        fout.write('---\n') 
        for line in fp :
            if (line[0:7] == 'typedef') : break
            if len(line.split()) > 0 :
                fout.write(line.split()[0])
                fout.write(' : ')
                for field in line.split()[1:] : fout.write(field)
                fout.write('\n')
        fout.close()
        holes=yaml.safe_load(open('tmp.yml','r'))
        for key in ['standard','science','sky'] :
            try: allplate['N'+key.upper()][i] = holes['napogee_'+key]
            except: allplate['N'+key.upper()][i] = holes['napogee_south_'+key]
        allplate['PLATEDESIGN_VERSION'][i] = holes['platedesign_version']
        try: allplate['RADIUS'][i] = holes['tilerad']
        except: allplate['RADIUS'][i] = 1.49
        os.remove('tmp.yml')

    if out is not None : 
        hdulist = fits.HDUList()
        hdulist.append(fits.BinTableHDU(allplate))
        hdulist[0].header['VERSION'] = (os.environ['APOGEE_VER'],'APOGEE software version APOGEE_VER')
        hdulist.writeto(out,overwrite=True)

    return allplate 

def doppler_rv(planfile,survey='apogee',telescope='apo25m',apred='r13',apstar_vers=None,obj=None,
               nobj=0,threads=8,maxvisit=500,snmin=3,nres=[5,4.25,3.5],rv_reject=10,vmedian=None,
               save=False,clobber=False,rvclobber=False,vcclobber=False,verbose=False,tweak=False,plot=False,windows=None) :
    """ Run DOPPLER RVs for a field
    """ 
  
    plan=yaml.safe_load(open(planfile,'r'))
    if plan['apogee_ver'] != os.environ['APOGEE_VER'] :
        print('apogee_ver {:s} does not match running version {:s}'.format(plan['apogee_ver'],os.environ['APOGEE_VER']))
        pdb.set_trace()

    apred=plan['apred_vers']
    if apstar_vers is None : apstar_vers=plan['apstar_vers'] if plan.get('apstar_vers') else 'stars'

    telescope=plan['telescope']
    field=plan['field']
    rv_reject = plan['rv_reject'] if plan.get('rv_reject') else rv_reject
    snmin = plan['snmin'] if plan.get('snmin') else snmin
    try : vmedian = plan['vmedian'] 
    except KeyError : pass

    if clobber :
        rvclobber= True 
        vcclobber= True 

    # get all the VisitSum files for this field and concatenate them using astropy Table
    files=glob.glob(os.environ['APOGEE_REDUX']+'/'+apred+'/visit/'+telescope+'/'+field+'/apVisitSum*')
    if len(files) == 0 :
        print('no apVisitSum files found for {:s}'.format(field))
        return
    visitsum=[]
    for file in files :
        visittab=Table.read(file)
        target.add_design(visittab)
        visitsum.append(visittab)
    allvisits=vstack(visitsum)

    # sort by JD
    allvisits.sort('JD')

    # strip spaces from APOGEE_ID
    for visit in allvisits : visit['APOGEE_ID'] = visit['APOGEE_ID'].strip()

    # change datatype of STARFLAG to 64-bit
    allvisits['STARFLAG'] = allvisits['STARFLAG'].astype(np.int64)

    # SDSS-V convention
    #for filt in ['J','H','K'] :
    #    allvisits[filt] = allvisits[filt+'MAG']
    #    allvisits[filt+'_ERR'] = allvisits[filt+'ERR']

    # select out gd visits to process
    starmask=bitmask.StarBitMask()
    gd=np.where(((allvisits['STARFLAG'] & starmask.badval()) == 0) & 
                 (allvisits['APOGEE_ID'] != b'') &
                 (allvisits['SNR'] > snmin) )[0]
    print(len(allvisits),len(gd))

    # output directory
    load=apload.ApLoad(apred=apred,telescope=telescope)
    outfield=load.filename('Field',field=field)
    if apstar_vers != 'stars' :
        outfield=outfield.replace('/stars/','/'+apstar_vers+'/')
    try : os.makedirs(os.path.dirname(outfield))
    except FileExistsError: pass
    outfieldvisits=load.filename('FieldVisits',field=field)
    if apstar_vers != 'stars' :
        outfieldvisits=outfieldvisits.replace('/stars/','/'+apstar_vers+'/')

    # get all unique (or requested) objects
    if obj is None :
        if nobj > 0 :
            allobj=set(allvisits['APOGEE_ID'][0:nobj])
        else :
            allobj=set(allvisits['APOGEE_ID'])
    else :
        allobj = obj

    # output apField structur
    fieldtype = np.dtype([('FILE','S64'),('APOGEE_ID','S30'),('TELESCOPE','S6'),('LOCATION_ID',np.int32),('FIELD','S20'),
                          ('ALT_ID','S30'),('RA',float),('DEC',float),('GLON',float),('GLAT',float),
                          ('J',np.float32),('J_ERR',np.float32),('H',np.float32),('H_ERR',np.float32),('K',np.float32),('K_ERR',np.float32),
                          ('SRC_H','S16'),('WASH_M',np.float32),('WASH_M_ERR',np.float32),('WASH_T2',np.float32),('WASH_T2_ERR',np.float32),
                          ('DDO51',np.float32),('DDO51_ERR',np.float32),('IRAC_3_6',np.float32),('IRAC_3_6_ERR',np.float32),
                          ('IRAC_4_5',np.float32),('IRAC_4_5_ERR',np.float32),('IRAC_5_8',np.float32),('IRAC_5_8_ERR',np.float32),
                          ('IRAC_8_0',np.float32),('IRAC_8_0_ERR',np.float32),
                          ('WISE_4_5',np.float32),('WISE_4_5_ERR',np.float32),('TARG_4_5',np.float32),('TARG_4_5_ERR',np.float32),
                          ('WASH_DDO51_GIANT_FLAG',np.int32),('WASH_DDO51_STAR_FLAG',np.int32),
                          ('TARG_PMRA',np.float32),('TARG_PMDEC',np.float32),('TARG_PM_SRC','S16'),
                          ('AK_TARG',np.float32),('AK_TARG_METHOD','S32'),
                          ('AK_WISE',np.float32),('SFD_EBV',np.float32),
                          ('APOGEE_TARGET1',np.int32),('APOGEE_TARGET2',np.int32),
                          ('APOGEE2_TARGET1',np.int32),('APOGEE2_TARGET2',np.int32),('APOGEE2_TARGET3',np.int32),('APOGEE2_TARGET4',np.int32),
                          ('TARGFLAGS','S132'),('SURVEY','S32'),('PROGRAMNAME','S32'),
                          ('NINST',np.int32),('NVISITS',np.int32),('COMBTYPE',np.int32),('COMMISS',np.int32),
                          ('SNR',np.float32),('SNREV',np.float32),
                          ('STARFLAG',np.int64),('STARFLAGS','S132'),('ANDFLAG',np.int64),('ANDFLAGS','S132'),
                          ('VHELIO_AVG',np.float32),('VSCATTER',np.float32),('VERR',np.float32),
                          ('RV_TEFF',np.float32),('RV_LOGG',np.float32),('RV_FEH',np.float32),('RV_ALPHA',np.float32),('RV_CARB',np.float32),
                          ('RV_CHI2',np.float32),('RV_CCFWHM',np.float32),('RV_AUTOFWHM',np.float32),('RV_FLAG',np.int32),
                          ('N_COMPONENTS',np.int32),('MEANFIB',np.float32),('SIGFIB',np.float32),
                          ('MIN_H',np.float32),('MAX_H',np.float32),('MIN_JK',np.float32),('MAX_JK',np.float32)
                         ])
    #default initialize to 0/blank
    allfield = np.zeros(len(allobj),dtype=fieldtype)
    # numpy structured arrays use encoded byte strings. Convert to astropy table allows these to be referenced as strings
    allfield = Table(allfield)

    allfield['TELESCOPE'] = telescope
    allfield['FIELD'] = field
    #initialize some tags to NaN in case RV fails
    for tag in ['VHELIO_AVG', 'RV_TEFF', 'RV_LOGG', 'RV_FEH', 'RV_CHI2' ] : allfield[tag] = np.nan
    allfield['N_COMPONENTS'] = -1

    allfiles=[]
    allv=[]
    nobj=0
    nvisit=0
    pixelmask=bitmask.PixelBitMask()
    rvmask=bitmask.RVBitMask()

    # loop over requested objects, building up allfiles list of 
    #  [(field,obj,clobber,verbose,tweak,plot,windows),filenames....] to pass to dorv()
    for iobj,star in enumerate(sorted(allobj)) :
        allfield['APOGEE_ID'][iobj] = star

        # copy basic information from first visit in case star fails
        visits=np.where(allvisits['APOGEE_ID'] == star)[0]
        visit0=visits[0]
        keys=['RA','DEC','GLON','GLAT','LOCATION_ID','ALT_ID','J','J_ERR','H','H_ERR','K','K_ERR',
              'SRC_H','WASH_M','WASH_M_ERR','WASH_T2','WASH_T2_ERR',
              'DDO51','DDO51_ERR','IRAC_3_6','IRAC_3_6_ERR',
              'IRAC_4_5','IRAC_4_5_ERR','IRAC_5_8','IRAC_5_8_ERR','IRAC_8_0','IRAC_8_0_ERR',
              'WISE_4_5','WISE_4_5_ERR','TARG_4_5','TARG_4_5_ERR',
              'WASH_DDO51_GIANT_FLAG','WASH_DDO51_STAR_FLAG',
              'AK_TARG','AK_TARG_METHOD','AK_WISE','SFD_EBV']
        for key in keys :
            try: allfield[key][iobj] = allvisits[key][visit0]
            except KeyError: pass
        # rename targeting proper motions and fill design columns in case of RV failure
        keys = ['PMRA','PMDEC','PM_SRC']
        for key in keys :
            try: allfield['TARG_'+key][iobj] = allvisits[key][visit0]
            except KeyError: pass
        allfield['MIN_H'] = allvisits['MIN_H'][visits].min()
        allfield['MAX_H'] = allvisits['MAX_H'][visits].max()
        allfield['MIN_JK'] = allvisits['MIN_JK'][visits].min()
        allfield['STARFLAG'] = np.bitwise_or.reduce(allvisits['STARFLAG'][visits])
        allfield['ANDFLAG'] = np.bitwise_and.reduce(allvisits['STARFLAG'][visits])

        # we will only consider good visits for RVs and combination
        visits=np.where(allvisits['APOGEE_ID'][gd] == star)[0]
        print('object: {:}  nvisits: {:d}'.format(star,len(visits)))
        nobj+=1
        nvisit+=len(visits)

        if len(visits) > 0 :
            allfiles.append([allvisits[gd[visits]],load,(field,star,rvclobber,verbose,tweak,plot,windows,apstar_vers,save)])
        else :
            allfield['STARFLAG'][iobj] |= starmask.getval('RV_FAIL')
            allfield['ANDFLAG'][iobj] |= starmask.getval('RV_FAIL')
            allfield['RV_FLAG'][iobj] = rvmask.getval('NO_GOOD_VISITS')
    print('total objects: ', nobj, ' total visits: ', nvisit) 

    # add GAIA information
    allfield=gaia.add_gaia(allfield)

    # now do the RVs, in parallel if requested
    if threads == 0 :
        output=[]
        for speclist in allfiles :
            print(speclist)
            output.append(dorv(speclist))
    else :
        pool = mp.Pool(threads)
        output = pool.map_async(dorv, allfiles).get()
        pool.close()
        pool.join()
    print('done pool')

    # load up the individual visit RV information

    # remove old visit RV tags and initialize new ones
    for col in ['VTYPE','VREL','VRELERR','VHELIO','BC','RV_TEFF','RV_LOGG','RV_FEH','RV_CARB','RV_ALPHA'] :
        try : allvisits.remove_column(col)
        except KeyError: pass
        if col != 'VTYPE' : allvisits[col] = np.nan
        #allvisits.rename_column(col,'EST'+col)
        #if col == 'VTYPE' : allvisits[col] = 0
        #else : allvisits[col] = np.nan
    for col in ['XCORR_VREL','XCORR_VRELERR','XCORR_VHELIO','BC','CCFWHM','AUTOFWHM','RV_CHI2'] :
        allvisits[col] = np.float32(np.nan)

    # add columns for RV components
    allvisits['N_COMPONENTS'] = np.int32(-1)
    rv_components = Column(name='RV_COMPONENTS',dtype=np.float32,shape=(3,),length=len(allvisits))
    allvisits.add_column(rv_components)
    rvtab = Column(name='RVTAB',dtype=Table,length=len(allvisits))
    allvisits.add_column(rvtab)
    allvisits['RV_FLAG'] = np.int32(0)

    # now load the new ones with the dorv() output
    allv=[]
    for out,files in zip(output,allfiles) :
        apogee_id=files[-1][1]
        if len(out) == 3 :
            visits=[]
            ncomponents=0
            mask = out[2]
            # for special targets, set all velocities to average
            if vmedian is not None :
                if type(vmedian) == float : vmed = vmedian
                else : vmed=np.median(out[0][1]['vhelio'])
                out[0][1]['vhelio'] = vmed
                out[0][1]['vrel'] = vmed-out[0][1]['bc']

            for i,(v,g) in enumerate(zip(out[0][1],out[1])) :
                # match by filename components in case there was an error reading in doppler
                name=os.path.basename(v['filename']).replace('.fits','').split('-')
                if telescope == 'apo1m' :
                    visit = np.where( np.char.strip(allvisits['FILE']).astype(str) == os.path.basename(v['filename'].strip()) )[0]
                    if len(visit) == 0 :
                        # special case for incremental release...yuck
                        visit = np.where( np.char.strip(allvisits['FILE']).astype(str) == 
                                    os.path.basename(v['filename'].strip()).replace('-r13-','-r12-') )[0]
                else :
                    visit = np.where( (np.char.strip(allvisits['PLATE']).astype(str) == name[-3]) &
                                      (allvisits['MJD'] == int(name[-2])) &
                                      (allvisits['FIBERID'] == int(name[-1])) )[0]
                if len(visit) > 0 : visit=visit[0]
                else : continue
                visits.append(visit)
                allvisits[visit]['VREL']=v['vrel']
                allvisits[visit]['VRELERR']=v['vrelerr']
                allvisits[visit]['VHELIO']=v['vhelio']
                allvisits[visit]['XCORR_VREL']=v['xcorr_vrel']
                allvisits[visit]['XCORR_VRELERR']=v['xcorr_vrelerr']
                allvisits[visit]['XCORR_VHELIO']=v['xcorr_vhelio']
                allvisits[visit]['BC']=v['bc']
                allvisits[visit]['CCFWHM']=v['ccpfwhm']
                allvisits[visit]['AUTOFWHM']=v['autofwhm']
                if allvisits[visit]['CCFWHM'] > 2.0*allvisits[visit]['AUTOFWHM']:
                    allvisits[visit]['STARFLAG'] |= starmask.getval('SUSPECT_ROTATION')
                if allvisits[visit]['AUTOFWHM'] > 300 :
                    allvisits[visit]['STARFLAG'] |= starmask.getval('SUSPECT_BROAD_LINES')
                allvisits[visit]['RV_TEFF']=v['teff']
                allvisits[visit]['RV_LOGG']=v['logg']
                allvisits[visit]['RV_FEH']=v['feh']
                allvisits[visit]['RV_CHI2']=v['chisq']
                if g is None or allvisits[visit]['SNR'] < 10 : allvisits[visit]['N_COMPONENTS']=0
                else : allvisits[visit]['N_COMPONENTS']=g['N_components']
                if allvisits[visit]['N_COMPONENTS'] > 1 :
                    allvisits[visit]['STARFLAG'] |= starmask.getval('MULTIPLE_SUSPECT')
                    n=len(g['best_fit_parameters'])//3
                    gd=np.where(np.array(g['best_fit_parameters'])[0:n] > 0)[0]
                    rv_comp = np.array(g['best_fit_parameters'])[2*n+gd]
                    n_rv_comp = np.min([3,len(rv_comp)])
                    allvisits[visit]['RV_COMPONENTS'][0:n_rv_comp] = rv_comp[0:n_rv_comp]
                allvisits[visit]['RVTAB'] = v
                # flag visits with suspect RVs
                if allvisits[visit]['RV_TEFF'] < 6000 : bd_diff = rv_reject
                else : bd_diff = 5*rv_reject
                if (np.abs(allvisits[visit]['VHELIO']-allvisits[visit]['XCORR_VHELIO']) > bd_diff) :
                    allvisits[visit]['STARFLAG'] |= starmask.getval('RV_REJECT')
                    allvisits[visit]['RV_FLAG'] |= rvmask.getval('RV_REJECT')
                elif not np.isclose(allvisits[visit]['VHELIO'],allvisits[visit]['XCORR_VHELIO']) :
                    allvisits[visit]['STARFLAG'] |= starmask.getval('RV_SUSPECT')
                    allvisits[visit]['RV_FLAG'] |= rvmask.getval('RV_SUSPECT')
                allvisits[visit]['STARFLAGS'] = starmask.getname(allvisits[visit]['STARFLAG'])
                allvisits[visit]['RV_FLAG'] |= mask

            if len(visits) > 0 :
                visits=np.array(visits)
                # set up visit combination, removing visits with suspect RVs
                gdrv = np.where((allvisits[visits]['STARFLAG'] & starmask.getval('RV_REJECT')) == 0)[0]
                if len(gdrv) > 0 : 
                    allv.append([allvisits[visits[gdrv]],load,(field,apogee_id,vcclobber,apstar_vers,nres,save)])
                else :
                    bd = np.where(allfield['APOGEE_ID'] == apogee_id)[0]
                    if len(bd) > 0 : 
                        allfield['STARFLAG'][bd] |= starmask.getval('RV_FAIL') 
                        allfield['ANDFLAG'][bd] |= starmask.getval('RV_FAIL') 
                        allfield['STARFLAG'][bd] |= starmask.getval('RV_REJECT') 
                        allfield['ANDFLAG'][bd] |= starmask.getval('RV_REJECT') 
                        allfield['RV_FLAG'][bd] |= mask
                        allfield['RV_FLAG'][bd] |= rvmask.getval('ALL_VISITS_REJECTED')
        else :
            bdvisits = np.where(allvisits['APOGEE_ID'] == apogee_id)[0]
            starflag,rvflag = np.int64(0),np.int32(0)
            andflag = allvisits['STARFLAG'][bdvisits[0]]
            for bd in bdvisits :
                allvisits['STARFLAG'][bd] |= starmask.getval('RV_FAIL') 
                allvisits['RV_FLAG'][bd] |= out[0]
                starflag |= allvisits['STARFLAG'][bd]
                andflag &= allvisits['STARFLAG'][bd]
                rvflag |= allvisits['RV_FLAG'][bd]
            bd = np.where(allfield['APOGEE_ID'] == apogee_id)[0]
            if len(bd) > 0 : 
                allfield['ANDFLAG'][bd] = andflag 
                allfield['STARFLAG'][bd] = starflag 
                allfield['RV_FLAG'][bd] = rvflag

    # do the visit combination, in parallel if requested
    if threads == 0 :
        output=[]
        for v in allv :
            output.append(dovisitcomb(v))
    else :
        pool = mp.Pool(threads)
        output = pool.map_async(dovisitcomb, allv).get()
        pool.close()
        pool.join()
    print('done visitcomb pool pool',len(output),len(allv))

    # now load the combined star information into allfield structure
    # note that dovisitcomb() returns an apstar structure, with header
    # information in FITS header, which limits card names to 8 characters
    # Some of these are renamed in allField structure to use different
    # (longer, more clear) names
    for apstar,v in zip(output,allv) :
        j = np.where(allfield['APOGEE_ID'] == v[-1][1])[0]
        # basic target information loaded above for all targets
        if allfield['APOGEE_ID'][j[0]] != apstar.header['OBJID'] or \
           allfield['APOGEE_ID'][j[0]] != v[-1][1] :
            print("IDs don't match: ",allfield['APOGEE_ID'][j],apstar.header['OBJID'],v[-1][1])
            pdb.set_trace()
        #try: allfield['APOGEE_ID'][j] = apstar.header['OBJID']
        #except: allfield['APOGEE_ID'][j] = v[-1][1]

        # targeting flags have different names and are combined from visits
        apogee_target1 = apstar.header['APTARG1']
        apogee_target2 = apstar.header['APTARG2']
        apogee2_target1 = apstar.header['AP2TARG1']
        apogee2_target2 = apstar.header['AP2TARG2']
        apogee2_target3 = apstar.header['AP2TARG3']
        apogee2_target4 = apstar.header['AP2TARG4']
        allfield['APOGEE_TARGET1'][j] = apogee_target1
        allfield['APOGEE_TARGET2'][j] = apogee_target2
        allfield['APOGEE2_TARGET1'][j] = apogee2_target1
        allfield['APOGEE2_TARGET2'][j] = apogee2_target2
        allfield['APOGEE2_TARGET3'][j] = apogee2_target3
        allfield['APOGEE2_TARGET4'][j] = apogee2_target4
        # add character string for target flags
        allfield['TARGFLAGS'][j] = (bitmask.targflags(apogee_target1,apogee_target2,0,0,survey='apogee')+
                                    bitmask.targflags(apogee2_target1,apogee2_target2,apogee2_target3,apogee2_target4,survey='apogee2'))
        # some modified names
        print(j,v[-1][1],allfield['APOGEE_ID'][j],apstar.header['N_COMP'])
        allfield['N_COMPONENTS'][j] = apstar.header['N_COMP']
        allfield['VHELIO_AVG'][j] = apstar.header['VHELIO']
        allfield['RV_CHI2'][j] = apstar.header['RV_CHI2']
        allfield['RV_CCFWHM'][j] = apstar.header['CCFWHM']
        allfield['RV_AUTOFWHM'][j] = apstar.header['AUTOFWHM']
        if allfield['RV_CCFWHM'][j] > 2.0*allfield['RV_AUTOFWHM'][j] :
            allfield['STARFLAG'][j] |= starmask.getval('SUSPECT_ROTATION')
        if allfield['RV_AUTOFWHM'][j] > 300 :
            allfield['STARFLAG'][j] |= starmask.getval('SUSPECT_BROAD_LINES')

        # mostly unmodified names
        for key in ['STARFLAG','ANDFLAG','SNR','SNREV','VSCATTER','VERR','RV_TEFF','RV_LOGG','RV_FEH','NVISITS','MEANFIB','SIGFIB',
                    'MIN_H','MAX_H','MIN_JK','MAX_JK' ] :
            allfield[key][j] = apstar.header[key]
        # add character string for star flags
        allfield['STARFLAGS'][j] = starmask.getname(allfield['STARFLAG'][j])
        allfield['ANDFLAGS'][j] = starmask.getname(allfield['ANDFLAG'][j])

        # tags that are not from apStar
        allfield['SURVEY'][j] =  ','.join(set(v[0]['SURVEY']))
        allfield['PROGRAMNAME'][j] = ','.join(set(v[0]['PROGRAMNAME']))

        # apStar file name
        outfile=load.filename('Star',field=apstar.header['FIELD'],obj=apstar.header['OBJID'])
        allfield['FILE'][j] =  os.path.basename(outfile)

    #populate character string flags for ALL targets, including failed ones
    for star in allfield :
        star['STARFLAGS'] = starmask.getname(star['STARFLAG'])
        star['ANDFLAGS'] = starmask.getname(star['ANDFLAG'])
        star['TARGFLAGS'] = (bitmask.targflags(star['APOGEE_TARGET1'],star['APOGEE_TARGET2'],0,0,survey='apogee')+
                             bitmask.targflags(star['APOGEE2_TARGET1'],star['APOGEE2_TARGET2'],star['APOGEE2_TARGET3'],star['APOGEE2_TARGET4'],survey='apogee2'))
    

    #output apField and apFieldVisits
    hdulist=fits.HDUList()
    hdulist.append(fits.table_to_hdu(Table(allfield)))
    hdulist[0].header['VERSION'] = (os.environ['APOGEE_VER'],'APOGEE software version APOGEE_VER')
    hdulist.writeto(outfield,overwrite=True)

    hdulist=fits.HDUList()
    allvisits.remove_column('RVTAB')
    hdulist.append(fits.table_to_hdu(allvisits))
    hdulist[0].header['VERSION'] = (os.environ['APOGEE_VER'],'APOGEE software version APOGEE_VER')
    hdulist.writeto(outfieldvisits,overwrite=True)

    # make web page
    if obj is not None : suffix='_obj'
    else : suffix=''
    if tweak: suffix=suffix+'_tweak'
    print('making HTML page ....')
    mkhtml(field,suffix=suffix,apred=apred,telescope=telescope,apstar_vers=apstar_vers)

    return allfield,allvisits

def dorv(visitfiles) :            
    """ do the rv jointfit from list of files
    """
    # last list elements has configuration variables in a tuple
    allvisit = visitfiles[0]
    load = visitfiles[1]
    field=visitfiles[-1][0]
    obj=visitfiles[-1][1]
    clobber=visitfiles[-1][2]
    verbose=visitfiles[-1][3]
    tweak=visitfiles[-1][4]
    plot=visitfiles[-1][5]
    windows=visitfiles[-1][6]
    apstar_vers=visitfiles[-1][7]
    save=visitfiles[-1][8]
    #rvrange=visitfiles[-1][7]
    if tweak: suffix='_tweak'
    else : suffix='_out'
    outdir = os.path.dirname(load.filename('Star',field=field,obj=obj))
    if apstar_vers != 'stars' :
        outdir=outdir.replace('/stars/','/'+apstar_vers+'/')

    # if we have saved result and not clobber, return previous result
    if os.path.exists(outdir+'/'+obj+suffix+'.pkl') and not clobber:
        print(obj,' already done')
        fp=open(outdir+'/'+obj+suffix+'.pkl','rb')
        try: 
            out=pickle.load(fp)
            fp.close()
            return out
        except: 
            print('error loading: ', obj+suffix+'.pkl')
            pass

    speclist=[]
    pixelmask=bitmask.PixelBitMask()
    badval=pixelmask.badval()|pixelmask.getval('SIG_SKYLINE')|pixelmask.getval('LITTROW_GHOST')
    rvmask=bitmask.RVBitMask()
    rvmaskval=0
   
    # if we have a significant number of low S/N visits, combine first using
    #    barycentric correction only, use that to get an estimate of systemic
    #    velocity, then do RV determination restricting RVs to within 50 km/s
    #    of estimate. This seems to help significant for faint visits
    lowsnr_visits=np.where(allvisit['SNR']<10)[0]
    if (len(lowsnr_visits) > 1) & (len(lowsnr_visits)/len(allvisit) > 0.1) :
        try :
            print('running BC combination and jointfit for :',obj)
            apstar_bc=visitcomb(allvisit,bconly=True,load=load,write=False,dorvfit=False,apstar_vers=apstar_vers) 
            apstar_bc.setmask(badval)
            spec=doppler.Spec1D(apstar_bc.flux[0,:],err=apstar_bc.err[0,:],bitmask=apstar_bc.bitmask[0,:],
                 mask=apstar_bc.mask[0,:],wave=apstar_bc.wave,lsfpars=np.array([0]),
                 lsfsigma=apstar_bc.wave/22500/2.354,instrument='APOGEE',
                 filename=apstar_bc.filename)
            out= doppler.rv.jointfit([spec],verbose=verbose,plot=plot,tweak=tweak,maxvel=[-500,500])
            rvrange=[out[1][0]['vrel']-50,out[1][0]['vrel']+50]
            rvmaskval |= rvmask.getval('RV_BCFIT' )
        except :
            print('  BC jointfit failed')
            rvrange=[-500,500]
            rvmaskval |= rvmask.getval('RV_BCFIT_FAIL' )
    elif allvisit['H'].max() > 13.5 : 
        # if it's faint, restrict to +/- 500 km/s
        rvrange=[-500,500]
        rvmaskval |= rvmask.getval('RV_FAINT_FIT')
    else :
        # otherwise, restrict to +/ 1000 km/s
        rvrange=[-1000,1000]

    for i in range(len(allvisit)) :

        # load all of the visits into doppler Spec1D objects
        if load.telescope == 'apo1m' :
            visitfile= load.filename('Visit',plate=allvisit['PLATE'][i],
                                 mjd=allvisit['MJD'][i],reduction=allvisit['APOGEE_ID'][i])
            # use FILE tag in case we have MJDFRAC
            visitfile=os.path.dirname(visitfile)+'/'+allvisit['FILE'][i]
        else :
            visitfile= load.filename('Visit',plate=int(allvisit['PLATE'][i]),
                                 mjd=allvisit['MJD'][i],fiber=allvisit['FIBERID'][i])
        spec=doppler.read(visitfile,badval=badval)

        if windows is not None :
            # if we have spectral windows to mask, do so here
            for ichip in range(3) :
                mask = np.full_like(spec.mask[:,ichip],True)
                gd = []
                for window in windows :
                    gd.extend(np.where((spec.wave[:,ichip] > window[0]) & (spec.wave[:,ichip] < window[1]))[0])
                mask[gd] = False
                spec.mask[:,ichip] |= mask
                rvmaskval |= rvmask.getval('RV_WINDOW_MASK')

        if spec is not None : speclist.append(spec)

    # now do the doppler jointfit to get RVs
    # dump empty pickle to stand in case of failure (to prevent redo if not clobber)
    try:
        # dump empty pickle to stand in case of failure (to prevent redo if not clobber)
        if save :
            fp=open(outdir+'/'+obj+suffix+'.pkl','wb')
            pickle.dump(None,fp)
            fp.close()
        print('running jointfit for : {:s}  rvrange:[{:.1f},{:.1f}]  nvisits: {:d}'.format(obj,*rvrange,len(speclist)))
        out= doppler.rv.jointfit(speclist,maxvel=rvrange,verbose=verbose,
                                 plot=plot,saveplot=plot,outdir=outdir+'/',tweak=tweak)
        print('running decomp for :',obj)
        gout = gauss_decomp(out[1],phase='two',filt=True)
        if save :
            fp=open(outdir+'/'+obj+suffix+'.pkl','wb')
            pickle.dump([out,gout,rvmaskval],fp)
            fp.close()
        print('running plots for :',obj,outdir)
        try : os.makedirs(outdir+'/plots/')
        except : pass
        dop_plot(outdir+'/plots/',obj,out,decomp=gout)
    except KeyboardInterrupt : 
        raise
    except ValueError as err:
        print('Exception raised in dorv for: ', field, obj)
        print("ValueError: {0}".format(err))
        rvmaskval |= rvmask.getval('RV_VALUE_ERROR')
        return [rvmaskval]
    except RuntimeError as err:
        print('Exception raised in dorv for: ', field, obj)
        print("Runtime error: {0}".format(err))
        rvmaskval |= rvmask.getval('RV_RUNTIME_ERROR')
        return [rvmaskval]
    except :
        raise
        print('Exception raised in dorv for: ', field, obj)
        rvmaskval |= rvmask.getval('RV_ERROR')
        return [rvmaskval]

    # return summary RV info, visit RV info, decomp info 
    return [out[0:2],gout,rvmaskval]

def dovisitcomb(allv) :
    """ Routine to do visit combination in parallel
    """
    allvisits = allv[0]
    load = allv[1]
    field = allv[2][0]
    apogee_id = allv[2][1]
    clobber = allv[2][2]
    apstar_vers = allv[2][3]
    nres = allv[2][4]
    save = allv[2][5]
    pixelmask=bitmask.PixelBitMask()

    # already done?
    outdir=os.path.dirname(load.filename('Field',field=field))
    if apstar_vers != 'stars' :
        outdir=outdir.replace('/stars/','/'+apstar_vers+'/')
    if os.path.exists(outdir+'/'+apogee_id+'.pkl') and not clobber:
        print(apogee_id,' already done visitcomb')
        fp=open(outdir+'/'+apogee_id+'.pkl','rb')
        try: 
            out=pickle.load(fp)
            fp.close()
            return out
        except: 
            print('error loading: ', apogee_id+'.pkl')
            pass

    # do the combination
    apstar=visitcomb(allvisits,load=load,plot=False,apstar_vers=apstar_vers,nres=nres)

    # dump
    if save: pickle.dump(apstar,open(outdir+'/'+apogee_id+'.pkl','wb'))

    return apstar

def gaussian(amp, fwhm, mean):
    """ Gaussian as defined by gausspy
    """
    return lambda x: amp * np.exp(-4. * np.log(2) * (x-mean)**2 / fwhm**2)

import gausspy.gp as gp

def gauss_decomp(out,phase='one',alpha1=0.5,alpha2=1.5,thresh=[4,4],plot=None,filt=False) :
    """ Do Gaussian decomposition of CCF using gausspy

        Parameters:
        out : list of dictionaries for each frame, giving x_ccf, ccf, and ccferr
        phase : gausspy paramater
        alpha1 : gausspy parameter
        alpha2 : gausspy parameter for second set of gaussians if phase=='two'
        thresh : gausspy parameter
        plot (str) : if not None, do plot and use as root file name for plot
        filt (bool) : if true, apply filtering to remove components judged to be insignificant
    """
    g = gp.GaussianDecomposer()
    g.set('phase',phase)
    g.set('SNR_thresh',thresh)
    g.set('alpha1',alpha1)
    g.set('alpha2',alpha2)
    gout=[]
    if plot is not None : fig,ax=plots.multi(1,len(out),hspace=0.001,figsize=(6,2+n))
    for i,final in enumerate(out) :
        gd=np.where(np.isfinite(final['x_ccf']))[0]
        x=final['x_ccf'][gd]
        y=final['ccf'][gd] 
        # high pass filter for better performance
        if filt : final['ccf'][gd]-= gaussian_filter(final['ccf'][gd],50,mode='nearest')
        try : 
            decomp=g.decompose(x,final['ccf'][gd],final['ccferr'][gd])
            n=decomp['N_components']
        except :
            print('Exception in Gaussian decomposition, setting to 0 components')
            n=0
            decomp=None
        if filt and n>0 :
            # remove components if they are within width of brighter component, or <0.25 peak ,
            #   or more than twice as wide, or if primary component is wide
            for j in range(1,n) :
                # decompose can return component with 0 amplitude, need to remove these first
                # also remove narrow peaks 
                pars_j = decomp['best_fit_parameters'][j::n]
                if np.isclose(pars_j[0],0.) or np.abs(pars_j[1])< 1.:
                    decomp['best_fit_parameters'][j] = 0.
                    decomp['N_components'] -= 1
            for j in range(1,n) :
                pars_j = decomp['best_fit_parameters'][j::n]
                for k in range(j) :
                    pars_k = decomp['best_fit_parameters'][k::n]
                    #reject likely spurious components
                    if (pars_j[0]>pars_k[0] and pars_k[0]>0 and 
                                #closer than width
                                (abs(pars_j[2]-pars_k[2])<abs(pars_j[1])  or 
                                #peak less than 0.5*primary peak
                                 pars_k[0]<0.25*pars_j[0] or 
                                #primary peak less than 0.15
                                 pars_j[0]<0.15 or 
                                #broad primary
                                 abs(pars_j[1])>100 or
                                #broad secondary
                                 np.abs(pars_k[1])>2*np.abs(pars_j[1]) ) ) :
                        decomp['best_fit_parameters'][k] = 0
                        decomp['N_components'] -= 1
                    elif (pars_k[0]>pars_j[0] and pars_j[0]>0 and
                                (abs(pars_j[2]-pars_k[2])<abs(pars_k[1]) or 
                                 pars_j[0]<0.25*pars_k[0] or 
                                 pars_k[0]<0.15 or
                                 abs(pars_k[1])>100 or
                                 np.abs(pars_j[1])>2*np.abs(pars_k[1]) ) )  :
                        decomp['best_fit_parameters'][j] = 0
                        pars_j = decomp['best_fit_parameters'][j::n]
                        decomp['N_components'] -= 1
                  
        gout.append(decomp)
        if plot is not None:
            plots.plotl(ax[i],final['x_ccf'],final['ccf'])
            ax[i].plot(final['x_ccf'],final['ccferr'],color='r')
            for j in range(n) :
                pars=gout[i]['best_fit_parameters'][j::n]
                ax[i].plot(x,gaussian(*pars)(x))
                if pars[0] > 0 : color='k'
                else : color='r'
                ax[i].text(0.1,0.8-j*0.1,'{:8.2f}{:8.1f}{:8.1f}'.format(*pars),transform=ax[i].transAxes,color=color)
            fig.savefig(plot+'_ccf.png')
    del g
    return gout

def dop_plot(outdir,obj,out,decomp=None) :
    """ RV diagnostic plots
    """
    matplotlib.use('Agg')
    n = len(out[2])
    #plot final spectra and final models
    # full spectrum
    fig,ax=plots.multi(1,n,hspace=0.001,figsize=(8,2+n))
    ax=np.atleast_1d(ax)
    # continuum
    figc,axc=plots.multi(1,n,hspace=0.001,figsize=(8,2+n))
    axc=np.atleast_1d(axc)
    # windows
    windows=[[15700,15780],[15850,16000],[16700,16930]]
    fig2,ax2=plots.multi(len(windows),n,hspace=0.001,wspace=0.001,figsize=(12,2+n))
    ax2=np.atleast_2d(ax2)

    # loop over visits
    for i,(mod,spec) in enumerate(zip(out[2],out[3])) :
        ax[i].plot(spec.wave,spec.flux,color='k')
        for iorder in range(3) :
            gd = np.where(~spec.mask[:,iorder])[0]
            ax[i].plot(spec.wave[gd,iorder],spec.flux[gd,iorder],color='g')
        ax[i].plot(mod.wave,mod.flux,color='r')
        ax[i].text(0.1,0.1,'{:d}'.format(spec.head['MJD5']),transform=ax[i].transAxes)
        for iwind,wind in enumerate(windows) :
            ax2[i,iwind].plot(spec.wave,spec.flux,color='k')
            for iorder in range(3) :
                gd = np.where(~spec.mask[:,iorder])[0]
                ax2[i,iwind].plot(spec.wave[gd,iorder],spec.flux[gd,iorder],color='g')
            ax2[i,iwind].plot(mod.wave,mod.flux,color='r')
            ax2[i,iwind].set_xlim(wind[0],wind[1])
            ax2[i,iwind].set_ylim(0.5,1.3)
            if iwind == 0 : ax2[i,iwind].text(0.1,0.1,'{:d}'.format(spec.head['MJD5']),transform=ax2[i,0].transAxes)
        axc[i].plot(spec.wave,spec.flux*spec.cont,color='k')
        axc[i].plot(spec.wave,spec.cont,color='g')
        axc[i].text(0.1,0.1,'{:d}'.format(spec.head['MJD5']),transform=axc[i].transAxes)
    fig.savefig(outdir+'/'+obj+'_spec.png')
    plt.close()
    fig2.savefig(outdir+'/'+obj+'_spec2.png')
    plt.close()
    figc.savefig(outdir+'/'+obj+'_cont.png')
    plt.close()

    # plot cross correlation functions with final model
    fig,ax=plots.multi(1,n,hspace=0.001,figsize=(6,2+n))
    ax=np.atleast_1d(ax)
    vmed=np.median(out[1]['vrel'])
    for i,(final,spec) in enumerate(zip(out[1],out[3])) :
        ax[i].plot(final['x_ccf'],final['ccf'],color='k')
        ax[i].plot(final['x_ccf'],final['ccferr'],color='r')
        ax[i].plot([final['vhelio'],final['vhelio']],ax[i].get_ylim(),color='g',label='fit vhelio')
        ax[i].plot([final['xcorr_vhelio'],final['xcorr_vhelio']],ax[i].get_ylim(),color='r',label='xcorr vhelio')
        ax[i].text(0.1,0.9,'{:d}'.format(spec.head['MJD5']),transform=ax[i].transAxes)
        ax[i].set_xlim(vmed-200,vmed+200)
        ax[i].legend()
        if decomp is not None :
            try: n=decomp[i]['N_components']
            except: n=0
            if n>0 : n=len(decomp[i]['best_fit_parameters'])//3
            x=final['x_ccf']
            for j in range(n) :
                pars=decomp[i]['best_fit_parameters'][j::n]
                ax[i].plot(x,gaussian(*pars)(x))
                if pars[0] > 0 : color='k'
                else : color='r'
                ax[i].text(0.1,0.8-j*0.1,'{:8.2f}{:8.1f}{:8.1f}'.format(*pars),transform=ax[i].transAxes,color=color)
    fig.savefig(outdir+'/'+obj+'_ccf.png')
    plt.close()

def mkhtml(field,suffix='',apred='r13',telescope='apo25m',apstar_vers='stars') :
    """ Make web pages with tables/plots of RV output
        c.f., Doppler vs DR16
    """

    starmask=bitmask.StarBitMask()
    rvmask=bitmask.RVBitMask()
    # get new RV results
    load=apload.ApLoad(apred=apred,telescope=telescope)
    #apf=load.apField(field)[1].data
    infile=load.filename('Field',field=field)
    if apstar_vers != 'stars' :
        infile=infile.replace('/stars/','/'+apstar_vers+'/')
    apf=fits.open(infile)[1].data

    infile=load.filename('FieldVisits',field=field)
    if apstar_vers != 'stars' :
        infile=infile.replace('/stars/','/'+apstar_vers+'/')
    #apfv=load.apFieldVisits(field)[1].data
    apfv=fits.open(infile)[1].data

    outdir=os.path.dirname(infile)
    try: os.makedirs(outdir+'/plots/')
    except: pass
 
    # get old DR16/IDL results
    dr16 = apload.ApLoad(apred='r12',telescope=telescope)
    try :
        apfieldvisits = dr16.apFieldVisits(field)[1].data
        apfield = dr16.apField(field)[1].data
        doapfield = True
    except : 
        print('No apField files found ...')
        doapfield = False

    # match
    if doapfield: i1,i2 = match.match(apfv['FILE'],apfieldvisits['FILE'])
    fig,ax=plots.multi(1,3,figsize=(12,6),hspace=0.5)
    ax[0].hist(apf['VHELIO_AVG'],bins=np.arange(-600,600,5),label='doppler',color='g',histtype='step')
    if doapfield: ax[0].hist(apfield['VHELIO_AVG'],bins=np.arange(-600,600,5),label='DR16',color='r',histtype='step')
    ax[0].legend()
    ax[0].set_xlabel('VHELIO_AVG')
    ax[1].hist(apf['VSCATTER'],bins=np.arange(0,5,0.02),label='doppler',color='g',histtype='step')
    if doapfield: ax[1].hist(apfield['VSCATTER'],bins=np.arange(0,5,0.02),label='DR16',color='r',histtype='step')
    ax[1].legend()
    ax[1].set_xlabel('VSCATTER')
    ax[2].hist(apf['SNR'],bins=np.arange(0,200,5),label='doppler',color='g',histtype='step')
    if doapfield: ax[2].hist(apfield['SNR'],bins=np.arange(0,200,5),label='DR16',color='r',histtype='step')
    ax[2].legend()
    ax[2].set_xlabel('SNR')
    fig.savefig(outdir+'/plots/'+field+'_rvhist.png')

    # create HTML and loop over objects
    fp=open(outdir+'/'+field+suffix+'.html','w')
    fp.write('<HTML>\n')
    fp.write('<HEAD><script type=text/javascript src=../../../html/sorttable.js></script></head>')
    fp.write('<BODY>\n')
    fp.write('<H2> Field: {:s}</H2><p>\n'.format(field))
    fp.write('<A HREF=plots/{:s}_rvhist.png> <IMG SRC=plots/{:s}_rvhist.png> </A>'.format(field,field))
   
    fp.write('<BR>Click on column headers to sort by column value<BR>\n') 
    fp.write('<TABLE BORDER=2 CLASS=sortable>\n')
    fp.write('<TR><TD>Obj<TD>H<TD>RV_TEFF<TD>RV_CHI2<TD>N_comp<TD>Combined spectrum<TD>RV plot'+
             '<TD>CCFs<TD>Combined CCF/autoccf <TD> Spectrum windows<TD>Spectrum<TD> continuum\n')
    for star in apf :
        obj=star['APOGEE_ID']
        print(obj)

        # get visits in Doppler allvisit table
        j=np.where(apfv['APOGEE_ID'] == obj)[0]
        if len(j) == 0 : 
            print('missing {:s} in apfv'.format(obj))
            continue

        # get object in apField
        try: k=np.where(apfield['APOGEE_ID'] == obj)[0][0]
        except: k=-1
        if doapfield :jj=np.where(apfieldvisits['APOGEE_ID'] == obj)[0] 

        # star information
        if star['TARGFLAGS'].find('TELLURIC') >=0 :
            fp.write('<TR><TD bgcolor=lightblue>')
        else :
            fp.write('<TR><TD>')
        fp.write('{:s}'.format(obj))
        fp.write('(<A HREF="http://simbad.cfa.harvard.edu/simbad/sim-basic?Ident={:12.5f}%09{:12.5f}++&submit=SIMBAD+search"> SIMBAD </A>)<BR>'.format
                   (star['RA'],star['DEC']))
        fp.write('H  = {:7.2f}<br>'.format(star['H']))
        fp.write('SNR  = {:7.2f}<br>'.format(star['SNR']))
        fp.write('RV_CHI2  = {:7.2f}<br>'.format(star['RV_CHI2']))
        fp.write('RV_CCFWHM  = {:7.2f}  RV_AUTOFWHM = {:7.2f} <br>'.format(star['RV_CCFWHM'],star['RV_AUTOFWHM']))
        fp.write('{:s}<br>'.format(star['TARGFLAGS']))
        fp.write('{:s}<br>'.format(star['STARFLAGS']))
        fp.write('{:s}<br>'.format(rvmask.getname(star['RV_FLAG'])))

        # average velocities
        fp.write('<TABLE BORDER=2>\n')
        fp.write('<TR><TD><TD>VHELIO_AVG<TD>VSCATTER<TD>TEFF<TD>LOGG<TD>[FE/H]\n')
        fp.write('<TR><TD>Doppler<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.0f}<TD>{:8.2f}<TD>{:8.2f}\n'.format(
                 star['VHELIO_AVG'],star['VSCATTER'],
                 star['RV_TEFF'],star['RV_LOGG'],star['RV_FEH']))
        fp.write('<TR><TD>Doppler Xcorr<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.0f}<TD>{:8.2f}<TD>{:8.2f}\n'.format(
                 np.median(apfv['XCORR_VHELIO'][j]),
                 apfv['XCORR_VHELIO'][j].std(ddof=1),
                 star['RV_TEFF'],star['RV_LOGG'],star['RV_FEH']))
        if k>=0 :
            gd = np.where(np.abs(apfieldvisits['VHELIO']) < 999)[0]
            fp.write('<TR><TD>DR16<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.0f}<TD>{:8.2f}<TD>{:8.2f}\n'.format(
                     apfield['VHELIO_AVG'][k],apfield['VSCATTER'][k],
                     apfield['RV_TEFF'][k],apfield['RV_LOGG'][k],apfield['RV_FEH'][k]))
        fp.write('</TABLE><br>')

        # flag bad RVs
        vhelio=apfv['VHELIO'][j]

        # individual visit velocities
        fp.write('<TABLE BORDER=2>')
        fp.write('<TR><TD>JD<TD>PLATE<TD>MJD<TD>FIBER<TD>S/N<TD>Doppler xcorr<TD> xcorr_err<TD>Doppler<TD>VERR<TD>DR16<TD>VERR<TD>Dop BC<TD>apS BC\n')
        for ind,i in enumerate(j) :
            try : 
                #ii = np.where(apfieldvisits['FILE'] == apfv['FILE'][i])[0][0]
                ii = np.where((apfieldvisits['APOGEE_ID'] == apfv['APOGEE_ID'][i]) &
                              (apfieldvisits['MJD'] == apfv['MJD'][i]) )[0][0]
                vhelio_idl =  apfieldvisits['VHELIO'][ii]
                vrelerr_idl =  apfieldvisits['VRELERR'][ii]
                bc_idl =  apfieldvisits['BC'][ii]
                vscatter_idl = apfield['VSCATTER'][k]
            except : 
                vhelio_idl,vrelerr_idl,bc_idl = -99999,-99999,-99999
                vscatter_idl = -99999
            if np.isfinite(apfv['VHELIO'][i]) == False :
                bgcolor='bgcolor=red'
            elif apfv['STARFLAG'][i] & starmask.getval('RV_REJECT') > 0 :
                bgcolor='bgcolor=lightpink'
            elif apfv['STARFLAG'][i] & starmask.getval('RV_SUSPECT') > 0 :
                bgcolor='bgcolor=#F4DEDE'
            else : bgcolor=''
            fp.write(('<TR {:s}> <TD> <A HREF={:s} TARGET="_obj"> {:12.3f}</A> <TD> {:s} <TD> {:5d} <TD> {:5d}'+
                     '<TD> {:8.1f} <TD> {:8.2f} <TD> {:8.2f} <TD> {:8.2f} <TD> {:8.2f} <TD> {:8.2f} <TD> {:8.2f}'+
                     '<TD> {:8.2f} <TD>{:8.2f}\n').format(
                      bgcolor,
                      apfv['FILE'][i].replace('.fits','_dopfit.png').replace('-r12-','-r13-'),
                      apfv['JD'][i],apfv['PLATE'][i],apfv['MJD'][i],apfv['FIBERID'][i],
                      apfv['SNR'][i],
                      apfv['XCORR_VHELIO'][i],apfv['XCORR_VRELERR'][i],
                      apfv['VHELIO'][i],apfv['VRELERR'][i],
                      vhelio_idl, vrelerr_idl,apfv['BC'][i],bc_idl))
        fp.write('</TABLE>\n')

        # vscatter difference with DR16
        fp.write('<TD> {:8.2f}\n'.format(star['H']))
        fp.write('<TD> {:8.2f}\n'.format(star['RV_TEFF']))
        fp.write('<TD> {:8.2f}\n'.format(star['RV_CHI2']))
        fp.write('<TD> {:d}\n'.format(star['N_COMPONENTS']))

        # plot visit RVs
        if doapfield : 
            vidl=apfieldvisits['VHELIO'][jj]
            gd = np.where(np.abs(vidl) < 999)[0]
            vmax=np.nanmax(np.append(vhelio,vidl[gd]))
            vmin=np.nanmin(np.append(vhelio,vidl[gd]))
        else :
            vmax=np.nanmax(vhelio)
            vmin=np.nanmin(vhelio)
        yr=[vmin-0.1*(vmax-vmin),vmax+0.1*(vmax-vmin)]
        try :
            fig,ax=plots.multi(1,1)
            gd_dop = np.where((apfv['STARFLAG'][j] & starmask.getval('RV_REJECT')) == 0)[0]
            if len(gd_dop) > 0 : 
                plots.plotp(ax,apfv['MJD'][j[gd_dop]],vhelio[gd_dop],size=15,color='g',yr=yr,label='Doppler')
            bd_dop = np.where((apfv['STARFLAG'][j] & starmask.getval('RV_REJECT')) > 0)[0]
            if len(bd_dop) > 0 : ax.scatter(apfv['MJD'][j[bd_dop]],vhelio[bd_dop],s=15,
                                            facecolors='none',edgecolors='g',label='rejected Doppler')
            ax.plot(ax.get_xlim(),[star['VHELIO_AVG'],star['VHELIO_AVG']],color='g')
            if doapfield : 
                plots.plotp(ax,apfieldvisits['MJD'][jj[gd]],vidl[gd],size=15,color='r',yr=yr,label='DR16')
                ax.plot(ax.get_xlim(),[apfield['VHELIO_AVG'][k],apfield['VHELIO_AVG'][k]],color='r')
            ax.legend()
            fig.savefig(outdir+'/plots/'+obj+'_rv.png')
            plt.close()
        except KeyboardInterrupt: raise
        except :
            print('Plotting error....')
            plt.close()
            pass

        # include plots
        fp.write('<TD><a HREF=plots/{:s}.png TARGET="_obj"> <IMG SRC=plots/{:s}.png WIDTH=600></A>\n'.format(obj,obj))
        fp.write('<TD><IMG SRC=plots/{:s}_rv.png TARGET="_obj">\n'.format(obj))
        fp.write('<TD><A HREF=plots/{:s}_ccf.png TARGET="_obj"> <IMG SRC=plots/{:s}_ccf.png></A>\n'.format(obj,obj))
        fp.write('<TD><A HREF=plots/{:s}_autoccf.png TARGET="_obj"> <IMG SRC=plots/{:s}_autoccf.png></A>\n'.format(obj,obj))
        fp.write('<TD><A HREF=plots/{:s}_spec2.png TARGET="_obj"> <IMG SRC=plots/{:s}_spec2.png></a>\n'.format(obj,obj))
        fp.write('<TD><A HREF=plots/{:s}_spec.png TARGET="_obj"> <IMG SRC=plots/{:s}_spec.png></a>\n'.format(obj,obj))
        fp.write('<TD><A HREF=plots/{:s}_cont.png TARGET="_obj"> <IMG SRC=plots/{:s}_cont.png></a>\n'.format(obj,obj))
    fp.close() 


def visitcomb(allvisit,load=None,nres=[5,4.25,3.5],bconly=False,
              plot=False,write=True,dorvfit=True,apstar_vers='stars') :
    """ Combine multiple visits with individual RVs to rest frame sum
    """

    if load is None : 
        print('need to supply load=')
        pdb.set_trace()

    cspeed = 2.99792458e5  # speed of light in km/s

    wnew=aspcap.apStarWave()  
    nwave=len(wnew)
    nvisit=len(allvisit)

    print('doing visitcomb for {:s}, nvisit: {:d}'.format(allvisit['APOGEE_ID'][0],nvisit))
    # initialize array for stack of interpolated spectra
    zeros = np.zeros([nvisit,nwave])
    izeros = np.zeros([nvisit,nwave],dtype=int)
    stack=apload.ApSpec(zeros,err=zeros.copy(),bitmask=izeros,cont=zeros.copy(),
                sky=zeros.copy(),skyerr=zeros.copy(),telluric=zeros.copy(),telerr=zeros.copy())
    stack_homogenized=apload.ApSpec(zeros.copy())
    nlsf=7
    stack_lsf = np.zeros([nvisit,nwave,2*nlsf+1])

    apogee_target1, apogee_target2, apogee_target3 = 0, 0, 0
    apogee2_target1, apogee2_target2, apogee2_target3, apogee2_target4 = 0, 0, 0, 0
    starflag,andflag = np.int64(0),np.int64(0)
    andflag = allvisit['STARFLAG'][0]
    starmask=bitmask.StarBitMask()

    # loop over each visit and interpolate to final wavelength grid
    if plot : fig,ax=plots.multi(1,2,hspace=0.001)
    pixlim=np.zeros([3,2],dtype=np.int32)
    pixlim_overlap=np.zeros([3,2],dtype=np.int32)
    for chip in range(3) :
        pixlim[chip,0] = nwave
        pixlim[chip,1] = 0
        pixlim_overlap[chip,0] = 0
        pixlim_overlap[chip,1] = nwave
    for i,visit in enumerate(allvisit) :

        if bconly : 
            #vrel = -visit['BC']
            if load.telescope == 'lco25m' :
                vrel=-bc.getbc(visit['RA'],visit['DEC'],visit['JD'],obs='LCO') / 1000.
            else :
                vrel=-bc.getbc(visit['RA'],visit['DEC'],visit['JD'],obs='APO') / 1000.
        else : vrel = visit['VREL']

        # skip if we don't have an RV
        if np.isfinite(vrel) is False : continue

        # load the visit
        if load.telescope == 'apo1m' :
            #apvisit=load.apVisit1m(visit['PLATE'],visit['MJD'],visit['APOGEE_ID'],load=True)
            visitfile=load.filename('Visit',plate=visit['PLATE'],mjd=visit['MJD'],reduction=visit['APOGEE_ID'])
            # use FILE tag in case we have MJDFRAC
            apvisit=load.apVisit1m(visit['PLATE'],visit['MJD'],visit['APOGEE_ID'],usefile=visit['FILE'],load=True)
        else :
            apvisit=load.apVisit(int(visit['PLATE']),visit['MJD'],visit['FIBERID'],load=True)
        pixelmask=bitmask.PixelBitMask()

        # rest-frame wavelengths transformed to this visit spectra
        w=aspcap.apStarWave()*(1.0+vrel/cspeed)

        # loop over the chips
        for chip in range(3) :

            # get the pixel values to interpolate to
            pix=wave.wave2pix(w,apvisit.wave[chip,:])
            gd=np.where(np.isfinite(pix))[0]

            pixlim[chip,0]=min([pixlim[chip,0],gd[0]])
            pixlim[chip,1]=max([pixlim[chip,1],gd[-1]])
            pixlim_overlap[chip,0]=max([pixlim_overlap[chip,0],gd[0]])
            pixlim_overlap[chip,1]=min([pixlim_overlap[chip,1],gd[-1]])

            # normalize if we have very large numbers (VESTA!)
            fluxnorm = np.median(apvisit.flux[chip,:])
            if fluxnorm > 1.e10 :
                apvisit.flux[chip,:] /= fluxnorm
                apvisit.err[chip,:] /= fluxnorm

            # get a smoothed, filtered spectrum to use as replacement for bad values when we do the sinc interpolation
            cont = gaussian_filter(median_filter(apvisit.flux[chip,:],[501],mode='reflect'),100)
            errcont = gaussian_filter(median_filter(apvisit.flux[chip,:],[501],mode='reflect'),100)
            bd = np.where(apvisit.bitmask[chip,:]&pixelmask.badval())[0]
            if len(bd) > 0 : 
                apvisit.flux[chip,bd] = cont[bd] 
                apvisit.err[chip,bd] = errcont[bd] 

            # load up quantity/error pairs for interpolation
            raw=[[apvisit.flux[chip,:],apvisit.err[chip,:]**2],
                 [apvisit.sky[chip,:],apvisit.skyerr[chip,:]**2],
                 [apvisit.telluric[chip,:],apvisit.telerr[chip,:]**2]]

            # load up individual mask bits
            for ibit,name in enumerate(pixelmask.name) :
                if name is not '' and len(np.where(apvisit.bitmask[chip,:]&2**ibit)[0]) > 0 :
                    raw.append([np.clip(apvisit.bitmask[chip,:]&2**ibit,None,1),None])

            # do the sinc interpolation
            out=sincint.sincint(pix[gd],nres[chip],raw)

            # from output flux, get continuum to remove, so that all spectra are
            #   on same scale. We'll later multiply in the median continuum
            flux = out[0][0]
            stack.cont[i,gd] = gaussian_filter(median_filter(flux,[501],mode='reflect'),100)

            # load interpolated spectra into output stack
            stack.flux[i,gd] = out[0][0] / stack.cont[i,gd]
            stack.err[i,gd] = out[0][1] / stack.cont[i,gd]
            stack.sky[i,gd] = out[1][0]
            stack.skyerr[i,gd] = out[1][1]
            stack.telluric[i,gd] = out[2][0]
            stack.telerr[i,gd] = out[2][1]
            # for mask, set bits where interpolated value is above some threshold
            #   defined for each mask bit
            iout=3
            for ibit,name in enumerate(pixelmask.name) :
                if name is not '' and len(np.where(apvisit.bitmask[chip,:]&2**ibit)[0]) > 0 :
                    j = np.where(np.abs(out[iout][0]) > pixelmask.maskcontrib[ibit])[0]
                    stack.bitmask[i,gd[j]] |= 2**ibit
                    iout+=1

        # ignore bad pixels
        bd = np.where((stack.bitmask[i,:]&pixelmask.badval()) > 0)[0]
        if len(bd) > 0 : stack.err[i,bd] = 1.e10

        # increase uncertainties for persistence pixels
        bd = np.where((stack.bitmask[i,:]&pixelmask.getval('PERSIST_HIGH')) > 0)[0]
        if len(bd) > 0 : stack.err[i,bd] *= np.sqrt(5)
        bd = np.where(((stack.bitmask[i,:]&pixelmask.getval('PERSIST_HIGH')) == 0) &
                      ((stack.bitmask[i,:]&pixelmask.getval('PERSIST_MED')) > 0) )[0]
        if len(bd) > 0 : stack.err[i,bd] *= np.sqrt(4)
        bd = np.where(((stack.bitmask[i,:]&pixelmask.getval('PERSIST_HIGH')) == 0) &
                      ((stack.bitmask[i,:]&pixelmask.getval('PERSIST_MED')) == 0) &
                      ((stack.bitmask[i,:]&pixelmask.getval('PERSIST_LOW')) > 0) )[0]
        if len(bd) > 0 : stack.err[i,bd] *= np.sqrt(3)
        bd = np.where((stack.bitmask[i,:]&pixelmask.getval('SIG_SKYLINE')) > 0)[0]
        if len(bd) > 0 : stack.err[i,bd] *= np.sqrt(100)

        # downweight spectrum if MTPFLUX_LT_50 bit set
        if visit['STARFLAG'] & starmask.getval('MTPFLUX_LT_50') :
            stack.err[i,:] *= 3.

        if plot :
            ax[0].plot(aspcap.apStarWave(),stack.flux[i,:])
            ax[1].plot(aspcap.apStarWave(),stack.flux[i,:]/stack.err[i,:])
            plt.draw()
            pdb.set_trace()

        # accumulate for header of combined frame. Turn off visit specific RV flags first
        visitflag = visit['STARFLAG'] & ~starmask.getval('RV_REJECT') & ~starmask.getval('RV_SUSPECT')
        starflag |= visitflag
        andflag &= visitflag
        if visit['SURVEY'] == 'apogee' :
            apogee_target1 |= visit['APOGEE_TARGET1'] 
            apogee_target2 |= visit['APOGEE_TARGET2'] 
        elif visit['SURVEY'].find('apogee2') >=0  :
            apogee2_target1 |= visit['APOGEE_TARGET1'] 
            apogee2_target2 |= visit['APOGEE_TARGET2'] 
            apogee2_target3 |= visit['APOGEE_TARGET3'] 
            try: apogee2_target4 |= visit['APOGEE_TARGET4'] 
            except: pass
        elif visit['SURVEY'] == 'apo1m' :
            apogee_target2 |= visit['APOGEE_TARGET2'] 
            apogee2_target2 |= visit['APOGEE_TARGET2'] 

        # LSF for this image in resampled apStar frame
        aplsf=True
        if aplsf :
            nlsf_homogenized=35
            x,ls=lsf.getlsf(apvisit.header['LSFID'],apvisit.header['WAVEID'],apred=load.apred,telescope=load.telescope,
                            nlsf2=nlsf_homogenized,fiber='all',highres=1,fill=True)

            # for apStarLSF, use trimmed LSF
            stack_lsf[i]=ls[visit['FIBERID']-1][:8575,-nlsf+nlsf_homogenized:nlsf+nlsf_homogenized+1]
            sum=np.sum(stack_lsf[i],axis=1)
            for ipix in range(8575) : stack_lsf[i,ipix,:]/=sum[i]

        # LSF homogenized spectrum: still in development!
        homogenize = False
        if homogenize :
            spec=copy.deepcopy(stack.flux[i,:])
            err=copy.deepcopy(stack.err[i,:])
            # set bad pixels to 1
            bd=np.where(np.isnan(spec) | np.isclose(spec,0.) | 
                        (stack.bitmask[i,:]&pixelmask.badval()) |
                        (stack.bitmask[i,:]&pixelmask.getval('SIG_SKYLINE')) )[0]
            spec[bd] = 1.
            err[bd] = 1.e5

            # taper the window
            lsf_homogenize = ls[visit['FIBERID']-1][:8575,:]*tukey(2*nlsf_homogenized+1, alpha=0.2, sym=True)
            # normalize LSF
            for ipix in range(nlsf_homogenized,8575-nlsf_homogenized-1) : 
                lsf_homogenize[ipix,:]/=lsf_homogenize[ipix,:].sum()

            # test delta functions
            test = False
            if test :
                spec*=0
                spec[::100] = 1.
                new=copy.copy(spec)
                for ipix in range(nlsf_homogenized,len(spec)-nlsf_homogenized) :
                    kernel=lsf_homogenize[ipix]
                    new[ipix] = np.sum(spec[ipix-nlsf_homogenized:ipix+nlsf_homogenized+1]*kernel[::-1])
                spec=copy.copy(new)

            lsq = False
            if lsq :
                # least-squares fit for underlying spectrum, still in development, not working! and slow
                rhs=np.zeros(8575)
                curv=np.zeros([8575,4*nlsf_homogenized+1])
                for ipix in range(8575) :
                    print(ipix)
                    i1 = np.max([0,ipix-nlsf_homogenized])
                    i2 = np.min([8574,ipix+nlsf_homogenized+1])
                    for j,jpix in enumerate(range(np.max([0,ipix-2*nlsf_homogenized]),np.min([8574,ipix+2*nlsf_homogenized+1]))) :
                        j1 = np.max([0,jpix-nlsf_homogenized])
                        j2 = np.min([8574,jpix+nlsf_homogenized+1])
                        istart = np.max([i1,j1])
                        iend = np.min([i2,j2])
                        if iend>istart :
                            curv[ipix,j] = (lsf_homogenize[ipix,istart-ipix+nlsf_homogenized:iend-ipix+nlsf_homogenized] *
                                            lsf_homogenize[jpix,istart-jpix+nlsf_homogenized:iend-jpix+nlsf_homogenized]).sum()
               
                    #print(ipix,jpix,curv[ipix,jpix])
                    rhs[ipix] = (lsf_homogenize[ipix,i1-ipix+nlsf_homogenized:i2-ipix+nlsf_homogenized] *
                                 spec[i1:i2]).sum()
                pdb.set_trace()
                # lsq
                curv[:,2*nlsf_homogenized]+=0.001
                out=scipy.linalg.solve_banded((2*nlsf_homogenized,2*nlsf_homogenized), curv.T, rhs) #/err**2)

            else :
                # deconvolution:
                out=scipy.linalg.solve_banded((nlsf_homogenized,nlsf_homogenized), lsf_homogenize.T, spec) 

                R_homogenized=15000.
                outsig=1/(aspcap.dlogw*np.log(10))/R_homogenized/2.354
                kernel=lsf.gauss(np.arange(-nlsf_homogenized,nlsf_homogenized+1),0,outsig)
                stack_homogenized.flux[i]=np.convolve(out,kernel,mode='same')
                stack_homogenized.flux[i,bd] = np.nan

            plot=True
            if plot :
                plt.clf()
                plt.plot(out,color='y')
                plt.plot(spec)
                plt.plot(stack_homogenized.flux[i])
                plt.ylim(-0.5,2)
                plt.draw()
                pdb.set_trace()

    # create final spectrum
    if nvisit > 1 :
        zeros = np.zeros([nvisit+2,nwave])
        izeros = np.zeros([nvisit+2,nwave],dtype=int)
    else :
        zeros = np.zeros([2,nwave])
        izeros = np.zeros([2,nwave],dtype=int)

    if not bconly :
        if len(allvisit) == 1 :
            rvtab=Table(np.squeeze(np.vstack(allvisit['RVTAB']),axis=0))
        else :
            rvtab=Table(np.squeeze(np.vstack(allvisit['RVTAB'])))
        gdccf=np.where(np.isfinite(rvtab['x_ccf'][0,:]))[0]
        rvtab['x_ccf']=rvtab['x_ccf'][:,gdccf]
        rvtab['ccf']=rvtab['ccf'][:,gdccf]
        rvtab['ccferr']=rvtab['ccferr'][:,gdccf]
        rvtab['autoccf']=rvtab['autoccf'][:,gdccf]
    else : rvtab = None

    apstar=apload.ApSpec(zeros,err=zeros.copy(),bitmask=izeros,wave=aspcap.apStarWave(),
                sky=zeros.copy(),skyerr=zeros.copy(),telluric=zeros.copy(),telerr=zeros.copy(),
                cont=zeros.copy(),template=zeros.copy(),rvtab=rvtab)
    apstar.header['CRVAL1'] = (aspcap.logw0,'Start log10(wavelength) in subsequent HDUs')
    apstar.header['CDELT1'] = (aspcap.dlogw,'Dispersion in log10(wave) in subsequent HDUs')
    apstar.header['CRPIX1'] = (1,'Pixel of starting wavelength in subsequent HDUs')
    apstar.header['CTYPE1'] = ('LOG-LINEAR','Logarithmic wavelength scale in subsequent HDUs')
    apstar.header['DC-FLAG'] = (1,'Logarithmic wavelength scale in subsequent HDUs')
    apstar.header['NWAVE'] = (nwave,'Number of wavelengths in subsequent HDUs')

    # pixel-by-pixel weighted average
    cont = np.sum(stack.cont/stack.err**2,axis=0)/np.sum(1./stack.err**2,axis=0)
    apstar.flux[0,:] = np.sum(stack.flux/stack.err**2,axis=0)/np.sum(1./stack.err**2,axis=0) * cont
    apstar.err[0,:] =  np.sqrt(1./np.sum(1./stack.err**2,axis=0)) * cont
    apstar.bitmask[0,:] = np.bitwise_and.reduce(stack.bitmask,0)
    apstar.cont[0,:] = cont

    # LSF homogenized combination
    if homogenize :
        apstar.flux[1,:] = np.sum(stack_homogenized.flux/stack.err**2,axis=0)/np.sum(1./stack.err**2,axis=0) * cont
        apstar.err[1,:] =  np.sqrt(1./np.sum(1./stack.err**2,axis=0)) * cont

    # apStarLSF
    apstarlsf = ( np.sum(stack_lsf/np.repeat(stack.err**2,2*nlsf+1).reshape(nvisit,8575,2*nlsf+1),axis=0)/
                  np.sum(1./np.repeat(stack.err**2,2*nlsf+1).reshape(nvisit,8575,2*nlsf+1),axis=0) )

    # global weighting and individual visits
    if nvisit > 1 :
        # "global" weighted average
        if not homogenize :
            newerr = median_filter(stack.err,[1,100],mode='reflect')
            bd = np.where((stack.bitmask&pixelmask.getval('SIG_SKYLINE')) > 0)[0]
            if len(bd) > 0 : newerr[bd[0],bd[1]] *= np.sqrt(100)
            apstar.flux[1,:] = np.sum(stack.flux/newerr**2,axis=0)/np.sum(1./newerr**2,axis=0) * cont
            apstar.err[1,:] =  np.sqrt(1./np.sum(1./newerr**2,axis=0)) * cont

        apstar.flux[2:,:] = stack.flux * stack.cont
        apstar.err[2:,:] = stack.err * stack.cont
        apstar.bitmask[2:,:] = stack.bitmask
        apstar.sky[2:,:] = stack.sky
        apstar.skyerr[2:,:] = stack.skyerr
        apstar.telluric[2:,:] = stack.telluric
        apstar.telerr[2:,:] = stack.telerr

    if bconly : return apstar

    # populate header
    apstar.header['FIELD'] = (allvisit['FIELD'][0], 'APOGEE field name')
    apstar.header['OBJID'] = (allvisit['APOGEE_ID'][0], 'APOGEE object name')
    try :apstar.header['SNR'] = (np.nanmedian(apstar.flux[0,:]/apstar.err[0,:]), 'Median S/N per apStar pixel')
    except :apstar.header['SNR'] = (0., 'Median S/N per apStar pixel')
    gdwave=np.where((wnew>15850) & (wnew<16440))[0]
    try :apstar.header['SNREV'] = (np.nanmedian(apstar.flux[0,gdwave]/apstar.err[0,gdwave]), 'Median S/N per apStar pixel in middle chip')
    except :apstar.header['SNREV'] = (0., 'Median S/N per apStar pixel in middle chip')
    apstar.header['RA'] = (allvisit['RA'].max(), 'right ascension, deg, J2000')
    apstar.header['DEC'] = (allvisit['DEC'].max(), 'declination, deg, J2000')
    apstar.header['GLON'] = (allvisit['GLON'].max(), 'Galactic longitude')
    apstar.header['GLAT'] = (allvisit['GLAT'].max(), 'Galactic latitude')
    apstar.header['J'] = (allvisit['J'].max(), '2MASS J magnitude')
    apstar.header['J_ERR'] = (allvisit['J_ERR'].max(), '2MASS J magnitude uncertainty')
    apstar.header['H'] = (allvisit['H'].max(), '2MASS H magnitude')
    apstar.header['H_ERR'] = (allvisit['H_ERR'].max(), '2MASS H magnitude uncertainty')
    apstar.header['K'] = (allvisit['K'].max(), '2MASS K magnitude')
    apstar.header['K_ERR'] = (allvisit['K_ERR'].max(), '2MASS K magnitude uncertainty')
    try: apstar.header['SRC_H'] = (allvisit[0]['SRC_H'], 'source of H magnitude')
    except KeyError: pass
    keys=[ 'WASH_M','WASH_T2', 'DDO51','IRAC_3_6',
           'IRAC_4_5','IRAC_5_8', 'WISE_4_5','TARG_4_5']
    for key in keys :
        try: apstar.header[key] = allvisit[key].max()
        except KeyError: pass

    apstar.header['AKTARG'] = (allvisit['AK_TARG'].max(), 'Extinction used for targeting')
    apstar.header['AKMETHOD'] = (allvisit[0]['AK_TARG_METHOD'],'Extinction method using for targeting')
    apstar.header['AKWISE'] = (allvisit['AK_WISE'].max(),'WISE all-sky extinction')
    apstar.header['SFD_EBV'] = (allvisit['SFD_EBV'].max(),'SFD E(B-V)')
    apstar.header['APTARG1'] = (apogee_target1, 'APOGEE_TARGET1 targeting flag')
    apstar.header['APTARG2'] = (apogee_target2, 'APOGEE_TARGET2 targeting flag')
    apstar.header['AP2TARG1'] = (apogee2_target1, 'APOGEE2_TARGET1 targeting flag')
    apstar.header['AP2TARG2'] = (apogee2_target2, 'APOGEE2_TARGET2 targeting flag')
    apstar.header['AP2TARG3'] = (apogee2_target3, 'APOGEE2_TARGET3 targeting flag')
    apstar.header['AP2TARG4'] = (apogee2_target4, 'APOGEE2_TARGET4 targeting flag')
    apstar.header['NVISITS'] = (len(allvisit), 'Number of visit spectra combined flag')
    apstar.header['STARFLAG'] = (starflag,'bitwise OR of individual visit starflags')
    apstar.header['ANDFLAG'] = (andflag,'bitwise AND of individual visit starflags')

    try: apstar.header['N_COMP'] = (allvisit['N_COMPONENTS'].max(),'Maximum number of components in RV CCFs')
    except: pass
    apstar.header['VHELIO'] = ((allvisit['VHELIO']*allvisit['SNR']).sum() / allvisit['SNR'].sum(),'S/N weighted mean barycentric RV')
    if len(allvisit) > 1 : apstar.header['VSCATTER'] = (allvisit['VHELIO'].std(ddof=1), 'standard deviation of visit RVs')
    else : apstar.header['VSCATTER'] = (0., 'standard deviation of visit RVs')
    verr = np.sqrt((allvisit['VRELERR']**2*allvisit['SNR']**2).sum() / (allvisit['SNR']**2).sum())
    apstar.header['VERR'] = (verr,'weighted error in VHELIO')
    apstar.header['VERR_MED'] = (np.median(allvisit['VRELERR']),'median error in VHELIO')

    apstar.header['RV_TEFF'] = (allvisit['RV_TEFF'].max(),'Effective temperature from RV fit')
    apstar.header['RV_LOGG'] = (allvisit['RV_LOGG'].max(),'Surface gravity from RV fit')
    apstar.header['RV_FEH'] = (allvisit['RV_FEH'].max(),'Metallicity from RV fit')
    # these are filled below, but set here to make sure cards exist even in case of failure
    apstar.header['RV_CHI2'] = (0.,' chisq from Doppler RV fit')
    apstar.header['CCFWHM'] = (0.,' FWHM of RV CCF of star with template (km/s)')
    apstar.header['AUTOFWHM'] = (0.,' FWHM of RV CCF of template with template (km/s)')

    if len(allvisit) > 0 : meanfib=(allvisit['FIBERID']*allvisit['SNR']).sum()/allvisit['SNR'].sum()
    else : meanfib = 999999.
    if len(allvisit) > 1 : sigfib=allvisit['FIBERID'].std(ddof=1)
    else : sigfib = 0.
    apstar.header['MEANFIB'] = (meanfib,'S/N weighted mean fiber number')
    apstar.header['SIGFIB'] = (sigfib,'standard deviation (unweighted) of fiber number')
    apstar.header['NRES'] = ('{:5.2f}{:5.2f}{:5.2f}'.format(*nres),'number of pixels/resolution used for sinc')
    apstar.header['MIN_H'] = (allvisit['MIN_H'].min(),'minimum H mag of cohort')
    apstar.header['MAX_H'] = (allvisit['MAX_H'].max(),'maximum H mag of cohort')
    apstar.header['MIN_JK'] = (allvisit['MIN_JK'].min(),'minimum J-K of cohort')
    apstar.header['MAX_JK'] = (allvisit['MAX_JK'].max(),'maximum J-K of cohort')
    apstar.header['RMIN'] = (pixlim[0,0],'Min pixel of red chip contrib, any frame')
    apstar.header['RMAX'] = (pixlim[0,1],'Maxmm pixel of red chip contrib, any frame')
    apstar.header['GMIN'] = (pixlim[1,0],'Min pixel of green chip contrib, any frame')
    apstar.header['GMAX'] = (pixlim[1,1],'Maxm pixel of green chip contrib, any frame')
    apstar.header['BMIN'] = (pixlim[2,0],'Min pixel of blue chip contrib, any frame')
    apstar.header['BMAX'] = (pixlim[2,1],'Maxm pixel of blue chip contrib, any frame')
    apstar.header['ROVERMIN'] = (pixlim_overlap[0,0],'Min pixel of red chip contrib, all frames')
    apstar.header['ROVERMAX'] = (pixlim_overlap[0,1],'Max pixel of red chip contrib, all frames')
    apstar.header['GOVERMIN'] = (pixlim_overlap[1,0],'Min pixel of green chip contrib, all frames')
    apstar.header['GOVERMAX'] = (pixlim_overlap[1,1],'Max pixel of green chip contrib, all frames')
    apstar.header['BOVERMIN'] = (pixlim_overlap[2,0],'Minimum pixel of blue chip contrib, all frames')
    apstar.header['BOVERMAX'] = (pixlim_overlap[2,1],'Maximum pixel of blue chip contrib, all frames')

    # individual visit information in header
    for i0,visit in enumerate(allvisit) :
        i=i0+1
        apstar.header['SFILE{:d}'.format(i)] = (visit['FILE'],' Visit #{:d} spectrum file'.format(i))
        apstar.header['DATE{:d}'.format(i)] = (visit['DATEOBS'], 'DATE-OBS of visit {:d}'.format(i))
        apstar.header['JD{:d}'.format(i)] = (visit['JD'], 'Julian date of visit {:d}'.format(i))
        # hjd = helio_jd(visitstr[i].jd-2400000.0,visitstr[i].ra,visitstr[i].dec)
        #apstar.header['HJD{:d}'.format(i)] = 
        apstar.header['FIBER{:d}'.format(i)] = (visit['FIBERID'],'Fiber, visit {:d}'.format(i))
        apstar.header['BC{:d}'.format(i)] = (visit['BC'],'Barycentric correction (km/s), visit {:d}'.format(i))
        apstar.header['VRAD{:d}'.format(i)] = (visit['VREL'],'Radial velocity (km/s) of visit {:d}'.format(i))
        apstar.header['VERR{:d}'.format(i)] =  (visit['VRELERR'],'Uncertainty in radial velocity (km/s)')
        apstar.header['VHELIO{:d}'.format(i)] = (visit['VHELIO'],'Barycentric velocity (km/s), visit {:d}'.format(i))
        #apstar.header['SNRVIS{:d}'.format(i)] = (visit['SNR'],'Signal/Noise ratio, visit {:d}'.format(i))
        try :apstar.header['SNRVIS{:d}'.format(i)] = (np.nanmedian(apstar.flux[i0+2,:]/apstar.err[i0+2,:]), 
                                                      'Median S/N per apStar pixel,, visit {:d}'.format(i))
        except :apstar.header['SNRVIS{:d}'.format(i)] = (0., 'Median S/N per apStar pixel, visit {:d}'.format(i))
        apstar.header['FLAG{:d}'.format(i)] = (visit['STARFLAG'],'STARFLAG for visit {:d}'.format(i))
        apstar.header.insert('SFILE{:d}'.format(i),('COMMENT','VISIT {:d} INFORMATION'.format(i)))

    # Do a RV fit just to get a template and normalized spectrum, for plotting, and FWHM and AUTOFWHM
    if dorvfit :
        try :
            apstar.setmask(pixelmask.badval())
            spec=doppler.Spec1D(apstar.flux[0,:],err=apstar.err[0,:],bitmask=apstar.bitmask[0,:],
                 mask=apstar.mask[0,:],wave=apstar.wave,lsfpars=np.array([0]),
                 lsfsigma=apstar.wave/22500/2.354,instrument='APOGEE',
                 filename=apstar.filename)
            out= doppler.rv.jointfit([spec],verbose=False,plot=False,tweak=False,maxvel=[-500,500])
            apstar.cont=out[3][0].flux
            apstar.template=out[2][0].flux
            apstar.header['RV_CHI2'] = (out[1]['chisq'][0],'chisq from Doppler RV fit')
            apstar.header['CCFWHM'] = (out[1]['ccpfwhm'][0],'FWHM of RV CCF of star with template (km/s)')
            apstar.header['AUTOFWHM'] = (out[1]['autofwhm'][0],'FWHM of RV CCF of template with template (km/s)')
            rvtab=Table(out[1])
            gdccf=np.where(np.isfinite(rvtab['x_ccf'][0,:]))[0]
            rvtab['x_ccf']=rvtab['x_ccf'][:,gdccf]
            rvtab['ccf']=rvtab['ccf'][:,gdccf]
            rvtab['ccferr']=rvtab['ccferr'][:,gdccf]
            rvtab['autoccf']=rvtab['autoccf'][:,gdccf]
            apstar.comb_rvtab=rvtab
        except ValueError as err:
            print('Exception raised in visitcomb RV for: ', apstar.header['FIELD'],apstar.header['OBJID'])
            print("ValueError: {0}".format(err))
        except RuntimeError as err:
            print('Exception raised in visitcomb RV for: ', apstar.header['FIELD'],apstar.header['OBJID'])
            print("Runtime error: {0}".format(err))
        except : 
            print('Exception raised in visitcomb RV fit for: ',apstar.header['FIELD'],apstar.header['OBJID'])

    if write :
        outfile=load.filename('Star',field=apstar.header['FIELD'],obj=apstar.header['OBJID'])
        outlsffile=load.filename('StarLSF',field=apstar.header['FIELD'],obj=apstar.header['OBJID'])
        if apstar_vers != 'stars' :
            outfile=outfile.replace('/stars/','/'+apstar_vers+'/')
            outlsffile=outlsffile.replace('/stars/','/'+apstar_vers+'/')
        outdir = os.path.dirname(outfile)
        try: os.makedirs(os.path.dirname(outfile))
        except : pass
        apstar.write(outfile)
        hdulist=fits.HDUList()
        hdulist.append(fits.ImageHDU(apstarlsf))
        hdulist.writeto(outlsffile,overwrite=True)
        
        # plot
        gd=np.where((apstar.bitmask[0,:] & (pixelmask.badval()|pixelmask.getval('SIG_SKYLINE'))) == 0) [0]
        fig,ax=plots.multi(1,3,hspace=0.001,figsize=(48,6))
        med=np.nanmedian(apstar.flux[0,:])
        plots.plotl(ax[0],aspcap.apStarWave(),apstar.flux[0,:],color='k',yr=[0,2*med])
        ax[0].plot(aspcap.apStarWave()[gd],apstar.flux[0,gd],color='g')
        plots.plotl(ax[0],aspcap.apStarWave(),apstar.flux[1,:],color='y',linewidth=1)
        ax[0].set_ylabel('Flux')
        try :
            ax[1].plot(aspcap.apStarWave()[gd],apstar.cont[gd],color='g')
            ax[1].set_ylabel('Normalized')
            ax[1].plot(aspcap.apStarWave(),apstar.template,color='r')
        except : pass
        plots.plotl(ax[2],aspcap.apStarWave(),apstar.flux[0,:]/apstar.err[0,:],yt='S/N')
        for i in range(3) : ax[i].set_xlim(15100,17000)
        ax[0].set_xlabel('Wavelength')
        fig.savefig(outdir+'/plots/'+apstar.header['OBJID']+'.png')
        plt.close()
       
        fig,ax=plots.multi(1,1) 
        try:
            ax.plot(rvtab['x_ccf'][0],rvtab['ccf'][0])
            ax.plot(rvtab['x_ccf'][0],rvtab['ccf'][0]/rvtab['ccf'][0].max())
            n=len(rvtab['x_ccf'][0])
            ax.plot(rvtab['x_ccf'][0],rvtab['autoccf'][0])
            ax.text(0.1,0.9,'ccpfwhm: {:8.2f}   autofwhm: {:8.2f}'.format(rvtab['ccpfwhm'][0],rvtab['autofwhm'][0]),transform=ax.transAxes)
        except: pass
        fig.savefig(outdir+'/plots/'+apstar.header['OBJID']+'_autoccf.png')
        plt.close()

    # plot
    if plot : 
        ax[0].plot(aspcap.apStarWave(),apstar.flux,color='k')
        ax[1].plot(aspcap.apStarWave(),apstar.flux/apstar.err,color='k')
        plt.draw()
        pdb.set_trace()

    return apstar
