import numpy as np
from esutil import htm
from astropy.table import Table
from astropy.io import fits
from apogee.utils import bitmask
from tools import match
import os
import pdb
 
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
              'phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag',
              'radial_velocity','radial_velocity_error','r_est','r_lo','r_hi')
    dtypes=('i8','f8','f8','f8','f8','f8','f8','f4','f4','f4','f8','f8','f8','f8','f8')
    out_names=[]
    for name in in_names: out_names.append(('gaia_'+name).upper())
    newcols=Table(np.zeros([len(tab),len(out_names)])-9999.,names=out_names,dtype=dtypes)
    # for source_id, default to 0, not -9999.
    newcols['GAIA_SOURCE_ID'] = 0
    # get rid of targetting proper motions to avoid confusion!
    tab.remove_columns(['PMRA','PMDEC','PM_SRC'])
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
        bd = np.where(np.isnan(tab[outname]))[0]
        tab[outname][bd] = -9999.

    return tab

def add_spec(data) :
    tab=Table(data)
    names=['TEFF_SPEC','LOGG_SPEC']
    newcols=Table(np.zeros([len(tab),len(names)])-9999.,names=names)
    tab.add_columns(newcols.columns.values())
    tab['TEFF_SPEC'] = tab['FPARAM'][:,0]
    tab['LOGG_SPEC'] = tab['FPARAM'][:,1]

    return(tab)

        
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
    for i in range(len(old)):
        if i%1000 == 0 : print(i)
        if ((new['FIELD'][i] != old['FIELD'][i]) or
            (new['TEFF'][i] != old['TEFF'][i]) or
            (new['FPARAM'][i,0] != old['FPARAM'][i,0]) or
            (new['FPARAM'][i,1] != old['FPARAM'][i,1]) or
            (new['STARFLAG'][i] != old['STARFLAG'][i]) or
            (new['TARGFLAGS'][i] != old['TARGFLAGS'][i]) or
            (new['FE_H'][i] != old['FE_H'][i]) or
            (new['MG_FE'][i] != old['MG_FE'][i]) or
            (new['SNR'][i] != old['SNR'][i]) or
            (new['TELESCOPE'][i] != old['TELESCOPE'][i]) or
            (new['H'][i] != old['H'][i]) or
            (new['AK_TARG'][i] != old['AK_TARG'][i]) or
            (new['AK_WISE'][i] != old['AK_WISE'][i]) or
            (new['APOGEE_ID'][i] != old['APOGEE_ID'][i]) or
            (new['MEANFIB'][i] != old['MEANFIB'][i]) 
           ) :
             print(i,new['APOGEE_ID'][i],old['APOGEE_ID'][i],new['H'][i],old['H'][i])
             nbad+=1

    print('nbad:',nbad)
 
