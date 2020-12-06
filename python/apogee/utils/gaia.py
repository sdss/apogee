import numpy as np
import tempfile
try: from esutil import htm
except: print('esutil not available!')
from astropy.table import Table, Column
from astropy.io import fits
from apogee.utils import bitmask
from tools import match
import os
import pdb
from astroquery.gaia import Gaia

def getdata(data) :
    """ Given input structure, get GAIA information from 2MASS matches
        and positional match
        Returns two tables
    """

    tab=Table()
    tab.add_column(Column(data['APOGEE_ID'],name='twomass'))
    tab.add_column(Column(data['RA'],name='apogee_ra'))
    tab.add_column(Column(data['DEC'],name='apogee_dec'))
    #if type(data['APOGEE_ID'][0]) is str or type(data['APOGEE_ID'][0]) is np.str_ : 
    try:
        j=np.where(np.core.defchararray.find(data['APOGEE_ID'],'2M') == 0)[0]
        out,ind=np.unique(np.core.defchararray.replace(data['APOGEE_ID'][j],'2M',''),return_index=True)
    except :
        j=np.where(np.core.defchararray.find(data['APOGEE_ID'],b'2M') == 0)[0]
        out,ind=np.unique(np.core.defchararray.replace(data['APOGEE_ID'][j],b'2M',b''),return_index=True)
    tab['twomass'][ind] = out
    #tab.add_column(Column(out,name='twomass'))
    #tab.add_column(Column(data['RA'][ind],name='apogee_ra'))
    #tab.add_column(Column(data['DEC'][ind],name='apogee_dec'))
    xmlfilename= tempfile.mktemp('.xml',dir=os.getcwd())
    tab.write(xmlfilename,format='votable',overwrite=True)
    try :
        job= Gaia.launch_job_async(
                """SELECT tmass_match.original_ext_source_id, g.source_id, g.ra, g.dec, g.parallax, g.parallax_error, 
                           g.pmra, g.pmra_error, g.pmdec, g.pmdec_error, 
                           g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag, 
                           g.radial_velocity, g.radial_velocity_error, g.a_g_val, g.e_bp_min_rp_val, 
                           dist.r_est, dist.r_lo, dist.r_hi
                   FROM gaiadr2.gaia_source AS g
                   INNER JOIN gaiadr2.tmass_best_neighbour AS tmass_match ON tmass_match.source_id = g.source_id
                   INNER JOIN tap_upload.my_table as ids on ids.twomass = tmass_match.original_ext_source_id
                   LEFT OUTER JOIN external.gaiadr2_geometric_distance as dist ON  g.source_id = dist.source_id""",
                   upload_resource=xmlfilename,upload_table_name='my_table')
        twomass_gaia = job.get_results()
    except:
        print("error with gaia 2mass search")
        twomass_gaia = None

    try: 
        job= Gaia.launch_job_async(
                """SELECT  g.source_id, g.ra, g.dec, g.parallax, g.parallax_error, 
                           g.pmra, g.pmra_error, g.pmdec, g.pmdec_error, 
                           g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag, 
                           g.radial_velocity, g.radial_velocity_error, g.a_g_val, g.e_bp_min_rp_val, 
                           dist.r_est, dist.r_lo, dist.r_hi,
                           distance(
                             point('', ids.apogee_ra, ids.apogee_dec),
                             point('', g.ra, g.dec)
                           ) * 3600 as dist_arcsec
                   FROM gaiadr2.gaia_source as g
                   JOIN tap_upload.my_table as ids on 1 = contains(
                     point('', ids.apogee_ra, ids.apogee_dec),
                     circle('', g.ra, g.dec, 3. / 3600)
                   )
                   LEFT OUTER JOIN external.gaiadr2_geometric_distance as dist ON  g.source_id = dist.source_id""",
                   upload_resource=xmlfilename,upload_table_name='my_table')
        posn_gaia = job.get_results()
    except: 
        print("error with gaia position search")
        posn_gaia = None
    finally: os.remove(xmlfilename)

    return twomass_gaia, posn_gaia


def add_gaia(data) :
    """ Add GAIA data to input structure, with 2MASS match and coordinate match to (cross-matched) GAIA reference file
    """

    # get the GAIA data from both matches
    gaia_twomass, gaia_posn = getdata(data)

    # add new columns
    tab=Table(data)
    in_names=('source_id','parallax','parallax_error','pmra','pmra_error','pmdec','pmdec_error',
              'phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag','a_g_val', 'e_bp_min_rp_val',
              'radial_velocity','radial_velocity_error', 'r_est','r_lo','r_hi')
    dtypes=('i8','f4','f4','f4','f4','f4','f4','f4','f4','f4','f4','f4','f4','f4','f4','f4','f4')
    out_names=[]
    for name in in_names: out_names.append(('gaia_'+name).upper())
    # initialize
    newcols=Table(np.zeros([len(tab),len(out_names)])-9999.,names=out_names,dtype=dtypes)
    # for source_id, default to 0, not -9999.
    newcols['GAIA_SOURCE_ID'] = 0

    # rename targetting proper motions to avoid confusion!
    try: tab.rename_column('PMRA','TARG_PMRA')
    except: pass
    try: tab.rename_column('PMDEC','TARG_PMDEC')
    except: pass
    try: tab.rename_column('PM_SRC','TARG_PM_SRC')
    except: pass
    # add unpopulated columns
    for col in newcols.columns.values() :
        try: tab.add_column(col)
        except ValueError: pass
    #tab.add_columns(newcols.columns.values())

    #if gaia_twomass is None or gaia_posn is None : return tab

    # remove dups in GAIA twomass in favor of brightest
    print('number in GAIA-2MASS xmatch catalog: ',len(gaia_twomass),len(set(gaia_twomass['original_ext_source_id'])))
    ind=[]
    for tid in set(gaia_twomass['original_ext_source_id']) :
        j=np.where(gaia_twomass['original_ext_source_id'] == tid)[0]
        if len(j)> 1:
            ii=np.argsort(gaia_twomass['phot_rp_mean_mag'][j])
            ind.append(j[ii[0]])
            print('duplicate 2MASS: ',gaia_twomass['phot_rp_mean_mag'][j[ii]])
        else : ind.append(j)

    # read gaia 2MASS matched file, match by 2MASS ID, and populate
    while len(gaia_twomass)>0 :
        # loop for matches since we may have repeats and want them all matched
        j=np.where(tab['GAIA_SOURCE_ID'] == 0)[0]
        print('Number missing gaia_source_id: ', len(j))
        if len(j) == 0 : break
        if type(tab['APOGEE_ID'][0]) is np.str_ : 
            m1,m2=match.match(np.core.defchararray.replace(tab['APOGEE_ID'][j],'2M',''),gaia_twomass['original_ext_source_id'])
        else :
            m1,m2=match.match(np.core.defchararray.replace(tab['APOGEE_ID'][j],b'2M',b''),gaia_twomass['original_ext_source_id'])
        print('Number matched by 2MASS: ', len(m1))
        if len(m1) == 0 : break
        for inname,outname in zip(in_names,out_names) :
            tab[outname][j[m1]] = gaia_twomass[inname][m2]

    j=np.where(tab['GAIA_SOURCE_ID'] > 0)[0]
    print('number of unique APOGEE_ID matches: ',len(set(tab['APOGEE_ID'][j])))

    j=np.where(tab['GAIA_SOURCE_ID'] == 0)[0]
    print('missing sources after 2MASS matches: ',len(j))
    h=htm.HTM()
    # now do a positional match, take the brightest object within 3 arcsec (which is the max from the GAIA crossmatch)
    maxrad=3./3600.
    m1,m2,rad=h.match(tab['RA'][j],tab['DEC'][j],gaia_posn['ra'],gaia_posn['dec'],maxrad,maxmatch=10)
    for m in set(m1) :
        jj=np.where(m1 == m)[0]
        ii=np.argsort(gaia_posn['phot_rp_mean_mag'][m2[jj]])
        for inname,outname in zip(in_names,out_names) :
            tab[outname][j[m]] = gaia_posn[inname][m2[jj[ii[0]]]]

    j=np.where(tab['GAIA_SOURCE_ID'] == 0)[0]
    print('missing sources after second match: ',len(j))

    # replace NaNs
    for name in out_names :
        bd = np.where(np.isnan(tab[name]))[0]
        tab[name][bd] = -9999.

    return tab

