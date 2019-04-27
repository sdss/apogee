import numpy as np
from esutil import htm
from astropy.table import Table
from astropy.io import fits
import pdb
 
def mkcoord(file='allStar-r12-l33-58358.fits') :
    """ Create coordinate CSV from allStar file, to use for GAIA cross-match
    """
    hdulist=fits.open(file)
    data=hdulist[1].data
    coords=np.vstack((data['RA'],data['DEC'])).T
    np.savetxt('coords.csv',coords,delimiter=',',header='ra,dec')

def add_gaia(data,gaiafile='gaia.fits.gz') :
    """ Add GAIA data to allStar file, with coordinate match to (cross-matched) GAIA reference file
    """
    tab=Table(data)
    in_names=('source_id','parallax','parallax_error','pmra','pmra_error','pmdec','pmdec_error')
    out_names=[]
    for name in in_names: out_names.append(('gaia_'+name).upper())
    newcols=Table(np.zeros([len(tab),len(out_names)])-9999.,names=out_names)
    # get rid of targetting proper motions to avoid confusion!
    tab.remove_columns(['PMRA','PMDEC','PM_SRC'])
    # add unpopulated columns
    tab.add_columns(newcols.columns.values())

    # read gaia file, match by coordinates, and populate
    gaia=fits.open(gaiafile)[1].data
    h=htm.HTM()
    maxrad=1./3600.
    m1,m2,rad=h.match(tab['RA'],tab['DEC'],gaia['RA'],gaia['DEC'],maxrad,maxmatch=1)
    for inname,outname in zip(in_names,out_names) :
        tab[outname][m1] = gaia[inname][m2]

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

def new(infile='allStar-r12-l33-58358.fits',new='allStar-r12-l33.fits',trim='allStarLite-r12-l33.fits') :
    """ take allStar file, add GAIA info and new _spec columns, outpu
        also output allStarLite version
    """
    hdulist=fits.open(infile)
    print('adding gaia....')
    tab=add_gaia(hdulist[1].data)
    print('adding _SPEC columns....')
    tab=add_spec(tab)

    # write out the modified file
    print('writing file ',new)
    out=fits.HDUList()
    out.append(fits.BinTableHDU(tab))
    out.append(hdulist[2])
    out.append(hdulist[3])
    out.writeto(new,overwrite=True)

    print('writing file ',trim)
    out=fits.HDUList()
    out.append(fits.BinTableHDU(trimfile(tab)))
    out.writeto(trim,overwrite=True)
