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

def addgaia(file='allStar-r12-l33-58358.fits',gaiafile='gaia.fits.gz',outfile='test.fits') :
    """ Add GAIA data to allStar file, with coordinate match to (cross-matched) GAIA reference file
    """
    hdulist=fits.open(file)
    data=hdulist[1].data
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

    # write out the modified file
    out=fits.HDUList()
    out.append(fits.BinTableHDU(tab))
    out.append(hdulist[2])
    out.append(hdulist[3])
    out.writeto(outfile,overwrite=True)
        
def trimfile(file='test.fits',outfile='trim.fits') :
    """ Write a 'lite' allStar file, removing some of the big space users
    """
    hdulist=fits.open(file)
    data=hdulist[1].data
    tab=Table(data)
    remove=['ALL_VISITS','VISITS','ALL_VISIT_PK','VISIT_PK','FPARAM_CLASS','CHI2_CLASS',
            'FPARAM','FELEM','FPARAM_COV','PARAM','PARAM_COV','ELEMFLAG','FELEM_ERR',
            'APSTAR_ID','TARGET_ID','ASPCAP_ID','FILE','LOCATION_ID']
    out=fits.HDUList()
    tab.remove_columns(remove)
    out.append(fits.BinTableHDU(tab))
    out.writeto(outfile,overwrite=True)

