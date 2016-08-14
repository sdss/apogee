# routines for reduction of ARCES data using PYRAF

import glob
from astropy.io import fits
from pyraf import iraf
from pyraf.iraf import noao,imred,twodspec,apextract,onedspec,echelle
from holtz.pyvista import imred
import pdb
import os

no='no'
yes='yes'
INDEF='INDEF'

def process(root,dir='./',verbose=True) :
    # setup
    imred.setup(root,dir=dir,idet=16)

    # create bias frames
    biases = imred.getfiles('zero',verbose=verbose)
    bias = imred.combine(biases,verbose=verbose)
    
    # create flat field
    flats = imred.getfiles('flat',filter='Open',verbose=verbose)
    flat = imred.combine(flats,bias=bias,trim=True,verbose=verbose)
    flats = imred.getfiles('flat',filter='Blue',verbose=verbose)
    flat += imred.combine(flats,bias=bias,trim=True,verbose=verbose)
    hdu=fits.PrimaryHDU(flat)
    hdu.writeto('flat.fits',clobber=True)
    # create normalized 1D flat field
    for file in ['flat_mag.fits','flat_mag.ec.fits','flat.ec.fits','norm_flat.ec.fits'] :
        if os.path.isfile(file) : os.remove(file)
    iraf.magnify(input='flat.fits',output='flat_mag.fits',xmag=1,ymag=4)
    iraf.hedit('flat_mag.fits',fields='CCDSEC',value='[200:1850,1:8189]',
           add=no,addonly=no,delete=no,verify=no,show=yes,update=yes)
    iraf.hedit('flat_mag.fits','dispaxis','1',add=yes,verify=no,show=yes,update=yes)
    if verbose: print( 'Modeling and extracting the superflat...')
    iraf.apall(input='flat_mag.fits',ref='echtrace130522',format='echelle',
           interactive=no,find=no,recenter=yes,resize=yes,edit=no,trace=yes,
           fittrace=no,extract=yes,extras=no,review=no,line=825,nsum=10,
           lower=-14.0,upper=14.0,b_function='chebyshev',b_order=2,
           b_niterate=3,b_naverage=-3,b_sample='-22:-15,15:22',width=18.0,
           radius=18.0,npeaks=INDEF,shift=no,ylevel=.05,t_nsum=5,t_step=1,
           t_function='legendre',t_order=5,t_naverage=3,t_niterate=3,
           t_low_reject=2.5,t_high_reject=2.5,t_nlost=3,t_sample='*',
           background='fit',weights='none',clean=no,saturation=40000.0)
    flat=fits.open('flat_mag.ec.fits')[0].data
    specflat = imred.specflat(flat,indiv=True)
    hdu=fits.PrimaryHDU(specflat)
    hdu.writeto('normflat.ec.fits',clobber=True)
    #iraf.sfit (input="flat_mag.ec.fits",output="norm_flat_mag.ec.fits",type="ratio",
    #       replace=no,wavesca=no,logscal=no,override=yes,interac=no,sample="*",
    #       naverag=1,funct="spline3",order=5,low_rej=2,high_rej=0,niterate=10,grow=1)

    pdb.set_trace()

    # arcs
    arcs = imred.getfiles('comp',verbose=verbose,listfile='arcs.lis')
    for arc in arcs :
        if verbose: print('arc: ', arc)
        data=imred.reduce(arc,bias=bias,trim=True)
        arcfile='arc.{:04d}.fits'.format(arc)
        outfile='arc.{:04d}.ec.fits'.format(arc)
        data.writeto(arcfile,clobber=True)
        if verbose: print('Resampling the arc by a factor of 4 in the y direction...')
        iraf.magnify(input=arcfile,output=arcfile,xmag=1,ymag=4)
        iraf.hedit(images=arcfile,fields='CCDSEC',value='[200:1850,1:8189]',add=no,
               addonly=no,delete=no,verify=no,show=yes,update=yes)
        iraf.hedit(arcfile,'dispaxis','1',add=yes,verify=no,show=yes,update=yes)
        if verbose: print('Applying model apertures to the arc and extracting spectra...')
        if os.path.isfile(outfile) : os.remove(outfile)
        iraf.apall(input=arcfile,output=outfile,reference='flat_mag',format='echelle',
               interactive=no,find=no,recenter=no,resize=no,edit=no,
               trace=no,fittrace=no,extract=yes,extras=no,review=no,
               shift=no,background='none',weights='none')
        iraf.ecreidentify(outfile,reference='arcnewref.ec',
                  shift=INDEF,cradius=2,threshold=45,refit=yes)

    pdb.set_trace()

    # loop over objects
    objs = imred.getfiles('object',verbose=verbose)
    for obj in objs :
        data=imred.reduce(obj,bias=bias,flat=specflat,trim=True)
    pdb.set_trace()

def test() :
    process('131102',dir='/home/holtz/raw/apo/nov13/131102/Echelle')
