#!/usr/bin/env python

"""APQUICK.PY - APOGEE Quick reduction of on-the-mountain data.

This performs a very quick reduction of an APOGEE exposure to
give feedback about the S/N.  This program can be run on an exposure
in progress or on that has already finished because it uses the apRaw
files (individual reads). 

Here are the main steps:
1) Use Fowler/CDS to collapse cube to a 2D image using very few reads.
    This is only done for the green chip.
2) Construct the noise model assuming we used up-the-ramp as is done
     by the full pipeline.
3) Boxcar extract all fibers but only for ~50 columns at the center
     of the image.
4) Fit a line to log(S/N) vs. Hmag and use it to find the S/N at the
     fiducial H magnitude.
5) Write out the results to a single multi-extension FITS file.

"""

from __future__ import print_function

__authors__ = 'David Nidever <dnidever@noao.edu>'
__version__ = '20180922'  # yyyymmdd                                                                                                                           

import os
import numpy as np
import warnings
from astropy.io import fits
from astropy.table import Table, Column
from glob import glob
from apogee.utils import yanny, apload
from sdss_access.path import path

# Ignore these warnings, it's a bug
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")


class Frame:
    """Object class for a single APOGEE read image.

    Parameters
    ----------
    imfile: str
          The filename of the FITS file.
    head : astropy header
         The header of the FITS file.
    im : numpy array
       The 2D raw image.
    num : int
       The read number of the image.

    Attributes
    ----------
    Same as the parameters.

    """

    def __init__(self,imfile,head,im,num):
        self.file = imfile
        self.head = head
        self.im = im
        self.num = num


class SpecFrame:
    """Object class for an extracted APOGEE spectrum.

    Parameters
    ----------
    flux : numpy array
        The 2D flux array [Nfibers,Ncolumns].
    err : numpy array
        The 2D error array [Nfibers,Ncolumns].
    head : astropy header
         The header of the FITS file (from the last read used).

    """

    def __init__(self,flux,err,head):
        self.flux = flux
        self.err = err
        self.head = head


# Load the frames
def loadframes(rawdir,framenum,nfowler=2,chip=2,lastread=None):
    """Loads apRaw reads at beginning and end of an exposure.

    This function loads apRaw reads for an exposure at the beginning
    and end to be used for Fowler/CDS sampling.
    
    Parameters
    ----------
    rawdis : str
           Directory for the apRaw files.
    framenum : str
          The 8 digit APOGEE frame number.
    nfowler : int
            Number of reads to load.  The default is 2.
    chip : int
            Which detector to use (1, 2, or 3).  The default is 2. 
    lastread : int
            The number of the last read to use.

    Returns
    -------
    bframes : list of Frame objects
            List of Frame objects with images and header information for the
            reads at the beginning of the exposure.  Note, the first read is
            always skipped.
    eframes : list of Frame objects
            List of Frame objects with images and header information for the
            reads at the end of the exposure.
    nreads : int
           The number of reads (files) used.

    Example
    -------

    .. code-block:: python

        bframes,eframes,nreads = loadframes('/data/apogee/spectro/raw/58382/',28200054,nfowler=2,chip=2)

    """

    # Get the file list
    files = glob(rawdir+"/apRaw-"+framenum+"-???.fits")
    nfiles = len(files)
    if nfiles==0:
        print("No files for "+framenum)
        return None
    # Sort the files
    files = np.sort(files)
    # Get the read numbers
    readnum = np.zeros(nfiles,dtype=int)
    for i in range(nfiles):
        base = os.path.basename(files[i])
        # apRaw-28190009-059.fits
        readnum[i] = np.int(base[15:18])
    # If readnum input then only use reads up to that number
    if lastread is not None:
        gdf, = np.where(readnum <= int(lastread))
        ngdf = len(gdf)
        if ngdf < 2:
            raise Exception("Not enough reads")
        # Only keep the files that we want to use
        files = files[gdf]
        readnum = readnum[gdf]
        nfiles = ngdf
    # What nfowler are we using
    nfowler_use = nfowler
    # Raise an error if we don't have enough reads for nfowler
    if nfiles<3:
        raise Exception("Not enough reads ("+str(nfiles)+")")
    if nfiles<(2*nfowler+1):
        nfowler_use = np.int(np.floor((nfiles-1)/2.))
        print("Not enough reads ("+str(nfiles)+") for Nfowler="+str(nfowler)+".  Using "+str(nfowler_use)+" instead")
    #if nfiles<nfowler:
    #    raise Exception("Not enough reads ("+str(nfiles)+") for Nfowler="+str(nfowler))
    # Load the begging set of frames
    #  skip the first one, bad
    bframes = []
    for i in range(nfowler_use):
        imfile = files[i+1]
        num = readnum[i+1]
        im,head = fits.getdata(imfile,header=True)
        if chip is not None:
            im = im[:,(chip-1)*2048:chip*2048]
        frame = Frame(imfile,head,im,num)
        bframes.append(frame)
    # Load the ending set of frames
    eframes= []
    for i in range(nfowler_use):
        imfile = files[nfiles-nfowler_use+i]
        num = readnum[nfiles-nfowler_use+i]
        im,head = fits.getdata(imfile,header=True)
        if chip is not None:
            im = im[:,(chip-1)*2048:chip*2048]
        frame = Frame(imfile,head,im,num)
        eframes.append(frame)
    # Return the lists
    return bframes,eframes,nfiles


# Perform Fowler sampling
def fowler(bframes,eframes):
    """Collapses exposure to 2D using Fowler/CDS sampling.

    This function performs Fowler/CDS sampling given a list of
    raw frames from the beginning and end of the exposure.
    
    Parameters
    ----------
    bframes : list of Frame objects
            The read images at the beginning of the exposure.
    eframes : list of Frame objects
            The read images at the end of the exposure.

    Returns
    -------
    im : numpy array
       The collapsed 2D image.


    Example
    -------

    .. code-block:: python

        im = fowler(bframes,eframes)

    """

    # Begining sample
    nbeg = len(bframes)
    nx,ny = bframes[0].im.shape
    if (nbeg==1):
        im_beg = np.array(bframes[0].im.copy(),dtype=float)
    else:
        im_beg = np.array(bframes[0].im.copy(),dtype=float)*0
        for i in range(nbeg): im_beg += np.array(bframes[i].im.copy(),dtype=float)
        im_beg /= np.float(nbeg)
    # End sample
    nend = len(eframes)
    if (nend==1):
        im_end = np.array(eframes[0].im.copy(),dtype=float)
    else:
        im_end = np.array(eframes[0].im.copy(),dtype=float)*0
        for i in range(nend): im_end += np.array(eframes[i].im.copy(),dtype=float)
        im_end /= np.float(nend)
    # The middle read will be used twice for 3 reads
    # Subtract beginning from end
    im = im_end - im_beg

    return im


# Construct the noise model images
def noisemodel(im,nreads,noise,gain):
    """Computes the noise model for an image.

    This function computes the noise image for the case where
    the data cube had been collapsed using UTR as the full
    reduction code does.
    
    Parameters
    ----------
    im : numpy array
       The 2D collapsed image.
    nreads : int
       The number of reads in the image.
    noise : float
       The read noise for a single read in ADU.
    gain : float
       The gain in electrons/ADU.

    Returns
    -------
    err : numpy array
        The noise model image in ADU.

    Example
    -------

    .. code-block:: python

        err = noisemodel(im,nreads,noise,gain)

    """

    # need nread, ngdreads, noise, gain
    ngdreads = np.float(nreads-1)

    ## READ NOISE
    #if n_elements(rdnoiseim) gt 0 then begin
    #  noise = median(rdnoiseim)
    #endif else begin
    #  noise = 12.0  ; default value
    #endelse

    # See Equation 1 in Rauscher et al.(2007), SPIE
    #  with m=1
    #  noise and image/flux should be in electrons, sample_noise is in electrons
    impos = im.copy()
    bad = (impos < 0)
    impos[bad] = 0
    sample_noise = np.sqrt( 12*(ngdreads-1)/(np.float(nreads)*(ngdreads+1))*noise**2 + 6.*(ngdreads**2+1)/(5*ngdreads*(ngdreads+1))*impos*gain )
    sample_noise /= gain  # convert to ADU

    # Start the variance array
    varim = im.copy()*0
    # Sample/read noise
    varim += sample_noise**2 

    #Now convert to ELECTRONS
    #if n_elements(detcorr) gt 0 and keyword_set(outelectrons) then begin
    #  varim *= gain^2
    #  im *= gain

    # Convert variance to error
    err = im.copy()*0+1
    good = (varim > 0)
    err[good] = np.sqrt(varim[good])

    return err


# Boxcar extract fibers
def boxextract(frame,tracestr,fibers=None,xlo=0,xhi=None):
    """This function performs boxcar extraction on a 2D image.
    
    Parameters
    ----------
    frame : Frame object
          The Frame object that includes the flux and error images.
    tracestr : numpy structured array
          Numpy structured array that gives information on the traces.
    fibers : array or list, optional
           List/array of fibers to extract.  By default all fibers
           in tracestr are extracted.
    xlo : int, optional
        The starting column to extract.  By default this is 0.
    xhi : int, optional
        The ending column to extract.  This is set to the final column
        of the 2D iamge.

    Returns
    -------
    spec : SpecFrame object
         The SpecFrame object that contains extracted flux and error arrays.

    Example
    -------

    .. code-block:: python

        spec = boxextract(frame,tracestr)

    """

    ntrace = len(tracestr)
    if fibers is None:
        fibers = np.arange(ntrace)
    if xhi is None:
        xhi = frame.im.shape[1]-1
    # transpose b/c x/y are flipped in python
    flux = frame.im.copy().T
    err = frame.err.copy().T
    nx,ny = flux.shape
    # Loop through the Fibers
    ncol = np.int(xhi-xlo+1)
    x = np.arange(xlo,xhi+1)
    nfibers = len(fibers)
    oflux = np.zeros((nfibers,ncol))
    oerr = np.zeros((nfibers,ncol))
    for i in range(nfibers):
        fwhm = tracestr['FWHM'][fibers[i]]
        coef = tracestr['COEF'][fibers[i]]
        ymid = np.polyval(coef[::-1],x)   # reverse order
        ylo = np.int(np.min(np.floor(ymid-fwhm)))
        if ylo<0: ylo=0
        yhi = np.int(np.max(np.ceil(ymid+fwhm)))
        if yhi>(ny-1): yhi=ny-1
        num = yhi-ylo+1

        # Make a MASK based on the trace and FWHM
        yy = np.resize(np.repeat(np.arange(ylo,yhi+1),ncol),(ncol,num))
        ymid2d = np.resize(np.repeat(ymid,num),(ncol,num))
        mask = np.array((yy >= np.floor(ymid2d-fwhm)) & (yy <= np.ceil(ymid2d+fwhm)),dtype=int)
        # Flux
        oflux[i,:] = np.sum( flux[xlo:xhi+1,ylo:yhi+1]*mask, axis=1)
        # Error
        #  add in quadrature
        oerr[i,:] = np.sqrt( np.sum( (err[xlo:xhi+1,ylo:yhi+1]**2)*mask, axis=1) )

    return SpecFrame(oflux,oerr,frame.head)


# Calculate mean sky spectrum
def skysub(spec,plugmap):
    """This subtracts the median sky spectrum from all of the fiber spectra.
    
    Parameters
    ----------
    spec : SpecFrame object
         The SpecFrame object that constrains the 1D extracted spectra.
    plugmap : numpy structured array
            The plugmap information for each fiber including which fiber contains
            sky or stars.

    Returns
    -------
    spec : SpecFrame object
         The same SpecFrame object but now with the sky spectrum subtracted.

    Example
    -------

    .. code-block:: python

        spec2 = skysub(spec,plugmap)

    """

    # Find the object and sky fibers
    pcat = plugmap['PLUGMAPOBJ']
    #ofibs, = np.where( (pcat['fiberId']>=0) & (pcat['holeType']=='OBJECT') & (pcat['spectrographId']==2) &
    #                   ( (pcat['objType']=='STAR') | (pcat['objType']=='HOT_STD') ) )
    sfibs, = np.where( (pcat['fiberId']>=0) & (pcat['holeType']=='OBJECT') & (pcat['spectrographId']==2) & (pcat['objType']=='SKY') )
    skyindex = 300-pcat[sfibs]['fiberId']
    # Calculate the median sky spectrum
    skyspec = np.median(spec.flux[skyindex,:],axis=0)
    # Subtract from all fibers
    nspec,ncol = spec.flux.shape
    for i in range(nspec): spec.flux[i,:] -= skyspec
    return spec


# Calculate median S/N per fiber
def snrcat(spec,plugmap):
    """This function calculates the S/N for each fiber.
    
    Parameters
    ----------
    spec : SpecFrame object
         The SpecFrame object that constrains the 1D extracted spectra.
    plugmap : numpy structured array
            The plugmap information for each fiber including which fiber contains
            sky or stars.

    Returns
    -------
    cat : numpy structured array
         A catalog containing information on each object in the fibers and the
         median S/N.

    Example
    -------

    .. code-block:: python

        cat = snrcat(spec,plugmap)

    """

    dtype = np.dtype([('apogee_id',np.str,30),('ra',np.float64),('dec',np.float64),('hmag',np.float),('objtype',np.str,30),
                      ('fiberid',np.int),('fiberindex',np.int),('flux',np.float),('err',np.float),('snr',np.float)])
    cat = np.zeros(300,dtype=dtype)

    # Load the spectral data
    cat['fiberindex'] = np.arange(300)
    cat['flux'] = np.median(spec.flux,axis=1)
    cat['err'] = np.median(spec.err,axis=1)
    err = cat['err']
    bad = (err <= 0.0)
    err[bad] = 1.0
    cat['snr'] = cat['flux']/err

    # Load the plugging data
    pcat = plugmap['PLUGMAPOBJ']
    fibs, = np.where( (pcat['fiberId']>=0) & (pcat['holeType']=='OBJECT') & (pcat['spectrographId']==2) )
    fiberindex = 300-pcat[fibs]['fiberId']
    cat['apogee_id'][fiberindex] = pcat[fibs]['tmass_style']
    cat['ra'][fiberindex] = pcat[fibs]['ra']
    cat['dec'][fiberindex] = pcat[fibs]['dec']
    cat['hmag'][fiberindex] = pcat[fibs]['mag'][:,1]
    cat['objtype'][fiberindex] = pcat[fibs]['objType']
    cat['fiberid'][fiberindex] = pcat[fibs]['fiberId']
    cat = Table(cat)

    return cat


# Linear fit to log(S/N) vs. H
def snrhmag(cat,nreads,nframes,hfid=12.2):
    """Returns the S/N for the fiducial H magnitude.

    This fits a line to log(S/N) vs. Hmag for the stars and
    predicts what the S/N should be at the end of the exposure.
    
    Parameters
    ----------
    cat : numpy structured array
         A catalog containing information on each object in the fibers and the
         median S/N.
    nreads: int
         The number of reads used in the current processing.
    nframes: int
         The total number of reads in the entire exposure.
    hfid : float, optional
         The fiducial Hmag.  The default is 12.2.

    Returns
    -------
    coefstr : numpy structured array
            An structure that contains all of the information computed:
            hmag_fid : The fiducial Hmag
            logsnr_hmag_coef : The linear coefficients of the fit of log(S/N) vs. Hmag
            snr_fid : The S/N at the fiducial Hmag using the linear fit.
            snr_predict : The predicted S/N at the fiducial Hmag at the end of the expsure.

    Example
    -------

    .. code-block:: python

        coefstr = snrhmag(cat,30,47)

    """

    gd, = np.where( (cat['objtype'] != 'SKY') & (cat['hmag'] > 4) & (cat['hmag'] < 20) & (cat['snr'] > 0) )
    coef = np.polyfit(cat[gd]['hmag'],np.log10(cat[gd]['snr']),1)
    snr_fid = 10**np.polyval(coef,hfid)
    # Predicted S/N at end of exposure
    #  (S/N)^2 should scale with time
    snr_predict = np.sqrt( snr_fid**2*np.float(nframes)/np.float(nreads) )

    dtype = np.dtype([('hmag_fid',np.float),('snr_fid',np.float),('logsnr_hmag_coef',(np.float,2)),('snr_predict',np.float)])
    coefstr = np.zeros(1,dtype=dtype)
    coefstr['hmag_fid'] = hfid
    coefstr['snr_fid'] = snr_fid
    coefstr['logsnr_hmag_coef'] = coef
    coefstr['snr_predict'] = snr_predict
    coefstr = Table(coefstr)

    return coefstr


# Run everything
def runquick(framenum,observatory='apo',lastread=None,hfid=12.2):
    """This runs all of the main steps of the quick reduction.
    
    Parameters
    ----------
    framenum : string
             The APOGEE frame number of the exposure to process.
    observatory : string
               The observatory.  Either "apo" or "lco".
    lastread : int or string
             The number for the last read to use.
    hfid : float, optional
         The fiducial Hmag.  The default is 12.2.

    Returns
    -------
    frame : Frame object
          The extracted 2D image and noise model.
    spec : SpecFrame object
         The extracted spectra and errors with the sky subtracted.
    cat : numpy structured array
        The catalog of information on each fiber including the S/N.
    coefstr : numpy structured array
            Information on the fit of log(S/N) vs. Hmag and the
            S/N at the fiducial Hmag.

    Example
    -------

    .. code-block:: python

        frame, spec, cat, coefstr = runquick('28200054','apo',10)

    """

    mjd = apload.cmjd(np.int(framenum))
    if lastread is None:
        print('Running APQUICK on '+str(framenum)+' MJD='+str(mjd)+" all reads ")
    else:
        print('Running APQUICK on '+str(framenum)+' MJD='+str(mjd)+" "+str(lastread)+" reads")

    # Get apRaw directory
    sdss_path = path.Path()
    rawdir = sdss_path.dir('apRaw',num=framenum,read=1,mjd=mjd)
    if observatory=='apo':
        prefix = 'ap'
    else:
        prefix = 'as'
    psfdir = sdss_path.dir('apPSF',apred='quickred',prefix=prefix,chip='b',num=1)
    plugdir = "/data-ql/plugmaps/"
    detdir = sdss_path.dir('apDetector',apred='quickred',prefix=prefix,chip='b',num=1)

    # Load the reads
    bframes,eframes,nreads = loadframes(rawdir,framenum,2,lastread=lastread)
    head = eframes[-1].head.copy()   # header of last read
    nframes = np.int(head['NFRAMES'])
    # Do Fowler/CDS collapse
    im = fowler(bframes,eframes)
    # Get rdnoise/gain from apDetector file
    detfiles = glob(detdir+'/apDetector-b-????????.fits')
    detfiles = np.sort(detfiles)
    print('Using '+detfiles[-1])
    rdnoiseim = fits.getdata(detfiles[-1],1)
    rdnoise = np.median(rdnoiseim)
    gainim = fits.getdata(detfiles[-1],2)
    gain = np.median(gainim)
    # Generate the noise image
    err = noisemodel(im,nreads,rdnoise,gain)
    frame = Frame("",head,im,0)
    frame.err = err
    frame.head['FRAMENUM'] = framenum
    # Load the trace information
    psffiles = np.sort(glob(psfdir+'/apPSF-b-????????.fits'))
    print('Using '+psffiles[-1])
    tracestr = Table.read(psffiles[-1],1)
    # Boxcar extract the fibers
    spec = boxextract(frame,tracestr,xlo=950,xhi=1000)
    # Load the plugmap file
    plugfile = plugdir+'plPlugMapA-'+head['name']+'.par'
    if os.path.exists(plugfile) is False:
        raise Exception(plugfile+' NOT FOUND')
    plugmap = yanny.yanny(plugfile,np=True)
    # Subtract the sky
    spec = skysub(spec,plugmap)
    # Create the S/N catalog
    cat = snrcat(spec,plugmap)
    # Fit line and get fiducial S/N
    coefstr = snrhmag(cat,nreads,nframes,hfid=hfid)
    # Print out the S/N values
    print('S/N (H=%5.2f) = %5.2f (%3d reads)' % (hfid,coefstr['snr_fid'],nreads))
    print('S/N (H=%5.2f) = %5.2f (prediction for %3d reads)' % (hfid,coefstr['snr_predict'],nframes))

    return frame, spec, cat, coefstr


# Write out results
def writeresults(outfile, frame, spec, cat, coefstr):
    """This writes out the results of the quick reduction to an output file.
    
    Parameters
    ----------
    outfile : str
          The output filename to write the results to.
    frame : Frame object
          The extracted 2D image and noise model.
    spec : SpecFrame object
         The extracted spectra and errors with the sky subtracted.
    cat : numpy structured array
        The catalog of information on each fiber including the S/N.
    coefstr : numpy structured array
            Information on the fit of log(S/N) vs. Hmag and the
            S/N at the fiducial Hmag.
             
    Returns
    -------
    The function doesn't return anything but writes the results to an
    output file.

    Example
    -------

    .. code-block:: python

        writeresults('apq-28200054.fits', frame, spec, cat, coefstr)

    """

    print('Final catalog = '+outfile)
    Table(coefstr).write(outfile,overwrite=True)

    # THE CODE BELOW DOES *NOT* CURRENTLY WORK AT APO
    # ASTROPY NEEDS TO BE UPDATED

    ## HDU0: header only
    #head = frame.head
    #head.add_history("APQUICK results for "+str(head['FRAMENUM']))
    #head.add_history("HDU0: Header only")
    #head.add_history("HDU1: coefstr, S/N for fiducial Hmag")
    #head.add_history("HDU2: cat, catalog of S/N values")
    #head.add_history("HDU3: image/err, numpy structured array")
    #head.add_history("HDU4: spec/err, numpy structured array")
    #fits.writeto(outfile,None,head,overwrite=True)
    #hdulist = fits.open(outfile)
    ## HDU1: Coef structure
    #hdu = fits.table_to_hdu(Table(coefstr))
    #hdu.header.add_history("APQUICK results for "+str(head['FRAMENUM']))
    #hdu.header.add_history("Final fiducial S/N values")
    #hdulist.append(hdu)
    ## HDU2: Catalog
    #hdu = fits.table_to_hdu(cat)
    #hdu.header.add_history("APQUICK results for "+str(head['FRAMENUM']))
    #hdu.header.add_history("Catalog of S/N values and other fiber information")
    #hdulist.append(hdu)
    ## HDU3: Spectra
    #dtype = np.dtype([('flux',(np.float32,spec_flux.shape)),('err',(np.float32,spec_err.shape))])
    #spec = np.zeros(1,dtype=dtype)
    #spec['flux'] = spec.flux
    #spec['err'] = spec.err
    #spec = Table(spec)
    #hdu = fits.table_to_hdu(spec)
    #hdu.header.add_history("APQUICK results for "+str(head['FRAMENUM']))
    #hdu.header.add_history("Extracted spectra and errors (ADU)")
    #hdulist.append(hdu)
    ## HDU4: Image
    #dtype = np.dtype([('im',(np.float32,im.shape)),('err',(np.float32,err.shape))])
    #image = np.zeros(1,dtype=dtype)
    #image['im'] = im
    #image['err'] = err
    #image = Table(image)
    #hdu = fits.table_to_hdu(image)
    #hdulist.append(hdu)
    #hdu.header.add_history("APQUICK results for "+str(head['FRAMENUM']))
    #hdu.header.add_history("Collapsed 2D image and noise model (ADU)")
    #hdulist.writeto(outfile,overwrite=True)
    #hdulist.close() 


# Main command-line program
#if __name__ == "__main__":
def main(args) :

    from argparse import ArgumentParser
    parser = ArgumentParser(description='The APOGEE Quick Reduction/Look software that runs on a single exposure.')
    parser.add_argument('framenum', nargs='?', help='APOGEE 8-digit frame/exposure number, e.g. 28200054')
    parser.add_argument('observatory', nargs='?', default='apo', help='Observatory: apo or lco')
    parser.add_argument('lastread', nargs='?', default=None, help='OPTIONAL. Final read to use')
    parser.add_argument('--outfile', type=str, default=None, help='Output filename')
    parser.add_argument('--hfid', type=int, default=12.2, help='Fiducial H-magnitude')
    args = parser.parse_args()

    # Run everything
    frame, spec, cat, coefstr = runquick(args.framenum,args.observatory,args.lastread,hfid=args.hfid)

    # Write to file
    if args.outfile is None:
        if args.lastread is None:
            outfile = ('apq-%8d.fits' % int(args.framenum))
        else:
            outfile = ('apq-%8d_%03d.fits' % (int(args.framenum),int(args.lastread)))
    else:
        outfile = args.outfile
    writeresults(outfile, frame, spec, cat, coefstr)
