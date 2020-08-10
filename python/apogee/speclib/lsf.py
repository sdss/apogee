###############################################################################
# apogee.spec.lsf: Utilities to work with APOGEE LSFs
# routines from github jobovy, modified by holtz to work in apogee product context
###############################################################################
import os, os.path
from functools import wraps
import warnings
import math
import numpy
from scipy import special, interpolate, sparse, ndimage
import scipy.sparse.linalg
import time
import sys
import pdb
from apogee.utils import spectra
from apogee.aspcap import aspcap
from astropy.io import fits
import matplotlib.pyplot as plt

from apogee.utils import apload
try: load=apload.ApLoad()
except: print('cant initialize load in lsf...')

_SQRTTWO= numpy.sqrt(2.)


def showtime(string) :
    """ Utiltiy routine to print a string and clock time, with flush to stdout

    Args:
        string (str) : string to print with clock time
    """
    print(string+' {:8.2f}'.format(time.time()))
    sys.stdout.flush()

def get(lsfid,waveid,fiber,highres=9,apred=None,telescope=None) :
    """  Return standard sparsified LSF

    Args:
        lsfid (int) : ID of apLSF file
        waveid (int) : ID of apWave file
        fiber (int) : fiber(s) of LSF
        highres (int) : number of subpixels for LSF calculation (default=9)
        apred (str) : apred version to get apLSF and apWave fromn (default=None --> uses default from apload)
    Returns :
        x (np.array) : array of relative pixel locations for LSF
        l () : output LSF

    """
    if apred is not None: load.apred = apred
    if telescope is not None: load.settelescope(telescope)
    x=numpy.arange(-15.,15.01,1./highres)
    x=numpy.arange(-7.,7.01,1./highres)
    l=eval(x,fiber=fiber,waveid=waveid,lsfid=lsfid)
    return x,l

def convolve(wav,spec,
             lsf=None,xlsf=None,dxlsf=None,fiber='combo',
             vmacro=6.,vrot=None, highout=False):
    """ convolve an input spectrum with APOGEE LSF and resample to APOGEE's apStar wavelength grid
    Args:
       wav - wavelength array (linear in wavelength in \AA)
       spec - spectrum on wav wavelength grid [nspec,nwave]
       lsf= (None) pre-calculated LSF array from apogee.spec.lsf.eval
       Either:
          xlsf= (None) 1/integer equally-spaced pixel offsets at which the lsf=lsf input is calculated
          dxlsf= (None) spacing of pixel offsets
       fiber= if lsf is None, the LSF is calculated for this fiber
       vmacro= (6.) Gaussian macroturbulence smoothing to apply as well (FWHM or a [sparse] matrix like lsf on the same x grid; can be computed with apogee.modelspec.vmacro)
    Returns :
       spectrum on apStar wavelength grid
    HISTORY:
       2015-03-14 - Written - Bovy (IAS)
    """
    # Parse LSF input
    if lsf is None:
        xlsf= numpy.linspace(-7.,7.,43)
        lsf= eval(xlsf,fiber=fiber)
    if dxlsf is None:
        dx= xlsf[1]-xlsf[0]
    else:
        dx= dxlsf
    hires= int(round(1./dx))
    waveout=aspcap.apStarWave()
    if wav[0]-waveout[0] > 10 :
        # if first wavelength is far from first apStar wavelength, minigrid: trim the output and LSF 
        gd = numpy.where((waveout > wav[0]) & (waveout < wav[-1]) )[0]
    else :
        gd=numpy.arange(len(waveout))
    l10wav= numpy.log10(waveout[gd])
    dowav= l10wav[1]-l10wav[0]
    tmpwav= 10.**numpy.arange(l10wav[0],l10wav[-1]+dowav/hires,dowav/hires)
    tmp= numpy.empty(len(l10wav)*hires)   
    # for minigrid, extract same length for LSF
    if wav[0]-waveout[0] > 10 :
        lsf = lsf[gd[0]*hires:gd[0]*hires+len(tmpwav),:]

    # sparsify, need to do after trimming lsf for minigrid
    if not isinstance(lsf,sparse.dia_matrix):
        lsf= sparsify(lsf)

    # Interpolate the input spectrum, starting from a polynomial baseline
    if len(spec.shape) == 1: spec= numpy.reshape(spec,(1,len(spec)))
    nspec= spec.shape[0]
    tmp= numpy.empty((nspec,len(tmpwav)))
    for ii in range(nspec):
        baseline= numpy.polynomial.Polynomial.fit(wav,spec[ii],4)
        ip= interpolate.InterpolatedUnivariateSpline(wav,
                                                     spec[ii]/baseline(wav),
                                                     k=3)
        tmp[ii]= baseline(tmpwav)*ip(tmpwav)
    # Add macroturbulence, allowing for different kernel for each spectrum
    if not vmacro is None :
        sigvm= vmacro/3./10.**5./numpy.log(10.)*hires/dowav/2./numpy.sqrt(2.*numpy.log(2.))
        if isinstance(vmacro,float): sigvm=numpy.tile(sigvm,(nspec,1))
        for ii in range(nspec) :
            tmp[ii,:]= ndimage.gaussian_filter1d(tmp[ii,:],sigvm[ii],mode='constant')
    # Add rotation
    if vrot is not None :
        deltav=dowav/hires*3.e5*numpy.log(10)
        print('rot: ',deltav,vrot)
        kernel=rotate(deltav,vrot,epsilon=0.25)
        kernel=sparsify(numpy.tile(kernel,(len(tmpwav),1)))
        print(kernel.shape)
        tmp=kernel.dot(tmp.T).T
    if not isinstance(tmp,sparse.csr_matrix):
        # Use sparse representations to quickly calculate the convolution
        tmp= sparse.csr_matrix(tmp)

    if highout :
        return lsf.dot(tmp.T).T.toarray(),waveout[gd]
    else :
        return lsf.dot(tmp.T).T.toarray()[:,::hires],waveout[gd]

def sparsify(lsf):
    """convert an LSF matrix calculated with eval [ncen,npixoff] to a sparse [ncen,ncen] matrix with the LSF on the diagonals (for quick convolution with the LSF)
    Args:
       lsf - lsf matrix [ncen,npixoff] calculated by eval
    Returns :
       sparse matrix with the lsf on the diagonals
    HISTORY:
       2015-03-14 - Written - Bovy (IAS)
    """
    nx= lsf.shape[1]
    diagonals= []
    offsets= []
    for ii in range(nx):
        offset= nx//2-ii
        offsets.append(offset)
        if offset < 0:
            diagonals.append(lsf[:offset,ii])
        else:
            diagonals.append(lsf[offset:,ii])
    return sparse.diags(diagonals,offsets)

def dummy(dx=1./3.,sparse=False):
    """ return a 'dummy' LSF that is a delta function
    Args:
       dx= (1/3) spacing between LSF centers in the apStar grid
       sparse= (False) if True, return a sparse representation that can be passed to apogee.spec.lsf.convolve for easy convolution
    Returns:
       LSF(x|pixel center);
       pixel centers are apStarWave if dx=1, and denser 1/integer versions if dx=1/integer
    HISTORY:
       2015-03-23 - Written - Bovy (IAS)
    """
    # Are the x unit pixels or a fraction 1/hires thereof?
    hires= int(round(1./dx))
    # Setup output
    wav= aspcap.apStarWave()
    l10wav= numpy.log10(wav)
    dowav= l10wav[1]-l10wav[0]
    # Hi-res wavelength for output
    hireswav= 10.**numpy.arange(l10wav[0],l10wav[-1]+dowav/hires,dowav/hires)
    out= numpy.ones((len(hireswav),1))
    if sparse: out= sparsify(out)
    return out

def eval(x,fiber='combo',lsfid=5440020,waveid=2420038,sparse=False):
    """ evaluate the LSF for a given fiber
    Args :
       x - Array of X values for which to compute the LSF, in pixel offset relative to pixel centers; the LSF is calculated at the x offsets for each pixel center; x need to be 1/integer equally-spaced pixel offsets
       fiber= ('combo') fiber number or 'combo' for an average LSF (uses the same one-based indexing as the APOGEE fibers [i.e., fibers range from 1 to 300])
       sparse= (False) if True, return a sparse representation that can be passed to apogee.spec.lsf.convolve for easy convolution
    Returns :
       LSF(x|pixel center);
       pixel centers are apStarWave if dx=1, and denser 1/integer versions if dx=1/integer
    HISTORY:
       2015-03-12 - Written based on Jon H's code (based on David N's code) - Bovy (IAS)
    """
    # Parse fiber input
    if (isinstance(fiber,str) or isinstance(fiber,unicode)) and fiber.lower() == 'combo':
        fiber= [50,100,150,200,250]
    elif isinstance(fiber,int):
        fiber= [fiber]
    elif not isinstance(fiber,list) and isinstance(fiber[0],int):
        raise ValueError('fiber input to apogee.spec.lsf.eval not understood ...')
    # Are the x unit pixels or a fraction 1/hires thereof?
    hires= int(round(1./(x[1]-x[0])))
    # Setup output
    wav= aspcap.apStarWave()
    l10wav= numpy.log10(wav)
    dowav= l10wav[1]-l10wav[0]
    # Hi-res wavelength for output
    hireswav= 10.**numpy.arange(l10wav[0],l10wav[-1]+dowav/hires,dowav/hires)
    out= numpy.zeros((len(hireswav),len(x)))
    lsfpars=load.apLSF(lsfid,hdu=0)[0]
    for chip in ['a','b','c']:
        # Get pixel array for this chip, use fiber[0] for consistency if >1 fib
        pix= wave2pix(hireswav,chip,fiber=300-fiber[0],waveid=waveid)
        dx= numpy.roll(pix,-hires,)-pix
        dx[-1]= dx[-1-hires]
        dx[-2]= dx[-2-hires]
        dx[-3]= dx[-3-hires]
        xs= numpy.tile(x,(len(hireswav),1))\
            *numpy.tile(dx,(len(x),1)).T # nwav,nx       
        gd= True^numpy.isnan(pix)
        # Read LSF file for this chip
        #lsfpars= apread.apLSF(chip,ext=0)
        # Loop through the fibers
        for fib in fiber:
            out[gd]+= raw(xs[gd],pix[gd],lsfpars[chip][:,300-fib])
    out[out<0.]= 0.
    out/= numpy.tile(numpy.sum(out,axis=1),(len(x),1)).T
    if sparse: out= sparsify(out)
    return out

def raw(x,xcenter,params,nowings=False):
    """ Evaluate the raw APOGEE LSF (on the native pixel scale)
    Args :
       x - Array of X values for which to compute the LSF (in pixel offset relative to xcenter; the LSF is calculated at the x offsets for each xcenter if x is 1D, otherwise x has to be [nxcenter,nx]))
       xcenter - Position of the LSF center (in pixel units)
       lsfarr - the parameter array (from the LSF HDUs in the APOGEE data products)
    Returns :
       LSF(x|xcenter))
    HISTORY:
       2015-02-26 - Written based on Nidever's code in apogeereduce - Bovy (IAS)
    """
    # Parse x
    if len(x.shape) == 1:
        x= numpy.tile(x,(len(xcenter),1))
    # Unpack the LSF parameters
    params= unpack_lsf_params(params)
    # Get the wing parameters at each x
    wingparams= numpy.empty((params['nWpar'],len(xcenter)))
    for ii in range(params['nWpar']):
        poly= numpy.polynomial.Polynomial(params['Wcoefs'][ii])       
        wingparams[ii]= poly(xcenter+params['Xoffset'])
    # Get the GH parameters at each x
    ghparams= numpy.empty((params['Horder']+2,len(xcenter)))

    # note that this is modified/corrected a bit from Bovy's routines based on comparison
    # with LSF from IDL routines, noticeable when wings are non-negligible
    for ii in range(params['Horder']+2):
        if ii == 1:
            ghparams[ii]= 1.
        else:
            poly= numpy.polynomial.Polynomial(params['GHcoefs'][ii-(ii > 1)])
            ghparams[ii]= poly(xcenter+params['Xoffset'])
        # normalization
        if ii > 0: 
            ghparams[ii]/= numpy.sqrt(2.*numpy.pi*math.factorial(ii-1))
            if not nowings: ghparams[ii] *= (1.-wingparams[0])
    # Calculate the GH part of the LSF
    out= _gausshermitebin(x,ghparams,params['binsize'])
    # Calculate the Wing part of the LSF
    if not nowings: out+= _wingsbin(x,wingparams,params['binsize'],params['Wproftype'])
    return out

def _gausshermitebin(x,params,binsize):
    """ Evaluate the integrated Gauss-Hermite function
    """
    ncenter= params.shape[1]
    out= numpy.empty((ncenter,x.shape[1]))
    integ= numpy.empty((params.shape[0]-1,x.shape[1]))
    for ii in range(ncenter):
        poly= numpy.polynomial.HermiteE(params[1:,ii])
        # Convert to regular polynomial basis for easy integration
        poly= poly.convert(kind=numpy.polynomial.Polynomial)
        # Integrate and add up
        w1= (x[ii]-0.5*binsize)/params[0,ii]
        w2= (x[ii]+0.5*binsize)/params[0,ii]
        eexp1= numpy.exp(-0.5*w1**2.)
        eexp2= numpy.exp(-0.5*w2**2.)
        integ[0]= numpy.sqrt(numpy.pi/2.)\
            *(special.erf(w2/_SQRTTWO)-special.erf(w1/_SQRTTWO))
        out[ii]= poly.coef[0]*integ[0]
        if params.shape[0] > 1:
            integ[1]= -eexp2+eexp1
            out[ii]+= poly.coef[1]*integ[1]
        for jj in range(2,params.shape[0]-1):
            integ[jj]= (-w2**(jj-1)*eexp2+w1**(jj-1)*eexp1)\
                +(jj-1)*integ[jj-2]
            out[ii]+= poly.coef[jj]*integ[jj]
    return out

def _wingsbin(x,params,binsize,Wproftype):
    """Evaluate the wings of the LSF
    """
    ncenter= params.shape[1]
    out= numpy.empty((ncenter,x.shape[1]))
    for ii in range(ncenter):
        if Wproftype == 1: # Gaussian
            w1=(x[ii]-0.5*binsize)/params[1,ii]
            w2=(x[ii]+0.5*binsize)/params[1,ii]
            out[ii]= params[0,ii]/2.*(special.erf(w2/_SQRTTWO)\
                                          -special.erf(w1/_SQRTTWO))
    return out

def unpack_lsf_params(lsfarr):
    """ Unpack the LSF parameter array into its constituents
    Args :
       lsfarr - the parameter array
    Returns :
       dictionary with unpacked parameters and parameter values:
          binsize: The width of a pixel in X-units
          Xoffset: An additive x-offset; used for GH parameters that vary globally
          Horder: The highest Hermite order
          Porder: Polynomial order array for global variation of each LSF parameter
          GHcoefs: Polynomial coefficients for sigma and the Horder Hermite parameters
          Wproftype: Wing profile type
          nWpar: Number of wing parameters
          WPorder: Polynomial order for the global variation of each wing parameter          
          Wcoefs: Polynomial coefficients for the wings parameters
    HISTORY:
       2015-02-15 - Written based on Nidever's code in apogeereduce - Bovy (IAS@KITP)
    """
    out= {}
    # binsize: The width of a pixel in X-units
    out['binsize']= lsfarr[0]
    # X0: An additive x-offset; used for GH parameters that vary globally
    out['Xoffset']= lsfarr[1]
    # Horder: The highest Hermite order
    out['Horder']= int(lsfarr[2])
    # Porder: Polynomial order array for global variation of each LSF parameter
    out['Porder']= lsfarr[3:out['Horder']+4]
    out['Porder']= out['Porder'].astype('int')
    nGHcoefs= numpy.sum(out['Porder']+1)
    # GHcoefs: Polynomial coefficients for sigma and the Horder Hermite parameters
    maxPorder= numpy.amax(out['Porder'])
    GHcoefs= numpy.zeros((out['Horder']+1,maxPorder+1))
    GHpar= lsfarr[out['Horder']+4:out['Horder']+4+nGHcoefs] #all coeffs
    CoeffStart= numpy.hstack((0,numpy.cumsum(out['Porder']+1)))
    for ii in range(out['Horder']+1):
        GHcoefs[ii,:out['Porder'][ii]+1]= GHpar[CoeffStart[ii]:CoeffStart[ii]+out['Porder'][ii]+1]
    out['GHcoefs']= GHcoefs
    # Wproftype: Wing profile type
    wingarr= lsfarr[3+out['Horder']+1+nGHcoefs:]
    out['Wproftype']= int(wingarr[0])
    # nWpar: Number of wing parameters
    out['nWpar']= int(wingarr[1])
    # WPorder: Polynomial order for the global variation of each wing parameter
    out['WPorder']= wingarr[2:2+out['nWpar']]
    out['WPorder']= out['WPorder'].astype('int')
    # Wcoefs: Polynomial coefficients for the wings parameters
    maxWPorder= numpy.amax(out['WPorder'])
    Wcoefs= numpy.zeros((out['nWpar'],maxWPorder+1))
    Wpar= wingarr[out['nWpar']+2:]
    WingcoeffStart= numpy.hstack((0,numpy.cumsum(out['WPorder']+1)))
    for ii in range(out['nWpar']):
        Wcoefs[ii,:out['WPorder'][ii]+1]= Wpar[WingcoeffStart[ii]:WingcoeffStart[ii]+out['WPorder'][ii]+1]
    out['Wcoefs']= Wcoefs
    return out

def scalarDecorator(func):
    """Decorator to return scalar outputs for wave2pix and pix2wave
    """
    @wraps(func)
    def scalar_wrapper(*args,**kwargs):
        if numpy.array(args[0]).shape == ():
            scalarOut= True
            newargs= (numpy.array([args[0]]),)
            for ii in range(1,len(args)):
                newargs= newargs+(args[ii],)
            args= newargs
        else:
            scalarOut= False
        result= func(*args,**kwargs)
        if scalarOut:
            return result[0]
        else:
            return result
    return scalar_wrapper

@scalarDecorator
def wave2pix(wave,chip,fiber=300,waveid=2420038):
    """ convert wavelength to pixel
    Args :
       wavelength - wavelength (\AA)
       chip - chip to use ('a', 'b', or 'c')
       fiber= (300) fiber to use the wavelength solution of
    Returns :
       pixel in the chip
    HISTORY:
        2015-02-27 - Written - Bovy (IAS)
    """
# Load wavelength solutions
    wave0 =load.apWave(waveid,hdu=2)[0][chip][300-fiber]
    pix0= numpy.arange(len(wave0))
    # Need to sort into ascending order
    sindx= numpy.argsort(wave0)
    wave0= wave0[sindx]
    pix0= pix0[sindx]
    # Start from a linear baseline
    baseline= numpy.polynomial.Polynomial.fit(wave0,pix0,1)
    ip= interpolate.InterpolatedUnivariateSpline(wave0,pix0/baseline(wave0),
                                                 k=3)
    out= baseline(wave)*ip(wave)
    # NaN for out of bounds
    out[wave > wave0[-1]]= numpy.nan
    out[wave < wave0[0]]= numpy.nan
    return out

@scalarDecorator
def pix2wave(pix,chip,fiber=300,waveid=2420038):
    """ convert pixel to wavelength
    Args :
       pix - pixel
       chip - chip to use ('a', 'b', or 'c')
       fiber= (300) fiber to use the wavelength solution of
    Returns :
       wavelength in \AA
    HISTORY:
        2015-02-27 - Written - Bovy (IAS)
    """
    wave0 =load.apWave(waveid,hdu=2)[0][chip][300-fiber]
    pix0= numpy.arange(len(wave0))
    # Need to sort into ascending order
    sindx= numpy.argsort(pix0)
    wave0= wave0[sindx]
    pix0= pix0[sindx]
    # Start from a linear baseline
    baseline= numpy.polynomial.Polynomial.fit(pix0,wave0,1)
    ip= interpolate.InterpolatedUnivariateSpline(pix0,wave0/baseline(pix0), k=3)
    out= baseline(pix)*ip(pix)
    # NaN for out of bounds
    out[pix < 0]= numpy.nan
    out[pix > 2047]= numpy.nan
    return out

def _load_precomp(dr=None,fiber='combo',sparse=True):
    """Load a precomputed LSF
    """
    if dr is None: dr= appath._default_dr()
    if dr is 'current':
        warnings.warn("Preloaded LSFs for current DR not yet available, falling back on DR12 files")
        dr= '12'
    fileDir= os.path.dirname(appath.apLSFPath('a',dr=dr))
    filePath= os.path.join(fileDir,'apogee-lsf-dr%s-%s.fits' % (dr,fiber))
    # Download the file if necessary
    if not os.path.exists(filePath):
        dlink= \
            filePath.replace(fileDir,'https://zenodo.org/record/16147/files')
        _download_file(dlink,filePath,dr)
    x= numpy.linspace(-7.,7.,43)
    elsf= fits.getdata(filePath)
    if sparse:
        return (x,sparsify(elsf))
    else:
        return (x,elsf)

def deconvolve(spec,specerr,
               lsf=None,eps=2500.,smooth=None):
    """ deconvolve the LSF
    Args :
       spec - spectrum (nwave)
       specerr - spectrum uncertainty array (nwave)
       lsf= (None) LSF to deconvolve, needs to be specified in non-sparse format
       eps= (2500.) smoothness parameter
       smooth= (None) if set to a resolution, smooth with a FWHM resolution of 'smooth' and return the spectrum on the apStar wavelength grid
    Returns :
       high-resolution deconvolved spectrum or smoothed deconvolved spectrum on apStar wavelength grid is smooth= is set
    HISTORY:
       2015-04-24 - Written - Bovy (IAS)
    """
    # Parse LSF input
    if lsf is None:
        raise ValueError("lsf= keyword with LSF in non-sparse format required for apogee.spec.lsf.deconvolve")
    if isinstance(lsf,sparse.dia_matrix):
        raise ValueError("lsf= keyword with LSF needs to be in non-sparse format")
    lsf[numpy.isnan(lsf)]= 0.
    # How much higher resolution is the LSF than the data?
    hires= int(round(lsf.shape[0]/8575.))
    # Setup output
    out= numpy.zeros(lsf.shape[0])
    # Loop through the detectors and analyze each one separately
    for sindx, eindx in zip([140,3450,6250],[3370,6200,8450]):
        # Get the LSF for this detector
        slsf= sparsify(lsf[hires*sindx:hires*eindx])
        # Parse the spectrum and its error for this detector, normalize
        tspec= numpy.ones((eindx-sindx)*hires)
        tinvspecerr= numpy.zeros((eindx-sindx)*hires)
        norm= numpy.nanmean(spec[sindx:eindx])
        tspec[::hires]= spec[sindx:eindx]/norm
        tinvspecerr[::hires]= norm/specerr[sindx:eindx]
        # Deal with NaNs
        tinvspecerr[numpy.isnan(tspec)]= 0.
        tspec[numpy.isnan(tspec)]= 1.
        # Set up the necessary sparse matrices
        Cinv= sparse.diags([tinvspecerr**2.],[0])
        CinvL= Cinv.dot(slsf)
        LTCinvL= (slsf.T).dot(CinvL)
        # P smoothness matrix
        diags1= -numpy.ones(slsf.shape[1])
        diags1[-1]= 0.
        diags2= numpy.ones(slsf.shape[1]-1)
        P= sparse.diags([diags1,diags2],[0,1])
        A= LTCinvL+eps*(P.T).dot(P)
        # b
        Cinvs= Cinv.dot(tspec)
        b= (slsf.T).dot(Cinvs)
        tmp= scipy.sparse.linalg.bicg(A,b)
        if tmp[1] == 0:
            tmp= tmp[0]
        else:
            raise RuntimeError("Deconvolution did not converge")
        out[sindx*hires:eindx*hires]= tmp*norm
    if not smooth is None:
        wav= aspcap.apStarWave()
        l10wav= numpy.log10(wav)
        dowav= l10wav[1]-l10wav[0]
        sigvm= hires/dowav/smooth/numpy.log(10.)\
            /2./numpy.sqrt(2.*numpy.log(2.))
        out= ndimage.gaussian_filter1d(out,sigvm,mode='constant')[::hires]
    return out


def test(highres=9.,plot=False):
    """ test convolution routine for apogee context
    """
    s=fits.open(os.environ['APOGEE_SPECLIB']+'/synth/turbospec/kurucz/giantisotopes/tgGK_150714/ap00cp00np00vp20.fits')[0].data
    sh=fits.open(os.environ['APOGEE_SPECLIB']+'/synth/turbospec/kurucz/giantisotopes/tgGK_150714/ap00cp00np00vp20.fits')[0].header
    ws=sh['CRVAL1']+numpy.arange(sh['NAXIS1'])*sh['CDELT1']
    ws=spectra.airtovac(ws)
 
    wa=aspcap.apStarWave()
    x=numpy.arange(-15.,15.,1./highres)
    l=eval(x,fiber='combo',waveid=2420038,lsfid=5440020)
    ls=sparsify(l)

    sm=fits.open(os.environ['APOGEE_SPECLIB']+'/synth/turbospec/kurucz/giantisotopes/tgGK_150714_lsfcombo5_l31c/ap00cp00np00vp20.fits')[0].data
    smh=fits.open(os.environ['APOGEE_SPECLIB']+'/synth/turbospec/kurucz/giantisotopes/tgGK_150714_lsfcombo5_l31c/ap00cp00np00vp20.fits')[0].header
    smw=numpy.exp(smh['CRVAL1']+numpy.arange(smh['NAXIS1'])*smh['CDELT1'])

    for iz in range(sh['NAXIS4']) :
      mh=sh['CRVAL4']+iz*sh['CDELT4']
      vmacro = 10.**(0.470794-0.254*mh)
      vmacro = vmacro if vmacro<15 else 15.
      for ig in range(sh['NAXIS3']) :
        for it in range(sh['NAXIS2'])  :
          print(iz,ig,it,mh,vmacro)
          z=convolve(ws,s[iz,ig,it,:],lsf=ls,xlsf=x,vmacro=vmacro)
 
          if plot :
              plt.clf()
              plt.plot(wa,z[0,:],color='r')
              plt.plot(smw,sm[iz,ig,it,:],color='g')
              plt.draw()
              pdb.set_trace()

def rotate(deltav,vsini,epsilon=0.6) :
    """ rotation kernel from IDL Users library routine
    
    Args :
        deltav (float) : velocity spacing of pixels
        vsini (float ) : v sin i for profile
        epsilon (float ) : parameter for rotation kernel (default=0.6)

    Returns :
        normalized rotation profile
    """
    e1 = 2.0*(1.0 - epsilon)
    e2 = numpy.pi*epsilon/2.0
    e3 = numpy.pi*(1.0 - epsilon/3.0)

    npts = numpy.ceil(2*vsini/deltav)
    if npts%2 == 0 : npts = npts +1
    nwid = int(npts/2.)
    x = (numpy.arange(0,npts)- nwid)
    x = x*deltav/vsini
    x1 = abs(1.0 - x**2)
    kernel=(e1*numpy.sqrt(x1) + e2*x1)/e3
    return kernel/kernel.sum()


