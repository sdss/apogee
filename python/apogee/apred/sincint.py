import pdb
import numpy as np

def sincint(x, nres, f, ferr) :
    """ Use sinc interpolation to get resampled values

        x : desired positions
        nres : number of pixels per resolution element (2=Nyquist)
        f : input function
        ferr : uncertainties on input function
    """

    dampfac = 3.25*nres/2.
    ksize = int(21*nres/2.)
    if ksize%2 == 0 : ksize +=1
    nhalf = ksize//2 

    #number of output and input pixels
    nx = len(x)
    nf = len(f)

    # integer and fractional pixel location of each output pixel
    ix = x.astype(int)
    fx = x-ix

    # output array
    out = np.full_like(x,0)

    for i in range(len(x)) :
        xkernel = np.arange(ksize)-nhalf - fx[i]
        # in units of Nyquist
        xkernel /= (nres/2.)
        u1 = xkernel/dampfac
        u2 = np.pi*xkernel
        sinc = np.exp(-(u1**2)) * np.sin(u2) / u2
        sinc /= (nres/2.)

        lobe = np.arange(ksize) - nhalf + ix[i]
        vals = np.zeros(ksize)
        vars = np.zeros(ksize)
        gd = np.where( (lobe>=0) & (lobe<nf) )[0]
        vals = f[lobe[gd]]
        out[i] = (sinc[gd]*vals).sum()

    return out
