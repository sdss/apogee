import pdb
import numpy as np

def sincint(x, nres, speclist) :
    """ Use sinc interpolation to get resampled values

        x : desired positions
        nres : number of pixels per resolution element (2=Nyquist)
        speclist : list of [quantity, variance] pairs (variance can be None)
    """

    dampfac = 3.25*nres/2.
    ksize = int(21*nres/2.)
    if ksize%2 == 0 : ksize +=1
    nhalf = ksize//2 

    #number of output and input pixels
    nx = len(x)
    nf = len(speclist[0][0])

    # integer and fractional pixel location of each output pixel
    ix = x.astype(int)
    fx = x-ix

    # outputs
    outlist=[]
    for spec in speclist :
        if spec[1] is None :
            outlist.append([np.full_like(x,0),None])
        else :
            outlist.append([np.full_like(x,0),np.full_like(x,0)])

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

        for spec,out in zip(speclist,outlist) :
            vals = spec[0][lobe[gd]]
            out[0][i] = (sinc[gd]*vals).sum()
            if spec[1] is not None : 
                var = spec[1][lobe[gd]]
                out[1][i] = (sinc[gd]**2*var).sum()

    for out in outlist :
       if out[1] is not None : out[1] = np.sqrt(out[1])
    
    return outlist
