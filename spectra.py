import matplotlib.pyplot as plt
import pdb
import scipy.signal
import numpy as np
from astropy.modeling import models, fitting
from astropy.io import ascii

def mash(hd,sp=None) :
    """
    Mash image into spectra using requested window
    """
    if sp is None :
        sp=[0,hd.data.shape[0]]
    return hd.data[sp[0]:sp[1]].sum(axis=0)

def wavecal(hd,file=None,wref=None,disp=None,wid=[3],rad=5,snr=3,degree=2):
    """
    Get wavelength solution for single 1D spectrum
    """
    sz=hd.data.shape
    spec=hd.data[sz[0]/2,:]
    fig=plt.figure()
    plt.subplot(211)
    plt.plot(spec)
    if wref is None :
        w0=hd.header['DISPWC']
        pix0=sz[1]/2 
    else :
        w0=wref[0]
        pix0=wref[1]
    if disp is None:
        disp=hd.header['DISPDW']
    f=open(file,'r')
    lines=[]
    for line in f :
        w=float(line.split()[0])
        name=line[10:].strip()
        pix=(w-w0)/disp+pix0
        print(pix, w, name)
        if pix > 0 and pix < sz[1] :
            plt.text(pix,0.,'{:7.1f}'.format(w),rotation='vertical',va='top',ha='center')
            lines.append(w)
    lines=np.array(lines)

    peaks=scipy.signal.find_peaks_cwt(spec,np.array(wid),min_snr=snr)
    cents=[]
    for peak in peaks :
        if peak > rad and peak < sz[1]-rad :
            cents.append((spec[peak-rad:peak+rad]*np.arange(peak-rad,peak+rad)).sum()/spec[peak-rad:peak+rad].sum())
    cents=np.array(cents)
    print('cents:', cents)
    waves=[]
    weight=[]
    for cent in cents :
        w=(cent-pix0)*disp+w0
        plt.plot([cent,cent],[0,10000],'k')
        print(cent, w, lines[np.abs(w-lines).argmin()])
        waves.append(lines[np.abs(w-lines).argmin()])
        weight.append(1.)
    waves=np.array(waves)
    weight=np.array(weight)
    
    done = False
    fit=fitting.LinearLSQFitter()
    mod=models.Polynomial1D(degree=degree)
    while not done :
        gd=np.where(weight>0.)[0]
        bd=np.where(weight<=0.)[0]
        p=fit(mod,cents[gd],waves[gd],weights=weight[gd])
        plt.subplot(212)
        plt.cla()
        plt.plot(cents[gd],p(cents[gd])-waves[gd],'go')
        plt.plot(cents[bd],p(cents[bd])-waves[bd],'ro')
        diff=p(cents[gd])-waves[gd]
        plt.ylim(diff.min()-1,diff.max()+1)
        for i in range(len(cents)) :
            print(cents[i],p(cents[i])-waves[i],'{:2d}'.format(i))
            plt.subplot(212)
            plt.text(cents[i],p(cents[i])-waves[i],'{:2d}'.format(i),va='top',ha='center')
            plt.subplot(211)
            if weight[i] > 0 :
              plt.plot([cents[i],cents[i]],[0,10000],'g')
            else :
              plt.plot([cents[i],cents[i]],[0,10000],'r')
        plt.draw()
        for i in range(len(cents)) :
            print(i, cents[i], p(cents[i]), waves[i], waves[i]-p(cents[i]))
        i = raw_input('enter ID of line to remove (-n for all lines<n, +n for all lines>n, return to continue): ')
        if i is '' :
            done = True
        elif '+' in i :
            weight[int(i):] = 0.
        elif '-' in i :
            weight[0:abs(int(i))] = 0.
        elif int(i) >= 0 :
            weight[int(i)] = 0.
        else :
            print('invalid input')

    return p

def fluxcal(obs,wobs,file=None) :
    """
    flux calibration
    """

    fluxdata=ascii.read(file)
    stan=np.interp(wobs,fluxdata['col1'],fluxdata['col2'])
    return stan/obs
    
