# some snippets of code to show how one might
# get initial trace and wavecal files for APO instruments

import copy
import numpy as np
import os
import pickle
import pdb
import yaml
import matplotlib.pyplot as plt
from pyvista import imred
from pyvista import image
from pyvista import spectra
from pyvista import tv
from astropy.io import fits
from astropy import units
from astropy.nddata import CCDData, StdDevUncertainty
import scipy.signal

ROOT = os.path.dirname(os.path.abspath(__file__)) + '/../../'

import readmultispec
import os
import pickle

def dis() :
    pdb.set_trace()
    #dred=imred.Reducer(inst='DIS',dir='UT191019/DIS/')
    dred=imred.Reducer(inst='DIS',dir='UT191027/DIS/')
    dcomb=imred.Combiner(reducer=dred)
    arc=dcomb.sum([1,2,3])

    # define a flat trace model
    def model(x) : return(np.zeros(len(x))+500)
    traces=spectra.Trace(rad=10,model=[model],sc0=1024)

    spec=traces.extract(arc[0])
    spec.data-=scipy.signal.medfilt(spec.data[0,:],kernel_size=101)

    wcal=spectra.WaveCal(type='chebyshev')
    #wcal.set_spectrum(spec)
    #f=open(os.environ['PYVISTA_DIR']+'/data/dis/dis_blue_lowres.pkl','rb')
    f=open(os.environ['PYVISTA_DIR']+'/data/DIS/DIS_hires_6400_waves.pkl','rb')
    wcal0=pickle.load(f)
    wcal0.orders=[1]
    wav=wcal0[0][0].wave(image=spec.shape)

    wcal.identify(spec,wav=wav,file=os.environ['PYVISTA_DIR']+'/data/lamps/henear.dat',rad=3)
    wcal.fit()
    w=wcal.wave(image=np.array(spec.data.shape))

    spec2d=traces.extract2d(arc[0],rows=range(200,900))
    pdb.set_trace()

    star=dred.reduce(26)
    return spec,w,wcal

def arces() :
    ered=imred.Reducer(inst='ARCES',dir='UT191020/ARCES/')
    ecomb=imred.Combiner(reducer=ered)
    flat=ecomb.sum([11,15])
    thar=ecomb.sum([19,20])

    t=tv.TV()

    apertures=np.loadtxt('newap')[:,1]
    #traces=trace(flat,apertures/4)
    traces=Trace()
    traces.trace(flat,apertures/4,sc0=1024,thresh=1000,plot=t)
    ec=traces.extract(thar,scat=True)

    wcal=spectra.WaveCal(type='chebyshev2D',orders=54+np.arange(107),spectrum=ec)
    wav=readmultispec.readmultispec('w131102.0004.ec.fits')['wavelen']
    new=ec.data*0.
    new[:,204:1855]=wav
    wcal.identify(ec,wav=new,file=os.environ['PYVISTA_DIR']+'/data/lamps/thar_arces',xmin=201,xmax=1849,thresh=200,rad=3,plot=t)
    wcal.fit()

    wcal.identify(ec,file=os.environ['PYVISTA_DIR']+'/data/lamps/thar_arces',xmin=150,xmax=1900,thresh=200,rad=3,plot=t)
    wcal.fit()

    w=wcal.wave(image=np.array(ec.data.shape))

    hr7950=ered.reduce(1,superflat=flat)
    hr7950ec=traces.extract(hr7950,scat=True)

    return flat,thar,traces,ec,wcal,w

def tspec() :

    tspec=imred.Reducer(inst='TSPEC',dir='UT191026/TSPEC',nfowler=8)   
    a=tspec.reduce(21)
    acr=tspec.reduce(21,crbox=[11,1])
    b=tspec.reduce(22)

    t=tv.TV()
    fig=plt.figure()

    rows=[[135,235],[295,395],[435,535],[560,660],[735,830]]
    apers=[155,316,454,591,761]    
    traces=[]
    waves=[]

    # wavelength arrays from TSpexTool reduction, note that flips wavelengths, and has K band first, not last
    wav=fits.open('tspec_wave.fits')[0].data
    #w0=[23147.049,17388.174,13927.875,11621.16,9972.335]
    #w0.reverse()
    #disp=[-2.8741264,-2.156124, -1.7261958,-1.4418244,-1.2400593]
    #disp.reverse()
    order=7
    for aper,row in zip(apers,rows) :
        trace=spectra.Trace(order=3,rows=row,lags=range(-75,75))
        trace.trace(a.subtract(b),aper,sc0=350,thresh=100,plot=t)
        traces.append(trace)
        trace.pix0 +=30
        out=trace.extract(acr,rad=20,plot=t)
        out.data=out.data - scipy.signal.medfilt(out.data,kernel_size=[1,201])
        wcal=spectra.WaveCal(type='chebyshev',orders=[order],degree=3)
        w=np.atleast_2d(wav[order-3,0,:][::-1])*1.e4
        bd=np.where(~np.isfinite(w))
        w[bd[0],bd[1]]=9000.
        #wcal.identify(out,file=os.environ['PYVISTA_DIR']+'/data/sky/OHll.dat',rad=10,wref=[wr,2048-1500],disp=d,plot=fig)
        wcal.identify(out,file=os.environ['PYVISTA_DIR']+'/data/sky/OHll.dat',rad=3,wav=w,plot=fig)
        wcal.fit()
        delattr(wcal,'ax')
        waves.append(wcal)
        order-=1

    #spec2d=traces.extract2d(a.subtract(b),rows=range(200,900))
    pickle.dump([traces],open('TSPEC_traces.pkl','wb'))
    pickle.dump([waves],open('TSPEC_waves.pkl','wb'))

    return traces,waves

def test():
    a=np.loadtxt('ecnewarc.ec')
    x=a[:,2]
    y=a[:,1]
    w=a[:,4]
    fitter=fitting.LinearLSQFitter()
    mod=models.Chebyshev2D(x_degree=3,y_degree=3,x_domain=[1,2000],y_domain=[55,160])
    fit=fitter(mod,x,y,w*y)
    plt.figure()
    plt.plot(w,w-fit(x,y)/y,'ro')

def testfit() :
    a=np.loadtxt('ecnewarc.ec')
    a[:,2]
    x=a[:,2]
    y=a[:,1]
    w=a[:,4]
    fitter=fitting.LinearLSQFitter()
    mod=models.Chebyshev2D(x_degree=3,y_degree=3,x_domain=[1,2000],y_domain=[55,160])
    fit=fitter(mod,x,y,w*y)
    plt.figure()
    plt.plot(w,w-fit(x,y)/y)
    plt.clf()
    plt.plot(w,w-fit(x,y)/y,'ro')


def echfit() :
    mod=models.custom_model(echelle_model)

def ripple(w,ord,amp=1,alpha=1,wc=5500) :
#def echelle_model(ind,amp,alpha,wc) :
#    w=ind[0,:]
#    ord=ind[1,:]
    x = ord * (1 - wc/w)
    out=amp*(np.sin(np.pi*alpha*x)/(np.pi*alpha*x))**2
    bd=np.where(np.pi*alpha*x == 0)[0]
    out[bd] = amp[bd]
    return out

def ripple2d(w,ord,amp0=1,amp1=0,amp2=0,alpha0=1,alpha1=0,alpha2=0,wc0=5500,wc1=0,wc2=0) :
#def echelle_model(ind,amp,alpha,wc) :
#    w=ind[0,:]
#    ord=ind[1,:]
    wc = wc0 + wc1*(ord-109) + wc2*(ord-109)**2
    amp = amp0 + amp1*(ord-109) + amp2*(ord-109)**2
    alpha = alpha0 + alpha1*(ord-109) + alpha2*(ord-109)**2
    #print(ord,amp,alpha,wc)
    x = ord * (1 - wc/w)
    out=amp*(np.sin(np.pi*alpha*x)/(np.pi*alpha*x))**2
    bd=np.where(np.pi*alpha*x == 0)[0]
    out[bd] = amp[bd]
    return out

def echelle_model_deriv(w,ord,amp=1,alpha=1,wc=5500) :
    x = ord * (1 - wc/w)
    f=amp*(np.sin(np.pi*alpha*x)/(np.pi*alpha*x))**2
    deriv=[]
    deriv.append(f/amp)
    deriv.append(2*amp*np.sin(np.pi*alpha*x)/(np.pi*alpha*x)**2*np.cos(np.pi*alpha*x)*np.pi*x
                 - 2*f/alpha )
    deriv.append(2*amp*np.sin(np.pi*alpha*x)/(np.pi*alpha*x)**2*np.cos(np.pi*alpha*x)*np.pi*alpha*-1*ord/w
                 + 2*f/x*ord/w)

    return np.array(deriv)

def mktrace(ymlfile) :

    # trace from scratch
    try :
        traces=group['traces']
        print('create trace')
        if group['traces']['frame'] == 'sflat' : 
            a=sflat
        else :
            a= red.reduce(group['traces']['frame'],superbias=sbias,superdark=sdark,superflat=sflat,scat=red.scat)
        try : 
            b= red.reduce(group['traces']['dark'],superbias=sbias,superdark=sdark,superflat=sflat,scat=red.scat)
            a=a.subtract(b)
        except: pass
        if group['inst'] == 'TSPEC' :
            apers=[155,316,454,591,761]
            apers.reverse()
            traces=spectra.Trace(inst=group['inst'])
            traces.trace(a,apers,sc0=350,plot=display)
            traces=[traces]
        elif group['inst'] == 'ARCES' :
            apers=np.loadtxt('newap')[:,1]/4
            traces=spectra.Trace(inst=group['inst'])
            traces.trace(a,apers,sc0=1024,plot=display)
            traces=[[traces]]
        elif group['inst'] == 'DIS' :
            apers=[662,562]
            traces=[]
            for i,(im,ap) in enumerate(zip(a,apers)) :
                trace=spectra.Trace(inst=group['inst'],channel=i)
                trace.trace(im,ap,sc0=1024,plot=display)
                traces.append(trace)
            traces=[[traces[0]],[traces[1]]]
        pickle.dump(traces,open(group['inst']+'_trace.pkl','wb'))
    except:
        print('no trace frames given')

    # existing trace template
    trace=pickle.load(open(group['inst']+'_trace.pkl','rb'))

