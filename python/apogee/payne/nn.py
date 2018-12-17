import numpy as np

import matplotlib.pyplot as plt
import multiprocessing as mp
from scipy.ndimage.filters import gaussian_filter
from scipy.optimize import curve_fit
from keras import models
from keras import layers
from keras import optimizers
from keras import regularizers
from astropy.io import fits
from astropy.table import Table, TableColumns, Column
import pickle
import copy
import sys
import pdb
import time
from apogee.utils import apload
from apogee.aspcap import aspcap
from apogee.aspcap import norm
from tools import plots

nepochs=10000
nodes=20
reg=0.0005
batch_size=1000
verbose=0
reg=0.
nepochs=25000

nodes=20
nepochs=50000

def train(file,plot=False,pixels=[1000,9000,1000],suffix='',fitfrac=1.0, order=0, threads=32,payne=False,
          teff=[0,10000],logg=[-1,6],mh=[-3,1],am=[-1,1],cm=[-2,2],nm=[-2,2],raw=False,rot=False,nolog=True,elem=False) :
    """ Train a neural net model on an input training set
    """
    global nfit, verbose, nepochs

    # Get input spectra and parameters
    pixels=np.arange(pixels[0],pixels[1],pixels[2])
    npix=len(pixels)
    print('npix: ', npix)
    if '.npz' in file :
        a = np.load(file)
        pars = a['labels'].T
        spec = a['spectra'] 
    else :
        pars=fits.open(file+'.fits')[0].data
        if raw:
            spec=fits.open(file+'.fits')[1].data
        else :
            spec=fits.open(file+'.fits')[2].data
    head = {}
    head['nin'] = pars.shape[0]
    head['npix'] = npix

    # limit parameter range
    if payne :
        gd = np.where( (pars[:,1] > (pars[:,0]-6000.)*(1.0/1000.)+1.99) &
                       (pars[:,1] < (pars[:,0]-5000.)*(-0.45/3000.)+5.01) )[0]
    else :
        gd=np.where((pars[:,0]>=teff[0]) & (pars[:,0]<=teff[1]) &
                (pars[:,1]>=logg[0]) & (pars[:,1]<=logg[1]) &
                (pars[:,2]>=mh[0]) & (pars[:,2]<=mh[1]) &
                (pars[:,3]>=am[0]) & (pars[:,3]<=am[1]) &
                (pars[:,4]>=cm[0]) & (pars[:,4]<=cm[1])  &
                (pars[:,5]>=nm[0]) & (pars[:,5]<=nm[1]) 
               )[0]
    spec=spec[gd,:]
    pars=pars[gd,:]
    head['ntot'] = pars.shape[0]
    head['teff'] = teff
    head['logg'] = logg
    head['mh'] = mh
    head['am'] = am
    head['cm'] = cm

    # limit parameters?
    if not elem :
        if rot :
            pars=pars[:,0:8]
        else :
            pars=pars[:,0:7]

    if nolog : pars[:,2] = 10.**pars[:,2]

    #normalize spectra
    print('normalizing...')
    x=np.arange(0,spec.shape[1])
    specerr = np.full_like(spec[0,:],1.)
    for i in range(spec.shape[0]) :
        gd = np.where(np.isfinite(spec[i,:]))[0]
        if len(gd) == 0 :
          print(i,pars[i,:])
        if order >= 0 : 
            cont = norm.cont(spec[i,:],specerr,poly=True,order=order,chips=True)
            spec[i,:] /= cont

    if plot :
        fig,ax=plots.multi(2,2)
        plots.plotc(ax[0,0],pars[:,0],spec[:,1000],pars[:,1])
        plots.plotc(ax[1,0],pars[:,0],spec[:,1000],pars[:,2])
    print(spec.shape,pars.shape)
    shape=pars.shape

    # shuffle them and get fit and validation set
    print('shuffling...')
    p=np.random.permutation(shape[0])
    spec=spec[p,:]
    pars=pars[p,:]
    nfit=int(len(p)*fitfrac)
    shape=pars.shape

    # scale parameters to zero mean and unit standard deviation, and save scaling parameters
    pmeans=[]
    pstds=[]
    normpars = copy.copy(pars)
    for i in range(shape[1]) :
      mn=pars[:,i].mean()
      std=pars[:,i].std()
      normpars[:,i] -= mn
      if std > 0. : normpars[:,i] /= std
      pmeans.append(mn)
      pstds.append(std)

    # replot to check
    if plot :
        plots.plotc(ax[0,1],normpars[:,0],spec[:,1000],normpars[:,1])
        plots.plotc(ax[1,1],normpars[:,0],spec[:,1000],normpars[:,2])
        plt.show()

    # loop over the requested pixels and normalize data to
    #   zero mean and unit standard deviation: save parameters
    weights=[]
    biases=[]
    means=[]
    stds=[]
    data=[]
    print('preparing to fit...')
    for ipix in pixels :
        pix=spec[:,ipix]
        mn=pix.mean()
        std=pix.std()
        if np.isfinite(mn) :
          pix-=mn
          pix /= std
          data.append((normpars,pix,ipix))
        means.append(mn)
        stds.append(std)

    # get the model in parallel for different pixels
    print('fitting: ',len(data))
    pool = mp.Pool(threads)
    output = pool.map_async(fit, data).get()
    pool.close()
    pool.join()
    print('done pool')

    if plot: 
        fig,ax=plots.multi(npix,7,wspace=0.001,hspace=0.001,figsize=(15,10),xtickrot=60)
        fig2,ax2=plots.multi(npix,2,wspace=0.001,hspace=0.5,figsize=(15,4),xtickrot=90)
    ifit=0
    for i,ipix in enumerate(pixels) :
      if np.isfinite(means[i]) :
        w,b,mod,loss,vloss=output[ifit]
        ifit+=1
        if plot :
              mod=mod*stds[i]+means[i]
              m=[]
              for ip in range(pars.shape[0]) : m.append(model(normpars[ip,:],means[i],stds[i],w,b)[0])
              pix=spec[:,ipix]*stds[i]+means[i]
              plots.plotc(ax[0,i],pars[0:nfit,0],pix[0:nfit]-mod[0:nfit,0],pars[0:nfit,1],xr=[3000,8000],zr=[0,5])
              plots.plotc(ax[1,i],pars[0:nfit,0],pix[0:nfit]-mod[0:nfit,0],pars[0:nfit,2],xr=[3000,8000],zr=[-2.0,0.75])
              plots.plotc(ax[2,i],pars[0:nfit,0],pix[0:nfit]-mod[0:nfit,0],pars[0:nfit,3],xr=[3000,8000],zr=[-1,1.])
              plots.plotc(ax[3,i],pars[0:nfit,0],pix[0:nfit]-mod[0:nfit,0],pars[0:nfit,4],xr=[3000,8000],zr=[-1,1.])
              plots.plotc(ax[4,i],pars[0:nfit,0],pix[0:nfit]-mod[0:nfit,0],pars[0:nfit,5],xr=[3000,8000],zr=[-1,1.])
              plots.plotc(ax[5,i],pars[0:nfit,0],pix[0:nfit]-mod[0:nfit,0],pars[0:nfit,6],xr=[3000,8000],zr=[0.5,3.])
              plots.plotc(ax[6,i],pars[0:nfit,0],pix[0:nfit]-mod[0:nfit,0],pars[0:nfit,7],xr=[3000,8000],zr=[0,50.])
              n=len(loss)
              plots.plotl(ax2[0,i],range(n),np.log10(loss),xr=[0,nepochs],yr=[-4,0],color='b')
              if fitfrac < 1.0 :
                  plots.plotl(ax2[0,i],range(n),np.log10(vloss),xr=[0,nepochs],yr=[-4,0],color='r')
              try : 
                  ax2[1,i].hist(np.abs(pix-mod[:,0]),bins=np.logspace(-7,3,50),histtype='step',normed=True,cumulative=True,color='k')
                  ax2[1,i].hist(np.abs(pix[0:nfit]-mod[0:nfit,0]),bins=np.logspace(-7,3,50),histtype='step',normed=True,cumulative=True,color='b')
              except: pass
              if fitfrac < 1.0 :
                  ax2[1,i].hist(np.abs(pix[nfit:]-mod[nfit:,0]),bins=np.logspace(-7,3,50),histtype='step',normed=True,cumulative=True,color='r')
              ax2[1,i].set_xlim(0,0.01)
              ax2[1,i].set_ylim(0,1.1)
              plt.draw()
              plt.show()
  
        weights.append(w)
        biases.append(b)
    if plot: fig.savefig(file+suffix+'_pixels.jpg')
    pdb.set_trace()

    # Save the model as a dictionary into a pickle file
    head['nodes'] = nodes
    head['reg'] = reg
    head['batch_size'] = batch_size
    head['nepochs'] = nepochs
    head['pmeans'] = pmeans
    head['pstds'] = pstds
    head['means'] = means
    head['stds'] = stds
    head['weights'] = weights
    head['biases'] = biases
    with open(file+suffix+'.pkl', 'w') as f:  # Python 3: open(..., 'wb')
        pickle.dump(head, f)
        #pickle.dump([head, pmeans, pstds, means, stds, weights, biases], f)

    return head

def fit(data) :
    """ Routine to do a single NN model fit given input data=(pars,pix)
    """
    pars=data[0]
    pix=data[1]

    showtime('fitting pixel: '+str(data[2]))

    net=models.Sequential()
    #net.add(layers.Dense(32, activation='sigmoid', input_shape=(pars.shape[1],),
    #        kernel_regularizer=regularizers.l2(reg)))
    #net.add(layers.Dense(64, activation='sigmoid',
    #        kernel_regularizer=regularizers.l2(reg)))
    #net.add(layers.Dense(128, activation='sigmoid',
    #        kernel_regularizer=regularizers.l2(reg)))
    #net.add(layers.Dense(256, activation='sigmoid',
    #        kernel_regularizer=regularizers.l2(reg)))
    net.add(layers.Dense(nodes, activation='sigmoid', input_shape=(pars.shape[1],),
            kernel_regularizer=regularizers.l2(reg)))
    net.add(layers.Dense(nodes, activation='sigmoid',kernel_regularizer=regularizers.l2(reg)))
    net.add(layers.Dense(1, activation='linear'))
    ##opt=optimizers.RMSprop(lr=0.01)
    opt=optimizers.Adam(lr=0.001)
    net.compile(optimizer=opt,loss='mse')
    if verbose > 0 : net.summary()

    history=net.fit(pars[0:nfit],pix[0:nfit],epochs=nepochs,batch_size=batch_size,verbose=verbose,validation_data=(pars[nfit:],pix[nfit:]))

    w=(net.get_weights()[0],net.get_weights()[2])
    b=(net.get_weights()[1],net.get_weights()[3])
    mod=net.predict(pars)

    showtime('done fitting pixel: '+str(data[2]))

    try :
        return w,b,mod,history.history['loss'],history.history['val_loss']
    except :
        return w,b,mod,history.history['loss'],0.

def merge(file,n=8) :
    """ Merge pieces of a model (e.g., run on different nodes for pixel subsets) into a single model
    """
    pm,ps,m,s,w,b=[],[],[],[],[],[]
    for i in range(n) :
      with open(file+'_{:d}.pkl'.format(i+1)) as f: 
        head = pickle.load(f)
        # pmeans, pstds same for all pixels
        m.extend(head['means'])
        s.extend(head['stds'])
        w.extend(head['weights'])
        b.extend(head['biases'])

    head['means'] = m
    head['stds'] = s
    head['weights'] = w
    head['biases'] = b
    with open(file+'.pkl', 'w') as f:  
        pickle.dump(head, f)


def sigmoid(z):
    """ sigmoid function
    """
    return 1.0/(1.0+np.exp(-z))
 
def model(pars, mn, std, weights, biases) :
    """ function to return single pixel model given normalized input parameters and pixel normalization
    """
    return mn + std * (np.dot( sigmoid((np.dot(weights[0].T,pars)+biases[0])).T, weights[1] ) +biases[1])

def spectrum(x,*pars) :
    """ Return full spectrum given input list of pixels, parameters
    """
    spec=np.full_like(x.astype(float),np.nan)
    pnorm= (pars-pmeans)/pstds
    for j,i in enumerate(x) :
        if  np.isfinite(means[i]) :
            spec[j]= model(pnorm, means[i], stds[i], weights[ifit[i]], biases[ifit[i]]) 

    return spec

def get(file) :
    """ Load a model pickle file into global variables
    """
    global head, pmeans, pstds, means, stds, weights, biases, ifit

    # Getting back the objects:
    with open(file+'.pkl') as f: 
        head = pickle.load(f)
        #pmeans, pstds, means, stds, weights, biases = pickle.load(f)
    pmeans = np.array(head['pmeans'])
    pstds = np.array(head['pstds'])
    means = head['means']
    stds = head['stds']
    weights = head['weights']
    biases = head['biases']
    #pmeans = np.array(pmeans)
    #pstds = np.array(pstds)

    # get correspondence of pixel number with weight/bias index (since NaNs are not fit)
    ifit = np.zeros(len(means)).astype(int)
    j=0
    for i in range(len(means)) :
        if np.isfinite(means[i]) : 
            ifit[i] = j
            j += 1

def test(pmn, pstd, mn, std, weights, biases,n=100, t0=[3750.,4500.], g0=2., mh0=0.) :
    """ Plots cross-sections of model for fit pixels
    """
    fig,ax=plots.multi(2,6,figsize=(8,12))

    xt=['Teff','logg','[M/H]','[alpha/M]','[C/M]','[N/M]']
    for i,ipar in enumerate([0,1,2,3,4,5]) : 
      for ipix in range(len(weights)) :
       for it0 in range(2) :
        pars=np.tile([t0[it0], g0, mh0, 0.0, 0., 0., 2.],(n,1))
        if ipar == 0 : pars[:,ipar]=np.linspace(3000.,8000.,n)
        elif ipar == 1 : pars[:,ipar]=np.linspace(-0.5,5.5,n)
        elif ipar == 2 : pars[:,ipar]=np.linspace(-2.5,1.,n)
        elif ipar == 3 : pars[:,ipar]=np.linspace(-0.5,1.0,n)
        elif ipar == 4 : pars[:,ipar]=np.linspace(-1.,1.,n)
        elif ipar == 5 : pars[:,ipar]=np.linspace(-0.5,2.,n)
        m=[]
        for ip in range(pars.shape[0]) : m.append(model((pars[ip,:]-pmn)/pstd,mn[ipix],std[ipix],weights[ipix],biases[ipix]))
        plots.plotl(ax[i,it0],pars[:,ipar],m,xt=xt[i])
        #m=[]
        #for ip in range(pars.shape[0]) : m.append(nets[ipix].predict((pars[ip,:].reshape(1,7)-pmn)/pstd)[0,0]*std[ipix]+mn[ipix])
        #plots.plotl(ax[i,it0],pars[:,ipar],m)
        if i == 0 : ax[i,it0].set_title('{:8.0f}{:7.2f}{:7.2f}'.format(t0[it0],g0,mh0))
    fig.tight_layout()

def fitinput(file,threads=8,nfit=8,dofit=True,order=4) :
    """ Solves for parameters using input spectra and NN model
    """
    p=fits.open(file+'.fits')[0].data
    s=fits.open(file+'.fits')[2].data
    p=p[:,0:8]
    if nfit == 0 : nfit = p.shape[0]
    get(file)

    specerr=np.full_like(s[0,:],0.005)
    if dofit :
        specs=[]
        for i in range(nfit) :
            cont = norm.cont(s[i,:],specerr,poly=True,order=order,chips=True)
            specs.append((s[i,:]/cont, specerr))

        pool = mp.Pool(threads)
        output = pool.map_async(solve, specs).get()
        pool.close()
        pool.join()

        # plot output minus input parameters 
        output=np.array(output)
        fig,ax=plots.multi(2,4,hspace=0.001)
        plots.plotc(ax[0,0],p[0:nfit,0],output[:,0]-p[0:nfit,0],p[0:nfit,2],yr=[-250,250],yt='Teff')
        plots.plotc(ax[1,0],p[0:nfit,0],output[:,1]-p[0:nfit,1],p[0:nfit,2],yr=[-1,1],yt='logg')
        plots.plotc(ax[2,0],p[0:nfit,0],output[:,2]-p[0:nfit,2],p[0:nfit,2],yr=[-0.5,0.5],yt='[M/H]')
        plots.plotc(ax[3,0],p[0:nfit,0],output[:,3]-p[0:nfit,3],p[0:nfit,2],yr=[-0.5,0.5],yt='[alpha/M]')
        plots.plotc(ax[0,1],p[0:nfit,0],output[:,4]-p[0:nfit,4],p[0:nfit,2],yr=[-0.5,0.5],yt='[C/M]')
        plots.plotc(ax[1,1],p[0:nfit,0],output[:,5]-p[0:nfit,5],p[0:nfit,2],yr=[-0.5,0.5],yt='[N/M]')
        plots.plotc(ax[2,1],p[0:nfit,0],output[:,6]-p[0:nfit,6],p[0:nfit,2],yr=[-0.5,0.5],yt='vmicro')
        fig.savefig(file+'_out.png')
        # write the spectra out
        hdu=fits.HDUList()
        hdu.append(fits.ImageHDU(output))
        hdu.writeto(file+'_out.fits',overwrite=True)
        hdu.close()

    # save model and fit spectra
    pix = np.arange(0,8575,1)
    model=[]
    for i in range(nfit) :
        snorm = s[i,:] / np.nanmean(s[i,:])
        spec=spectrum(pix, *p[i,:])    
        model.append(spec)

    hdu=fits.HDUList()
    hdu.append(fits.ImageHDU(np.array(model)))
    hdu.writeto(file+'_model.fits',overwrite=True)   
    hdu.close()


def dclip(d,lim=[-0.5,0.5]) :
    d[np.where(d < lim[0])[0]]=lim[0]
    d[np.where(d > lim[1])[0]]=lim[1]
    return d

def comp(file,order=4,z=2)  :
    """ Plot results of testinput vs true parameters
    """
    p=fits.open(file+'.fits')[0].data
    s=fits.open(file+'.fits')[2].data
    out=fits.open(file+'_out.fits')[0].data
    fit=fits.open(file+'_model.fits')[0].data
    specerr=np.full_like(s[0,:],0.005)

    fig,ax=plots.multi(2,7,hspace=0.001,wspace=0.5)
    plots.plotc(ax[0,0],p[:,0],out[:,0]-p[:,0],p[:,z],xt='Teff',yt=r'$\Delta$Teff') #,yr=[-200,200])
    plots.plotc(ax[1,0],p[:,0],out[:,1]-p[:,1],p[:,z],xt='Teff',yt=r'$\Delta$logg') #,yr=[-0.5,0.5])
    plots.plotc(ax[2,0],p[:,0],out[:,2]-p[:,2],p[:,z],xt='Teff',yt=r'$\Delta$[M/H]') #,yr=[-0.5,0.5])
    plots.plotc(ax[3,0],p[:,0],out[:,3]-p[:,3],p[:,z],xt='Teff',yt=r'$\Delta$[a/M]') #,yr=[-0.5,0.5])
    plots.plotc(ax[4,0],p[:,0],out[:,4]-p[:,4],p[:,z],xt='Teff',yt=r'$\Delta$[C/M]') #,yr=[-0.5,0.5])
    plots.plotc(ax[5,0],p[:,0],out[:,5]-p[:,5],p[:,z],xt='Teff',yt=r'$\Delta$[N/M]') #,yr=[-0.5,0.5])
    plots.plotc(ax[6,0],p[:,0],out[:,6]-p[:,6],p[:,z],xt='Teff',yt=r'$\Delta$vmicro') #,yr=[-0.5,0.5])
    ax[0,1].hist(dclip(out[:,0]-p[:,0],lim=[-200,200]),bins=np.arange(-200,201,10),histtype='step')
    ax[1,1].hist(dclip(out[:,1]-p[:,1]),bins=np.arange(-0.5,0.51,0.01),histtype='step')
    ax[2,1].hist(dclip(out[:,2]-p[:,2]),bins=np.arange(-0.5,0.51,0.01),histtype='step')
    ax[3,1].hist(dclip(out[:,3]-p[:,3]),bins=np.arange(-0.5,0.51,0.01),histtype='step')
    ax[4,1].hist(dclip(out[:,4]-p[:,4]),bins=np.arange(-0.5,0.51,0.01),histtype='step')
    ax[5,1].hist(dclip(out[:,5]-p[:,5]),bins=np.arange(-0.5,0.51,0.01),histtype='step')
    ax[6,1].hist(dclip(out[:,6]-p[:,6]),bins=np.arange(-0.5,0.51,0.01),histtype='step')
    fig.suptitle(file)
    pdb.set_trace()

    for i in range(s.shape[0]) :
        cont = norm.cont(s[i,:],specerr,poly=True,order=order,chips=True)
        print('{:8.1f}{:7.2f}{:7.2f}{:7.2f}{:7.2f}{:7.2f}{:7.2f}{:7.2f}'.format(
               p[i,0],p[i,1],p[i,2],p[i,3],p[i,4],p[i,5],p[i,6],p[i,7]))
        print('{:8.1f}{:7.2f}{:7.2f}{:7.2f}{:7.2f}{:7.2f}{:7.2f}{:7.2f}'.format(
               out[i,0],out[i,1],out[i,2],out[i,3],out[i,4],out[i,5],out[i,6],out[i,7]))

        plt.clf()
        plt.plot(s[i,:]/cont,color='b')
        plt.plot(fit[i,:],color='r')
        plt.plot(fit[i,:]/(s[i,:]/cont)+0.1,color='g')
        #if dofit :
        #    print(output[i])
        #    fit=spectrum(pix, *output[i])
        #    gd=np.where(np.isfinite(snorm))[0]
        #    print(np.sum((spec[gd]-snorm[gd])**2),np.sum((fit[gd]-snorm[gd])**2))
        #    plt.plot(fit,color='b')
        plt.draw()
        pdb.set_trace()

def solve(spec) :
    """ Solve for parameters for a single input spectrum
    """
    s=spec[0]
    serr=spec[1]
    pix = np.arange(0,len(s),1)
    init = np.array([4000.,2.5,0.,0.,0.,0.,1.5,0.])
    bounds = (np.array([3000.,-0.5,-3.,-1.,-1.,-1.,0.5,0.]),
              np.array([8000., 5.5, 1., 1., 1., 1.,4.5,100.]))
    gd = np.where(np.isfinite(s))[0]
    fpars,fcov = curve_fit(spectrum,pix[gd],s[gd],sigma=serr[gd],p0=init,bounds=bounds)
    return fpars

def fitfield(model,field,stars=None,nfit=0,order=4,threads=8,plot=False,write=True) :
    """ Fit observed spectra in an input field, given a model
    """

    # get model and list of stars
    get(model)
    load=apload.ApLoad(dr='dr14')
    apfield=load.apField(field)[1].data
    aspcap_param=load.aspcapField(field)[1].data
    aspcap_spec=load.aspcapField(field)[2].data
    if stars is None :
        stars=apfield['apogee_id']
        if nfit == 0 : stars = stars[0:nfit]
 
    # load up normalized spectra and uncertainties 
    specs=[]
    if plot : fig,ax=plots.multi(1,2,hspace=0.001,figsize=(15,3))
    pix = np.arange(0,8575,1)
    for star in stars :
        apstar=dr14load.apStar(field,star)
        spec = apstar[1].data[0,:].squeeze()
        specerr = apstar[2].data[0,:].squeeze()
        cont = norm.cont(spec,specerr,poly=True,order=order,chips=True)
        specs.append((spec/cont,specerr/cont))
        if plot :
            ax[0].cla()
            ax[1].cla()
            ax[0].plot(spec)
            ax[0].plot(cont)
            ax[1].plot(spec/cont)
            j=np.where(aspcap_param['APOGEE_ID'] == star)[0][0]
            aspec=aspcap.aspcap2apStar(aspcap_spec[j]['spec'])
            ax[1].plot(aspec,color='r')
            plt.draw()
            plt.show()
            pdb.set_trace()

    # do the fits in parallel
    pool = mp.Pool(threads)
    output = pool.map_async(solve, specs).get()
    pool.close()
    pool.join()
    print('done pool')

    # output FITS table
    output=np.array(output)
    out=Table()
    out['APOGEE_ID']=stars
    length=len(out)
    out.add_column(Column(name='FPARAM',data=output))
    spec=[]
    err=[]
    bestfit=[]
    chi2=[]
    for i,star in enumerate(stars) :
        spec.append(specs[i][0])
        err.append(specs[i][1])
        fit=spectrum(pix, *output[i])
        bestfit.append(fit)
        chi2.append(np.nansum((specs[i][0]-fit)**2/specs[i][1]**2))
    out.add_column(Column(name='SPEC',data=np.array(spec)))
    out.add_column(Column(name='ERR',data=np.array(err)))
    out.add_column(Column(name='SPEC_BESTFIT',data=np.array(bestfit)))
    out.add_column(Column(name='CHI2',data=np.array(chi2)))
    if write : out.write('nn-'+field+'.fits',format='fits',overwrite=True)
    return out

def aspcap_comp(model,fields,plot=True,save=None,loggmax=99) :
    """ Compare NN results with ASPCAP 
    """

    get(model)
    load=apload.ApLoad(dr='dr14')

    apars_all=[]
    npars_all=[]
    if plot :    fig,ax=plots.multi(1,3,hspace=0.001,figsize=(15,3),sharex=True)
    for field in fields :
      try :
        out=fits.open('nn-'+field+'.fits')[1].data
        nfit = len(out)
        print(field,nfit)
        apfield=load.apField(field)[1].data
        aspcap_param=load.aspcapField(field)[1].data
        aspcap_spec=load.aspcapField(field)[2].data
        stars=apfield['apogee_id']

        # look at results for each spectrum
        pix = np.arange(0,8575,1)

        for i in range(nfit) :
            print(i,stars[i])
            j=np.where((aspcap_param['APOGEE_ID'] == stars[i]) & (aspcap_param['FPARAM'][:,1] < loggmax))[0]
            if len(j) > 0 :
              j=j[0]
              fp=aspcap_param[j]['FPARAM'].squeeze()
              apars=np.array([fp[0],fp[1],fp[3],fp[6],fp[4],fp[5],10.**fp[2]])
              apars_all.append(apars)
              npars_all.append(out['FPARAM'][i,:])
              if plot and np.abs(out['FPARAM'][i,0]-fp[0]) > 1000:
                pprint(out[i]['FPARAM'])
                ax[0].cla()
                ax[0].plot(out[i]['SPEC'],color='k')
                ax[0].plot(out[i]['ERR'],color='k',ls='dotted')
                ax[0].set_ylim(0.5,1.5)
                ax[1].cla()
                ax[1].plot(out[i]['SPEC'],color='k')
                ax[1].plot(out[i]['SPEC_BESTFIT'],color='b')
                ax[1].set_ylim(0.5,1.5)
                print('nn chi2: ',np.nansum((out[i]['SPEC']-out[i]['SPEC_BESTFIT'])**2/out[i]['ERR']**2))
                ax[2].plot((out[i]['SPEC']-out[i]['SPEC_BESTFIT'])**2/out[i]['ERR']**2,color='b')
                # using ASPCAP parameters and NN model
                # plot ASPCAP normalized spectrum
                aspec=aspcap.aspcap2apStar(aspcap_spec[j]['spec'])
                aerr=aspcap.aspcap2apStar(aspcap_spec[j]['err'])
                print('nn chi2 with aerr: ',np.nansum((out[i]['SPEC']-out[i]['SPEC_BESTFIT'])**2/aerr**2))
                ax[0].plot(aspec,color='r')
                ax[0].plot(aerr,color='r',ls='dotted')
                pprint(apars)
                print('rot: ',fp[7])
                print(aspcap_param[j]['FPARAM_CLASS'][0:3],aspcap_param[j]['CHI2_CLASS'][0:3])
                # NN model with ASPCAP params
                fit=spectrum(pix, *apars)
                print('nn(ASPCAP) chi2',np.nansum((out[i]['SPEC']-fit)**2/out[i]['ERR']**2))
                print('nn(ASPCAP) chi2 with aerr',np.nansum((out[i]['SPEC']-fit)**2/aerr**2))
                ax[1].plot(fit,color='g')
                ax[2].plot((out[i]['SPEC']-fit)**2/out[i]['ERR']**2,color='g')
                # ASPCAP model
                aspec=aspcap.aspcap2apStar(aspcap_spec[j]['spec_bestfit'])
                ax[1].plot(aspec,color='r')
                plt.draw()
                plt.show()
                pdb.set_trace()
      except : pass

    apars_all=np.array(apars_all)
    npars_all=np.array(npars_all)
    fig,ax=plots.multi(2,1,hspace=0.001,wspace=0.001)
    plots.plotc(ax[0],npars_all[:,0],npars_all[:,1],npars_all[:,2],xr=[8000,3000],yr=[5,0],zr=[-2,0.5])
    plots.plotc(ax[1],apars_all[:,0],apars_all[:,1],apars_all[:,2],xr=[8000,3000],yr=[5,0],zr=[-2,0.5])
    if save is not None :
        fig.savefig(save+'_hr.png')

    fig,ax=plots.multi(2,7,hspace=0.001,wspace=0.001)
    yt=['Teff','logg','[M/H]','[alpha/M]','[C/M]','[N/M]','vmicro']
    for i in range(7) :
      plots.plotc(ax[i,0],apars_all[:,0],npars_all[:,i]-apars_all[:,i],apars_all[:,2],xr=[8100,2900],zr=[-2,0.5],size=20,yt=yt[i])
      plots.plotc(ax[i,1],apars_all[:,i],npars_all[:,i]-apars_all[:,i],apars_all[:,0],zr=[2900,8100],size=20)
    if save is not None :
        fig.savefig(save+'.png')

def pprint(pars) :
    print('{:8.1f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}'.format(*pars))


def showtime(string) :
    """ Utiltiy routine to print a string and clock time
    """
    print(string+' {:8.2f}'.format(time.time()))
    sys.stdout.flush()

