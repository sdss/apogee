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
import pickle
import pdb
from tools import plots

nepochs=25000
nodes=20
reg=0.0005
threads=24
batch_size=1000

def train(file,plot=False,pixels=[1000,9000,1000],suffix='') :
    """ Train a neural net model on an input training set
    """

    # Get input spectra and parameters
    pixels=np.arange(pixels[0],pixels[1],pixels[2])
    npix=len(pixels)
    print('npix: ', npix)
    pars=fits.open(file+'.fits')[0].data
    spec=fits.open(file+'.fits')[2].data

    # limit parameter range
    #gd=np.where((pars[:,2]>-1.5)&(abs(pars[:,3]) < 0.35) & (abs(pars[:,4]) < 0.35) & (abs(pars[:,5]) < 0.35) )[0]
    #spec=spec[gd,:]
    #pars=pars[gd,:]
    #print(len(gd))

    # limit parameters?
    pars=pars[:,0:7]
    #for ipar in range(2,6) : pars[:,ipar] = 10.**pars[:,ipar]

    #normalize spectra
    for i in range(spec.shape[0]) :
        spec[i,:] /= np.nanmean(spec[i,:])

    if plot :
        fig,ax=plots.multi(2,2)
        plots.plotc(ax[0,0],pars[:,0],spec[:,1000],pars[:,1])
        plots.plotc(ax[1,0],pars[:,0],spec[:,1000],pars[:,2])
    print(spec.shape,pars.shape)
    shape=pars.shape

    # shuffle them and get fit and validation set
    p=np.random.permutation(shape[0])
    spec=spec[p,:]
    pars=pars[p,:]
    nfit=int(len(p)*7/8.)
    shape=pars.shape

    # scale parameters to zero mean and unit standard deviation, and save scaling parameters
    pmeans=[]
    pstds=[]
    for i in range(shape[1]) :
      mn=pars[:,i].mean()
      std=pars[:,i].std()
      pars[:,i] -= mn
      if std > 0. : pars[:,i] /= std
      pmeans.append(mn)
      pstds.append(std)

    # replot to check
    if plot :
        plots.plotc(ax[0,1],pars[:,0],spec[:,1000],pars[:,1])
        plots.plotc(ax[1,1],pars[:,0],spec[:,1000],pars[:,2])
        plt.show()

    # loop over the requested pixels and normalize data to
    #   zero mean and unit standard deviation: save parameters
    weights=[]
    biases=[]
    means=[]
    stds=[]
    data=[]
    for ipix in pixels :
        pix=spec[:,ipix]
        mn=pix.mean()
        std=pix.std()
        if np.isfinite(mn) :
          pix-=mn
          pix /= std
          data.append((pars,pix,ipix))
        means.append(mn)
        stds.append(std)

    # get the model in parallel for different pixels
    print('len: ',len(data))
    pool = mp.Pool(threads)
    output = pool.map_async(fit, data).get()
    pool.close()
    pool.join()
    print('done pool')

    if plot: fig,ax=plots.multi(npix,3,wspace=0.001,hspace=0.001,figsize=(20,8))
    ifit=0
    for i,ipix in enumerate(pixels) :
      print('ipix: ', ipix,plot)
      if np.isfinite(means[i]) :
        w,b,mod,loss=output[ifit]
        ifit+=1
        if plot :
              mod=mod*stds[i]+means[i]
              m=[]
              for ip in range(pars.shape[0]) : m.append(model(pars[ip,:],means[i],stds[i],w,b)[0])
              pix=spec[:,ipix]*stds[i]+means[i]
              ax[0,i].plot(pars[:,0],pix,'ko')
              ax[0,i].plot(pars[:,0],mod,'ro')
              ax[0,i].set_xlim(-2.5,2.5)
              ax[0,i].set_ylim(0.5,1.5)
              n=len(loss)
              plots.plotl(ax[1,i],range(n),np.log10(loss),xr=[0,nepochs],yr=[-4,0])
              ax[2,i].hist(pix-mod[:,0],bins=np.logspace(-7,3,50),histtype='step',normed=True,cumulative=True,color='k')
              ax[2,i].hist(pix[0:nfit]-mod[0:nfit,0],bins=np.logspace(-7,3,50),histtype='step',normed=True,cumulative=True,color='b')
              ax[2,i].hist(pix[nfit:]-mod[nfit:,0],bins=np.logspace(-7,3,50),histtype='step',normed=True,cumulative=True,color='r')
              ax[2,i].set_xlim(0,0.01)
              ax[2,i].set_ylim(0,1.1)
              plt.draw()
              plt.show()
  
        weights.append(w)
        biases.append(b)
    if plot: fig.savefig(file+suffix+'.jpg')


    # Saving the objects:
    with open(file+suffix+'.pkl', 'w') as f:  # Python 3: open(..., 'wb')
        pickle.dump([pmeans, pstds, means, stds, weights, biases], f)

    return pmeans, pstds, means, stds, weights, biases

def report(data) :
    print("done: ",data[2])

def fit(data) :
    """ Routine to do a single NN model fit given input data=(pars,pix)
    """
    pars=data[0]
    pix=data[1]

    net=models.Sequential()
    net.add(layers.Dense(nodes, activation='sigmoid', input_shape=(pars.shape[1],),
            kernel_regularizer=regularizers.l2(reg)))
    #net.add(layers.Dense(1, activation='sigmoid'))
    net.add(layers.Dense(1, activation='linear'))
    #opt=optimizers.RMSprop(lr=0.01)
    opt=optimizers.Adam(lr=0.001)
    net.compile(optimizer=opt,loss='mse')
    net.summary()

    history=net.fit(pars,pix,epochs=nepochs,batch_size=batch_size,verbose=1)

    w=(net.get_weights()[0],net.get_weights()[2])
    b=(net.get_weights()[1],net.get_weights()[3])
    mod=net.predict(pars)

    return w,b,mod,history.history['loss']

# define sigmoid function
def sigmoid(z):
    return 1.0/(1.0+np.exp(-z))
 
def model(pars, mn, std, weights, biases) :
    return mn + std * (np.dot( sigmoid((np.dot(weights[0].T,pars)+biases[0])).T, weights[1] ) +biases[1])

def test(pmn, pstd, mn, std, weights, biases,n=100, t0=[3750.,4500.], g0=2., mh0=0.) :
    fig,ax=plots.multi(2,6,figsize=(8,12))

    xt=['Teff','logg','[M/H]','[alpha/M]','[C/M]','[N/M]']
    for i,ipar in enumerate([0,1,2,3,4,5]) : 
      for ipix in range(len(weights)) :
       for it0 in range(2) :
        pars=np.tile([t0[it0], g0, mh0, 0.0, 0., 0., 2.],(n,1))
        if ipar == 0 : pars[:,ipar]=np.linspace(3500.,5000.,n)
        elif ipar == 1 : pars[:,ipar]=np.linspace(0.,5.,n)
        elif ipar == 2 : pars[:,ipar]=np.linspace(-2.5,1.,n)
        elif ipar == 3 : pars[:,ipar]=np.linspace(-0.5,0.75,n)
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


def get(file) :
    global pmeans, pstds, means, stds, weights, biases, ifit

    # Getting back the objects:
    with open(file+'.pkl') as f: 
        pmeans, pstds, means, stds, weights, biases = pickle.load(f)

    ifit = np.zeros(len(means)).astype(int)
    j=0
    for i in range(len(means)) :
        if np.isfinite(means[i]) : 
            ifit[i] = j
            j += 1
    pmeans = np.array(pmeans)
    pstds = np.array(pstds)

    #return np.array(pmeans), np.array(pstds), means, stds, weights, biases 

def spectrum(x,*pars) :
    """ Return spectrum given input list of pixels, parameters
    """
    spec=np.zeros(len(x))
    pnorm= (pars-pmeans)/pstds
    for j,i in enumerate(x) :
        if  np.isfinite(means[i]) :
            spec[j]= model(pnorm, means[i], stds[i], weights[ifit[i]], biases[ifit[i]]) 

    return spec

def comp(file,threads=8,nfit=8) :
    p=fits.open(file+'.fits')[0].data
    s=fits.open(file+'.fits')[2].data
    p=p[:,0:9]
    if nfit == 0 : nfit = p.shape[0]
    get(file)

    specs=[]
    #for i in range(s.shape[0]) :
    for i in range(nfit) :
        specs.append(s[i,:]/np.nanmean(s[i,:]))

    pool = mp.Pool(threads)
    output = pool.map_async(solve, specs).get()
    pool.close()
    pool.join()
 
    pix = np.arange(0,8575,1)
    for i in range(nfit) :
        print(p[i,:])
        print(output[i])
        snorm = s[i,:] / np.nanmean(s[i,:])
        spec=spectrum(pix, *p[i,:])
        fit=spectrum(pix, *output[i])
        gd=np.where(np.isfinite(snorm))[0]
        print(np.sum((spec[gd]-snorm[gd])**2),np.sum((fit[gd]-snorm[gd])**2))
        plt.clf()
        plt.plot(snorm)
        plt.plot(spec)
        plt.plot(fit)
        plt.draw()

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

    pdb.set_trace()

def solve(spec) :
    """ Solve for parameters for a single input spectrum
    """
    pix = np.arange(0,8575,1)
    sig = np.ones(len(pix))
    init = np.array([4000.,2.5,0.,0.,0.,0.,1.5,0.,0.])
    bounds = (np.array([3000.,-0.5,-3.,-1.,-1.,-1.,0.5,-0.001,-0.001]),
              np.array([6000., 5.5, 1., 1., 1., 1.,4.5, 0.001, 0.001]))
    gd = np.where(np.isfinite(spec))[0]
    fpars,fcov = curve_fit(spectrum,pix[gd],spec[gd],sigma=sig[gd],p0=init,bounds=bounds)
    return fpars

