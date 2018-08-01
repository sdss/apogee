import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.io import fits
from astropy.modeling import models,fitting
from astropy.convolution import convolve, Box2DKernel
from apogee.utils import apload
from pyvista import image
from tools import plots
import os
import pdb

def maps() :
    apload.apred='current'
    apload.instrument='apogee-s'
    modelmap(lsfid=22670019,waveid=22600042)
    apload.instrument='apogee-n'
    modelmap(lsfid=13400033,waveid=13400053)

def modelmap(lsfid=22670019,waveid=22600042,cols=np.arange(5,2000,40),fibers=np.arange(00,300,6),apStar=False) :
    '''
    Make LSF map from a model a[ps]LSF file
    '''

    lsfmap=np.zeros([len(fibers),len(cols),3])
    w=np.zeros([len(fibers),len(cols),3])
    dw=np.zeros([len(fibers),len(cols),3])
    r=np.zeros([len(fibers),len(cols),3])

    for ichip,chip in enumerate(['a','b','c']) :
      print('chip ',chip)
      lsf=apload.apLSF(lsfid)[chip][1].data
      wave=apload.apWave(waveid)[chip][2].data

      nx=lsf.shape[0]
      x=np.arange(nx)
      col=1000
      fit=fitting.LevMarLSQFitter()
      g_init=models.Gaussian1D(mean=nx/2,stddev=2.5/2.354,amplitude=0.3)

      for i,col in enumerate(cols) :
        for j,fiber in enumerate(fibers) :
          y=lsf[:,fiber,col]
          g=fit(g_init,x,y)
          w[j,i,ichip] = wave[fiber,col]
          dw[j,i,ichip] = wave[fiber,col]-wave[fiber,col-1]
          print(fiber,col,g.stddev*2.354,g.stddev*2.354/(6e-6*wave[fiber,col]*np.log(10.))*np.abs(dw[j,i,ichip]))
          lsfmap[j,i,ichip] = g.stddev*2.354
          r[j,i,ichip] = w[j,i,ichip]/(g.stddev*2.354*np.abs(dw[j,i,ichip]))
          if apStar :
              lsfmap[j,i,ichip]=lsfmap[j,i,ichip]/(6e-6*wave[fiber,col]*np.log(10.))*np.abs(dw[j,i,ichip])
      lsfmap[:,:,ichip] = convolve(np.squeeze(lsfmap[:,:,ichip]), Box2DKernel(11),boundary='extend')
      r[:,:,ichip] = convolve(np.squeeze(r[:,:,ichip]), Box2DKernel(11),boundary='extend')
      plt.imshow(lsfmap[:,:,ichip],vmin=1.5,vmax=3.5,interpolation='nearest')    
      plt.draw()
      plt.show()

    lsf=apload.apLSF(lsfid)
    fig,ax=plots.multi(1,1,figsize=(12,6))
    im=ax.imshow(np.reshape(lsfmap,(len(fibers),3*len(cols)),order='F'),vmin=1.5,vmax=3.5,interpolation='nearest')    
    ax.set_xlabel('Wavelength')
    ax.set_ylabel('Fiber')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    plt.colorbar(im,ax=ax,orientation='horizontal')
    plt.show()
    fig.savefig(lsf['filename'].strip('.fits')+'_fwhm.png')
    fig.savefig(lsf['filename'].strip('.fits')+'_fwhm.eps')

    plt.close()
    fig,ax=plots.multi(1,1,figsize=(12,6))
    im=ax.imshow(np.reshape(r,(len(fibers),3*len(cols)),order='F'),vmin=18000,vmax=25000,interpolation='nearest',cmap='jet_r')    
    ax.set_xlabel('Wavelength')
    ax.set_ylabel('Fiber')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    plt.colorbar(im,ax=ax,orientation='horizontal')
    plt.show()
    fig.savefig(lsf['filename'].strip('.fits')+'_R.png')
    fig.savefig(lsf['filename'].strip('.fits')+'_R.eps')
    plt.close()
    pdb.set_trace()

def datamap(frameid=22600042,waveid=22600042,psfid=22600030,lamp='UNe',fibers=np.arange(50,300,50)) :

    lines=ascii.read(os.environ['APOGEEREDUCE_DIR']+'/lib/linelists/'+lamp+'.vac.apogee')
    pdb.set_trace()
    bright=np.where((lines['FLUX'] > 500) & (lines['USEWAVE'] ==1))[0]
    wbright=lines['WAVE'][bright]
    data=apload.ap2D(frameid)
    psf=apload.apEPSF(psfid)
    wave=apload.apWave(waveid)

    fig,ax=plots.multi(1,3)
    lsfmap=np.zeros([300,2048,3])
    for ichip,chip in enumerate(['a','b','c']) :
      x=[]
      y=[]
      z=[]
      for w in wbright :
        for fiber in fibers :
            col=np.abs(wave[chip][2].data[fiber,:]-w).argmin()
            if (col>50) and (col<2000) :
                row=psf[chip][fiber].data['CENT'][0,col]
                g=image.gfit(data[chip][1].data,col,row,sub=False,pafixed=True,size=3)
                print(chip,w,fiber,col,row,g[0].x_stddev.value*2.354,g[0].y_stddev.value*2.354)
                lsfmap[fiber,col,ichip] = g[0].x_stddev.value*2.354
                x.append(col)
                y.append(row)
                z.append(g[0].x_stddev.value*2.354)
                pdb.set_trace()
      x=np.array(x)
      y=np.array(y)
      z=np.array(z)
      plots.plotc(ax[ichip],x,y,z,zr=[2,4],size=50) 
      pdb.set_trace()
