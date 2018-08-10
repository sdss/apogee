from __future__ import print_function
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

def maps(out='R.dat') :
    fp=open(out,'w')
    apload.apred='current'
    apload.instrument='apogee-s'
    f,r=modelmap(lsfid=22670019,waveid=22600042)
    fp.write("apogee-s 22670019\n")
    stats(r,fp)
    apload.instrument='apogee-n'
    f,r=modelmap(lsfid=13400033,waveid=13400053)
    fp.write("apogee-n 13400033\n")
    stats(r,fp)
    f,r=modelmap(lsfid=5440020,waveid=13400053)
    fp.write("apogee-n 05440020\n")
    stats(r,fp)
    fp.close()

def stats(r,fp) :
    for fib in [6,50,94] :
        fp.write("{:3d}".format(3*fib))
        for ichip in range(3) :
          fp.write("{:8.1f}{:8.1f}{:8.1f}   ".format(np.median(r[fib,:,ichip],axis=0),np.min(r[fib,:,ichip],axis=0),np.max(r[fib,:,ichip],axis=0)))
        fp.write("\n")

def modelmap(lsfid=22670019,waveid=22600042,cols=np.arange(5,2000,20),fibers=np.arange(1,300,3),apStar=False,smooth=5) :
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
          y=lsf[:,300-fiber,col]
          g=fit(g_init,x,y)
          w[j,i,ichip] = wave[300-fiber,col]
          dw[j,i,ichip] = wave[300-fiber,col]-wave[300-fiber,col-1]
          print(300-fiber,col,g.stddev*2.354,g.stddev*2.354/(6e-6*wave[300-fiber,col]*np.log(10.))*np.abs(dw[j,i,ichip]))
          lsfmap[j,i,ichip] = g.stddev*2.354
          r[j,i,ichip] = w[j,i,ichip]/(g.stddev*2.354*np.abs(dw[j,i,ichip]))
          if apStar :
              lsfmap[j,i,ichip]=lsfmap[j,i,ichip]/(6e-6*wave[300-fiber,col]*np.log(10.))*np.abs(dw[j,i,ichip])
      lsfmap[:,:,ichip] = convolve(np.squeeze(lsfmap[:,:,ichip]), Box2DKernel(smooth),boundary='extend')
      r[:,:,ichip] = convolve(np.squeeze(r[:,:,ichip]), Box2DKernel(smooth),boundary='extend')
      plt.imshow(lsfmap[:,:,ichip],vmin=1.5,vmax=3.5,interpolation='nearest')    
      plt.draw()
      plt.show()

    plt.close()
    lsf=apload.apLSF(lsfid)

    for idata,data in enumerate([lsfmap, r ]) :
      if idata == 0 : 
          zr=[1.5,3.5]
          suffix='_fwhm'
          cmap='jet_r'
      else : 
          zr=[18000,26000]
          suffix='_R'
          cmap='jet'
      fig,ax=plots.multi(3,1,figsize=(12,6),wspace=0.1)
      #im=ax.imshow(np.reshape(r,(len(fibers),3*len(cols)),order='F'),vmin=18000,vmax=25000,interpolation='nearest',cmap='jet_r')    
      for ichip,chip in enumerate(['a','b','c']) :
        im=ax[ichip].imshow(data[:,:,ichip],vmin=zr[0],vmax=zr[1],interpolation='nearest',cmap=cmap,
                            extent=(cols[0],cols[-1],fibers[-1],fibers[0]),origin='upper',aspect=2048./300.)    
        ax[ichip].set_xlabel('Wavelength')
        wave=apload.apWave(waveid)[chip][2].data
        xpix=[]
        xlab=[]
        wlab=[15200,15500,15750,15900,16150,16400,16550,16750,16900] 
        for w in wlab :
          ipix=abs(wave[150,:]-w).argmin()
          print(w,ipix)
          if ipix > 0 and ipix<2000 :
            xpix.append(ipix)
            xlab.append(str(w))
          ax[ichip].set_xticks(xpix)
          ax[ichip].set_xticklabels(xlab)
        if ichip == 0 : 
          ax[ichip].set_ylabel('Fiber')
        else : 
          #ax.set_xticklabels([])
          ax[ichip].set_yticklabels([])
      cbar_ax = fig.add_axes([0.25, 0.05, 0.5, 0.03])
      fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
      plt.show()
      fig.savefig(lsf['filename'].strip('.fits')+suffix+'.png')
      fig.savefig(lsf['filename'].strip('.fits')+suffix+'.eps')
      pdb.set_trace()
      plt.close()

    return lsfmap, r

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
