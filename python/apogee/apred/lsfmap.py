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
from tools import html
import os
import pdb
from sklearn.cluster import KMeans

def maps(out='R.dat') :
    fp=open(out,'w')
    load.apred='current'
    load.instrument='apogee-s'
    lsf=load.apLSF(22670019)
    wave=load.apLSF(22600042)
    f,r=modelmap(lsf,wave)
    fp.write("apogee-s 22670019\n")
    stats(r,fp)

    load.instrument='apogee-n'
    lsf=load.apLSF(13400033)
    wave=load.apLSF(13400053)
    f,r=modelmap(lsf,wave)
    fp.write("apogee-n 13400033\n")
    stats(r,fp)

    lsf=load.apLSF(5440020)
    wave=load.apLSF(13400053)
    f,r=modelmap(lsf,wave)
    fp.write("apogee-n 05440020\n")
    stats(r,fp)
    fp.close()

def stats(r,fp) :
    for fib in [6,50,94] :
        fp.write("{:3d}".format(3*fib))
        for ichip in range(3) :
          fp.write("{:8.1f}{:8.1f}{:8.1f}   ".format(np.median(r[fib,:,ichip],axis=0),np.min(r[fib,:,ichip],axis=0),np.max(r[fib,:,ichip],axis=0)))
        fp.write("\n")

def modelmap(lsfframe,waveframe=None,cols=np.arange(5,2000,20),fibers=np.arange(1,300,3),apStar=False,smooth=5,hard=False) :
    '''
    Make LSF map from a model a[ps]LSF file
    '''

    lsfmap=np.zeros([len(fibers),len(cols),3])
    w=np.zeros([len(fibers),len(cols),3])
    dw=np.zeros([len(fibers),len(cols),3])
    r=np.zeros([len(fibers),len(cols),3])

    for ichip,chip in enumerate(['a','b','c']) :
      lsf=lsfframe[chip][1].data
      if waveframe is not None : wave=waveframe[chip][2].data

      nx=lsf.shape[0]
      x=np.arange(nx)
      col=1000
      fit=fitting.LevMarLSQFitter()
      g_init=models.Gaussian1D(mean=nx/2,stddev=2.5/2.354,amplitude=0.3)

      for i,col in enumerate(cols) :
        for j,fiber in enumerate(fibers) :
          y=lsf[:,300-fiber,col]
          g=fit(g_init,x,y)
          if waveframe is not None :
              w[j,i,ichip] = wave[300-fiber,col]
              dw[j,i,ichip] = wave[300-fiber,col]-wave[300-fiber,col-1]
              print(300-fiber,col,g.stddev*2.354,g.stddev*2.354/(6e-6*wave[300-fiber,col]*np.log(10.))*np.abs(dw[j,i,ichip]))
          lsfmap[j,i,ichip] = g.stddev*2.354
          r[j,i,ichip] = w[j,i,ichip]/(g.stddev*2.354*np.abs(dw[j,i,ichip]))
          if apStar and waveframe is not None :
              lsfmap[j,i,ichip]=lsfmap[j,i,ichip]/(6e-6*wave[300-fiber,col]*np.log(10.))*np.abs(dw[j,i,ichip])
      if smooth> 0: 
          lsfmap[:,:,ichip] = convolve(np.squeeze(lsfmap[:,:,ichip]), Box2DKernel(smooth),boundary='extend')
          r[:,:,ichip] = convolve(np.squeeze(r[:,:,ichip]), Box2DKernel(smooth),boundary='extend')
      plt.imshow(lsfmap[:,:,ichip],vmin=1.5,vmax=3.5,interpolation='nearest')    
      plt.draw()
      plt.show()

    plt.close()

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
        xpix=[]
        xlab=[]
        wlab=[15200,15500,15750,15900,16150,16400,16550,16750,16900] 
        if waveframe is not None :
          wave=waveframe[chip][2].data
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
      if hard :
          fig.savefig(lsfframe['filename'].strip('.fits')+suffix+'.png')
          fig.savefig(lsfframe['filename'].strip('.fits')+suffix+'.eps')
          plt.close()
      #pdb.set_trace()

    return lsfmap, r

def datamap(frameid=22600042,waveid=22600042,psfid=22600030,lamp='UNe',fibers=np.arange(50,300,50)) :

    lines=ascii.read(os.environ['APOGEEREDUCE_DIR']+'/lib/linelists/'+lamp+'.vac.apogee')
    pdb.set_trace()
    bright=np.where((lines['FLUX'] > 500) & (lines['USEWAVE'] ==1))[0]
    wbright=lines['WAVE'][bright]
    data=load.ap2D(frameid)
    psf=load.apEPSF(psfid)
    wave=load.apWave(waveid)

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

def group(lsf,wave=None,hard=None,groups=None) :
    """ Plot the FWHM from gaussian fits to the LSF of the 3 chips at column 1024
    """
    fibers=np.arange(1,301)
    fwhm,r=modelmap(lsf,waveframe=wave,cols=[512,1024,1536],fibers=fibers,smooth=0)
    plt.close()
    plt.close()
    fig,ax=plots.multi(2,2,figsize=(8,6))
    plots.plotc(ax[0,0],fwhm[:,1,0],fwhm[:,1,1],fibers,xr=[2.3,3.8],yr=[1.8,3.0],xt='red fwhm, col 1024',yt='green fwhm, col 1024')
    ax[0,0].set_title('Color coded by fiber')
    plots.plotc(ax[0,1],fwhm[:,1,1],fwhm[:,1,2],fibers,xr=[2.0,3.5],yr=[1.8,3.0],xt='green fwhm, col 1024',yt='blue fwhm, col 1024')
    ax[0,1].set_title('Color coded by fiber')
    ax[0,0].grid()
    ax[0,1].grid()
    ax[1,0].grid()
    ax[1,1].grid()
    rfig,rax=plots.multi(2,2,figsize=(8,6))
    plots.plotc(rax[0,0],r[:,1,0],r[:,1,1],fibers,xr=[18000,28000],yr=[18000,28000],xt='red R, col 1024',yt='green R, col 1024')
    rax[0,0].set_title('Color coded by fiber')
    plots.plotc(rax[0,1],r[:,1,1],r[:,1,2],fibers,xr=[18000,28000],yr=[18000,28000],xt='green R, col 1024',yt='blue R, col 1024')
    rax[0,1].set_title('Color coded by fiber')
    rax[0,0].grid()
    rax[0,1].grid()
    rax[1,0].grid()
    rax[1,1].grid()

    # do a Kmeans analysis to split into LSF groups and plot
    X = np.squeeze(fwhm[:,1,:])
    km = KMeans(n_clusters=4)
    labels = km.fit_predict(X)
    for j in range(3) :
      plots.plotc(ax[1,0],fwhm[:,j,0],fwhm[:,j,1],labels,xr=[2.3,3.8],yr=[1.8,3.0],xt='red fwhm, col 1024',yt='green fwhm, col 1024')
      plots.plotc(ax[1,1],fwhm[:,j,1],fwhm[:,j,2],labels,xr=[2.0,3.5],yr=[1.8,3.0],xt='green fwhm, col 1024',yt='blue fwhm, col 1024')
      plots.plotc(rax[1,0],r[:,j,0],r[:,j,1],labels,xr=[18000,28000],yr=[18000,28000],xt='red R, col 1024',yt='green R, col 1024')
      plots.plotc(rax[1,1],r[:,j,1],r[:,j,2],labels,xr=[18000,28000],yr=[18000,28000],xt='green R, col 1024',yt='blue R, col 1024')
    ax[1,0].set_title('Color coded by group')
    ax[1,1].set_title('Color coded by group')
    rax[1,0].set_title('Color coded by group')
    rax[1,1].set_title('Color coded by group')
    gfig,gax=plots.multi(1,1)
    plots.plotp(gax,fibers,labels)
    if groups is not None :
        # show default groups if input
        for start in groups : gax.plot([start-0.5,start-0.5],[-1,5],color='k')

    if hard is not None : 
        fig.tight_layout()
        fig.savefig(hard+'.png')
        rfig.tight_layout()
        rfig.savefig(hard+'_r.png')
        gfig.savefig(hard+'_group.png')
        plt.close()
        plt.close()
        plt.close()

    return fwhm

def parplot(lsf,hard=None) :
    """ Plot the LSF parameters of an LSF fit
    """
    fig,ax=plots.multi(4,3)

    colors=['r','g','b']
    for ichip,chip in enumerate(['a','b','c']) :
      ii=0
      for ipar in range(9,18)+range(24,26) :
        iy=ii%3
        ix=ii//3
        ax[iy,ix].plot(lsf[chip][0].data[ipar,:],color=colors[ichip])
        ii+=1
    if hard is not None : 
        fig.savefig(hard+'.png')
        plt.close()

def sum(apred='r11',telescope='apo25m',lsfs=[3430016,7510018,11130063,14600018,18430026,22330043,25560065],waveid=None,out='apogee-n' ,verbose=False,groups=None) :
    """ Make plots for a series of LSFs and a summary web page
    """
    load=apload.ApLoad(apred=apred,telescope=telescope,verbose=verbose)
    if telescope == 'apo25m' : prefix = 'ap'
    else : prefix = 'as'
    if waveid is not None : wave=load.apWave(waveid)
    else : wave=None

    grid=[]
    ytit=[]
    for lsfid in lsfs :
        lsf=load.apLSF(lsfid)
        name1='pars_{:08d}'.format(lsfid)
        parplot(lsf,hard=name1) 
        name2='fwhm_{:08d}'.format(lsfid)
        group(lsf,wave=wave,hard=name2,groups=groups) 
        grid.append([name1+'.png',name2+'.png',name2+'_r.png',name2+'_group.png'])
        ytit.append('<A HREF={:s}LSF-{:08d}.html>{:08d}</A>'.format(prefix,lsfid,lsfid))

    xt=['LSF parameters','LSF FWHM','LSF R','LSF groups']

    html.htmltab(grid,xtitle=xt,ytitle=ytit,file=out+'.html')

def fibergroups(groups) :
    for i in range(len(groups)-1) :
        n=groups[i+1]-groups[i]
        for j in range(5) : print(int(groups[i]+j*n/5.+n/10.),end=" ")
        print("")

def dr16() :
    """ Make summary LSF pages/plots for DR16 LSFs
    """
    groups=[1,50,146,246,301]
    fibergroups(groups)
    sum(apred='r11',telescope='apo25m',lsfs=[3430016,7510018,11130063,14600018,18430026,22330043,25560065],waveid=24040000,out='apogee-n',groups=groups)
    groups=[1,32,89,151,301]
    fibergroups(groups)
    sum(apred='r11',telescope='lco25m',lsfs=[22940020,26990075],out='apogee-s', waveid=24040000,groups=groups)

