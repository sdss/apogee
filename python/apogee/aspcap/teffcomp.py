# routines for comparing gravities with photometric` sample

from apogee.utils import apload
from apogee.utils import apselect
from astropy.io import fits
from astropy.io import ascii
from tools import match
from tools import plots
from tools import fit
from apogee.utils import bitmask
from apogee.aspcap import err
import pdb
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib

def bindata(xdata,ydata,bins,median=True) :
    """
    Given xdata, ydata, and bins in x, returns mean of ydata in each of the bins
    """
    mean=bins*0.
    for i in range(len(bins)-1) :
      j=np.where((xdata>bins[i]) & (xdata<bins[i+1]))[0]
      if median :
          mean[i]=np.median(ydata[j])
      else :
          mean[i]=ydata[j].mean() 
    return mean

def ghb(allstar,glatmin=30.,ebvmax=0.03,trange=[3750,5500],loggrange=[-1,6],mhrange=[-2.5,0.75],alpha=False,out='teffcomp',yr=[-500,500],
        calib=False,dr13=False,grid=None) :
    """
    Compares allstar ASPCPAP Teff with photometric Teff from GHB for sample of stars with GLAT>glatmin and SFD_EBV<ebvmax,
    does fits

    Args:
        allstar   : allStar structure

    Keyword args:
        glatmin (float) : minimum GLAT for sample (default=30)
        ebvmax (float)  : maximum SFD_EBV for sample (default=0.03)
        dwarf (bool)    : use dwarfs and dwarf GHB  (default = False)
    """

    # select data to use
    if 'TARGFLAGS' in allstar.columns.names : badtarg = ['EMBEDDED','EXTENDED']
    else : badtarg = None
    gd=apselect.select(allstar,badval=['STAR_BAD'],badtarg=badtarg,teff=trange,mh=mhrange,logg=loggrange,raw=True)
    allstar=allstar[gd]
    #if dr13 :
    #  j=np.where((abs(allstar['GLAT'])>glatmin)&(allstar['SFD_EBV']<ebvmax))[0]
    #else :
    j=np.where((abs(allstar['GLAT'])>glatmin)&(allstar['SFD_EBV']>-0.01)&(allstar['SFD_EBV']<ebvmax)&(abs(allstar['J'])<90)&(abs(allstar['K'])<90))[0]

    # remove second gen GC stars
    #if not dr13 :
    gcstars = ascii.read(os.environ['APOGEE_DIR']+'/data/calib/gc_szabolcs.dat')
    bd=np.where(gcstars['pop'] != 1)[0]
    j = [x for x in j if allstar[x]['APOGEE_ID'] not in gcstars['id'][bd]]

    allstar=allstar[j]

    ghb,dtdjk=cte_ghb(allstar['J']-allstar['K'],allstar['FPARAM'][:,3],dwarf=False)
    ghb_dwarf,dtdjk_dwarf=cte_ghb(allstar['J']-allstar['K'],allstar['FPARAM'][:,3],dwarf=True)
    # use dwarf relation for dwarfs
    dw=np.where(allstar['FPARAM'][:,1] > 3.8)[0]
    ghb[dw]=ghb_dwarf[dw]
    dtdjk[dw]=dtdjk_dwarf[dw]
    gd=np.where(abs(allstar['FPARAM'][:,0]-ghb) < 500)[0]
    ghb=ghb[gd]
    dtdjk=dtdjk[gd]
    allstar=allstar[gd]
    print('Teff calibration, number of stars: ', len(allstar))

    if calib : 
        param='PARAM'
        teff=allstar[param][:,0]
        logg=allstar[param][:,1]
        mh=allstar[param][:,3]
        am=allstar[param][:,6]
    elif grid is None :
        param='FPARAM'
        teff=allstar[param][:,0]
        logg=allstar[param][:,1]
        mh=allstar[param][:,3]
        am=allstar[param][:,6]
    else :
        param='FPARAM_CLASS'
        teff=allstar[param][:,grid,0]
        logg=allstar[param][:,grid,1]
        mh=allstar[param][:,grid,3]
        am=allstar[param][:,grid,6]
        out=out+'_grid{:1d}'.format(grid)

    # plot Teff difference against metallicity, color-code by temperature
    fig,ax=plots.multi(1,1,hspace=0.001,wspace=0.001,figsize=(12,6))
    xr=[-3.0,1.0]
    zr=trange
    if dr13: zr=[3500,5500]
    binsize=0.25
    bins=np.arange(-2.5,0.75,binsize)
    # diff color-coded by gravity as f([M/H])

    if alpha :
        plots.plotc(ax,mh,teff-ghb,am,zr=[-0.1,0.4],xr=xr,yr=yr,xt='[M/H]',yt='ASPCAP-photometric Teff',colorbar=True,zt=r'[$\alpha$/M]',rasterized=True)
    else :
        plots.plotc(ax,mh,teff-ghb,teff,xr=xr,yr=yr,xt='[M/H]',yt='ASPCAP-photometric Teff',colorbar=True,zt='$T_{eff}$',rasterized=True,zr=trange)
    mean=bindata(mh,teff-ghb,bins,median=False)
    if not dr13: plots.plotp(ax,bins+binsize/2.,mean,marker='o',size=40)
    mean=bindata(mh,teff-ghb,bins,median=True)
    if not dr13: plots.plotp(ax,bins+binsize/2.,mean,marker='o',size=40,color='b')
    ax.text(0.1,0.9,'E(B-V)<{:6.2f}'.format(ebvmax),transform=ax.transAxes)
    gd=np.where(np.isfinite(mean))[0]
    tefit = fit.fit1d(bins[gd]+binsize/2.,mean[gd],degree=2,reject=0)
    # 1D quadratic fit as a function of metallicity
    allfit = fit.fit1d(mh,teff-ghb,ydata=teff,degree=2,reject=0)
    fig2,ax2=plots.multi(1,1)
    tefit2 = fit.fit2d(mh,teff,teff-ghb,reject=0,plot=ax2,zr=[-500,200],xt='[M/H]',yt=['Teff'],zt='$\Delta Teff$')
    #pfit = fit.fit2d(allstar[param][:,3],allstar[param][:,0],allstar[param][:,0]-ghb,plot=ax[0,0],zr=[-500,200],xt='[M/H]',yt=['Teff'],zt='$\Delta Teff$')
    #ejk=np.clip(np.sqrt(allstar['J_ERR']**2+allstar['K_ERR']**2),0.,0.02)
    #errpar = err.errfit(teff,allstar['SNR'],mh,teff-tefit(mh)-ghb,title='Teff',out=out+'_phot',zr=[0,250],meanerr=abs(dtdjk)*ejk)
    errpar = err.errfit(teff,allstar['SNR'],mh,teff-tefit(mh)-ghb,title='Teff',out=out,zr=[0,150])

    x=np.linspace(-3,1,200)
    rms = (teff-tefit(mh)-ghb).std()
    if dr13: 
      plots.plotl(ax,x,-36.17+95.97*x-15.09*x**2,color='k')
      print(allfit)
    else :
      plots.plotl(ax,x,tefit(x),color='k')
      ax.text(0.98,0.9,'rms: {:6.1f}'.format(rms),transform=ax.transAxes,ha='right')

      cmap = matplotlib.cm.get_cmap('rainbow')
      for t in np.arange(trange[0],trange[1],500.) :
          rgba=cmap((t-trange[0])/(trange[1]-trange[0]))
          y=x*0.+t
          plots.plotl(ax,x,tefit2(x,y),color=rgba)

    plots._data_x = mh
    plots._data_y = teff-ghb
    plots._data = allstar
    plots.event(fig)

    # separate fits for low/hi alpha/M if requested
    if alpha :
        gdlo=apselect.select(allstar,badval=['STAR_BAD'],teff=trange,mh=mhrange,logg=[0,3.8],alpha=[-0.1,0.1],raw=True)
        mean=bindata(mh[gdlo],teff[gdlo]-ghb[gdlo],bins)
        plots.plotp(ax,bins,mean,marker='o',size=40,color='g')
        tmpfit = fit.fit1d(mh[gdlo],teff[gdlo]-ghb[gdlo],ydata=teff[gdlo],degree=2)
        plots.plotl(ax,x,tmpfit(x))
        print('low alpha: ', len(gdlo))

        gdhi=apselect.select(allstar,badval=['STAR_BAD'],teff=trange,mh=mhrange,logg=[0,3.8],alpha=[0.1,0.5],raw=True)
        mean=bindata(mh[gdhi],teff[gdhi]-ghb[gdhi],bins)
        plots.plotp(ax,bins,mean,marker='o',size=40,color='b')
        tmpfit = fit.fit1d(mh[gdhi],teff[gdhi]-ghb[gdhi],ydata=teff[gdhi],degree=2)
        plots.plotl(ax,x,tmpfit(x))
        print('hi alpha: ', len(gdhi))

    fig.tight_layout()
    fig.savefig(out+'.png')
    plt.close()
    plt.rc('font',size=14)
    plt.rc('axes',titlesize=14)
    plt.rc('axes',labelsize=14)
    fig.savefig(out+'.pdf')
    plt.close()

    # auxiliary plots with different color-codings
    try:
        meanfib=allstar['MEANFIB']
    except:
        meanfib=teff*0.
    fig,ax=plots.multi(2,2,hspace=0.001,wspace=0.001)
    plots.plotc(ax[0,0],mh,teff-ghb,logg,zr=[0,5],xr=xr,yr=yr,xt='[M/H]',yt='ASPCAP-photometric Teff',colorbar=True,zt='log g')
    plots.plotc(ax[0,1],mh,teff-ghb,meanfib,zr=[0,300],xr=xr,yr=yr,xt='[M/H]',yt='ASPCAP-photometric Teff',colorbar=True,zt='mean fiber')
    pfit = fit.fit1d(mh,teff-ghb,ydata=teff,plot=ax[1,0],zr=[-500,200],xt='[M/H]',yt='$\Delta Teff$',xr=[-2.7,0.9],yr=[3500,5000],colorbar=True,zt='Teff')
    pfit = fit.fit1d(teff,teff-ghb,ydata=mh,plot=ax[1,1],zr=[-500,200],xt='Teff',xr=trange,yr=[-2.5,0.5],colorbar=True,zt='[M/H]')
    fig.tight_layout()
    fig.savefig(out+'_b.png')
    plt.close()
   
    # do some test 2D and 1D fits and plots 
    #fig,ax=plots.multi(2,2,hspace=0.5,wspace=0.001)
    #ax[0,1].xaxis.set_visible(False)
    #ax[0,1].yaxis.set_visible(False)
    #pfit = fit.fit2d(allstar[param][:,3],allstar[param][:,0],allstar[param][:,0]-ghb,plot=ax[0,0],zr=[-500,200],xt='[M/H]',yt=['Teff'],zt='$\Delta Teff$')
    #pfit = fit.fit1d(allstar[param][:,3],allstar[param][:,0]-ghb,ydata=allstar[param][:,0],plot=ax[1,0],zr=[-500,200],xt='[M/H]',yt='$\Delta Teff$',xr=[-2.7,0.9],yr=[3500,5000])
    #pfit = fit.fit1d(allstar[param][:,0],allstar[param][:,0]-ghb,ydata=allstar[param][:,3],plot=ax[1,1],zr=[-500,200],xt='Teff',xr=[3900,5100],yr=[-2.5,0.5])
    plt.draw()
    return {'caltemin': 3000., 'caltemax': 10000., 'temin' : trange[0], 'temax': trange[1], 
            'mhmin': mhrange[0], 'mhmax' : mhrange[1],
            'par': tefit.parameters, 'rms' :rms, 'par2d': tefit2.parameters, 'errpar' : errpar}


def irfm(allstar,trange=[4000,5000],mhrange=[-2.5,0.75],out='dteff') :
    '''
    Compares allstar ASPCPAP Teff with various photometric Teff from JAJ compilation (SAGA, CL, TH, SFD)
    Does fits 

    Args:
        allstar   : allStar structure

    '''

    # select stars
    gd=apselect.select(allstar,badval=['STAR_BAD'],teff=trange,mh=mhrange,raw=True)
    allstar=allstar[gd]

    # get IRFM data
    irfm=fits.open(os.environ['APOGEE_DIR']+'/data/calib/irfm_temp.fits')[1].data

    # get the subsamples and match. Note that we have to do this separately for each subsample because some
    #   stars appear in more than one subsample
    saga=np.where(irfm['SOURCE'] == 'SAGA')[0]
    saga1,saga2=match.match(np.chararray.strip(allstar['APOGEE_ID']),np.chararray.strip(irfm['2MASS ID'][saga]))
    cl=np.where(irfm['SOURCE'] == 'CL')[0]
    cl1,cl2=match.match(np.chararray.strip(allstar['APOGEE_ID']),np.chararray.strip(irfm['2MASS ID'][cl]))
    th=np.where(irfm['SOURCE'] == 'TH')[0]
    th1,th2=match.match(np.chararray.strip(allstar['APOGEE_ID']),np.chararray.strip(irfm['2MASS ID'][th]))
    sfd=np.where(irfm['SOURCE'] == 'SFD')[0]
    sfd1,sfd2=match.match(np.chararray.strip(allstar['APOGEE_ID']),np.chararray.strip(irfm['2MASS ID'][sfd]))

    # plot diff color-coded by gravity as f([M/H])
    fig,ax=plots.multi(2,2,hspace=0.001,wspace=0.001)
    xr=[-3.0,1.0]
    yr=[-400,300]
    zr=[3500,6000]
    bins=np.arange(-2.5,0.75,0.25)

    # SAGA
    plots.plotc(ax[0,0],allstar['FPARAM'][saga1,3],allstar['FPARAM'][saga1,0]-irfm['IRFM TEFF'][saga[saga2]],allstar['FPARAM'][saga1,0],zr=zr,xr=xr,yr=yr,xt='[M/H]',yt='ASPCAP-photometric Teff')
    mean=bindata(allstar['FPARAM'][saga1,3],allstar['FPARAM'][saga1,0]-irfm['IRFM TEFF'][saga[saga2]],bins)
    plots.plotp(ax[0,0],bins,mean,marker='o',size=40)
    ax[0,0].text(0.1,0.9,'SAGA',transform=ax[0,0].transAxes)

    # CL
    plots.plotc(ax[0,1],allstar['FPARAM'][cl1,3],allstar['FPARAM'][cl1,0]-irfm['IRFM TEFF'][cl[cl2]],allstar['FPARAM'][cl1,0],zr=zr,xr=xr,yr=yr,xt='[M/H]')
    mean=bindata(allstar['FPARAM'][cl1,3],(allstar['FPARAM'][cl1,0]-irfm['IRFM TEFF'][cl[cl2]]),bins)
    plots.plotp(ax[0,1],bins,mean,marker='o',size=40)
    ax[0,1].text(0.1,0.9,'CL',transform=ax[0,1].transAxes)

    # TH
    plots.plotc(ax[1,0],allstar['FPARAM'][th1,3],allstar['FPARAM'][th1,0]-irfm['IRFM TEFF'][th[th2]],allstar['FPARAM'][th1,0],zr=zr,xr=xr,yr=yr,xt='[M/H]',yt='ASPCAP-photometric Teff')
    mean=bindata(allstar['FPARAM'][th1,3],(allstar['FPARAM'][th1,0]-irfm['IRFM TEFF'][th[th2]]),bins)
    plots.plotp(ax[1,0],bins,mean,marker='o',size=40)
    ax[1,0].text(0.1,0.9,'TH',transform=ax[1,0].transAxes)

    # SFD
    plots.plotc(ax[1,1],allstar['FPARAM'][sfd1,3],allstar['FPARAM'][sfd1,0]-irfm['IRFM TEFF'][sfd[sfd2]],allstar['FPARAM'][sfd1,0],zr=zr,xr=xr,yr=yr,xt='[M/H]')
    mean=bindata(allstar['FPARAM'][sfd1,3],(allstar['FPARAM'][sfd1,0]-irfm['IRFM TEFF'][sfd[sfd2]]),bins)
    plots.plotp(ax[1,1],bins,mean,marker='o',size=40)
    ax[1,1].text(0.1,0.9,'SFD',transform=ax[1,1].transAxes)

    fig.savefig(out+'_mh.png')

    # plot diff color-coded by gravity as f([M/H])
    fig,ax=plots.multi(2,2,hspace=0.001,wspace=0.001)
    zr=[-2.0,0.5]
    yr=[-400,300]
    xr=[6000,3500]
    bins=np.arange(3500,5500,250)

    # SAGA
    plots.plotc(ax[0,0],allstar['FPARAM'][saga1,0],allstar['FPARAM'][saga1,0]-irfm['IRFM TEFF'][saga[saga2]],allstar['FPARAM'][saga1,3],zr=zr,xr=xr,yr=yr,xt='Teff',yt='ASPCAP-photometric Teff')
    mean=bindata(allstar['FPARAM'][saga1,0],(allstar['FPARAM'][saga1,0]-irfm['IRFM TEFF'][saga[saga2]]),bins)
    plots.plotp(ax[0,0],bins,mean,marker='o',size=40)
    ax[0,0].text(0.1,0.9,'SAGA',transform=ax[0,0].transAxes)

    # CL
    plots.plotc(ax[0,1],allstar['FPARAM'][cl1,0],allstar['FPARAM'][cl1,0]-irfm['IRFM TEFF'][cl[cl2]],allstar['FPARAM'][cl1,3],zr=zr,xr=xr,yr=yr,xt='Teff')
    mean=bindata(allstar['FPARAM'][cl1,0],(allstar['FPARAM'][cl1,0]-irfm['IRFM TEFF'][cl[cl2]]),bins)
    plots.plotp(ax[0,1],bins,mean,marker='o',size=40)
    ax[0,1].text(0.1,0.9,'CL',transform=ax[0,1].transAxes)

    # TH
    plots.plotc(ax[1,0],allstar['FPARAM'][th1,0],allstar['FPARAM'][th1,0]-irfm['IRFM TEFF'][th[th2]],allstar['FPARAM'][th1,3],zr=zr,xr=xr,yr=yr,xt='Teff',yt='ASPCAP-photometric Teff')
    mean=bindata(allstar['FPARAM'][th1,0],(allstar['FPARAM'][th1,0]-irfm['IRFM TEFF'][th[th2]]),bins)
    plots.plotp(ax[1,0],bins,mean,marker='o',size=40)
    ax[1,0].text(0.1,0.9,'TH',transform=ax[1,0].transAxes)

    # SFD
    plots.plotc(ax[1,1],allstar['FPARAM'][sfd1,0],allstar['FPARAM'][sfd1,0]-irfm['IRFM TEFF'][sfd[sfd2]],allstar['FPARAM'][sfd1,3],zr=zr,xr=xr,yr=yr,xt='Teff')
    mean=bindata(allstar['FPARAM'][sfd1,0],(allstar['FPARAM'][sfd1,0]-irfm['IRFM TEFF'][sfd[sfd2]]),bins)
    plots.plotp(ax[1,1],bins,mean,marker='o',size=40)
    ax[1,1].text(0.1,0.9,'SFD',transform=ax[1,1].transAxes)

    fig.savefig(out+'_teff.png')

    # do 2D fits with Teff and [M/H], and 1D fits with each

    fig,ax=plots.multi(2,2,hspace=0.5,wspace=0.001)
    ax[0,1].xaxis.set_visible(False)
    ax[0,1].yaxis.set_visible(False)
    pfit = fit.fit2d(ax[0,0],allstar['FPARAM'][sfd1,3],allstar['FPARAM'][sfd1,0],allstar['FPARAM'][sfd1,0]-irfm['IRFM TEFF'][sfd[sfd2]],plot=True,zr=[-500,200],xt='[M/H]',yt=['Teff'],zt='$\Delta Teff$')
    pfit = fit.fit1d(ax[1,0],allstar['FPARAM'][sfd1,3],allstar['FPARAM'][sfd1,0]-irfm['IRFM TEFF'][sfd[sfd2]],ydata=allstar['FPARAM'][sfd1,0],plot=True,zr=[-500,200],xt='[M/H]',yt='$\Delta Teff$',xr=[-2.7,0.9],yr=[3500,5000])
    pfit = fit.fit1d(ax[1,1],allstar['FPARAM'][sfd1,0],allstar['FPARAM'][sfd1,0]-irfm['IRFM TEFF'][sfd[sfd2]],ydata=allstar['FPARAM'][sfd1,3],plot=True,zr=[-500,200],xt='Teff',xr=[3900,5100],yr=[-2.5,0.5])

    pdb.set_trace()

    return pfit


def dr13dr12() :
    '''
    compare dr13 dr12 Teff
    '''

    dr12load=apload.ApLoad(dr='dr12')
    dr12=dr12load.allStar()[1].data
    dr13load=apload.ApLoad(dr='dr13')
    dr13=dr13load.allStar()[1].data
    i1,i2 = match.match(dr12['APOGEE_ID'],dr13['APOGEE_ID'])
    dr12=dr12[i1]
    dr13=dr13[i2]

    fig,ax=plots.multi(1,2,hspace=0.001,wspace=0.001)
    plots.plotc(ax[0],dr13['M_H'],dr13['TEFF']-dr12['TEFF'],dr13['TEFF'],xr=[-2.5,0.75],yr=[-300,300],zr=[3500,5000])

    plots.plotc(ax[1],dr13['TEFF'],dr13['TEFF']-dr12['TEFF'],dr13['M_H'],xr=[6500,3000],yr=[-300,300],zr=[-2,0.5])

def cte_ghb(jk0,feh,dwarf=False) :
    """
    Color-temperature relation from Gonzalez Hernandez & Bonifacio (2009):  (J-K)_0 - Teff
    """
    if dwarf :
      b0=0.6524 ; b1=0.5813 ; b2=0.1225 ; b3=-0.0646 ; b4=0.0370 ; b5=0.0016 # dwarfs
    else :
      b0=0.6517 ; b1=0.6312 ; b2=0.0168 ; b3=-0.0381 ; b4=0.0256 ; b5=0.0013 # giants
    theta=b0+b1*jk0+b2*jk0**2+b3*jk0*feh+b4*feh+b5*feh**2
    dtheta_djk = b1+2*b2*jk0+b3*feh
    dt_djk= -5040./theta**2*dtheta_djk

    return 5040./theta, dt_djk

