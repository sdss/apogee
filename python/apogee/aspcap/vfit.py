from __future__ import print_function
from astropy.io import fits
from astropy.io import ascii
from astropy.modeling import models, fitting
import astropy.constants as const
import matplotlib.pyplot as plt
import numpy as np
import math
import pdb
import matplotlib as mpl
import sys
import argparse
import os
from tools import fit
from tools import plots
from tools import match
from apogee.utils import bitmask

def fit_vmicro(data,teffrange=[3550,5500],mhrange=[-2.5,1],loggrange=[-0.5,4.75],vrange=[0,4],vmrange=[0,6],maxerr=0.1, degree=1,reject=0,apokasc='APOKASC_cat_v3.6.0.fits',nopersist=False,out='vmicro',func=None) :
    """ 
    Fit microturbulence relation  with 1D f(log g) and 2D f(Teff, logg) fits, plots

    Args:
        data : calib structure

    Keyword args :
        degree : degree of fit (default=1)
        mhrange : range of [M/H] (default=[-1,1])
        loggrange : range of log g (default=[-1,3.8])
        maxerr : maximum uncertainty in vmicro
        vrange  : scaling range for vmacro plot, vrange[1] sets maximum good vmacro

    Returns:
        fit1d, fit2d : 1D and 2D polynomial fits
    """
    vmicro = data['FPARAM'][:,2]
    vmacro = data['FPARAM'][:,7]
    # fix locked vmacro by hand (bad!)
    #j=np.where(np.isclose(vmacro,1.))[0]
    #vmacro[j] = 0.6
    teff = data['FPARAM'][:,0]
    logg = data['FPARAM'][:,1]
    mh = data['FPARAM'][:,3]
    try :
        meanfib = data['MEANFIB']
    except :
        meanfib = 0.*vmicro
    try :
        ninst= data['NINST'][:,1]-data['NINST'][:,2]
    except :
        ninst= data['NVISITS']*0

    apokasc=fits.open(os.environ['APOGEE_DIR']+'/data/apokasc/'+apokasc)[1].data
    i1,i2=match.match(data['APOGEE_ID'],apokasc['2MASS_ID'])
    rgb=np.where(apokasc['CONS_EVSTATES'][i2] == 'RGB')[0]
    rc=np.where(apokasc['CONS_EVSTATES'][i2] == 'RC')[0]
    type=mh*0.
    type[i1[rgb]]=1
    type[i1[rc]]=-1

    print('mhrange: ', mhrange)
    print('teffrange: ', teffrange)
    print('loggrange: ', loggrange)
    print('vmacrorange: ', vmrange)
    print('vmicrorange: ', vrange)
    print('maxerr: ', maxerr)
    print('nopersist: ', nopersist)
    print('reject: ', reject)
    gd = np.where((mh>mhrange[0]) & (mh<mhrange[1]) & (logg > loggrange[0]) & (logg < loggrange[1]) & 
                  (teff>teffrange[0]) & (teff<teffrange[1]) & (10.**vmacro>vmrange[0]) & (10.**vmacro<vmrange[1])  &
                  (np.sqrt(data['FPARAM_COV'][:,2,2]) < maxerr) & (10.**vmicro > vrange[0]) & (10.**vmicro < vrange[1]))[0]
    if nopersist :
        j = np.where ((data['STARFLAG'][gd] & bitmask.persist()) == 0)[0]
        gd=gd[j]

    # remove non-1st generation GC stars
    gcstars = ascii.read(os.environ['APOGEE_DIR']+'/data/calib/gc_szabolcs.dat')
    bd=np.where(gcstars['pop'] != 1)[0]
    gd = [x for x in gd if data[x]['APOGEE_ID'] not in gcstars['id'][bd]]


    # 1D plots a f(log g)
    #fig,ax = plots.multi(2,3,figsize=(15,15))
    #fit1d = fit.fit1d(logg[gd], vmicro[gd],degree=degree,reject=reject,plot=ax[0,0],ydata=mh[gd],log=True,xt='log g',yt='vmicro ([M/H])',yr=[-2.5,0.5],colorbar=True,zt='[M/H]',xr=[0,4.5])
    ## plot ALL points (even outside of fit range)
    #plots.plotc(ax[0,0],logg,10.**vmicro,mh,zr=[-2.5,0.5],xr=[-1,5],size=1)
#
#    junk = fit.fit1d(logg[gd], vmicro[gd],degree=degree,reject=reject,plot=ax[0,1],ydata=teff[gd],log=True,xt='log g',yt='vmicro',yr=[3500,5500],pfit=fit1d,colorbar=True,zt='Teff',xr=[0,4.5])
#    plots.plotc(ax[0,1],logg,10.**vmicro,teff,zr=[3500,5500],xr=[-1,5],size=1)
#    junk = fit.fit1d(logg[gd], vmicro[gd],degree=degree,reject=reject,plot=ax[1,0],ydata=meanfib[gd],log=True,xt='log g',yt='vmicro',yr=[0,300],pfit=fit1d,colorbar=True,zt='mean fiber',xr=[0,4.5])
#    plots.plotc(ax[1,0],logg,10.**vmicro,meanfib,zr=[0,300],xr=[-1,5],size=1)
#    junk = fit.fit1d(logg[gd], vmicro[gd],degree=degree,reject=reject,plot=ax[1,1],ydata=10.**vmacro[gd],log=True,xt='log g',yt='vmicro',yr=[0,15],pfit=fit1d,colorbar=True,zt='vmacro',xr=[0,4.5])
#    plots.plotc(ax[1,1],logg,10.**vmicro,10.**vmacro,zr=[0,15],xr=[-1,5],size=1)
#
#    junk = fit.fit1d(logg[gd], vmicro[gd],degree=degree,reject=reject,plot=ax[2,0],ydata=ninst[gd],log=True,xt='log g',yt='vmicro',yr=[-1,1],pfit=fit1d,colorbar=True,zt='ninst1-ninst2',xr=[0,4.5])
#    plots.plotc(ax[2,0],logg,10.**vmicro,ninst,zr=[-1,1],xr=[-1,5],size=1)
#    #plots.plotc(ax[3,1],logg[gd1],10.**vmicro[gd1],mh[gd1],zr=[-2.5,0.5],xr=[-1,5])
#    # 2d plot
#    #junk = fit.fit1d(logg, vmicro,degree=degree,reject=reject,plot=ax[2,0],plot2d=True,ydata=teff,log=True,yt='Teff',xt='log g',xr=[5,-0.5],yr=[6000,3000],pfit=fit1d,zr=[0,4])
#    junk = fit.fit1d(logg[gd], vmicro[gd],degree=degree,reject=reject,plot=ax[2,1],ydata=type[gd],log=True,xt='log g',yt='vmicro',yr=[-1,1],pfit=fit1d,colorbar=True,zt='RGB/RC',xr=[0,4.5])
#    plots.plotc(ax[2,0],logg,10.**vmicro,ninst,zr=[-1,1],xr=[-1,5],size=1)
#
#
#    # plot with DR13 relation
#    #dr13fit=models.Polynomial1D(degree=3)
#    #dr13fit.parameters=[0.226,-0.0228,0.0297,-0.013]
#    #junk = fit.fit1d(logg[gd], vmicro[gd],degree=degree,reject=reject,plot=ax[2,1],ydata=mh[gd],pfit=dr13fit,log=True,xt='log g',yt='vmicro',colorbar=True,zt='[M/H]',yr=[-2.5,0.5])
#    #plots.plotc(ax[2,0],logg,10.**vmicro,mh,zr=[-2.5,0.5],xr=[-1,5],size=1)
##
#    fig.tight_layout()
#    fig.savefig(out+'.png')
#
#    fig,ax=plots.multi(1,3,figsize=(6,10))
#    plots.plotc(ax[0],teff[gd],10.**vmicro[gd]-10.**fit1d(logg[gd]),logg[gd],xt='Teff',yt=r'$\Delta vmicro$',xr=[3500,6000],yr=[-1,1],zr=[0,5],zt='log g',size=5)
#    plots.plotc(ax[0],teff,10.**vmicro-10.**fit1d(logg),logg,xt='Teff',yt=r'$\Delta vmicro$',xr=[3500,6000],yr=[-1,1],zr=[0,5],zt='log g',colorbar=True,size=1)
#    plots.plotc(ax[1],mh[gd],10.**vmicro[gd]-10.**fit1d(logg[gd]),teff[gd],xt='[M/H]',yt=r'$\Delta vmicro$',xr=[-2.5,0.5],yr=[-1,1],zr=[3500,5500],zt='Teff',size=5)
#    plots.plotc(ax[1],mh,10.**vmicro-10.**fit1d(logg),teff,xt='[M/H]',yt=r'$\Delta vmicro$',xr=[-2.5,0.5],yr=[-1,1],zr=[3500,5500],zt='Teff',colorbar=True,size=1)
#    plots.plotc(ax[2],meanfib[gd],10.**vmicro[gd]-10.**fit1d(logg[gd]),logg[gd],xt='mean fiber',yt=r'$\Delta vmicro$',xr=[0,300],yr=[-1,1],zr=[0,5],zt='log g',size=5)
#    plots.plotc(ax[2],meanfib,10.**vmicro-10.**fit1d(logg),logg,xt='mean fiber',yt=r'$\Delta vmicro$',xr=[0,300],yr=[-1,1],zr=[0,5],zt='log g',colorbar=True,size=1)
#    fig.tight_layout()
#    fig.savefig(out+'_res.png')
#
    #fig,ax=plots.multi(1,2)
    #xr=[
    #junk = fit.fit1d(logg[gd], vmicro[gd],degree=degree,plot=ax[0,0],plot2d=True,yr=[0,5],xr=[0,5],xt='Teff',yt='log g')
    #y, x = np.mgrid[yr[1]:yr[0]:200j, xr[1]:xr[0]:200j]
 
    # 2D plots a f(teff, logg)
    fig,ax=plots.multi(1,5,figsize=(6,12))
    #fit2d = fit.fit2d(logg[gd], mh[gd], vmicro[gd],degree=degree,plot=ax[0],yr=[-2.5,1],xr=[0,5],xt='logg',yt='[M/H]',log=True)
    if func == vm1t :
      print('gd:',len(gd))
      x1=logg
      xr=[-0.5,5.5]
      xt='logg'
      x2=mh
      yr=[-2,0.5]
      yt='[M/H]'
      #x3=logg
      #x3_0=1.0
      #yr=[-0.5,5.5]
      #yt='log g'
      des=func(x1[gd],x2[gd])
      soln,inv=fit.linear(vmicro[gd],des)
      y, x = np.mgrid[yr[0]:yr[1]:200j, xr[0]:xr[1]:200j]
      ax[0].imshow(10.**func(x,y,soln=soln),
                extent=[xr[0],xr[1],yr[0],yr[1]],aspect='auto',vmin=0.,vmax=3., origin='lower')
      plots.plotc(ax[0],x1,x2,10.**vmicro,xr=xr,yr=yr,zr=[0.,3.],
                xt=xt,yt=yt,zt='vmicro',size=5,linewidth=0)
      plots.plotc(ax[0],x1[gd],x2[gd],10.**vmicro[gd],xr=xr,yr=yr,zr=[0.,3.],
                xt=xt,yt=yt,zt='vmicro',colorbar=True,size=15,linewidth=0.2)
      cs=ax[0].contour(x,y,10**func(x,y,soln=soln),colors='k')
      ax[0].clabel(cs)

    else :
      des=func(logg[gd],mh[gd])
      soln,inv=fit.linear(vmicro[gd],des)
      y, x = np.mgrid[-2.5:1.:200j, 0:5.:200j]
      ax[0].imshow(10.**func(x,y,soln=soln),
                extent=[0,5,-2.5,1.],aspect='auto',vmin=0.,vmax=3., origin='lower')
      plots.plotc(ax[0],logg[gd],mh[gd],10.**vmicro[gd],xr=[-0.5,5],yr=[-2.5,1],zr=[0.,3.],
                xt='log g',yt='[M/H]',zt='vmicro',colorbar=True,size=15,linewidth=1)
      cs=ax[0].contour(x,y,10**func(x,y,soln=soln),colors='k')
      ax[0].clabel(cs,color='k')
        # create independent variable grid for model and display
    plots.plotc(ax[1],logg[gd],10.**vmicro[gd],mh[gd],xt='logg',yt=r'$vmicro$',xr=[-0.5,5.],yr=[0,4],zr=[-2.,0.5],zt='[M/H]',size=5)
    plots.plotc(ax[1],logg,10.**vmicro,mh,xt='log g',yt=r'$vmicro$',xr=[-0.5,5.],yr=[0,4],zr=[-2.,0.5],zt='[M/H]',colorbar=True,size=1)
    plots.plotc(ax[2],logg[gd],10.**vmicro[gd],teff[gd],xt='logg',yt=r'$vmicro$',xr=[-0.5,5.],yr=[0,4],zr=[3000,5500],zt='Teff',size=5)
    plots.plotc(ax[2],logg,10.**vmicro,teff,xt='log g',yt=r'$vmicro$',xr=[-0.5,5.],yr=[0,4],zr=[3000,5500],zt='Teff',colorbar=True,size=1)
    plots.plotc(ax[3],logg[gd],10.**vmicro[gd]-10.**func(x1[gd],x2[gd],soln=soln),teff[gd],xt='logg',yt=r'$\Delta vmicro$',
                xr=[-0.5,5.],yr=[-1,1],zr=[3000,5500],zt='Teff',size=5)
    plots.plotc(ax[3],logg,10.**vmicro-10.**func(x1,x2,soln=soln),teff,xt='log g',yt=r'$\Delta vmicro$',
                xr=[-0.5,5.],yr=[-1,1],zr=[3000,5500],zt='Teff',colorbar=True,size=1)
    plots.plotc(ax[4],logg[gd],10.**vmicro[gd]-10.**func(x1[gd],x2[gd],soln=soln),mh[gd],xt='logg',yt=r'$\Delta vmicro$',
                xr=[-0.5,5.],yr=[-1,1],zr=[-2.0,0.5],zt='[M/H]',size=5)
    plots.plotc(ax[4],logg,10.**vmicro-10.**func(x1,x2,soln=soln),mh,xt='log g',yt=r'$\Delta vmicro$',
                xr=[-0.5,5.],yr=[-1,1],zr=[-2.0,0.5],zt='[M/H]',colorbar=True,size=1)

    fig.tight_layout()
    fig.savefig(out+'_res.png')

    #summary plots
    #plot(teff, logg, mh, meanfib, vmicro, vrange, fit1d, fit2d, vt='vmicro')

    print('{',end="")
    for i in range(len(soln)) :
        print('{:12.8f}'.format(soln[i]),end="")
    print('}')

#    return fit1d
    return

def vm3_1x(logg,mh,soln=None) :
    des=np.zeros([8,len(logg.flatten())])
    des[0,:]=1.
    des[1,:]=logg.flatten()
    des[2,:]=logg.flatten()**2
    des[3,:]=logg.flatten()**3
    des[4,:]=mh.flatten()
    des[5,:]=logg.flatten()*mh.flatten()
    des[6,:]=logg.flatten()**2*mh.flatten()
    des[7,:]=logg.flatten()**3*mh.flatten()
    if soln is not None:
        return np.dot(soln,des).reshape(logg.shape)
    else :
        return des

def vm3_1(logg,mh,soln=None) :
    des=np.zeros([5,len(logg.flatten())])
    des[0,:]=1.
    des[1,:]=logg.flatten()
    des[2,:]=logg.flatten()**2
    des[3,:]=logg.flatten()**3
    des[4,:]=mh.flatten()
    if soln is not None:
        return np.dot(soln,des).reshape(logg.shape)
    else :
        return des

def vm3(logg,mh,soln=None) :
    des=np.zeros([4,len(logg.flatten())])
    des[0,:]=1.
    des[1,:]=logg.flatten()
    des[2,:]=logg.flatten()**2
    des[3,:]=logg.flatten()**3
    if soln is not None:
        return np.dot(soln,des).reshape(logg.shape)
    else :
        return des

def vm1t(x1,x2,x3=None,soln=None) :
    if x3 is None : des=np.zeros([3,len(x1.flatten())])
    else : des=np.zeros([4,len(x1.flatten())])
    des[0,:]=1.
    des[1,:]=x1.flatten()
    des[2,:]=x2.flatten()
    if x3 is not None : des[3,:]=x3.flatten()
    if soln is not None:
        return np.dot(soln,des).reshape(x1.shape)
    else :
        return des

def fit_vmacro(data,teffrange=[3550,5500],mhrange=[-2.5,1.], loggrange=[-1,3.8], vrange=[1,15],degree=1,maxerr=0.1,apokasc='APOKASC_cat_v3.6.0.fits',nopersist=False, reject=0,out='vmacro') :
    """ 
    Fit macroturbulence relation  with 1D f(log g) and 2D f(Teff, logg) fits, plots

    Args:
        data : calib structure

    Keyword args :
        degree : degree of fit (default=1)
        mhrange : range of [M/H] (default=[1.,1])
        loggrange : range of log g (default=[1.,3.8])
        vrange  : scaling range for vmacro plot, vrange[1] sets maximum good vmacro

    Returns:
        fit1d, fit2d : 1D and 2D polynomial fits
    """

    teff = data['FPARAM'][:,0]
    logg = data['FPARAM'][:,1]
    vmicro = data['FPARAM'][:,2]
    mh = data['FPARAM'][:,3]
    vmacro = data['FPARAM'][:,7]
    verr = np.sqrt(data['FPARAM_COV'][:,7,7])
    try :
       meanfib = data['MEANFIB']
    except :
       meanfib = 0.*vmicro
    print('mhrange: ', mhrange)
    print('loggrange: ', loggrange)
    print('vrange: ', vrange)
    print('maxerr: ', maxerr)
    print('nopersist: ', nopersist)
    print('reject: ', reject)
    gd = np.where((mh>mhrange[0]) & (mh<mhrange[1]) & (logg > loggrange[0]) & (logg < loggrange[1]) & (10.**vmacro < vrange[1]) &
                  (teff>teffrange[0]) & (teff<teffrange[1]) & (verr<maxerr) )[0]
    if nopersist :
        j = np.where ((data['STARFLAG'][gd] & bitmask.persist()) == 0)[0]
        gd=gd[j]

    print('len(gd)', len(gd))

    # remove non-1st generation GC stars
    gcstars = ascii.read(os.environ['IDLWRAP_DIR']+'/data/gc_szabolcs.dat')
    bd=np.where(gcstars['pop'] != 1)[0]
    gd = [x for x in gd if data[x]['APOGEE_ID'] not in gcstars['id'][bd]]

    apokasc=fits.open(os.environ['IDLWRAP_DIR']+'/data/'+apokasc)[1].data
    i1,i2=match.match(data['APOGEE_ID'],apokasc['2MASS_ID'])
    rgb=np.where(apokasc['CONS_EVSTATES'][i2] == 'RGB')[0]
    rc=np.where(apokasc['CONS_EVSTATES'][i2] == 'RC')[0]
    type=mh*0.
    type[i1[rgb]]=1
    type[i1[rc]]=-1

    #fig,ax=plots.multi(2,2,figsize=(12,8))
    #fit1d = fit.fit1d(logg[gd], vmacro[gd],degree=degree,plot=ax[0,0],xt='log g', yt='vmacro',ydata=mh[gd],colorbar=True,zt='[M/H]',reject=reject,log=True)
    #fit1d = fit.fit1d(logg[gd], vmacro[gd],degree=degree,plot=ax[0,1],xt='log g', yt='vmacro',ydata=teff[gd],colorbar=True,zt='Teff',reject=reject,log=True)
    #fit1d = fit.fit1d(logg[gd], vmacro[gd],degree=degree,plot=ax[1,0],xt='log g', yt='vmacro',ydata=meanfib[gd],colorbar=True,zt='mean fib',reject=reject,log=True)
    #fit1d = fit.fit1d(logg[gd], vmacro[gd],degree=degree,plot=ax[1,1],xt='log g', yt='vmacro',ydata=type[gd],colorbar=True,zt='RGB/RC',reject=reject,log=True)
    #fig.tight_layout()
    #fig.savefig(out+'.png')

    #fig,ax=plots.multi(1,2)
    fig,ax=plots.multi(2,3,figsize=(12,8))
    print('2D logg, mh')
    fit2d = fit.fit2d(logg[gd], mh[gd], vmacro[gd],degree=degree,plot=ax[0,0],xt='log g', yt='[M/H]',zt='log(vmacro)',reject=reject,zr=[0,15],log=True)
    plots.plotc(ax[0,1],teff[gd],10.**vmacro[gd]-10.**fit2d(logg[gd],mh[gd]),mh[gd],xt='Teff',yt=r'$\Delta vmacro$',zt='[M/H]',xr=[5500,3500],yr=[-10,20],colorbar=True,yerr=verr[gd]*vmacro[gd]*math.log(10))
    parprint([fit2d.parameters[0],0.,fit2d.parameters[1],fit2d.parameters[2]])

    #print('2D teff, logg')
    fit2d_b = fit.fit2d(teff[gd], logg[gd], vmacro[gd],degree=degree,plot=ax[1,0],xt='Teff', yt='log g',zt='log(vmacro)',reject=reject,xr=[6000,3000],yr=[5,-1],zr=[0,15],log=True)
    plots.plotc(ax[1,1],teff[gd],10.**vmacro[gd]-10.**fit2d_b(teff[gd],logg[gd]),mh[gd],xt='Teff',yt=r'$\Delta vmacro$',zt='[M/H]',xr=[5500,3500],yr=[-10,20],colorbar=True,yerr=verr[gd]*vmacro[gd]*math.log(10))
    parprint([fit2d_b.parameters[0],fit2d_b.parameters[1],fit2d_b.parameters[2],0.])
    #print('2D teff, [M/H]')
    #fit2d_c = fit.fit2d(teff[gd], mh[gd], vmacro[gd],degree=degree,plot=ax[0,1],xt='Teff', yt='[M/H]',zt='log(vmacro)',reject=reject,zr=[0,15],log=True)
    #parprint([fit2d_c.parameters[0],fit2d_c.parameters[1],0.,fit2d_c.parameters[2]])

    fit1d = fit.fit1d(mh[gd], vmacro[gd],degree=degree,plot=ax[2,0],xt='[M/H]', yt='vmacro',colorbar=True,reject=reject,log=True,ydata=logg[gd],yr=[0,5],zt='log g')
    parprint([fit1d.parameters[0],0.,0.,fit1d.parameters[1]])
    plots.plotc(ax[2,1],teff[gd],10.**vmacro[gd]-10.**fit1d(mh[gd]),mh[gd],xt='Teff',yt=r'$\Delta vmacro$',zt='[M/H]',xr=[5500,3500],yr=[-10,20],colorbar=True,yerr=verr[gd]*vmacro[gd]*math.log(10),zr=[-2,0.5])

    #plots.plotc(ax[1,1],teff[gd],logg[gd],10.**vmacro[gd],xr=[6000,3000],yr=[5,-1],xt='Teff',yt='log g', zr=[0,15])

    # DR13 fit 
    #dr13fit=models.Polynomial2D(degree=1)
    #dr13fit.parameters=[0.741,-0.0998,-0.225]
    #junk = fit.fit2d(logg[gd], mh[gd], vmacro[gd],degree=degree,plot=ax[1],xt='log g', yt='[M/H]',zt='log(vmacro)',pfit=dr13fit)
    fig.tight_layout()
    fig.savefig(out+'_2d.png')

    # single plot with 1D fit
    fig,ax=plots.multi(1,1,figsize=(12,4))
    fit1d = fit.fit1d(mh[gd], vmacro[gd],degree=degree,plot=ax,xt='[M/H]', yt='vmacro',colorbar=True,reject=reject,log=True,ydata=logg[gd],yr=[0,5],zt='log g')
    fig.tight_layout()
    fig.savefig(out+'_1d.pdf')

    #pdb.set_trace()
    #plot(teff, logg, mh, meanfib, vmacro, vrange, fit1d, fit2d, vt='vmacro')
    return
    #return fit1d, fit2d

def parprint(par) :

    print('{',end="")
    print('{:10.6f}'.format(par[0]),end="")
    print('{:10.6f}'.format(par[1]),end="")
    print('{:10.6f}'.format(par[2]),end="")
    print('{:10.6f}'.format(par[3]),end="")
    print('}')


def litplot(mass=1, feh=0, out='litvmacro') : 
    '''
    Plots literature vmacro relations
    '''
#    parser = argparse.ArgumentParser()
#    parser.add_argument("--mass",default=1.)
#    parser.add_argument("--feh",default=0.)
#    args=parser.parse_args()
#    print('Using mass: ', args.mass)
#    mass=float(args.mass)
#    feh=float(args.feh)
#    print('Using feh: ', feh)
#    mass = 0.9

    # setup plots
    mpl.rcParams['xtick.labelsize'] = 8
    mpl.rcParams['ytick.labelsize'] = 8

    p,ax=plots.multi(2,2)
    lw=0
    ax[0,0].set_title('Massarotti vmacro')
    ax[1,0].set_title('Shetrone vmacro')
    ax[0,1].set_title('Bergemann vmacro')
    ax[1,1].set_title('Bergemann vmacro new')

    # L = 4*pi*R^2*sigma*Teff^4
    # R^2 = L / (4 pi sigma Teff^4)
    # log R^2 = log L - log (4*pi*sigma) - 4 logte
    # logg = log G + log M - log L + log(4*pi*sigma) + 4 log te
    sigma=const.sigma_sb.cgs.value
    G=const.G.cgs.value
    lsun=const.L_sun.cgs.value
    msun=const.M_sun.cgs.value
    m=mass*msun
    # get isochrone data
    for iso in [0.,-1.,-2.] :
        if iso >=0 : 
            c='p' 
        else: 
            c='m'
        cfeh='z'+c+'{:02d}'.format(int(round(abs(iso)*10)))
        print(cfeh)
        a=ascii.read(os.getenv('ISOCHRONE_DIR')+'/'+cfeh+'.dat')
        age=a['col2']
        j=np.where((age == 9.6) | (age == 9.0) | (age == 10.0))
        logl=a['col5'][j]
        iso_te=10.**a['col6'][j]
        iso_logg=a['col7'][j]
        iso_feh=iso_logg*0.+iso
        iso_logl=4*np.log10(iso_te)-iso_logg-np.log10(lsun)+np.log10(4*np.pi*sigma*G)+np.log10(m)

        vm=vmacro_massarotti(iso_te,iso_logl,iso_feh)
        plots.plotc(ax[0,0],iso_te,iso_logg,vm,xr=[6000,3000],yr=[5,-1.],zr=[0,15],linewidth=lw,size=15)
        vm=vmacro_shetrone(iso_te,iso_logl,iso_feh)
        plots.plotc(ax[1,0],iso_te,iso_logg,vm,xr=[6000,3000],yr=[5,-1.],zr=[0,15],linewidth=lw,size=15)
        vm=vmacro_bergemann(iso_te,iso_logg,iso_feh)
        plots.plotc(ax[0,1],iso_te,iso_logg,vm,xr=[6000,3000],yr=[5,-1.],zr=[0,15],linewidth=lw,size=15)
        vm=vmacro_bergemann_new(iso_te,iso_logg,iso_feh)
        plots.plotc(ax[1,1],iso_te,iso_logg,vm,xr=[6000,3000],yr=[5,-1.],zr=[0,15],linewidth=lw,size=15)


    # setup grids for Teff and logg
    o=np.ones([61,31])
    teff=o*(np.arange(31)*100.)+3000
    logg=np.transpose(np.transpose(o)*(np.arange(61)*.1))-1.
    feh=o*0.+feh
    logl=4*np.log10(teff)-logg-np.log10(lsun)+np.log10(4*np.pi*sigma*G)+np.log10(m)

    a=vmacro_massarotti(teff,logl,feh)
    #axim=ax[0,0].imshow(np.fliplr(a),vmin=0,vmax=10,extent=[7500,2500.,5.,-1.],aspect='auto',origin='upper')
    #p.colorbar(axim)

    a=vmacro_shetrone(teff,logl,feh)
    #axim=ax[1,0].imshow(np.fliplr(a),vmin=0,vmax=10,extent=[7500,2500.,5.,-1.],aspect='auto',origin='upper')
    #p.colorbar(axim)
    
    a=vmacro_bergemann(teff,logg,feh)
    #axim=ax[0,1].imshow(np.fliplr(a),vmin=0,vmax=10,extent=[7500.,2500.,5.,-1.],aspect='auto',origin='upper')
    #p.colorbar(axim)

    a=vmacro_bergemann_new(teff,logg,feh)
    #axim=ax[1,1].imshow(np.fliplr(a),vmin=0,vmax=10,extent=[7500.,2500.,5.,-1.],aspect='auto',origin='upper')
    #p.colorbar(axim)

    p.tight_layout()
    p.savefig(out+'_'+str(mass)+'.png')

def vmacro_massarotti(teff,logl,feh) :
    '''
    Returns Massarotti vmacro as f(teff,logl,feh)
   
    Args:
        teff : input Teff
        logg : input logg
        feh : input feh

    Returns:
        vmacro
    '''
    logte=np.log10(teff)
    vm=10**(3.5*logte+0.25*logl-12.67)
    return vm

def vmacro_shetrone(teff,logl,feh) :
    '''
    Returns Shetrone vmacro as f(teff, logl, feh)
   
    Args:
        teff : input Teff
        logg : input logg
        feh : input feh

    Returns:
        vmacro
    '''
    logte=np.log10(teff)
    return 10**(1.85*logte+0.22*logl-6.49)

def vmacro_bergemann(teff,logg,feh) :
    '''
    Returns Bergemann vmacro as f(teff, logl, feh)
   
    Args:
        teff : input Teff
        logg : input logg
        feh : input feh

    Returns:
        vmacro
    '''
    tref=5250. ; to = 5500. ; go = 4.0
    t=teff ; g=logg

    # set up output, which we will fill in pieces for different regimes
    vt=teff*0.

    # t >= Tref and logg >=3.5
    j= np.where((teff>tref) & (logg >= 3.5))
    a1 =  1.15 ; b1 =  7.e-4; c1 = 1.2e-6
    b2 = -0.13 ; c2 =  0.13
    b3 = -0.37 ; c3 = -0.07
    vt[j] = 3*(a1 + b1*(t[j]-to) + c1*(t[j]-to)**2 + b2*(g[j]-go) + c2*(g[j]-go)**2 + b3*feh[j] + c3*feh[j]**2)

    # t < Tref and logg >=3.5
    j= np.where((teff<tref) & (logg >= 3.5))
    a1 =  1.15 ; b1 =  2.e-4; c1 = 3.95e-7
    b2 = -0.13 ; c2 =  0.13
    b3 =  0.0  ; c3 =  0.0
    vt[j] = 3*(a1 + b1*(t[j]-to) + c1*(t[j]-to)**2 + b2*(g[j]-go) + c2*(g[j]-go)**2 + b3*feh[j] + c3*feh[j]**2)

    # logg < 3.5
    j= np.where(logg < 3.5)
    a1 =  1.15 ; b1 =  2.2e-5 ; c1 = -0.5e-7
    b2 = -0.1 ; c2 =  0.04 
    b3 = -0.37  ; c3 = -0.07
    vt[j] = 3*(a1 + b1*(t[j]-to) + c1*(t[j]-to)**2 + b2*(g[j]-go) + c2*(g[j]-go)**2 + b3*feh[j] + c3*feh[j]**2)
    return vt

def vmacro_bergemann_new(teff,logg,feh) :
    '''
    Returns newer Bergemann vmacro as f(teff, logl, feh)
   
    Args:
        teff : input Teff
        logg : input logg
        feh : input feh

    Returns:
        vmacro
    '''

    t0  = 5500
    g0  = 4.0
    vt=teff*0.

    # MS and RGB, Teff >=5250
    j= np.where((teff>5250.) & (logg >= 3.5))
    Rxa = ( 1.15 + 7e-4*(teff[j]-t0) + 1.2e-6*(teff[j]-t0)**2 
            - 0.13*(logg[j]-g0) +     0.13*(logg[j]-g0)**2 
            -    0.37*feh[j]     -    0.07*feh[j]**2 )   
    vt[j] = 3*Rxa

    # MS, Teff <=5250
    j= np.where((teff<5250.) & (logg >= 3.5))
    Rxa = (1.15 + 2e-4*(teff[j]-t0) + 3.95e-7*(teff[j]-t0)**2 
            - 0.13*(logg[j]-g0) +     0.13*(logg[j]-g0)**2 )
    vt[j] = 3*Rxa

    # RGB/AGB
    j= np.where(logg < 3.5)
    Rxa = ( 1.15 + 2.2e-5*(teff[j]-t0) - 0.5e-7*(teff[j]-t0)**2 
            -    0.1*(logg[j]-g0) +   0.04*(logg[j]-g0)**2  
            -    0.37*feh[j]     -    0.07*feh[j]**2 )
    vt[j] = 3*Rxa
    return vt

def vmicro_bergemann(teff,logg,feh) :
    '''
    Returns Bergemann vmicro as f(teff, logl, feh)
   
    Args:
        teff : input Teff
        logg : input logg
        feh : input feh

    Returns:
        vmicro
    '''

    # reference variables
    tref=5500. ; to = 5500. ; go = 4.0
    t=teff ; g=logg

    # set up output variable, which we will fill in pieces
    vt=teff*0.

    # t >= Tref and logg >=3.5
    j= np.where((teff>=tref) & (logg >= 3.5))
    a1 =  1.1 ; b1 =  2.e-4; c1 = 7.95e-7
    b2 = -0.13 ; c2 =  0.13
    b3 = -0.27 ; c3 = -0.1
    vt[j] = (a1 + b1*(t[j]-to) + c1*(t[j]-to)**2 + b2*(g[j]-go) + c2*(g[j]-go)**2 + b3*feh[j] + c3*feh[j]**2)

    # t < Tref and logg >=3.5
    j= np.where((teff<tref) & (logg >= 3.5))
    a1 =  1.1 ; b1 =  2.e-4; c1 = 3.95e-7
    b2 = -0.13 ; c2 =  0.13
    b3 =  -0.6  ; c3 =  -0.2
    vt[j] = (a1 + b1*(t[j]-to) + c1*(t[j]-to)**2 + b2*(g[j]-go) + c2*(g[j]-go)**2 + b3*feh[j] + c3*feh[j]**2)

    # logg < 3.5
    j= np.where(logg < 3.5)
    a1 =  1.3 ; b1 =  2.2e-5 ; c1 = -0.5e-7
    b2 = -0.1 ; c2 =  0.04 
    b3 =  0.0  ; c3 = 0.0
    vt[j] = (a1 + b1*(t[j]-to) + c1*(t[j]-to)**2 + b2*(g[j]-go) + c2*(g[j]-go)**2 + b3*feh[j] + c3*feh[j]**2)
    return vt

def plot(teff, logg, mh, meanfib, v, vr, fit1d, fit2d, 
         xr=[5500,3500], yr=[5,0], vt='vmicro') :
    '''
    auxiliary plotting routine for vmicro, vmacro plots
    '''
    fig = plt.figure(figsize=(14,8))
    y, x = np.mgrid[yr[1]:yr[0]:200j, xr[1]:xr[0]:200j]
    ax=fig.add_subplot(3,2,1)
    plots.plotc(ax,teff,logg,10.**v,xr=xr,yr=yr,zr=vr,
                xt='Teff',yt='log g',zt=vt,colorbar=True,size=15)
    ax.imshow(10.**fit2d(x,y),extent=[xr[1],xr[0],yr[1],yr[0]],
                aspect='auto',vmin=0,vmax=vr[1], origin='lower')

    ax=fig.add_subplot(3,2,2)
    plots.plotc(ax,teff,logg,10.**v,xr=xr,yr=yr,zr=vr,
                xt='Teff',yt='log g',zt=vt,colorbar=True,size=15)
    ax.imshow(10.**fit1d(y),extent=[xr[1],xr[0],yr[1],yr[0]],
              aspect='auto',vmin=0,vmax=vr[1], origin='lower')

    ax=fig.add_subplot(3,2,3)
    plots.plotc(ax, logg, 10.**fit2d(teff,logg), mh, xr=[0,5], yr=vr, zr=[-2.5,0.5],
                size=10,colorbar=True,xt='log g', yt=vt, zt='[M/H]')

    ax=fig.add_subplot(3,2,4)
    plots.plotc(ax, logg, 10.**fit1d(logg), mh, xr=[0,5], yr=vr, zr=[-2.5,0.5],
                size=10,colorbar=True,xt='log g', yt=vt, zt='[M/H]')
    x=np.arange(0,5,0.01)
    plots.plotl(ax,x, 10.**(0.226 - 0.0228 *x + 0.0297 *x**2 - 0.0113 *x**3))

    ax=fig.add_subplot(3,2,5)
    plots.plotc(ax, logg, 10.**v, mh, xr=[0,5], yr=vr, zr=[-2.5,0.5],
                size=10,colorbar=True,xt='log g', yt=vt, zt='[M/H]')

    ax=fig.add_subplot(3,2,6)
    plots.plotc(ax, logg, 10.**v, meanfib, xr=[0,5], yr=vr, zr=[0,300],
                size=10,colorbar=True,xt='log g', yt=vt, zt='<fiber>')

    plt.show()

def dr14_vmacro(teff, logg, mh, meanfib, v, vr, fit1d, fit2d, 
         xr=[5500,3500], yr=[5,0], vt='vmicro') :
    ax=fig.add_subplot(3,2,5)
    plots.plotc(ax, logg, 10.**v, mh, xr=[0,5], yr=vr, zr=[-2.5,0.5],
                size=10,colorbar=True,xt='log g', yt=vt, zt='[M/H]')
