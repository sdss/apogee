#!/usr/bin/env python

from tools import plots
from tools import html
from astropy.io import fits
from astropy.io import ascii
import numpy as np
import math
import pdb
import argparse
import os
import matplotlib.pyplot as plt

def throughplot(instrument='apogee-s',outfile=None,inter=False) :
    '''
    Routine to make zeropoint/throughput plots from apogeeSci summary files
    with information including FWHM, GDRMS, CART
    '''

    # instrument specific
    if instrument == 'apogee-s' :
        gain=3.
        carts=[20,25]
        fiber_rad=0.65
        telescope='lco25m'
    else :
        gain=1.9
        carts=[0,10]
        fiber_rad=1.
        telescope='apo25m'

    # read summary data made by mkmonitor
    a=fits.open(instrument+'Sci.fits')[1].data
    gd = np.where(a['NREADS'] >= 47)[0]
    a=a[gd]

    # use weather information if we can
    clouds=np.zeros(len(a)).astype(int)
    nmiss=0
    nhave=0
    try :
        c=ascii.read(os.environ['APOGEEREDUCEPLAN_DIR']+'/data/'+telescope+'/clouds.txt')
        try:
            for i,p in enumerate(a['PLATE']) :
                j=np.where((c['plate'] == p) & (c['MJD'] == a['MJD'][i]) )[0]
                if len(j)>0 : 
                    if len(j)>1 : 
                        print('double cloud match',p,a['MJD'][i])
                        pdb.set_trace()
                    clouds[i] = c['clouds_level'][j[0]]
                    nhave+=1
                else : 
                    nmiss+=1
                    print('no clouds match found for',a['MJD'][i],p)
        except :
            print('error!',i,p,j)
            pdb.set_trace()

        gd=np.where(clouds <= 1)[0]
        a=a[gd]
    except :
        print('cant open clouds file')


    # seeing correction factor
    sigma = a['FWHM']/2.354
    sigma = a['SEEING']/2.354
    ee = 1. - np.exp(-(fiber_rad**2)/(2*sigma**2))
    corr = a['ZERONORM']-2.5*np.log10(ee)

    gd = np.where(np.isfinite(corr))[0]
    a=a[gd]
    ee=ee[gd]
    corr=corr[gd]

    # run number for LCO
    run = ((a['MJD']-57850)/29.+0.5).astype(int)

    # rough throughput calculation
    h=6.63e-27
    c=3.e10
    lam=1.6e-4
    dlam=0.3
    dt=10.6
    area=math.pi*(125.**2-50.**2)
    fvega=11.38e-11
    through=10**(0.4*a['ZERONORM'])*h*c/lam/dlam/dt*gain/area/fvega/ee
   
    # straight DHA
    dha=a['HA']-a['DESIGN_HA'][:,0]
    #dha=np.abs(dha)
    # "normalized" DHA
    j=np.where(a['HA']<a['DESIGN_HA'][:,0])[0]
    dha[j]/=(a['DESIGN_HA'][j,0]-a['DESIGN_HA'][j,1])
    j=np.where(a['HA']>=a['DESIGN_HA'][:,0])[0]
    dha[j]/=(a['DESIGN_HA'][j,2]-a['DESIGN_HA'][j,0])

    #plots with MJD
    files=[]

    out='monitor/'+instrument+'/'+instrument

    # point size by FWHM
    psize=a['FWHM']/1.*40
    j=np.where(psize == 0.)[0]
    psize[j] = 10

    # histograms by run
    fig,ax=plots.multi(2,3,figsize=(8,12))
    file=out+'zero_hist.png'
    runs=list(set(run))
    runs.append(999)
    for r in runs :
        gd = np.where(run == r)[0]
        if r == 999 :
             gd = np.where(run < 999)[0]
        if r >= 8 : lw=2
        else : lw=1
        print(r,len(gd))
        try:
            n,b,p=plt.hist(a['GDRMS'][gd],histtype='step',bins=np.arange(0,1,0.05),label='{:3d}'.format(r),linewidth=lw,normed=False)
            if r == 999 : n/=2
            ax[0,0].plot(b[0:-1]+(b[1]-b[0])/2.,n,linewidth=lw,label='{:2d}'.format(r))
            ax[0,0].set_xlabel('GDRMS')
        except : pass
        try:
            n,b,p=plt.hist(a['ZERONORM'][gd],histtype='step',bins=np.arange(12,15.5,0.1),linewidth=lw,normed=False,label='{:2d}'.format(r))
            if r == 999 : n/=2
            ax[0,1].plot(b[0:-1]+(b[1]-b[0])/2.,n,linewidth=lw,label='{:2d}'.format(r))
            ax[0,1].set_xlabel('ZERONORM')
            n,b,p=plt.hist(corr[gd],histtype='step',bins=np.arange(12,16,0.1),linewidth=lw,normed=False,label='{:3d}'.format(r))
            if r == 999 : n/=2
            ax[1,0].plot(b[0:-1]+(b[1]-b[0])/2.,n,linewidth=lw,label='{:3d}'.format(r))
            ax[1,0].set_xlabel('ZERONORM (adjusted)')
            n,b,p=plt.hist(a['ZERORMS'][gd],histtype='step',bins=np.arange(0,1,0.05),linewidth=lw,normed=False,label='{:3d}'.format(r))
            if r == 999 : n/=2
            ax[1,1].plot(b[0:-1]+(b[1]-b[0])/2.,n,linewidth=lw,label='{:3d}'.format(r))
            ax[1,1].set_xlabel('ZERORMS')
            n,b,p=plt.hist(through[gd],histtype='step',bins=np.arange(0,0.34,0.02),linewidth=lw,normed=False,label='{:3d}'.format(r))
            if r == 999 : n/=2
            ax[2,0].plot(b[0:-1]+(b[1]-b[0])/2.,n,linewidth=lw,label='{:3d}'.format(r))
            ax[2,0].set_xlabel('THROUGHPUT (adjusted)')
        except : pass
    
    if instrument == 'apogee-s' :
        ax[0,0].legend(fontsize=6,loc=1,title='Run')
        ax[0,1].legend(fontsize=6,loc=2,title='Run')
        ax[1,0].legend(fontsize=6,loc=2,title='Run')
        ax[1,1].legend(fontsize=6,loc=1,title='Run')
    ax[2,1].remove()
    fig.tight_layout()
    fig.savefig(file)
    files.append([os.path.basename(file)])

    ctype = [a['FWHM'],a['SEEING'],a['GDRMS'],dha,a['CART']]
    name = ['zero_fwhm','zero_seeing','zero_gdrms','zero_dha','zero_cart']
    zr=[[0.5,2.],[0.5,2.],[0,0.8],[-2,2],carts]
    zt=['FWHM','SEEING','GDRMS','DHA','CART']
    for j,c in enumerate(ctype) :
      fig,ax=plots.multi(1,4,hspace=0.001,sharex=True,figsize=(24,6))
      file=out+name[j]+'.png'
      plots.plotc(ax[0],a['MJD'],a['ZERONORM'],c,yr=[12,15.5],zr=zr[j],size=psize,colorbar=True,xt='MJD',yt='ZERONORM',zt=zt[j])
      plots.plotc(ax[1],a['MJD'],corr,c,yr=[12,15.5],zr=zr[j],size=psize,colorbar=True,xt='MJD',yt='ZERONORM (adjusted)',zt=zt[j])
      plots.plotc(ax[2],a['MJD'],a['ZERORMS'],c,yr=[0,1],zr=zr[j],size=psize,colorbar=True,xt='MJD',yt='ZERORMS',zt=zt[j])
      plots.plotc(ax[3],a['MJD'],through,c,yr=[0,0.3],zr=zr[j],size=psize,colorbar=True,xt='MJD',yt='throughput',zt=zt[j])
      fig.savefig(file)
      files.append([os.path.basename(file)])

    fig,ax=plots.multi(1,1)
    plots.plotc(ax,a['SEEING'],a['ZERONORM'],a['GDRMS'],xr=[0.,3.0],yr=[13,15.5],zr=[0.2,1.2],xt='Seeing',yt='ZERONORM',zt='GDRMS',colorbar=True,size=1)
    #plots.plotc(ax[1],a['SEEING'],corr,a['GDRMS'],xr=[0.,3.0],yr=[13,15.5],zr=[0.2,1.2],xt='Seeing',yt='seeing-corrected ZERONORM',zt='GDRMS',colorbar=True,size=1)
    file=out+'_seeing.png'
    fig.savefig(file)
    files.append([os.path.basename(file)])

    out='monitor/'+instrument+'/'+instrument
    html.htmltab(files,file=out+'zero.html')
    if inter :
        pdb.set_trace()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make throughput plots",
                                     usage="through --instrument apogee-s")
    parser.add_argument("-i", "--instrument", type=str,
                        required=True,
                        help="instrument to plot",
                        choices=['apogee-s', 'apogee-n'])
    args = parser.parse_args()
    throughplot(instrument=args.instrument)

