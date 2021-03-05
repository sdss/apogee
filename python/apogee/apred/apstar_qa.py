# routines for assessing RVs from pipeline

import os
import copy
import glob
import pdb
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import esutil
import pickle
import yaml
from astropy.io import fits
from apogee.utils import apload, applot, bitmask, spectra
from apogee.aspcap import norm
from apogee.speclib import lsf
from tools import plots, html, match, struct
from sdss import yanny
from scipy import interpolate
from scipy.signal import correlate
from scipy.signal.windows import tukey
import scipy.linalg
from scipy.ndimage.filters import median_filter, gaussian_filter

import doppler 
import multiprocessing as mp
from astropy.table import Table, Column, vstack
from apogee.apred import bc, target

colors=['r','g','b','c','m','y','k']
chips=['a','b','c']

def vscat(a,fig=None,ls=None,marker='o',nmin=2,mhmin=-3,density=False,out=None,teff='RV_TEFF') :
    """ Make histograms of VSCATTER for different bins of Teff H], given min NVISITS, and min [M/H]
    """
    if fig == None : fig,ax=plots.multi(4,6,hspace=0.001,wspace=0.4,figsize=(12,8))
    else : fig,ax=fig
    tbins=[3000,3500,4000,4500,5500,8000,30000] 
    hbins=[8,11,12,13,15]
    try: snr = a['SNREV']
    except: snr=a['SNR']
    j=np.where(snr > 300) [0]
    snr[j] = 300
    for i in range(len(tbins)-1) :
        ax[i,0].text(0.9,0.9,'{:d}<={:s}<{:d}'.format(tbins[i],teff, tbins[i+1]),ha='right',transform=ax[i,0].transAxes,fontsize=8)
        for j in range(len(hbins)-1) :
            ax[0,j].set_title('{:d}<=H<{:d}'.format(hbins[j],hbins[j+1]))
            gd = np.where((a[teff]>=tbins[i]) & (a[teff]<tbins[i+1]) &
                          (a['H']>=hbins[j]) & (a['H']<hbins[j+1]) &
                           (a['NVISITS']>nmin) & (a['RV_FEH']>mhmin) & (a['VSCATTER'] > 0)) [0]
            print(tbins[i],tbins[i+1],hbins[j],hbins[j+1],nmin,len(gd))
            try :
                #plots.plotc(ax[i,2],snr[gd],a['VSCATTER'][gd],a['RV_FEH'][gd],marker=marker,xr=[0,310],yr=[0,1],xt='S/N',yt='VSCATTER')
                ax[i,j].hist(a['VSCATTER'][gd],bins=np.arange(0,1,0.01),ls=ls,histtype='step',color=colors[j],normed=density)
                ax[i,j].set_xlabel('VSCATTER (km/s)')
                ax[i,j].plot([0.1,0.1],ax[i,j].get_ylim())
                #ax[i,1].hist(a['VSCATTER'][gd],bins=np.arange(0,1,0.01),histtype='step',cumulative=True,normed=True,ls=ls,color=colors[j])
                #ax[i,1].set_xlabel('VSCATTER')
            except : pass

    if out is not None : 
        fig.savefig(out+'.png')
        plt.close()

    fig.suptitle('NVISITS>{:d} [M/H]>{:6.2f}'.format(nmin,mhmin))
    return fig,ax

def apolco(a,minfeh=-3,out=None) :
    """  VSCATTER histograms for APO vs LCO
    """
    apo=np.where((a['TELESCOPE'] == 'apo25m') & (a['RV_FEH']>minfeh) )[0]
    fig=vscat(a[apo],marker='o',density=True)
    lco=np.where((a['TELESCOPE'] == 'lco25m') & (a['RV_FEH']>minfeh) )[0]
    vscat(a[lco],fig=fig,ls=':',marker='+',density=True)
    if out is not None : 
        fig[0].savefig(out+'_1.png')
        plt.close()
    i1,i2=match.match(a['APOGEE_ID'][apo],a['APOGEE_ID'][lco])
    print('matched {:d} stars'.format(len(i1)))
    fig,ax=plots.multi(1,2)
    #plots.plotp(ax[0,0],a['SNR'][apo[i1]],a['VHELIO_AVG'][apo[i1]]-a['VHELIO_AVG'][lco[i2]],yr=[-3,3],yt=r'$\Delta$ VHELIO_AVG',xt='S/N')
    #plots.plotp(ax[0,1],a['SNR'][apo[i1]],a['VHELIO_AVG'][apo[i1]]-a['VHELIO_AVG'][lco[i2]],yr=[-50,50],yt=r'$\Delta$ VHELIO_AVG',xt='S/N')
    #plots.plotp(ax[1,0],a['SNR'][apo[i1]],a['VSCATTER'][apo[i1]]-a['VSCATTER'][lco[i2]],yr=[-0.5,0.5],yt=r'$\Delta$ VSCATTER',xt='S/N')
    #plots.plotp(ax[1,1],a['SNR'][apo[i1]],a['VSCATTER'][apo[i1]]-a['VSCATTER'][lco[i2]],yr=[-5,5],yt=r'$\Delta$ VSCATTER',xt='S/N')
    ax[0].hist(a['VHELIO_AVG'][apo[i1]]-a['VHELIO_AVG'][lco[i2]],bins=np.arange(-0.5,0.5,0.02),histtype='step')
    ax[0].set_xlabel(r'$\Delta$ VHELIO_AVG')
    ax[1].hist(a['VSCATTER'][apo[i1]]-a['VSCATTER'][lco[i2]],bins=np.arange(-0.25,0.25,0.01),histtype='step')
    ax[1].set_xlabel(r'$\Delta$ VSCATTER')
    if out is not None : 
        fig.savefig(out+'_2.png')
        plt.close()


def comp(a,b,av=None,bv=None,domatch=True,out=None) :
    """ VSCATTER comparison of two different data sets
    """
    if domatch :
        i1,i2=match.match(a['APOGEE_ID'],b['APOGEE_ID'])
        gd = np.where(a['NVISITS'][i1] == b['NVISITS'][i2])[0]
        a=a[i1[gd]]
        b=b[i2[gd]]

    fig = vscat(a)
    vscat(b,fig=fig,ls=':')
    if out is not None : 
        fig[0].savefig(out+'_1.png')
        plt.close()

    if domatch :
        fig,ax=plots.multi(1,2)
        #plots.plotp(ax[0,0],a['SNR'],a['VHELIO_AVG']-b['VHELIO_AVG'],yr=[-3,3],yt=r'$\Delta$ VHELIO_AVG')
        #plots.plotp(ax[0,1],a['SNR'],a['VHELIO_AVG']-b['VHELIO_AVG'],yr=[-50,50],yt=r'$\Delta$ VHELIO_AVG')
        #plots.plotp(ax[1,0],a['SNR'],a['VSCATTER']-b['VSCATTER'],yr=[-0.5,0.5],yt=r'$\Delta$ VSCATTER')
        #plots.plotp(ax[1,1],a['SNR'],a['VSCATTER']-b['VSCATTER'],yr=[-5,5],yt=r'$\Delta$ VSCATTER')
        ax[0].hist(a['VHELIO_AVG']-b['VHELIO_AVG'],bins=np.arange(-0.5,0.5,0.02),histtype='step')
        ax[0].set_xlabel(r'$\Delta$ VHELIO_AVG')
        ax[1].hist(a['VSCATTER']-b['VSCATTER'],bins=np.arange(-0.5,0.5,0.02),histtype='step')
        ax[1].set_xlabel(r'$\Delta$ VSCATTER')
        if out is not None : 
            fig.savefig(out+'_2.png')
            plt.close()

    return a,b

def visitsum_tel(all) :

    j=np.where(all['TELESCOPE'] == 'apo25m')[0]
    apo= all[j]
    j=np.where(all['TELESCOPE'] == 'lco25m')[0]
    lco= all[j]
    apoobjs = np.array(list(set(apo['APOGEE_ID'])))
    lcoobjs = np.array(list(set(lco['APOGEE_ID'])))
    i1,i2=match.match(apoobjs,lcoobjs)
    vhelio = []
    vscat = []
    verr = []
    sigfiber = []
    vdiff = []
    n = []
    mjd = []
    tel = []
    for i in i1 :
        j=np.where(all['APOGEE_ID'] == apoobjs[i])[0]
        vhelio.append(all['VHELIO'][j].mean())
        vscat.append(all['VHELIO'][j].std())
        verr.append(all['VRELERR'][j].max())
        sigfiber.append(all['FIBERID'][j].std())
        vdiff.extend(all['VHELIO'][j]-all['VHELIO'][j].mean())
        mjd.extend(all['MJD'][j])
        tel.extend(all['TELESCOPE'][j])
        n.append(len(j))

    fig,ax=plots.multi(1,2)
    mjd=np.array(mjd)
    tel=np.array(tel)
    vdiff=np.array(vdiff)
    plots.plotp(ax[0],mjd,vdiff,typeref=tel,types=['apo25m','lco25m'],color=['b','g'],yr=[-1,1])
    j=np.where(tel == 'apo25m')[0]
    ax[1].hist(vdiff[j],color='b',bins=np.arange(-1,1,0.01),histtype='step')
    mjds= [55800, 56130, 56512, 56876, 57230, 57600, 57966, 58360]
    for i in range(len(mjds)-1) :
        j=np.where((tel == 'apo25m') & (mjd >mjds[i]) & (mjd<mjds[i+1]) )[0]
        print(mjds[i],len(j))
        ax[1].hist(vdiff[j],bins=np.arange(-1,1,0.03),histtype='step')

    j=np.where(tel == 'lco25m')[0]
    ax[1].hist(vdiff[j],color='g',bins=np.arange(-1,1,0.01),histtype='step')
    mjds= [57829, 57966, 58360]
    plt.show()
    pdb.set_trace()


def visitsum(all,out=None,minvisit=1) :
    objs = set(all['APOGEE_ID'])
    if out is None :
        vhelio = []
        vscat = []
        verr = []
        sigfiber = []
        vdiff = []
        n = []
        print('n objects: ', len(objs))
        for iobj,obj in enumerate(objs) :
            j = np.where(all['APOGEE_ID'] == obj)[0]
            print(iobj,len(j))
            vhelio.append(all['VHELIO'][j].mean())
            vscat.append(all['VHELIO'][j].std())
            verr.append(all['VRELERR'][j].max())
            sigfiber.append(all['FIBERID'][j].std())
            vdiff.extend(all['VHELIO'][j]-all['VHELIO'][j].mean())
            n.append(len(j))
        vhelio=np.array(vhelio)
        vscat=np.array(vscat)
        verr=np.array(verr)
        sigfiber=np.array(sigfiber)
        vdiff=np.array(vdiff)
        n=np.array(n)
    else :
        vhelio,vscat,verr,sigfiber,vdiff,n = out

    vdiff=np.array(vdiff)
    fig,ax=plots.multi(2,3)
    gd = np.where(n>minvisit)[0]
    ax[0,0].hist(vscat[gd],bins=np.arange(0.01,1,0.01),histtype='step',cumulative=True,normed=True,color='b')
    ax[2,0].hist(vdiff,color='b',bins=np.arange(-1.,1,0.01),histtype='step')
    gd=np.where((n>minvisit) & (verr < 0.2))[0]
    ax[0,0].hist(vscat[gd],bins=np.arange(0.01,1,0.01),histtype='step',cumulative=True,normed=True,color='g')
    ax[2,1].hist(vscat[gd],bins=np.arange(0.01,1,0.01),histtype='step',color='g')
    fig.tight_layout()
    plt.show()
    return vhelio,vscat,verr,sigfiber,vdiff,n


def field(name,dr14=False,dir='./',minvisit=1) :
    """ look at a single field
    """
    all=struct.concat([dir+'/apVisitSum*.fits'])
    if name == 'M67' : locid=[os.environ['APOGEE_REDUX']+'/r8/fields/apo25m/4162//apVisitSum*']
    elif name == 'N188' : locid=[os.environ['APOGEE_REDUX']+'/r8/fields/apo25m/4217//apVisitSum*', 
                                 os.environ['APOGEE_REDUX']+'/r8/fields/apo25m/5067//apVisitSum*']
    alldr14=struct.concat(locid)
    objs = set(all['APOGEE_ID'])
    vhelio = []
    vscat = []
    verr = []
    sigfiber = []
    vdiff = []
    n = []
    dr14vhelio = []
    dr14vscat = []
    dr14sigfiber = []
    dr14n = []
    dr14vdiff = []
    for obj in objs :
        j = np.where(all['APOGEE_ID'] == obj)[0]
        vhelio.append(all['VHELIO'][j].mean())
        vscat.append(all['VHELIO'][j].std())
        verr.append(all['VRELERR'][j].max())
        sigfiber.append(all['FIBERID'][j].std())
        vdiff.extend(all['VHELIO'][j]-all['VHELIO'][j].mean())
        n.append(len(j))
        #print(all['MJD'][j],all['VHELIO'][j])
        j = np.where(alldr14['APOGEE_ID'] == obj)[0]
        dr14vhelio.append(alldr14['VHELIO'][j].mean())
        dr14vscat.append(alldr14['VHELIO'][j].std())
        dr14sigfiber.append(alldr14['FIBERID'][j].std())
        dr14n.append(len(j))
        dr14vdiff.extend(alldr14['VHELIO'][j]-alldr14['VHELIO'][j].mean())
        #print(all['MJD'][j],all['VHELIO'][j],all['VRELERR'][j])
        #print(alldr14['MJD'][j],alldr14['VHELIO'][j],alldr14['VRELERR'][j])
        #pdb.set_trace()
    vhelio=np.array(vhelio)
    vscat=np.array(vscat)
    verr=np.array(verr)
    sigfiber=np.array(sigfiber)
    n=np.array(n)
    dr14vhelio=np.array(dr14vhelio)
    dr14vscat=np.array(dr14vscat)
    dr14sigfiber=np.array(dr14sigfiber)
    dr14n=np.array(dr14n)
    fig,ax=plots.multi(2,3)
    gd =np.where(n > minvisit)[0]
    ax[0,0].hist(vscat[gd],bins=np.arange(0.01,1,0.01),histtype='step',cumulative=True,normed=True,color='b')
    ax[0,0].hist(dr14vscat[gd],bins=np.arange(0.01,1,0.01),histtype='step',cumulative=True,normed=True,color='r')
    gd=np.where((verr < 0.2) & (n>minvisit))[0]
    ax[0,0].hist(vscat[gd],bins=np.arange(0.01,1,0.01),histtype='step',cumulative=True,normed=True,color='g')
    ax[0,0].hist(dr14vscat[gd],bins=np.arange(0.01,1,0.01),histtype='step',cumulative=True,normed=True,color='m')
    ax[2,1].hist(vscat[gd],bins=np.arange(0.01,1,0.01),histtype='step',color='k')
    ax[2,1].hist(dr14vscat[gd],bins=np.arange(0.01,1,0.01),histtype='step',color='r')

    plots.plotc(ax[1,0],vhelio-dr14vhelio,vscat-dr14vscat,verr,xr=[-0.5,0.5],yr=[-0.3,0.3],zr=[0,0.15])
    plots.plotc(ax[0,1],sigfiber,vscat-dr14vscat,verr,zr=[0,0.15],yr=[-0.3,0.3])
    plots.plotc(ax[1,1],vscat,vscat-dr14vscat,verr,zr=[0,0.15],yr=[-0.3,0.3],xr=[0,0.5])
    ax[2,0].hist(vdiff,color='b',bins=np.arange(-1.,1,0.01),histtype='step')
    ax[2,0].hist(dr14vdiff,color='r',bins=np.arange(-1.,1,0.01),histtype='step')
    fig.tight_layout()
    plt.show()


def visitcomp(plate,mjd,indiv=False,apred='test') :
    """ Compare RVs for plate/mjd with DR14 RVs
    """
    #plt.close('all')
    load=apload.ApLoad(apred=apred)
    a=load.apVisitSum(plate,mjd)[1].data

    #dr14=apload.ApLoad(dr='dr14')
    p=yanny.yanny(os.environ['PLATELIST_DIR']+'/platePlans.par',np=True)
    j=np.where(p['PLATEPLANS']['plateid'] == plate)[0][0]
    locid=p['PLATEPLANS']['locationid'][j]
    b=fits.open(os.environ['APOGEE_REDUX']+'/r8/fields/apo25m/{:04d}/apVisitSum-{:04d}-{:05d}.fits'.format(
                   locid,plate,mjd))[1].data

    fig,ax=plots.multi(2,2)
    i1,i2=match.match(a['FIBERID'],b['FIBERID'])
    plots.plotc(ax[0,0],a['FIBERID'][i1],a['VHELIO'][i1]-b['VHELIO'][i2],a['RV_TEFF'][i1],zr=[3500,5500],xt='Fiber',yt=r'$\Delta$ VHELIO')
    plots.plotc(ax[0,1],a['FIBERID'][i1],a['VHELIO'][i1]-b['VHELIO'][i2],a['RV_TEFF'][i1],zr=[3500,5500],yr=[-2.,2.],xt='Fiber',yt=r'$\Delta$ VHELIO')
    plots.plotc(ax[1,0],a['FIBERID'][i1],a['RV_TEFF'][i1]-b['RV_TEFF'][i2],a['RV_TEFF'][i1],zr=[3500,5500],xt='Fiber',yt=r'$\Delta$ RV_TEFF')
    plots.plotc(ax[1,1],a['RV_TEFF'][i1]-b['RV_TEFF'][i2],a['VHELIO'][i1]-b['VHELIO'][i2],a['RV_TEFF'][i1],zr=[3500,5500],xt=r'$\Delta$ RV_TEFF',yt=r'$\Delta$ VHELIO')
    out=load.filename('Plate',chips=True,plate=plate,mjd=mjd)
    outdir=os.path.dirname(out)
    outname=os.path.basename(out).replace('-a','').replace('.fits','_dr14comp.png')
    fig.tight_layout()
    pdb.set_trace()
    print(outdir+'/plots/'+outname)
    fig.savefig(outdir+'/plots/'+outname)

    if indiv :
        va=load.apPlate(plate,mjd)
        vb={}
        for chip in chips :
            tmp=fits.open(os.environ['APOGEE_REDUX']+'/r8/apo25m/{:04d}/{:05d}/apPlate-{:s}-{:04d}-{:05d}.fits'.format(
                   plate,mjd,chip,plate,mjd))
            vb[chip] = tmp

        fig,ax=plots.multi(1,3,hspace=0.3) 
        pfig,pax=plots.multi(1,3,hspace=0.3) 
        wfig,wax=plots.multi(1,3,hspace=0.3) 
        for i in range(len(i1)) :
            fiber = a['FIBERID'][i1[i]]
            if (a['VHELIO'][i1[i]]-b['VHELIO'][i2[i]]) > 0.5 :
              print(fiber,a['VHELIO'][i1[i]],b['VHELIO'][i2[i]],a['RV_TEFF'][i1[i]],b['RV_TEFF'][i2[i]])
              applot.chip(va,ax=ax,row=300-fiber,color='r')
              #applot.chip(va,ax=pax,row=300-fiber,color='r',pixel=True)
              applot.chip(vb,ax=ax,row=300-fiber,color='b')
              #applot.chip(vb,ax=pax,row=300-fiber,color='b',pixel=True)
              for ichip,chip in enumerate(chips) :
                  pax[ichip].plot(va[chip][1].data[300-fiber,:]/vb[chip][1].data[300-fiber,:])
                  wax[ichip].plot(va[chip][4].data[300-fiber,:]-vb[chip][4].data[300-fiber,:])
              plt.show()

              pdb.set_trace()
              for ichip in range(3) :
                  ax[ichip].cla()
                  pax[ichip].cla()
                  wax[ichip].cla()
    plt.close()

def dr14comp(a,b,av,bv):
    """ compare multiple field RVs from, e.g. allField file with DR14
    """
    load=apload.ApLoad(apred='r11')
    dr14=apload.ApLoad(dr='dr14')

    i1,i2=match.match(a['APOGEE_ID'],b['APOGEE_ID'])
    gd = np.where((a['NVISITS'][i1] == b['NVISITS'][i2]) & (a['SNR'][i1]>75) )[0]
    a=a[i1[gd]]
    b=b[i2[gd]]
   
    j=np.argsort(a['VHELIO_AVG']-b['VHELIO_AVG'])
 
    fig,ax=plots.multi(1,3,hspace=0.3) 
    pfig,pax=plots.multi(1,3,hspace=0.3) 
    wfig,wax=plots.multi(1,3,hspace=0.3) 
    chips=['a','b','c']
    for jj in j :
       j1=np.where(av['APOGEE_ID'] == a['APOGEE_ID'][jj])[0]
       j2=np.where(bv['APOGEE_ID'] == a['APOGEE_ID'][jj])[0]
       print(a['APOGEE_ID'][jj],a['RV_TEFF'][jj],b['RV_TEFF'][jj],a['SNR'][jj],b['SNR'][jj])
       for jjj,kkk in zip(j1,j2) : 
           print(av['MJD'][jjj],av['PLATE'][jjj],av['FIELD'][jjj],av['SNR'][jjj],av['FIBERID'][jjj],av['VHELIO'][jjj],av['ESTVHELIO'][jjj])
           print(bv['MJD'][kkk],bv['PLATE'][kkk],bv['FIELD'][kkk],bv['SNR'][kkk],bv['FIBERID'][kkk],bv['VHELIO'][kkk],bv['ESTVHELIO'][kkk])
           va=load.apPlate(int(av['PLATE'][jjj]),av['MJD'][jjj])
           vsum=load.apVisitSum(int(av['PLATE'][jjj]),av['MJD'][jjj])[1].data
           f=np.where(vsum['FIBERID'] == av['FIBERID'][jjj])[0]
           print(vsum['RV_TEFF'][f])
           applot.chip(va,ax=ax,row=300-av['FIBERID'][jjj],color='r')
           applot.chip(va,ax=pax,row=300-av['FIBERID'][jjj],color='r',pixel=True)
           vb={}
           for chip in chips :
             tmp=fits.open(os.environ['APOGEE_REDUX']+'/r8/apo25m/{:04d}/{:05d}/apPlate-{:s}-{:04d}-{:05d}.fits'.format(
               int(bv['PLATE'][kkk]),bv['MJD'][kkk],chip,int(bv['PLATE'][kkk]),bv['MJD'][kkk]))
             vb[chip] = tmp
           vsum=fits.open(os.environ['APOGEE_REDUX']+'/r8/fields/apo25m/{:04d}/apVisitSum-{:04d}-{:05d}.fits'.format(
               int(bv['LOCATION_ID'][kkk]),int(bv['PLATE'][kkk]),bv['MJD'][kkk]))[1].data
           f=np.where(vsum['FIBERID'] == bv['FIBERID'][kkk])[0]
           print(vsum['RV_TEFF'][f])
           applot.chip(vb,ax=ax,row=300-bv['FIBERID'][kkk],color='b')
           applot.chip(vb,ax=pax,row=300-bv['FIBERID'][kkk],color='b',pixel=True)
           for ichip,chip in enumerate(chips) :
               wax[ichip].plot(va[chip][4].data[300-av['FIBERID'][jjj],:]-vb[chip][4].data[300-bv['FIBERID'][kkk],:])
           plt.show()
           pdb.set_trace()

       for ichip in range(3) :
           ax[ichip].cla()
           pax[ichip].cla()
           wax[ichip].cla()

def standards(a,out=None) :
    """ Compare RVs to standards
    """
    stan = fits.open(os.environ['APOGEE_DIR']+'/data/rv/rvstandards.fits')[1].data
    h=esutil.htm.HTM()
    m1,m2,rad=h.match(a['ra'],a['dec'],stan['ra'],stan['dec'],1./3600.,maxmatch=500)
    fig,ax=plots.multi(1,1)
    ax.hist(a['VHELIO_AVG'][m1]-stan['RV'][m2],histtype='step',bins=np.arange(-1,1,0.1))
    ax.set_xlabel('RV(APOGEE) - RV(lit)')
    if out is not None :
        fig.savefig(out+'.png')
        plt.close()

   
def all(a,name='DR16',dr='dr14') :
    """ Do a series of RV comparisons for input data and previous DR
    """
    grid=[]
    xtit=[]
    load=apload.ApLoad(dr=dr)
    b=load.allStar()[1].data
   
    # vscatter of new RVs
    vscat(a,out='plots/vscat')
    vscat(a,out='plots/vscat5',nmin=5)
    grid.append(['vscat.png','vscat5.png'])
    xtit.append(name+' : VSCATTER')

    # APO/LCO comparison
    apolco(a,out='plots/apolco')
    grid.append(['apolco_1.png','apolco_2.png'])
    xtit.append(name+' : APO (solid) and LCO (dotted), same stars')
    #apolco(a,out='plots/apolco_nolowz',minfeh=-0.6)
    #grid.append(['apolco_nolowz_1.png','apolco_nolowz_2'])
    #xtit.append(name+', no low [Fe/H]: APO (solid) and LCO (dotted)')

    # RV standards
    standards(a,out='plots/rvstan')
    grid.append(['rvstan.png',''])
    xtit.append(name+', comparison with literature RVs')

    # comparison with previous DR
    comp(a,b,domatch=True,out='plots/drcomp')
    grid.append(['drcomp_1.png','drcomp_2.png'])
    xtit.append(name+', comparison with '+dr+' : same stars, same NVISITS, new(solid) '+dr+'(dotted)')

    html.htmltab(grid,ytitle=xtit,file='plots/rv.html')


def repeatspec(a) :
    stars=set(a['APOGEE_ID'])
    fig,ax=plots.multi(1,2,hspace=0.001,sharex=True)
    for star in stars :
        j=np.where(a['APOGEE_ID'] == star)[0]
        if len(j) > 1 :
            for i in j :
                print(a['TELESCOPE'][i],a['FIELD'][i],a['NVISITS'][i])
                spec=fits.open(a['TELESCOPE'][i]+'/'+a['FIELD'][i]+'/apStar-r12-'+a['APOGEE_ID'][i]+'.fits')
                if i == j[0] : spec0=copy.copy(spec[1].data[0,:])
                plots.plotl(ax[0],spectra.fits2vector(spec[1].header,1),spec[1].data[0,:])
                plots.plotl(ax[1],spectra.fits2vector(spec[1].header,1),spec[1].data[0,:]/spec0,yr=[0.5,1.5])
            plt.show()
            pdb.set_trace()
            ax[0].cla()
            ax[1].cla()
  

from scipy.signal import convolve
def dop_comp(field) :
    """ Compare RVs from different data releases
    """
    dop = fits.open(field+'/'+field+'_rv.fits')
    r13 = apload.ApLoad(apred='r13')
    old = r13.apField(field)

    i1,i2 = match.match(dop[1].data['APOGEE_ID'],old[1].data['APOGEE_ID'])
    print(len(dop[1].data),len(old[1].data),len(i1))

    fig,ax=plots.multi(1,1)
    plots.plotc(ax,dop[1].data['RV_TEFF'][i1],dop[1].data['VHELIO_AVG'][i1]-old[1].data['VHELIO_AVG'][i2],dop[1].data['VSCATTER'][i1])

    j=np.argsort(np.abs(dop[1].data['VHELIO_AVG'][i1]-old[1].data['VHELIO_AVG'][i2],dop[1].data['VSCATTER'][i1]))

    plots._data = dop[1].data
    plots._id_cols=['APOGEE_ID']
    plots.event(fig)
    key=' '
    sf,sax=plots.multi(1,2,sharex=True,hspace=0.001)
    while key != 'e' :
        x,y,key,index = plots.mark(fig,index=True)
        obj = dop[1].data['APOGEE_ID'][i1[index]]
        #jv = np.where(dop[2].data['APOGEE_ID'] == dop[1].data['APOGEE_ID'][i1])[0]
        out=pickle.load(open(field+'/'+obj+'_out.pkl','rb'))
        print(obj,old[1].data['APOGEE_ID'][i2[index]])
        print(out[0])
        sax[0].cla()
        spec=old[2].data['SPEC'][i2[index]]
        plots.plotl(sax[0],old[3].data['WAVE'][0,:],spec/convolve(spec,np.ones(500)/500,mode='same'),xr=[15000,17000],yr=[0.5,1.5])
        for mod,obs in zip(out[2],out[3]) :
            sax[1].cla()
            for chip in range(3) :
                plots.plotl(sax[1],obs.wave[:,chip],obs.flux[:,chip],color='k',yr=[0.5,1.5])
                gd = np.where(obs.mask[:,chip] == False)[0]
                plots.plotl(sax[1],obs.wave[gd,chip],obs.flux[gd,chip],color='g')
                plots.plotl(sax[1],mod.wave[:,chip],mod.flux[:,chip],color='r')
            plt.draw()
            input('hit a key: ')

def overlap(fields) :
    """ compare RVs from different fields for overlapping stars
    """

    r13=apload.ApLoad(apred='r13')
    f=[]
    a=[]
    for field in fields :
        f.append(fits.open(field+'/'+field+'_rv.fits'))
        a.append( r13.apFieldVisits(field))

    outdir=fields[0]+'_'+fields[1]
    try: os.makedirs(outdir)
    except: pass

    fp=open(outdir+'/'+outdir+'.html','w')
    fp.write('<HTML>\n')
    fp.write('<HEAD><script type=text/javascript src=../html/sorttable.js></script></head>')
    fp.write('<BODY>\n')
    fp.write('<TABLE BORDER=2>\n')

    matplotlib.use('Agg')
    i1,i2=match.match(f[0][1].data['APOGEE_ID'],f[1][1].data['APOGEE_ID'])
    colors=['g','r','b','m']
    for star in f[0][1].data['APOGEE_ID'][i1] :
        print(star)
        fp.write('<TR><TD>{:s}<BR>\n'.format(star))
        fp.write('<TABLE BORDER=2>\n')
        fig,ax=plots.multi(1,1)
        for i,field in enumerate(f) :
            j=np.where(field[2].data['APOGEE_ID'] == star)[0]
            plots.plotp(ax,field[2].data['MJD'][j],field[2].data['VHELIO'][j],color=colors[i],size=10)
            j=np.where(field[1].data['APOGEE_ID'] == star)[0][0]
            fp.write('<TR><TD>Doppler<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.0f}<TD>{:8.2f}<TD>{:8.2f}\n'.format(
                 field[1].data['VHELIO_AVG'][j],field[1].data['VSCATTER'][j],
                 field[1].data['RV_TEFF'][j],field[1].data['RV_LOGG'][j],field[1].data['RV_FEH'][j]))
        for i,field in enumerate(a) :
            j=np.where(field[1].data['APOGEE_ID'] == star)[0]
            gd=np.where(np.abs(field[1].data['VHELIO'][j]) < 999)[0]
            plots.plotp(ax,field[1].data['MJD'][j[gd]],field[1].data['VHELIO'][j[gd]],color=colors[i+2],size=10)
            #fp.write('<TR><TD>DR16<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.0f}<TD>{:8.2f}<TD>{:8.2f}\n'.format(
            #     field[1].data['VHELIO_AVG'],field[1].data['VSCATTER'],
            #     field[1].data['RV_TEFF'],field[1].data['RV_LOGG'],field[1].data['RV_FEH']))
        plt.draw()
        plt.close()
        fig.savefig(outdir+'/'+star+'.png')
        fp.write('</TABLE>\n')
        fp.write('<TD><a HREF={:s}.png> <IMG SRC={:s}.png> </a>\n'.format(star,star))
    fp.write('</TABLE>')
    fp.close()

def comp_apstar(field,apred='r13',telescope='apo25m') :

    load=apload.ApLoad(apred=apred,telescope=telescope)
    apfield=load.apField(field)

    fig,ax=plots.multi(1,3,hspace=0.001)
    w=aspcap.apStarWave() 
    for i,star in enumerate(apfield[1].data['APOGEE_ID']) : 
        old=load.apStar(field,star)
        plots.plotl(ax[0],w,old[1].data[0,:]*1.e-17)
        for j in [2] :plots.plotl(ax[1],w,old[2].data[j,:]*1.e-17)
        plots.plotl(ax[2],w,old[3].data[0,:])
        #plots.plotl(ax[0],old[3].data['WAVE'][0,:],old[2].data['SPEC'][i,:]) 
        #plots.plotl(ax[1],old[3].data['WAVE'][0,:],old[2].data['ERR'][i,:]) 
        #plots.plotl(ax[2],old[3].data['WAVE'][0,:],old[2].data['MASK'][i,:]) 

        new=fits.open('apo25m/K01_082+17/apStar-r13-'+star+'.fits') 
        plots.plotl(ax[0],w,new[1].data[0,:]*1.e-17) 
        for j in [2] :plots.plotl(ax[1],w,new[2].data[j,:]*1.e-17)
        ax[1].set_ylim(0,1.e-16)
        plots.plotl(ax[2],w,new[3].data[0,:])
        plt.draw() 
        pdb.set_trace() 
        ax[0].cla()    
        ax[1].cla()    
        ax[2].cla()    

def star(apfv, objs, htmlfile='rv.html',visittable=True, apred='dr17', apstar='stars') :
    """ Make web page with RV info for specified stars
    """

    fp = open(htmlfile,'w')
    fp.write('<HTML><BODY>')
    fp.write('<TABLE BORDER=2>')

    load=apload.ApLoad(apred=apred,apstar=apstar)
    # individual visit velocities
    starmask = bitmask.StarBitMask()
    for obj in objs :
        fp.write('<TR bgcolor=lightblue><TD colspan=5>'+obj)
        j = np.where(apfv['APOGEE_ID'] == obj)[0]
        fields = set(apfv['FIELD'][j])
        for field in fields :
          j = np.where((apfv['APOGEE_ID'] == obj) & (apfv['FIELD'] == field))[0]

          telescope=apfv['TELESCOPE'][j[0]]
          load.settelescope(telescope)
          # try to get apStar file for STARFLAGS and make combined ccf plot
          try:
              aps=load.apStar(field,obj)
              rvtab=aps[10].data
              starflags=starmask.getname(aps[0].header['STARFLAG'])
              snr=aps[0].header['SNR']
              nccf= len(rvtab['x_ccf'])
              fig,ax=plots.multi(1,nccf,hspace=0.001)
              ax=np.atleast_1d(ax)
              for i in range(nccf) :
                  ax[i].plot(rvtab['x_ccf'][i],rvtab['ccf'][i])
                  ax[i].plot(rvtab['x_ccf'][i],rvtab['ccf'][i]/rvtab['ccf'][i].max())
                  n=len(rvtab['x_ccf'][i])
                  ax[i].plot(rvtab['x_ccf'][i],rvtab['autoccf'][i][:n])
                  ax[i].text(0.1,0.9,'ccpfwhm: {:8.2f}   autofwhm: {:8.2f}'.format(rvtab['ccpfwhm'][i],rvtab['autofwhm'][i]),transform=ax[i].transAxes)
              plt.savefig('plots/{:s}_autoccf.png'.format(obj))
              plt.close()
          except : 
              starflags=''
              snr=-999.

          fp.write('<TR><TD>')
          fp.write('APOGEE_ID: {:s}<BR>'.format(obj))
          fp.write('FIELD: {:s}<BR>'.format(field))
          star=apfv[j[0]]
          fp.write('H: {:8.2f}<BR>'.format(star['H']))
          fp.write('RV_TEFF: {:8.2f}<BR>'.format(star['RV_TEFF']))
          fp.write('STARFLAGS: {:s}<BR>'.format(starflags))
          fp.write('SNR: {:8.2f}<BR>'.format(snr))
          try: fp.write('CCFWHM: {:8.2f}  AUTOFWHM: {:8.2f} <BR>'.format(star['RV_CCFWHM'],star['RV_AUTOFWHM']))
          except: pass
          if visittable :
              fp.write('<TABLE BORDER=2>')
              fp.write('<TR><TD>JD<TD>PLATE<TD>MJD<TD>FIBER<TD>S/N<TD>Doppler xcorr<TD> xcorr_err<TD>Doppler<TD>VERR<TD>BC<TD>N_COMPONENTS<TD>CCWFHM<TD>AUTOFWHM\n')
              for i in j : 
                if np.isfinite(apfv['VHELIO'][i]) == False :
                    bgcolor='bgcolor=red'
                elif apfv['STARFLAG'][i] & starmask.getval('RV_REJECT') > 0 :
                    bgcolor='bgcolor=lightpink'
                elif apfv['STARFLAG'][i] & starmask.getval('RV_SUSPECT') > 0 :
                    bgcolor='bgcolor=#F4DEDE'
                else : bgcolor=''
                fp.write(('<TR {:s}> <TD> <A HREF={:s} TARGET="_obj"> {:12.3f}</A> <TD> {:s} <TD> {:5d} <TD> {:5d}'+
                     '<TD> {:8.1f} <TD> {:8.2f} <TD> {:8.2f} <TD> {:8.2f} ' +
                     '<TD> {:8.2f} <TD>{:8.2f} <TD>{:d} <TD>{:8.2f} <TD> {:8.2f}\n').format(
                      bgcolor,
                      apfv['FILE'][i].replace('.fits','_dopfit.png').replace('-r12-','-r13-'),
                      apfv['JD'][i],apfv['PLATE'][i],apfv['MJD'][i],apfv['FIBERID'][i],
                      apfv['SNR'][i],
                      apfv['XCORR_VHELIO'][i],apfv['XCORR_VRELERR'][i],
                      apfv['VHELIO'][i],apfv['VRELERR'][i],apfv['BC'][i],apfv['N_COMPONENTS'][i],apfv['CCFWHM'][i],apfv['AUTOFWHM'][i]))
              fp.write('</TABLE>\n')

          plotdir = apstar+'/{:s}/{:s}/plots/'.format(telescope,field)
          #fp.write('<TD><a HREF={:s}/{:s}.png TARGET="_obj"> <IMG SRC={:s}/{:s}.png WIDTH=600></A>\n'.format(plotdir,obj,plotdir,obj))
          #fp.write('<TD><IMG SRC={:s}/{:s}_rv.png TARGET="_obj">\n'.format(plotdir,obj))
          fp.write('<TD><A HREF={:s}/{:s}_ccf.png TARGET="_obj"> <IMG SRC={:s}/{:s}_ccf.png></A>\n'.format(plotdir,obj,plotdir,obj))
          #fp.write('<TD><A HREF={:s}/{:s}_spec.png TARGET="_obj"> <IMG SRC={:s}/{:s}_spec.png></a>\n'.format(plotdir,obj,plotdir,obj))
          fp.write('<TD><A HREF={:s}/{:s}_spec2.png TARGET="_obj"> <IMG SRC={:s}/{:s}_spec2.png></a>\n'.format(plotdir,obj,plotdir,obj))
          #fp.write('<TD><A HREF={:s}/{:s}_cont.png TARGET="_obj"> <IMG SRC={:s}/{:s}_cont.png></a>\n'.format(plotdir,obj,plotdir,obj))
          fp.write('<TD><A HREF=plots/{:s}_autoccf.png TARGET="_obj"> <IMG SRC=plots/{:s}_autoccf.png></a>\n'.format(obj,obj))



    fp.write('</TABLE></BODY></HTML>')
    fp.close() 
