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
from astropy.io import fits
from apogee.utils import apload
from apogee.utils import applot
from apogee.utils import bitmask
from apogee.utils import spectra
from apogee.aspcap import norm
from tools import plots
from tools import html
from tools import match
from tools import struct
from sdss import yanny
from scipy import interpolate
from scipy.signal import correlate
from scipy.ndimage.filters import median_filter, gaussian_filter

colors=['r','g','b','c','m','y','k']
chips=['a','b','c']


def allField(files=['apo*/*/apField-*.fits','apo*/*/apFieldC-*.fits','lco*/*/apField-*.fits'],out='allField.fits',verbose=False) :
    """ Concatenate set of apField files 
    """
    # concatenate the structures
    all=struct.concat(files,verbose=verbose)

    # write out the file
    if out is not None:
        print('writing',out)
        struct.wrfits(all,out)

    return all

def allFieldVisits(files=['apo*/*/apFieldVisits-*.fits','apo*/*/apFieldC-*.fits','lco*/*/apFieldVisits-*.fits'],out='allFieldVisits.fits',verbose=False) :
    """ Concatenate set of apFieldVisit files
    """
    # concatenate the structures
    all=struct.concat(files,verbose=verbose)

    # write out the file
    if out is not None:
        print('writing',out)
        struct.wrfits(all,out)

    return all


def vscat(a,fig=None,ls=None,marker='o',nmin=2,mhmin=-3,density=False,out=None) :
    """ Make histograms of VSCATTER for different bins of Teff H], given min NVISITS, and min [M/H]
    """
    if fig == None : fig,ax=plots.multi(4,6,hspace=0.001,wspace=0.4,figsize=(12,8))
    else : fig,ax=fig
    tbins=[3000,3500,4000,4500,5500,8000,30000] 
    hbins=[8,11,12,13,15]
    snr = a['SNREV']
    j=np.where(snr > 300) [0]
    snr[j] = 300
    for i in range(len(tbins)-1) :
        ax[i,0].text(0.9,0.9,'{:d}<=RV_TEFF<{:d}'.format(tbins[i],tbins[i+1]),ha='right',transform=ax[i,0].transAxes,fontsize=8)
        for j in range(len(hbins)-1) :
            ax[0,j].set_title('{:d}<=H<{:d}'.format(hbins[j],hbins[j+1]))
            gd = np.where((a['RV_TEFF']>=tbins[i]) & (a['RV_TEFF']<tbins[i+1]) &
                          (a['H']>=hbins[j]) & (a['H']<hbins[j+1]) &
                           (a['NVISITS']>nmin) & (a['RV_FEH']>mhmin) ) [0]
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


def visitspec(load,plate,mjd,fiber,gridfile='apg_rvsynthgrid',apstar=False) :
    """ Crude beginnings of an RV routine
    """
    grid = fits.open(os.environ['APOGEE_DIR']+'/data/synthgrid/'+gridfile+'.fits')
    if gridfile == 'apg_rvsynthgrid' : hdu=1
    elif gridfile == 'apg_rvsynthgrid_v2': hdu=0
    elif apstar : hdu=2
    else : hdu=1
    gridspec=grid[hdu].data
    gridwave = 10.**spectra.fits2vector(grid[hdu].header,2)
    griderr = np.ones(gridspec.shape[0])
    #for ispec in range(gridspec.shape[1]) :
    #    cont = norm.cont(gridspec[:,ispec],griderr)
    #    gridspec[:,ispec] /= cont

    data = load.apVisit(plate,mjd,fiber)

    # compare with DR14 
    comp(a,b,domatch=False,out='plots/dr14all')
    grid.append(['dr14all_1.png',''])
    xtit.append('all stars: DR14 (dotted) and test DR16 (solid)')

    comp(a,b,domatch=True,out='plots/dr14match')
    grid.append(['dr14match_1.png','dr14match_2.png'])
    xtit.append('same stars: DR14 (dotted) and test DR16 (solid)')
    # set bad pixels to nan
    shape=data[1].data.shape
    spec = copy.copy(data[1].data).flatten()
    specerr = copy.copy(data[2].data)
    specwave=data[4].data
    pixmask=bitmask.PixelBitMask()
    bd = np.where( ((data[3].data & pixmask.badval()) > 0) | 
                   ((data[3].data & pixmask.getval('SIG_SKYLINE')) > 0) ) [0]
    spec[bd] = np.nan
    spec = spec.reshape(shape)

    # continuum normalize and sample to grid
    outspec = np.full(len(gridwave),np.nan)
    if not apstar :
        # apVisit wavelengths are reversed
        spec=np.flip(spec)
        specwave=np.flip(specwave)
        specerr=np.flip(specerr)
        for ichip in range(3) :
            cont = norm.cont(spec[ichip,:],specerr[ichip,:])
            spec[ichip,:] /= cont
            gd=np.where(np.isfinite(spec[ichip,:]))[0]
            ip= interpolate.InterpolatedUnivariateSpline(specwave[ichip,gd],spec[ichip,gd],k=3)
            out = ip(gridwave)
            gd = np.where( (gridwave > specwave[ichip,0]) & (gridwave < specwave[ichip,-1]) )[0]
            outspec[gd] = out[gd]
            plt.plot(specwave[ichip,:],spec[ichip,:])
            plt.plot(gridwave[gd],out[gd])
            plt.show()

    for ispec in range(gridspec.shape[1]) :
        print(ispec)
        bd=np.where(np.isnan(outspec))
        outspec[bd]=1.
        out=correlate(outspec,gridspec[:,ispec])
        pdb.set_trace() 

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
  
import doppler 
import multiprocessing as mp
from astropy.table import Table
from apogee.apred import bc

def doppler_rv(field,telescope='apo25m',apred='r13',obj=None,nobj=0,threads=8,maxvisit=500,snmin=3,
               clobber=False,verbose=False,tweak=False,plot=False,windows=None) :
    """ Run DOPPLER RVs for a field
    """ 
    
    # get all the VisitSum files for this field and concatenate them
    files=glob.glob(os.environ['APOGEE_REDUX']+'/'+apred+'/visit/'+telescope+'/'+field+'/apVisitSum*')
    allvisits=struct.concat(files)
    starmask=bitmask.StarBitMask()
    gd=np.where(((allvisits['STARFLAG'] & starmask.badval()) == 0) & (allvisits['SNR'] > snmin) )[0]
    print(len(allvisits),len(gd))
    allvisits=Table(allvisits)

    # output directory
    try: os.mkdir(field)
    except FileExistsError: pass

    # get all unique (or requested) objects
    if obj is None :
        if nobj > 0 :
            allobj=set(allvisits['APOGEE_ID'][0:nobj])
        else :
            allobj=set(allvisits['APOGEE_ID'])
    else :
        allobj = obj

    # output apField structure
    fieldtype = np.dtype([('FILE','S64'),('APOGEE_ID','S20'),('TELESCOPE','S6'),('LOCATION_ID',int),('FIELD','S20'),
                          ('J',float),('J_ERR',float),('H',float),('H_ERR',float),('K',float),('K_ERR',float),
                          ('RA',float),('DEC',float),('GLON',float),('GLAT',float),
                          ('AK_TARG',float),('AK_TARG_METHOD','S32'),
                          ('AK_WISE',float),('SFD_EBV',float),
                          ('APOGEE_TARGET1',int),('APOGEE_TARGET2',int),('APOGEE_TARGET3',int),
                          ('APOGEE2_TARGET1',int),('APOGEE2_TARGET2',int),('APOGEE2_TARGET3',int),('APOGEE2_TARGET4',int),
                          ('TARGFLAGS','S132'),('SURVEY','S16'),('PROGRAMNAME','S32'),
                          ('NINST',int),('NVISITS',int),('COMBTYPE',int),('COMMISS',int),
                          ('SNR',float),('STARFLAG',float),('STARFLAGS','S132'),('ANDFLAG',float),('ANDFLAGS','S132'),
                          ('VHELIO_AVG',float),('VSCATTER',float),('VERR',float),
                          ('RV_TEFF',float),('RV_LOGG',float),('RV_FEH',float),('RV_ALPHA',float),('RV_CARB',float),
                          ('RV_CCPFWHM',float),('RV_AUTOFWHM',float),
                          ('N_COMPONENTS',int)
                         ])
    allfield = np.zeros(len(allobj),dtype=fieldtype)
    allfield['TELESCOPE'] = telescope
    allfield['FIELD'] = field

    allfiles=[]
    allv=[]
    load=apload.ApLoad(apred=apred,telescope=telescope)
    nobj=0
    nvisit=0
    pixelmask=bitmask.PixelBitMask()

    # loop over requested objects, building up allfiles list of 
    #  [(field,obj,clobber,verbose,tweak,plot,windows),filenames....] to pass to dorv()
    for iobj,star in enumerate(sorted(allobj)) :
        if type(star) is str : star=star.encode()
        allfield['APOGEE_ID'][iobj] = star
        # we will only consider good visits
        visits=np.where(allvisits['APOGEE_ID'][gd] == star)[0]
        print('object: {:}  nvisits: {:d}'.format(star,len(visits)))
        nobj+=1
        nvisit+=len(visits)

        if len(visits) > 0 :
            allfiles.append([allvisits[gd[visits]],load,(field,star,clobber,verbose,tweak,plot,windows)])

    print('total objects: ', nobj, ' total visits: ', nvisit) 
    if threads == 0 :
        output=[]
        for speclist in allfiles :
            print(speclist)
            output.append(dorv(speclist))
    else :
        pool = mp.Pool(threads)
        output = pool.map_async(dorv, allfiles).get()
        pool.close()
        pool.join()
    print('done pool')

    # rename old visit RV tags and initialize new ones
    for col in ['VTYPE','VREL','VRELERR','VHELIO','BC'] :
        allvisits.rename_column(col,'EST'+col)
        if col == 'VTYPE' : allvisits[col] = 0
        else : allvisits[col] = np.nan
    for col in ['XCORR_VREL','XCORR_VRELERR','XCORR_VHELIO','BC'] :
        allvisits[col] = np.nan
    allvisits['N_COMPONENTS'] = -1

    # load up the individual visit RV information
    allv=[]
    for out,files in zip(output,allfiles) :
        starflag, andflag = 0, 0
        apogee_target1, apogee_target2, apogee_target3 = 0, 0, 0
        apogee2_target1, apogee2_target2, apogee2_target3, apogee2_target4 = 0, 0, 0, 0
        try :
          if out is not None :
            visits=[]
            ncomponents=0
            for i,v in enumerate(out[0][1]) :
                # match by filename components in case there was an error reading in doppler
                name=os.path.basename(v['filename']).replace('.fits','').split('-')
                visit = np.where( (np.char.strip(allvisits['PLATE']).astype(str) == name[-3]) &
                                  (allvisits['MJD'] == int(name[-2])) &
                                  (allvisits['FIBERID'] == int(name[-1])) )[0][0]
                visits.append(visit)
                allvisits[visit]['VREL']=v['vrel']
                allvisits[visit]['VRELERR']=v['vrelerr']
                allvisits[visit]['VHELIO']=v['vhelio']
                allvisits[visit]['XCORR_VREL']=v['xcorr_vrel']
                allvisits[visit]['XCORR_VRELERR']=v['xcorr_vrelerr']
                allvisits[visit]['XCORR_VHELIO']=v['xcorr_vhelio']
                allvisits[visit]['BC']=v['bc']
                starflag |= allvisits[visit]['STARFLAG']
                andflag &= allvisits[visit]['STARFLAG']
                if allvisits[visit]['SURVEY'] == 'apogee' :
                    apogee_target1 |= allvisits[visit]['APOGEE_TARGET1'] 
                    apogee_target2 |= allvisits[visit]['APOGEE_TARGET2'] 
                    apogee_target3 |= allvisits[visit]['APOGEE_TARGET3'] 
                elif allvisits[visit]['SURVEY'].find('apogee2') >=0  :
                    apogee2_target1 |= allvisits[visit]['APOGEE_TARGET1'] 
                    apogee2_target2 |= allvisits[visit]['APOGEE_TARGET2'] 
                    apogee2_target3 |= allvisits[visit]['APOGEE_TARGET3'] 
                    try: apogee2_target4 |= allvisits[visit]['APOGEE_TARGET4'] 
                    except: pass
                allvisits[visit]['N_COMPONENTS']=out[1][i]['N_components']
                ncomponents=np.max([ncomponents,out[1][i]['N_components']])

            visits=np.array(visits)

            j = np.where(allfield['APOGEE_ID'] == files[-1][1])[0]
            allfield['RA'][j] = allvisits['RA'][visits[0]]
            allfield['DEC'][j] = allvisits['DEC'][visits[0]]
            for key in ['J','J_ERR','H','H_ERR','K','K_ERR'] :
                gd=np.where(abs(allvisits[key][visits])< 99.)[0]
                if len(gd) > 0 : allfield[key][j] = allvisits[key][visits[gd]].max()
            allfield['APOGEE_TARGET1'][j] = apogee_target1
            allfield['APOGEE_TARGET2'][j] = apogee_target2
            allfield['APOGEE_TARGET3'][j] = apogee_target3
            allfield['APOGEE2_TARGET1'][j] = apogee2_target1
            allfield['APOGEE2_TARGET2'][j] = apogee2_target2
            allfield['APOGEE2_TARGET3'][j] = apogee2_target3
            allfield['APOGEE2_TARGET4'][j] = apogee2_target4
            allfield['TARGFLAGS'][j] = (bitmask.targflags(apogee_target1,apogee_target2,apogee_target3,survey='apogee')+
                                        bitmask.targflags(apogee2_target1,apogee2_target2,apogee2_target3,survey='apogee2'))
            allfield['STARFLAG'][j] = starflag
            allfield['STARFLAGS'][j] = starmask.getname(starflag)
            allfield['ANDFLAG'][j] = andflag
            allfield['ANDFLAGS'][j] = starmask.getname(andflag)
            allfield['SNR'][j] = out[0][0]['medsnr']
            allfield['VHELIO_AVG'][j] = out[0][0]['vhelio']
            allfield['VSCATTER'][j] = out[0][0]['vscatter']
            allfield['VERR'][j] = out[0][0]['verr']
            allfield['RV_TEFF'][j] = out[0][0]['teff']
            allfield['RV_LOGG'][j] = out[0][0]['logg']
            allfield['RV_FEH'][j] = out[0][0]['feh']
            allfield['N_COMPONENTS'][j] = ncomponents

            # set up visit combination, removing visits with suspect RVs
            apogee_id=files[-1][1].decode() 
            gd = np.where(np.abs(allvisits[visits]['VHELIO']-allvisits[visits]['XCORR_VHELIO']) < 5) [0]
            if len(gd) < 1 : continue

            allv.append([allvisits[visits[gd]],load,(field,apogee_id,clobber)])

        except : pdb.set_trace()

    # do the visit combination
    if threads == 0 :
        output=[]
        for v in allv :
            output.append(dovisitcomb(v))
    else :
        pool = mp.Pool(threads)
        output = pool.map_async(dovisitcomb, allv).get()
        pool.close()
        pool.join()
    print('done visitcomb pool pool')

    #output file with allfield and allvisits
    hdulist=fits.HDUList()
    hdulist.append(fits.table_to_hdu(Table(allfield)))
    hdulist.append(fits.table_to_hdu(allvisits))
    hdulist.writeto(field+'/'+field+'_rv.fits',overwrite=True)

    # make web page
    if obj is not None : suffix='_obj'
    else : suffix=''
    if tweak: suffix=suffix+'_tweak'
    mkhtml(field,suffix=suffix)

    return allfield,allvisits

def dorv(visitfiles) :            
    """ do the rv jointfit from list of files
    """
    # last list elements has configuration variables in a tuple
    allvisit = visitfiles[0]
    load = visitfiles[1]
    field=visitfiles[-1][0]
    obj=visitfiles[-1][1].decode('UTF-8')
    clobber=visitfiles[-1][2]
    verbose=visitfiles[-1][3]
    tweak=visitfiles[-1][4]
    plot=visitfiles[-1][5]
    windows=visitfiles[-1][6]
    #rvrange=visitfiles[-1][7]
    if tweak: suffix='_tweak'
    else : suffix='_out'
    if os.path.exists(field+'/'+obj+suffix+'.pkl') and not clobber:
        print(obj,' already done')
        fp=open(field+'/'+obj+suffix+'.pkl','rb')
        try: 
            out=pickle.load(fp)
            fp.close()
            #if True : #len(out) == 5 :
            #    gout = gauss_decomp(out,phase='two',filt=True)
            #    dop_plot(field,obj,out[0],decomp=out[1])
            #    fp=open(field+'/'+obj+suffix+'.pkl','wb')
            #    pickle.dump([out,gout],fp)
            #    fp.close()
            #    return [out,gout]
            return out
        except: 
            print('error loading: ', obj+suffix+'.pkl')
            pass

    speclist=[]
    #for visitfile in visitfiles[0:-1] :
    lowsnr_visits=np.where(allvisit['SNR']<10)[0]
    if len(lowsnr_visits)/len(allvisit) > 0.1 :
        try :
            pixelmask=bitmask.PixelBitMask()
            apstar_bc=visitcomb(allvisit,bconly=True,load=load) 
            apstar_bc.mask(pixelmask.badval()|pixelmask.getval('SIG_SKYLINE'))
            spec=doppler.Spec1D(apstar_bc.flux,err=apstar_bc.err,bitmask=apstar_bc.bitmask,
                 mask=apstar_bc.mask,wave=apstar_bc.wave,lsfpars=np.array([0]),lsfsigma=apstar_bc.wave/22500/2.354,instrument='APOGEE',
                 filename=apstar_bc.filename)
            out= doppler.rv.jointfit([spec],verbose=verbose,plot=plot,tweak=tweak,maxvel=[-500,500])
            rvrange=[out[1][0]['vrel']-50,out[1][0]['vrel']+50]
        except :
            print('  BC jointfit failed')
            rvrange=[-500,500]
    else: rvrange=[-500,500]

    for i in range(len(allvisit)) :

        visitfile= load.allfile('Visit',plate=int(allvisit['PLATE'][i]),
                                 mjd=allvisit['MJD'][i],fiber=allvisit['FIBERID'][i])
        spec=doppler.read(visitfile)
        if windows is not None :
            for ichip in range(3) :
                mask = np.full_like(spec.mask[:,ichip],True)
                gd = []
                for window in windows :
                    gd.extend(np.where((spec.wave[:,ichip] > window[0]) & (spec.wave[:,ichip] < window[1]))[0])
                mask[gd] = False
                spec.mask[:,ichip] |= mask
                 
        if spec is not None : speclist.append(spec)
    try:
        print('running jointfit for :',obj)
        print('nvisits: ', len(speclist))
        out= doppler.rv.jointfit(speclist,verbose=verbose,plot=plot,saveplot=True,outdir=field+'/',tweak=tweak,maxvel=rvrange)
        print('running decomp for :',obj)
        gout = gauss_decomp(out[1],phase='two',filt=True)
        fp=open(field+'/'+obj+suffix+'.pkl','wb')
        pickle.dump([out,gout],fp)
        fp.close()
        print('running plots for :',obj)
        dop_plot(field,obj,out,decomp=gout)
    except KeyboardInterrupt : 
        raise
    except ValueError as err:
        print('Exception raised for: ', field, obj)
        print("ValueError: {0}".format(err))
        return
    except RuntimeError as err:
        print('Exception raised for: ', field, obj)
        print("Runtime error: {0}".format(err))
        return
    except :
        print('Exception raised for: ', field, obj)
        return

    return [out[0:2],gout]

def dovisitcomb(allv) :
    """ Routine to do visit combination in parallel
    """
    allvisits = allv[0]
    load = allv[1]
    field = allv[2][0]
    apogee_id = allv[2][1]
    clobber = allv[2][2]
    pixelmask=bitmask.PixelBitMask()

    # already done?
    if os.path.exists(field+'/'+apogee_id+'.pkl') and not clobber:
        print(apogee_id,' already done visitcomb')
        fp=open(field+'/'+apogee_id+'.pkl','rb')
        try: 
            out=pickle.load(fp)
            fp.close()
            return out
        except: 
            print('error loading: ', apogee_id+'.pkl')
            pass

    # do the combination
    apstar=visitcomb(allvisits,load=load,plot=False)

    # dump
    pickle.dump(apstar,open(field+'/'+apogee_id+'.pkl','wb'))

    # plot
    gd=np.where((apstar.bitmask & (pixelmask.badval()|pixelmask.getval('SIG_SKYLINE'))) == 0) [0]
    fig,ax=plots.multi(1,2,hspace=0.001,figsize=(24,6))
    med=np.nanmedian(apstar.flux)
    plots.plotl(ax[0],aspcap.apStarWave(),apstar.flux,color='k',yr=[0,2*med])
    ax[0].plot(aspcap.apStarWave()[gd],apstar.flux[gd],color='g')
    plots.plotl(ax[1],aspcap.apStarWave(),apstar.flux/apstar.err)
    fig.savefig(field+'/'+apogee_id+'.png')
    return apstar

def gaussian(amp, fwhm, mean):
    """ Gaussian as defined by gausspy
    """
    return lambda x: amp * np.exp(-4. * np.log(2) * (x-mean)**2 / fwhm**2)

import gausspy.gp as gp

def gauss_decomp(out,phase='one',alpha1=0.5,alpha2=1.5,thresh=[4,4],plot=None,filt=False) :
    """ Do Gaussian decomposition of CCF using gausspy

        Parameters:
        out : list of dictionaries for each frame, giving x_ccf, ccf, and ccferr
        phase : gausspy paramater
        alpha1 : gausspy parameter
        alpha2 : gausspy parameter for second set of gaussians if phase=='two'
        thresh : gausspy parameter
        plot (str) : if not None, do plot and use as root file name for plot
        filt (bool) : if true, apply filtering to remove components judged to be insignificant
    """
    g = gp.GaussianDecomposer()
    g.set('phase',phase)
    g.set('SNR_thresh',thresh)
    g.set('alpha1',alpha1)
    g.set('alpha2',alpha2)
    gout=[]
    if plot is not None : fig,ax=plots.multi(1,n,hspace=0.001,figsize=(6,2+n))
    for i,final in enumerate(out) :
        gd=np.where(np.isfinite(final['x_ccf']))[0]
        x=final['x_ccf'][gd]
        y=final['ccf'][gd]
        decomp=g.decompose(x,y,final['ccferr'][gd])
        n=decomp['N_components']
        if filt and n>0 :
            # remove components if they are within width of brighter component, or <0.25 peak 
            for j in range(1,n) :
                pars_j = decomp['best_fit_parameters'][j::n]
                for k in range(j) :
                    pars_k = decomp['best_fit_parameters'][k::n]
                    if pars_j[0]>pars_k[0] and (abs(pars_j[2]-pars_k[2])<pars_j[1]  or pars_k[0]<0.25*pars_j[0]): 
                        decomp['best_fit_parameters'][k] = 0
                        decomp['N_components'] -= 1
                    elif pars_k[0]>pars_j[0] and (abs(pars_j[2]-pars_k[2])<pars_k[1] or pars_j[0]<0.25*pars_k[0]) : 
                        decomp['best_fit_parameters'][j] = 0
                        decomp['N_components'] -= 1
                  
        gout.append(decomp)
        if plot is not None:
            plots.plotl(ax[i],final['x_ccf'],final['ccf'])
            ax[i].plot(final['x_ccf'],final['ccferr'],color='r')
            for j in range(n) :
                pars=gout[i]['best_fit_parameters'][j::n]
                ax[i].plot(x,gaussian(*pars)(x))
                print('pars: {:8.1f}{:8.1f}{:8.1f}'.format(*pars))
                if pars[0] > 0 : color='k'
                else : color='r'
                ax[i].text(0.1,0.8-j*0.1,'{:8.1f}{:8.1f}{:8.1f}'.format(*pars),transform=ax[i].transAxes,color=color)
            fig.savefig(plot+'_ccf.png')
    del g
    return gout

def dop_plot(field,obj,out,decomp=None) :
    """ RV diagnostic plots
    """
    matplotlib.use('Agg')
    n = len(out[2])
    #plot final spectra and final models
    # full spectrum
    fig,ax=plots.multi(1,n,hspace=0.001,figsize=(8,2+n))
    # continuum
    figc,axc=plots.multi(1,n,hspace=0.001,figsize=(8,2+n))
    # windows
    windows=[[15700,15780],[15850,16000],[16700,16930]]
    fig2,ax2=plots.multi(len(windows),n,hspace=0.001,wspace=0.001,figsize=(12,2+n))

    # loop over visitis
    for i,(mod,spec) in enumerate(zip(out[2],out[3])) :
        ax[i].plot(spec.wave,spec.flux,color='k')
        for iorder in range(3) :
            gd = np.where(~spec.mask[:,iorder])[0]
            ax[i].plot(spec.wave[gd,iorder],spec.flux[gd,iorder],color='g')
        ax[i].plot(mod.wave,mod.flux,color='r')
        ax[i].text(0.1,0.1,'{:d}'.format(spec.head['MJD5']),transform=ax[i].transAxes)
        for iwind,wind in enumerate(windows) :
            ax2[i,iwind].plot(spec.wave,spec.flux,color='k')
            for iorder in range(3) :
                gd = np.where(~spec.mask[:,iorder])[0]
                ax2[i,iwind].plot(spec.wave[gd,iorder],spec.flux[gd,iorder],color='g')
            ax2[i,iwind].plot(mod.wave,mod.flux,color='r')
            ax2[i,iwind].set_xlim(wind[0],wind[1])
            ax2[i,iwind].set_ylim(0.5,1.3)
            if iwind == 0 : ax2[i,iwind].text(0.1,0.1,'{:d}'.format(spec.head['MJD5']),transform=ax2[i,0].transAxes)
        axc[i].plot(spec.wave,spec.flux*spec.cont,color='k')
        axc[i].plot(spec.wave,spec.cont,color='g')
        axc[i].text(0.1,0.1,'{:d}'.format(spec.head['MJD5']),transform=axc[i].transAxes)
    fig.savefig(field+'/'+obj+'_spec.png')
    plt.close()
    fig2.savefig(field+'/'+obj+'_spec2.png')
    plt.close()
    figc.savefig(field+'/'+obj+'_cont.png')
    plt.close()

    # plot cross correlation functions with final model
    fig,ax=plots.multi(1,n,hspace=0.001,figsize=(6,2+n))
    vmed=np.median(out[1]['vrel'])
    for i,(final,spec) in enumerate(zip(out[1],out[3])) :
        ax[i].plot(final['x_ccf'],final['ccf'],color='k')
        ax[i].plot(final['x_ccf'],final['ccferr'],color='r')
        ax[i].plot([final['vrel'],final['vrel']],ax[i].get_ylim(),color='g',label='fit RV')
        ax[i].plot([final['xcorr_vrel'],final['xcorr_vrel']],ax[i].get_ylim(),color='r',label='xcorr RV')
        ax[i].text(0.1,0.9,'{:d}'.format(spec.head['MJD5']),transform=ax[i].transAxes)
        ax[i].set_xlim(vmed-200,vmed+200)
        ax[i].legend()
        if decomp is not None :
            n=decomp[i]['N_components']
            if n>0 : n=len(decomp[i]['best_fit_parameters'])//3
            x=final['x_ccf']
            for j in range(n) :
                pars=decomp[i]['best_fit_parameters'][j::n]
                ax[i].plot(x,gaussian(*pars)(x))
                print('pars: {:8.2f}{:8.1f}{:8.1f}'.format(*pars))
                if pars[0] > 0 : color='k'
                else : color='r'
                ax[i].text(0.1,0.8-j*0.1,'{:8.1f}{:8.1f}{:8.1f}'.format(*pars),transform=ax[i].transAxes,color=color)
    fig.savefig(field+'/'+obj+'_ccf.png')
    plt.close()

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

def mkhtml(field,suffix='') :
    """ Make web pages with tables/plots of RV output
        c.f., Doppler vs IDL
    """

    # get Doppler results
    dop = fits.open(field+'/'+field+'_rv.fits')
 
    # get IDL results
    r13 = apload.ApLoad(apred='r13')
    apfieldvisits = r13.apFieldVisits(field)[1].data
    apfield = r13.apField(field)[1].data

    # match
    i1,i2 = match.match(dop[2].data['FILE'],apfieldvisits['FILE'])
    print(len(dop[2].data),len(apfieldvisits),len(i1))
    fig,ax=plots.multi(1,2)
    ax[0].hist(dop[1].data['VHELIO_AVG'],bins=np.arange(-300,300,5),label='doppler',color='g',histtype='step')
    ax[0].hist(apfield['VHELIO_AVG'],bins=np.arange(-300,300,5),label='IDL',color='r',histtype='step')
    ax[0].legend()
    ax[1].hist(dop[1].data['VSCATTER'],bins=np.arange(0,1,0.02),label='doppler',color='g',histtype='step')
    ax[1].hist(apfield['VSCATTER'],bins=np.arange(0,1,0.02),label='IDL',color='r',histtype='step')
    ax[1].legend()
    fig.savefig(field+'/'+field+'_rvhist.png')

    # create HTML and loop over objects
    fp=open(field+'/'+field+suffix+'.html','w')
    fp.write('<HTML>\n')
    fp.write('<HEAD><script type=text/javascript src=../html/sorttable.js></script></head>')
    fp.write('<BODY>\n')
    fp.write('<H2> Field: {:s}</H2><p>\n'.format(field))
    fp.write('<A HREF={:s}_rvhist.png> <IMG SRC={:s}_rvhist.png> </A>'.format(field,field))
    
    fp.write('<TABLE BORDER=2 CLASS=sortable>\n')
    fp.write('<TR><TD>Obj<TD>Delta(VSCATTER)<TD>H<TD>Doppler RV_TEFF<TD>N_components<TD>RV plot<TD>Spectrum<TD>Spectrum windows<TD> continuum\n')
    for star in dop[1].data :
        obj=star['APOGEE_ID']
        print(obj)

        # get visits in Doppler allvisit table
        j=np.where(dop[2].data['APOGEE_ID'] == obj)[0]
        if len(j) == 0 : 
            print('missing {:s} in dop[2].data'.format(obj))
            continue

        # get object in apField
        try: k=np.where(apfield['APOGEE_ID'] == obj)[0][0]
        except: k=-1
        jj=np.where(apfieldvisits['APOGEE_ID'] == obj)[0]

        # star information
        if star['TARGFLAGS'].find('TELLURIC') >=0 :
            fp.write('<TR><TD bgcolor=lightblue>')
        else :
            fp.write('<TR><TD>')
        fp.write('{:s}'.format(obj))
        fp.write('(<A HREF="http://simbad.cfa.harvard.edu/simbad/sim-basic?Ident={:12.5f}%09{:12.5f}++&submit=SIMBAD+search"> SIMBAD </A>)<BR>'.format
                   (apfield['RA'][k],apfield['DEC'][k]))
        fp.write('H  = {:7.2f}<br>'.format(star['H']))
        fp.write('{:s}<br>'.format(star['TARGFLAGS']))
        fp.write('{:s}<br>'.format(star['STARFLAGS']))

        # average velocities
        fp.write('<TABLE BORDER=2>\n')
        fp.write('<TR><TD><TD>VHELIO_AVG<TD>VSCATTER<TD>TEFF<TD>LOGG<TD>[FE/H]\n')
        fp.write('<TR><TD>Doppler<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.0f}<TD>{:8.2f}<TD>{:8.2f}\n'.format(
                 star['VHELIO_AVG'],star['VSCATTER'],
                 star['RV_TEFF'],star['RV_LOGG'],star['RV_FEH']))
        fp.write('<TR><TD>Doppler Xcorr<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.0f}<TD>{:8.2f}<TD>{:8.2f}\n'.format(
                 np.median(dop[2].data['XCORR_VHELIO'][j]),
                 dop[2].data['XCORR_VHELIO'][j].std(ddof=1),
                 star['RV_TEFF'],star['RV_LOGG'],star['RV_FEH']))
        if k>=0 :
            gd = np.where(np.abs(apfieldvisits['VHELIO']) < 999)[0]
            fp.write('<TR><TD>IDL<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.0f}<TD>{:8.2f}<TD>{:8.2f}\n'.format(
                     apfield['VHELIO_AVG'][k],apfield['VSCATTER'][k],
                     apfield['RV_TEFF'][k],apfield['RV_LOGG'][k],apfield['RV_FEH'][k]))
        fp.write('</TABLE><br>')

        # individual visit velocities
        fp.write('<TABLE BORDER=2>')
        fp.write('<TR><TD>JD<TD>PLATE<TD>MJD<TD>FIBER<TD>S/N<TD>Doppler xcorr<TD> xcorr_err<TD>Doppler<TD>VERR<TD>IDL<TD>VERR,<TD>BC<TD>BC\n')
        for i in j :
            ii = np.where(apfieldvisits['FILE'] == dop[2].data['FILE'][i])[0][0]
            if np.isfinite(dop[2].data['VHELIO'][i]) == False :
                bgcolor='bgcolor=red'
            elif abs(dop[2].data['VHELIO'][i]-dop[2].data['XCORR_VHELIO'][i]) > 5 : 
                bgcolor='bgcolor=lightpink'
            elif abs(dop[2].data['VHELIO'][i]-dop[2].data['XCORR_VHELIO'][i]) > 0 : 
                bgcolor='bgcolor=#F4DEDE'
            else : bgcolor=''
            fp.write(('<TR {:s}> <TD> <A HREF={:s} TARGET="_obj"> {:12.3f}</A> <TD> {:s} <TD> {:5d} <TD> {:5d}'+
                     '<TD> {:8.1f} <TD> {:8.2f} <TD> {:8.2f} <TD> {:8.2f} <TD> {:8.2f} <TD> {:8.2f} <TD> {:8.2f}'+
                     '<TD> {:8.2f} <TD> {:8.2f} <TD>{:8.2f}\n').format(
                      bgcolor,
                      dop[2].data['FILE'][i].replace('.fits','_dopfit.png').replace('-r12-','-r13-'),
                      dop[2].data['JD'][i],dop[2].data['PLATE'][i],dop[2].data['MJD'][i],dop[2].data['FIBERID'][i],
                      dop[2].data['SNR'][i],
                      dop[2].data['XCORR_VHELIO'][i],dop[2].data['XCORR_VRELERR'][i],
                      dop[2].data['VHELIO'][i],dop[2].data['VRELERR'][i],
                      apfieldvisits['VHELIO'][ii],apfieldvisits['VRELERR'][ii],dop[2].data['ESTBC'][i],dop[2].data['BC'][i],apfieldvisits['BC'][ii]))
        fp.write('</TABLE>')

        # vscatter difference with IDL
        fp.write('<TD> {:8.2f}'.format(star['VSCATTER']-apfield['VSCATTER'][k]))
        fp.write('<TD> {:8.2f}'.format(star['H']))
        fp.write('<TD> {:8.2f}'.format(star['RV_TEFF']))
        fp.write('<TD> {:d}'.format(star['N_COMPONENTS']))

        # plot visit RVs
        vhelio=dop[2].data['VHELIO'][j]
        gd_dop = np.where(np.abs(vhelio - dop[2].data['XCORR_VHELIO'][j]) <= 5)[0]
        bd_dop = np.where(np.abs(vhelio - dop[2].data['XCORR_VHELIO'][j]) > 5)[0]
        vidl=apfieldvisits['VHELIO'][jj]
        gd = np.where(np.abs(vidl) < 999)[0]
        vmax=np.nanmax(np.append(vhelio,vidl[gd]))
        vmin=np.nanmin(np.append(vhelio,vidl[gd]))
        yr=[vmin-0.1*(vmax-vmin),vmax+0.1*(vmax-vmin)]
        try :
            fig,ax=plots.multi(1,1)
            if len(gd_dop) > 0 : plots.plotp(ax,dop[2].data['MJD'][j[gd_dop]],vhelio[gd_dop],size=15,color='g',yr=yr,label='Doppler')
            if len(bd_dop) > 0 : ax.scatter(dop[2].data['MJD'][j[bd_dop]],vhelio[bd_dop],s=15,
                                            facecolors='none',edgecolors='g',label='rejected Doppler')
            plots.plotp(ax,apfieldvisits['MJD'][jj[gd]],vidl[gd],size=15,color='r',yr=yr,label='IDL')
            ax.plot(ax.get_xlim(),[star['VHELIO_AVG'],star['VHELIO_AVG']],color='g')
            ax.plot(ax.get_xlim(),[apfield['VHELIO_AVG'][k],apfield['VHELIO_AVG'][k]],color='r')
            ax.legend()
            fig.savefig(field+'/'+obj+'_rv.png')
            plt.close()
        except KeyboardInterrupt: raise
        except :
             plt.close()
             pass

        # include plots
        fp.write('<TD><IMG SRC={:s}.png>\n'.format(obj))
        fp.write('<TD><IMG SRC={:s}_rv.png>\n'.format(obj))
        fp.write('<TD><A HREF={:s}_ccf.png> <IMG SRC={:s}_ccf.png></A>\n'.format(obj,obj))
        fp.write('<TD><A HREF={:s}_spec.png> <IMG SRC={:s}_spec.png></a>\n'.format(obj,obj))
        fp.write('<TD><A HREF={:s}_spec2.png> <IMG SRC={:s}_spec2.png></a>\n'.format(obj,obj))
        fp.write('<TD><A HREF={:s}_cont.png> <IMG SRC={:s}_cont.png></a>\n'.format(obj,obj))
    fp.close() 

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
            #fp.write('<TR><TD>IDL<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.0f}<TD>{:8.2f}<TD>{:8.2f}\n'.format(
            #     field[1].data['VHELIO_AVG'],field[1].data['VSCATTER'],
            #     field[1].data['RV_TEFF'],field[1].data['RV_LOGG'],field[1].data['RV_FEH']))
        plt.draw()
        plt.close()
        fig.savefig(outdir+'/'+star+'.png')
        fp.write('</TABLE>\n')
        fp.write('<TD><a HREF={:s}.png> <IMG SRC={:s}.png> </a>\n'.format(star,star))
    fp.write('</TABLE>')
    fp.close()

from apogee.aspcap import aspcap
from apogee.apred import wave
from apogee.apred import sincint

def visitcomb(allvisit,load=None, apred='r13',telescope='apo25m',nres=[5,4.25,3.5],bconly=False,plot=False) :
    """ Combine multiple visits with individual RVs to rest frame sum
    """

    if load is None : load = apload.ApLoad(apred=apred,telescope=telescope)
    cspeed = 2.99792458e5  # speed of light in km/s

    wnew=aspcap.apStarWave()  
    nwave=len(wnew)
    nvisit=len(allvisit)

    # initialize array for stack of interpolated spectra
    zeros = np.zeros([nvisit,nwave])
    izeros = np.zeros([nvisit,nwave],dtype=int)
    stack=apload.ApSpec(zeros,err=zeros.copy(),bitmask=izeros,cont=zeros.copy(),
                sky=zeros.copy(),skyerr=zeros.copy(),telluric=zeros.copy(),telerr=zeros.copy())

    # loop over each visit and interpolate to final wavelength grid
    if plot : fig,ax=plots.multi(1,2,hspace=0.001)
    for i,visit in enumerate(allvisit) :

        if bconly : vrel = -visit['BC']
        else : vrel = visit['VREL']

        # skip if we don't have an RV
        if np.isfinite(vrel) is False : continue

        # load the visit
        apvisit=load.apVisit(int(visit['PLATE']),visit['MJD'],visit['FIBERID'],load=True)
        pixelmask=bitmask.PixelBitMask()

        # rest-frame wavelengths transformed to this visit spectra
        w=aspcap.apStarWave()*(1.0+vrel/cspeed)
        print(vrel)

        # loop over the chips
        for chip in range(3) :

            # get the pixel values to interpolate to
            pix=wave.wave2pix(w,apvisit.wave[chip,:])
            gd=np.where(np.isfinite(pix))[0]

            # get a smoothed, filtered spectrum to use as replacement for bad values
            cont = gaussian_filter(median_filter(apvisit.flux[chip,:],[501],mode='reflect'),100)
            errcont = gaussian_filter(median_filter(apvisit.flux[chip,:],[501],mode='reflect'),100)
            bd = np.where(apvisit.bitmask[chip,:]&pixelmask.badval())[0]
            if len(bd) > 0 : 
                apvisit.flux[chip,bd] = cont[bd] 
                apvisit.err[chip,bd] = errcont[bd] 

            # load up quantity/error pairs for interpolation
            raw=[[apvisit.flux[chip,:],apvisit.err[chip,:]**2],
                 [apvisit.sky[chip,:],apvisit.skyerr[chip,:]**2],
                 [apvisit.telluric[chip,:],apvisit.telerr[chip,:]**2]]

            # load up individual mask bits
            for ibit,name in enumerate(pixelmask.name) :
                if name is not '' and len(np.where(apvisit.bitmask[chip,:]&2**ibit)[0]) > 0 :
                    raw.append([np.clip(apvisit.bitmask[chip,:]&2**ibit,None,1),None])

            # do the sinc interpolation
            out=sincint.sincint(pix[gd],nres[chip],raw)

            # from output flux, get continuum to remove, so that all spectra are
            #   on same scale. We'll later multiply in the median continuum
            flux = out[0][0]
            stack.cont[i,gd] = gaussian_filter(median_filter(flux,[501],mode='reflect'),100)

            # load interpolated spectra into output stack
            stack.flux[i,gd] = out[0][0] / stack.cont[i,gd]
            stack.err[i,gd] = out[0][1] / stack.cont[i,gd]
            stack.sky[i,gd] = out[1][0]
            stack.skyerr[i,gd] = out[1][1]
            stack.telluric[i,gd] = out[2][0]
            stack.telerr[i,gd] = out[2][1]
            # for mask, set bits where interpolated value is above some threshold
            #   defined for each mask bit
            iout=3
            for ibit,name in enumerate(pixelmask.name) :
                if name is not '' and len(np.where(apvisit.bitmask[chip,:]&2**ibit)[0]) > 0 :
                    j = np.where(np.abs(out[iout][0]) > pixelmask.maskcontrib[ibit])[0]
                    stack.bitmask[i,gd[j]] |= 2**ibit
                    iout+=1

        # increase uncertainties for persistence pixels
        bd = np.where((stack.bitmask[i,:]&pixelmask.getval('PERSIST_HIGH')) > 0)[0]
        if len(bd) > 0 : stack.err[i,bd] *= np.sqrt(5)
        bd = np.where(((stack.bitmask[i,:]&pixelmask.getval('PERSIST_HIGH')) == 0) &
                      ((stack.bitmask[i,:]&pixelmask.getval('PERSIST_MED')) > 0) )[0]
        if len(bd) > 0 : stack.err[i,bd] *= np.sqrt(4)
        bd = np.where(((stack.bitmask[i,:]&pixelmask.getval('PERSIST_HIGH')) == 0) &
                      ((stack.bitmask[i,:]&pixelmask.getval('PERSIST_MED')) == 0) &
                      ((stack.bitmask[i,:]&pixelmask.getval('PERSIST_LOW')) == 0) )[0]
        if len(bd) > 0 : stack.err[i,bd] *= np.sqrt(3)
        bd = np.where((stack.bitmask[i,:]&pixelmask.getval('SIG_SKYLINE')) > 0)[0]
        if len(bd) > 0 : stack.err[i,bd] *= np.sqrt(100)

        if plot :
            ax[0].plot(aspcap.apStarWave(),stack.flux[i,:])
            ax[1].plot(aspcap.apStarWave(),stack.flux[i,:]/stack.err[i,:])
            plt.draw()
            pdb.set_trace()

    # create final spectrum
    zeros = np.zeros(nwave)
    izeros = np.zeros(nwave,dtype=int)
    apstar=apload.ApSpec(zeros,err=zeros.copy(),bitmask=izeros,wave=aspcap.apStarWave(),
                sky=zeros.copy(),skyerr=zeros.copy(),telluric=zeros.copy(),telerr=zeros.copy())

    # pixel-by-pixel weighted average
    cont = np.median(stack.cont,axis=0)
    apstar.flux = np.sum(stack.flux/stack.err**2,axis=0)/np.sum(1./stack.err**2,axis=0) * cont
    apstar.err =  np.sqrt(1./np.sum(1./stack.err**2,axis=0)) * cont
    apstar.bitmask = np.bitwise_and.reduce(stack.bitmask,0)

    if plot : 
        ax[0].plot(aspcap.apStarWave(),apstar.flux,color='k')
        ax[1].plot(aspcap.apStarWave(),apstar.flux/apstar.err,color='k')
        plt.draw()
        pdb.set_trace()

    return apstar

def emission(spec) :
    """ Try to flag emission lines
    """

    for chip in range(3)  :
        cont = gaussian_filter(median_filter(spec.flux[chip,:],[501],mode='reflect'),100)
        emiss = np.where(((spec.flux[chip,:]/cont) > 1.10)  &
                         (spec.bitmask[chip,:]&(pixelmask.badval()|pixelmask.getval('SIG_SKYLINE'))) ) [0]
