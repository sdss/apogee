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

def doppler_rv(field,telescope='apo25m',apred='r13',obj=None,threads=8,maxvisit=500,maxobj=9999,snmin=5,
               clobber=False,verbose=False,tweak=False,plot=False) :
    """ Run DOPPLER RVs for a field
    """ 
    
    # get all the VisitSum files for this field and concatenate them
    files=glob.glob(os.environ['APOGEE_REDUX']+'/'+apred+'/visit/'+telescope+'/'+field+'/apVisitSum*')
    allvisits=struct.concat(files)
    starmask=bitmask.StarBitMask()
    gd=np.where(((allvisits['STARFLAG'] & starmask.badval()) == 0) & (allvisits['SNR'] > snmin) )[0]
    print(len(allvisits),len(gd))
    allvisits=Table(allvisits[gd])

    # output directory
    try: os.mkdir(field)
    except FileExistsError: pass

    # get all unique (or requested) objects
    if obj is None :
        allobj=set(allvisits['APOGEE_ID'][0:maxobj])
    else :
        allobj = obj
    # loop over requested objects
    """
    a={file: file_basename(starfile), apogee_id: objname, telescope: apstr.telescope, $
       location_id: long(locid), field: apstr.field, $
       j: apstr.j, j_err: apstr.j_err, h: apstr.h, h_err: apstr.h_err, k: apstr.k, k_err: apstr.k_err, $
       ra: apstr.ra, dec: apstr.dec, glon: apstr.glon, glat: apstr.glat, $
       ak_targ: float(apstr.ak_targ), ak_targ_method: apstr.ak_targ_method, $
       ak_wise: float(apstr.ak_wise), sfd_ebv: float(apstr.sfd_ebv),$
       apogee_target1: apogee_target1, apogee_target2: apogee_target2, apogee_target3: apogee_target3, $
       apogee2_target1: apogee2_target1, apogee2_target2: apogee2_target2, apogee2_target3: apogee2_target3, $
       targflags: targflags, survey: survey, programname: allvisits[obj[j[0]]].programname, $
       ninst: [n1,n2,n3],$
       nvisits: apstr.nvisits, combtype:apstr.combtype, commiss: commiss, $
       snr: float(apstr.snr), $
       starflag: apstr.starflag, starflags: starflag(apstr.starflag), $
       andflag: apstr.andflag, andflags: starflag(apstr.andflag), $
       vhelio_avg: float(apstr.vhelio), vscatter: float(apstr.vscatter),$
       verr: float(apstr.verr), verr_med: float(apstr.verr_med),$
       obsvhelio_avg: float(apstr.obsvhelio),  obsvscatter: float(apstr.obsvscatter),$
       obsverr: float(apstr.obsverr), obsverr_med: float(apstr.obsverr_med),$
       synthvhelio_avg: float(apstr.synthvhelio),  synthvscatter: float(apstr.synthvscatter),$
       synthverr: float(apstr.synthverr), synthverr_med: float(apstr.synthverr_med),$
       rv_teff: float(apstr.rv_teff), rv_logg: float(apstr.rv_logg), rv_feh: float(apstr.rv_feh),$
       rv_alpha: float(apstr.rv_alpha), rv_carb: float(apstr.rv_carb), $
       rv_ccfwhm: float(apstr.ccpfwhm), rv_autofwhm: float(apstr.autofwhm), synthscatter: float(apstr.synthscatter),$
       ;binary: apstr.rv.binary, $
       stablerv_chi2: apstr.rv.stablerv_chi2, $
       stablerv_rchi2: apstr.rv.stablerv_rchi2, chi2_threshold: apstr.rv.chi2_threshold, $
       stablerv_chi2_prob: apstr.rv.stablerv_chi2_prob}
    b={spec: apstr.spec[*,0], err: apstr.err[*,0], mask: apstr.mask[*,0]}
    c={wave: apstr.wavelength[*,0]}
    """
    fieldtype = np.dtype([('FILE','S64'),('APOGEE_ID','S20'),('TELESCOPE','S6'),('LOCATION_ID',int),('FIELD','S20'),
                          ('J',float),('J_ERR',float),('H',float),('H_ERR',float),('K',float),('K_ERR',float),
                          ('RA',float),('DEC',float),('GLON',float),('GLAT',float),
                          ('APOGEE_TARGET1',int),('APOGEE_TARGET2',int),('APOGEE_TARGET3',int),
                          ('APOGEE2_TARGET1',int),('APOGEE2_TARGET2',int),('APOGEE2_TARGET3',int),('APOGEE2_TARGET4',int),
                          ('TARGFLAGS','S132'),
                          ('SNR',float),('STARFLAG',float),('STARFLAGS','S132'),('ANDFLAG',float),('ANDFLAGS','S132'),
                          ('VHELIO_AVG',float),('VSCATTER',float),('VERR',float),
                          ('RV_TEFF',float),('RV_LOGG',float),('RV_FEH',float),
                          ('N_COMPONENTS',int)
                         ])
   
    allfield = np.zeros(len(allobj),dtype=fieldtype)
    allfield['TELESCOPE'] = telescope
    allfield['FIELD'] = field
    allfiles=[]
    allv=[]
    load=apload.ApLoad(apred=apred)
    nobj=0
    nvisit=0
    #allvisits['FWHM'] = np.zeros([len(allvisits),3])

    # loop over all objects, building up list of [(field,obj,clobber),filenames....] to 
    #   pass to dorv()
    for iobj,obj in enumerate(sorted(allobj)) :
        if type(obj) is str : obj=obj.encode()
        allfield['APOGEE_ID'][iobj] = obj
        visits=np.where(allvisits['APOGEE_ID'] == obj)[0]
        print('object: {:}  nvisits: {:d}'.format(obj,len(visits)))
        nobj+=1
        nvisit+=len(visits)

        specfiles=[(field,obj,clobber,verbose,tweak,plot)]
        for i,visit in enumerate(visits) :
            if i < maxvisit :
                # accumulate list of files
                visitfile= load.allfile('Visit',plate=int(allvisits['PLATE'][visit]),
                                        mjd=allvisits['MJD'][visit],fiber=allvisits['FIBERID'][visit])
                specfiles.append(visitfile)
    #            spec=doppler.read(visitfile)
    #            fwhm=[]
    #            for chip in range(3) : 
    #                try: fwhm.append(spec.lsf.fwhm(order=chip).min())
    #                except: fwhm.append(0.)
    #            allvisits['FWHM'][visit] = np.array(fwhm)
        allfiles.append(specfiles)
        allv.append(visits)
       
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
    for col in ['VTYPE','VREL','VRELERR','VHELIO'] :
        allvisits.rename_column(col,'est'+col)
        if col == 'VTYPE' : allvisits[col] = 0
        else : allvisits[col] = np.nan
    for col in ['XCORR_VREL','XCORR_VRELERR','XCORR_VHELIO'] :
        allvisits[col] = np.nan
    allvisits['N_COMPONENTS'] = -1

    # load up the individual visit RV information
    for out,files in zip(output,allfiles) :
        starflag, andflag = 0, 0
        apogee_target1, apogee_target2, apogee_target3 = 0, 0, 0
        apogee2_target1, apogee2_target2, apogee2_target3, apogee2_target4 = 0, 0, 0, 0
        if out is not None :
            allv=[]
            ncomponents=0
            for i,v in enumerate(out[0][1]) :
                # match by filename components in case there was an error reading in doppler
                name=os.path.basename(v['filename']).replace('.fits','').split('-')
                visit = np.where( (np.char.strip(allvisits['PLATE']).astype(str) == name[-3]) &
                                  (allvisits['MJD'] == int(name[-2])) &
                                  (allvisits['FIBERID'] == int(name[-1])) )[0][0]
                allv.append(visit)
                allvisits[visit]['VREL']=v['vrel']
                allvisits[visit]['VRELERR']=v['vrelerr']
                allvisits[visit]['VHELIO']=v['vhelio']
                allvisits[visit]['XCORR_VREL']=v['xcorr_vrel']
                allvisits[visit]['XCORR_VRELERR']=v['xcorr_vrelerr']
                allvisits[visit]['XCORR_VHELIO']=v['xcorr_vhelio']
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

            allv=np.array(allv)
            j = np.where(allfield['APOGEE_ID'] == files[0][1])[0]
            allfield['RA'][j] = allvisits['RA'][allv[0]]
            allfield['DEC'][j] = allvisits['DEC'][allv[0]]
            for key in ['J','J_ERR','H','H_ERR','K','K_ERR'] :
                gd=np.where(abs(allvisits[key][allv])< 99.)[0]
                if len(gd) > 0 : allfield[key][j] = allvisits[key][allv[gd]].max()
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

    #output file with allfield and allvisits
    hdulist=fits.HDUList()
    hdulist.append(fits.table_to_hdu(Table(allfield)))
    hdulist.append(fits.table_to_hdu(allvisits))
    hdulist.writeto(field+'/'+field+'_rv.fits',overwrite=True)

    # make web page
    mkhtml(field)

    return allfield,allvisits

def dorv(visitfiles) :            
    """ do the rv jointfit from list of files
    """
    field=visitfiles[0][0]
    obj=visitfiles[0][1].decode('UTF-8')
    clobber=visitfiles[0][2]
    verbose=visitfiles[0][3]
    tweak=visitfiles[0][4]
    plot=visitfiles[0][5]
    if os.path.exists(field+'/'+obj+'_out.pkl') and not clobber:
        print(obj,' already done')
        fp=open(field+'/'+obj+'_out.pkl','rb')
        try: 
            out=pickle.load(fp)
            fp.close()
            if len(out) == 5 :
                gout = gauss_decomp(out)
                dop_plot(field,obj,out,decomp=gout)
                fp=open(field+'/'+obj+'_out.pkl','wb')
                pickle.dump([out,gout],fp)
                fp.close()
                return [out,gout]
            return out
        except: 
            print('error loading: ', obj+'_out.pkl')
            pass

    speclist=[]
    for visitfile in visitfiles[1:] :
        spec=doppler.read(visitfile)
        if spec is not None : speclist.append(spec)
    try:
        print('running jointfit for :',obj)
        print('nvisits: ', len(speclist))
        out= doppler.rv.jointfit(speclist,verbose=verbose,plot=plot,saveplot=True,outdir=field+'/',tweak=tweak)
        gout = gauss_decomp(out)
        fp = open(field+'/'+obj+'_rv.txt','w')
        fp.write('{:s}  {:d} {:8.1f}'.format(obj,len(speclist),out[4]))
        fp.close()
        fp=open(field+'/'+obj+'_out.pkl','wb')
        pickle.dump([out,gout],fp)
        fp.close()
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

    return out

def gaussian(amp, fwhm, mean):
    return lambda x: amp * np.exp(-4. * np.log(2) * (x-mean)**2 / fwhm**2)

import gausspy.gp as gp

def gauss_decomp(out) :
    """ Do Gaussian decompoistion of CCF
    """
    g = gp.GaussianDecomposer()
    g.set('phase','one')
    g.set('SNR_thresh',[4,4])
    g.set('alpha1',0.5)
    g.set('alpha2',1.5)
    gout=[]
    for i,(final,spec) in enumerate(zip(out[1],out[3])) :
        x=final['x_ccf']
        y=final['ccf']
        gout.append(g.decompose(x,y,final['ccferr']))
    del g
    return gout

def dop_plot(field,obj,out,decomp=None) :
    """ RV diagnostic plots
    """
    n = len(out[2])
    #plot final spectra and final models
    # full spectrum
    fig,ax=plots.multi(1,n,hspace=0.001,figsize=(8,2+n))
    # two windows
    fig2,ax2=plots.multi(2,n,hspace=0.001,wspace=0.001,figsize=(8,2+n))
    for i,(mod,spec) in enumerate(zip(out[2],out[3])) :
        ax[i].plot(spec.wave,spec.flux,color='k')
        ax[i].plot(mod.wave,mod.flux,color='r')
        ax[i].text(0.1,0.1,'{:d}'.format(spec.head['MJD5']),transform=ax[i].transAxes)
        ax2[i,0].plot(spec.wave,spec.flux,color='k')
        ax2[i,0].plot(mod.wave,mod.flux,color='r')
        ax2[i,1].plot(spec.wave,spec.flux,color='k')
        ax2[i,1].plot(mod.wave,mod.flux,color='r')
        ax2[i,0].set_xlim(15700,15800)
        ax2[i,1].set_xlim(16700,16930)
        ax2[i,0].text(0.1,0.1,'{:d}'.format(spec.head['MJD5']),transform=ax[i].transAxes)
    fig.savefig(field+'/'+obj+'_spec.png')
    fig2.savefig(field+'/'+obj+'_spec2.png')
    plt.close()
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
        ax[i].set_xlim(vmed-100,vmed+100)
        ax[i].legend()
        if decomp is not None :
            n=decomp[i]['N_components']
            x=final['x_ccf']
            for j in range(n) :
                pars=decomp[i]['best_fit_parameters'][j::n]
                ax[i].plot(x,gaussian(*pars)(x))
                print('pars: {:8.1f}{:8.1f}{:8.1f}'.format(*pars))
                ax[i].text(0.1,0.8-j*0.1,'{:8.1f}{:8.1f}{:8.1f}'.format(*pars),transform=ax[i].transAxes)
    fig.savefig(field+'/'+obj+'_ccf.png')
    plt.close()

from scipy.signal import convolve
def dop_comp(field) :
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

def mkhtml(field) :
    """ Make web pages with tables/plots of RV output
        c.f., Doppler vs IDL
    """

    matplotlib.use('Agg')
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
    fp=open(field+'/'+field+'.html','w')
    fp.write('<HTML>\n')
    fp.write('<HEAD><script type=text/javascript src=../html/sorttable.js></script></head>')
    fp.write('<BODY>\n')
    fp.write('<H2> Field: {:s}</H2><p>\n'.format(field))
    fp.write('<A HREF={:s}_rvhist.png> <IMG SRC={:s}_rvhist.png> </A>'.format(field,field))
    
    fp.write('<TABLE BORDER=2 CLASS=sortable>\n')
    fp.write('<TR><TD>Obj<TD>Delta(VSCATTER)<TD>H<TD>Doppler RV_TEFF<TD>RV plot\n')
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

        # average veloicities
        fp.write('<TABLE BORDER=2>\n')
        fp.write('<TR><TD><TD>VHELIO_AVG<TD>VSCATTER<TD>VSIGMA<TD>TEFF<TD>LOGG<TD>[FE/H]\n')
        vsig=dop[2].data['VHELIO'][j].std(ddof=1)
        fp.write('<TR><TD>Doppler<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.0f}<TD>{:8.2f}<TD>{:8.2f}\n'.format(
                 star['VHELIO_AVG'],star['VSCATTER'],vsig,
                 star['RV_TEFF'],star['RV_LOGG'],star['RV_FEH']))
        fp.write('<TR><TD>Doppler Xcorr<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.0f}<TD>{:8.2f}<TD>{:8.2f}\n'.format(
                 np.median(dop[2].data['XCORR_VHELIO'][j]),
                 dop[2].data['XCORR_VHELIO'][j].std(ddof=1),dop[2].data['XCORR_VHELIO'][j].std(ddof=1),
                 star['RV_TEFF'],star['RV_LOGG'],star['RV_FEH']))
        if k>=0 :
            gd = np.where(np.abs(apfieldvisits['VHELIO']) < 999)[0]
            fp.write('<TR><TD>IDL<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.2f}<TD>{:8.0f}<TD>{:8.2f}<TD>{:8.2f}\n'.format(
                     apfield['VHELIO_AVG'][k],apfield['VSCATTER'][k],apfieldvisits['VHELIO'][gd].std(ddof=1),
                     apfield['RV_TEFF'][k],apfield['RV_LOGG'][k],apfield['RV_FEH'][k]))
        fp.write('</TABLE><br>')

        # individual visit velocities
        fp.write('<TABLE BORDER=2>')
        fp.write('<TR><TD>JD<TD>PLATE<TD>MJD<TD>FIBER<TD>S/N<TD>Doppler xcorr<TD> xcorr_err<TD>Doppler<TD>VERR<TD>IDL<TD>VERR\n')
        for i in j :
            ii = np.where(apfieldvisits['FILE'] == dop[2].data['FILE'][i])[0][0]
            fp.write(('<TR> <TD> <A HREF={:s} TARGET="_obj"> {:12.3f}</A> <TD> {:s} <TD> {:5d} <TD> {:5d}'+
                     '<TD> {:8.1f} <TD> {:8.2f} <TD> {:8.2f} <TD> {:8.2f} <TD> {:8.2f} <TD> {:8.2f} <TD> {:8.2f}\n').format(
                      dop[2].data['FILE'][i].replace('.fits','_dopfit.png').replace('-r12-','-r13-'),
                      dop[2].data['JD'][i],dop[2].data['PLATE'][i],dop[2].data['MJD'][i],dop[2].data['FIBERID'][i],
                      dop[2].data['SNR'][i],
                      dop[2].data['XCORR_VHELIO'][i],dop[2].data['XCORR_VRELERR'][i],
                      dop[2].data['VHELIO'][i],dop[2].data['VRELERR'][i],
                      apfieldvisits['VHELIO'][ii],apfieldvisits['VRELERR'][ii]))
        fp.write('</TABLE>')

        # vscatter difference with IDL
        fp.write('<TD> {:8.2f}'.format(star['VSCATTER']-apfield['VSCATTER'][k]))
        fp.write('<TD> {:8.2f}'.format(star['H']))
        fp.write('<TD> {:8.2f}'.format(star['RV_TEFF']))
        fp.write('<TD> {:d}'.format(star['N_COMPONENTS']))

        # plot visit RVs
        vhelio=dop[2].data['VHELIO'][j]
        vidl=apfieldvisits['VHELIO'][jj]
        gd = np.where(np.abs(vidl) < 999)[0]
        vmax=np.nanmax(np.append(vhelio,vidl[gd]))
        vmin=np.nanmin(np.append(vhelio,vidl[gd]))
        yr=[vmin-0.1*(vmax-vmin),vmax+0.1*(vmax-vmin)]
        try :
            fig,ax=plots.multi(1,1)
            plots.plotp(ax,dop[2].data['MJD'][j],vhelio,size=15,color='g',yr=yr,label='Doppler')
            plots.plotp(ax,apfieldvisits['MJD'][jj[gd]],vidl[gd],size=15,color='r',yr=yr,label='IDL')
            ax.plot(ax.get_xlim(),[star['VHELIO_AVG'],star['VHELIO_AVG']],color='g')
            ax.plot(ax.get_xlim(),[apfield['VHELIO_AVG'][k],apfield['VHELIO_AVG'][k]],color='r')
            ax.legend()
            fig.savefig(field+'/'+obj+'.png')
            plt.close()
        except KeyboardInterrupt: raise
        except :
             pdb.set_trace()
             pass

        # include plots
        fp.write('<TD><IMG SRC={:s}.png>\n'.format(obj))
        fp.write('<TD><A HREF={:s}_ccf.png> <IMG SRC={:s}_ccf.png></A>\n'.format(obj,obj))
        fp.write('<TD><A HREF={:s}_spec.png> <IMG SRC={:s}_spec.png></a>\n'.format(obj,obj))
        fp.write('<TD><A HREF={:s}_spec2.png> <IMG SRC={:s}_spec2.png></a>\n'.format(obj,obj))
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
