# routines for assessing RVs from pipeline

import os
import copy
import glob
import pdb
import numpy as np
import matplotlib.pyplot as plt
import esutil
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

def doppler_rv(field,telescope='apo25m',apred='r13',obj=None,threads=8,maxvisit=500,maxobj=9999) :
    """ Run DOPPLER RVs for a field
    """ 

    # get all the VisitSum files for this field and concatenate them
    files=glob.glob(os.environ['APOGEE_REDUX']+'/'+apred+'/visit/'+telescope+'/'+field+'/apVisitSum*')
    allvisits=struct.concat(files)

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
    fieldtype = np.dtype([('file','S64'),('apogee_id','S20'),('telescope',str),('location_id',int),('field',str),
                          ('snr',float),
                          ('vhelio_avg',float),('vscatter',float),('verr',float),
                          ('rv_teff',float),('rv_logg',float),('rv_feh',float) 
                         ])
   
                 #     ('medsnr',float),('totsnr',float),('vhelio',float),('vscatter',float),('verr',float),
                 #     ('teff',float),('tefferr',float),('logg',float),('loggerr',float),('feh',float),
                 #     ('feherr',float),('chisq',float)])
    allfield = np.zeros(len(allobj),dtype=fieldtype)
    allfield['telescope'] = telescope
    allfield['field'] = field
    allfiles=[]
    allv=[]
    load=apload.ApLoad(apred=apred)
    nobj=0
    nvisit=0
    for iobj,obj in enumerate(sorted(allobj)) :
        if type(obj) is str : obj=obj.encode()
        allfield['apogee_id'][iobj] = obj
        visits=np.where(allvisits['APOGEE_ID'] == obj)[0]
        print('object: {:}  nvisits: {:d}'.format(obj,len(visits)))
        nobj+=1
        nvisit+=len(visits)

        specfiles=[obj]
        for i,visit in enumerate(visits) :
            if i < maxvisit :
                # accumulate list of files
                visitfile= load.allfile('Visit',plate=int(allvisits['PLATE'][visit]),mjd=allvisits['MJD'][visit],fiber=allvisits['FIBERID'][visit])
                specfiles.append(visitfile)
        allfiles.append(specfiles)
        allv.append(visits)
       
    print('total objects: ', nobj, ' total visits: ', nvisit) 
    if threads == 0 :
        output=[]
        for speclist in allfiles :
            output.append(dorv(speclist))
    else :
        pool = mp.Pool(threads)
        output = pool.map_async(dorv, allfiles).get()
        pool.close()
        pool.join()
    print('done pool')

    allvisits=Table(allvisits)
    # rename old visit RV tags and initialize new ones
    for col in ['VTYPE','VREL','VRELERR','VHELIO'] :
        allvisits.rename_column(col,'est'+col)
        if col == 'VTYPE' : allvisits[col] = 0
        else : allvisits[col] = np.nan

    # load up the individual visit RV information
    for out,files in zip(output,allfiles) :
        for v in out[1] :
            # match by filename components in case there was an error reading in doppler
            name=os.path.basename(v['filename']).replace('.fits','').split('-')
            visit = np.where( (np.char.strip(allvisits['PLATE']).astype(str) == name[-3]) &
                              (allvisits['MJD'] == int(name[-2])) &
                              (allvisits['FIBERID'] == int(name[-1])) )[0][0]
            allvisits[visit]['VREL']=v['vrel']
            allvisits[visit]['VRELERR']=v['vrelerr']
            allvisits[visit]['VHELIO']=v['vhelio']
        j = np.where(allfield['apogee_id'] == files[0])[0]
        allfield['snr'][j] = out[0]['medsnr']
        allfield['vhelio_avg'][j] = out[0]['vhelio']
        allfield['vscatter'][j] = out[0]['vscatter']
        allfield['verr'][j] = out[0]['verr']
        allfield['rv_teff'][j] = out[0]['teff']
        allfield['rv_logg'][j] = out[0]['logg']
        allfield['rv_feh'][j] = out[0]['feh']

        #('medsnr',float),('totsnr',float),('vhelio',float),('vscatter',float),('verr',float),
        #('teff',float),('tefferr',float),('logg',float),('loggerr',float),('feh',float),
        #('feherr',float),('chisq',float)])

    return allfield,allvisits

def dorv(visitfiles) :            
    """ do the rv jointfit from list of files
    """
    speclist=[]
    print(visitfiles)
    for visitfile in visitfiles[1:] :
        spec=doppler.read(visitfile)
        if spec is not None : speclist.append(spec)
    out= doppler.rv.jointfit(speclist,verbose=True)
    obj=visitfiles[0].decode('UTF-8')
    fp = open(obj+'_rv.txt','w')
    fp.write('{:s}  {:d} {:8.1f}'.format(obj,len(speclist),out[4]))
    fp.close()

    return out

#        # do the RVs
#        try :
#            out=doppler.rv.jointfit(speclist,verbose=True)
#            allout.append(out)
#        except KeyboardInterrupt:
#            pass
#        except:
#            print('error in jointfit')
    return allout

