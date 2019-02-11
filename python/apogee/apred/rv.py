import os
import copy
import glob
import pdb
import numpy as np
import matplotlib.pyplot as plt
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
    """ Concatenate set of apField files
    """
    # concatenate the structures
    all=struct.concat(files,verbose=verbose)

    # write out the file
    if out is not None:
        print('writing',out)
        struct.wrfits(all,out)

    return all


def vscat(a,fig=None,ls=None,marker='o') :
    """ Make histograms of VSCATTER for different Teff and min NVISITS
    """
    if fig == None : fig,ax=plots.multi(3,6,hspace=0.001,wspace=0.4)
    else : fig,ax=fig
    tbins=[3000,3500,4000,4500,5500,8000,30000] 
    snr = a['SNR']
    j=np.where(snr > 300) [0]
    snr[j] = 300
    for i in range(len(tbins)-1) :
        ax[i,0].text(0.9,0.9,'{:d}<=RV_TEFF<{:d}'.format(tbins[i],tbins[i+1]),ha='right',transform=ax[i,0].transAxes,fontsize=8)
        for j,nmin in enumerate([2,3,5,10]) :
            gd = np.where((a['RV_TEFF']>=tbins[i]) & (a['RV_TEFF']<tbins[i+1]) &
                           (a['NVISITS']>nmin) ) [0]
            print(tbins[i],tbins[i+1],nmin,len(gd))
            try :
                plots.plotc(ax[i,2],snr[gd],a['VSCATTER'][gd],a['RV_FEH'][gd],marker=marker,xr=[0,310],yr=[0,1],xt='S/N',yt='VSCATTER')
                ax[i,0].hist(a['VSCATTER'][gd],bins=np.arange(0,1,0.01),ls=ls,histtype='step',color=colors[j])
                ax[i,0].set_xlabel('VSCATTER')
                ax[i,1].hist(a['VSCATTER'][gd],bins=np.arange(0,1,0.01),histtype='step',cumulative=True,normed=True,ls=ls,color=colors[j])
                ax[i,1].set_xlabel('VSCATTER')
            except : pass

    return fig,ax

def apolco(a,minfeh=-3,out=None) :
    """  VSCATTER histograms for APO vs LCO
    """
    gd=np.where((a['TELESCOPE'] == 'apo25m') & (a['RV_FEH']>minfeh) )[0]
    fig=vscat(a[gd],marker='o')
    gd=np.where((a['TELESCOPE'] == 'lco25m') & (a['RV_FEH']>minfeh) )[0]
    vscat(a[gd],fig=fig,ls=':',marker='+')
    if out is not None : 
        fig[0].savefig(out+'.png')
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
        fig,ax=plots.multi(2,2,hspace=0.001,wspace=0.4)
        plots.plotp(ax[0,0],a['SNR'],a['VHELIO_AVG']-b['VHELIO_AVG'],yr=[-3,3],yt=r'$\Delta$ VHELIO_AVG')
        plots.plotp(ax[0,1],a['SNR'],a['VHELIO_AVG']-b['VHELIO_AVG'],yr=[-50,50],yt=r'$\Delta$ VHELIO_AVG')
        plots.plotp(ax[1,0],a['SNR'],a['VSCATTER']-b['VSCATTER'],yr=[-0.5,0.5],yt=r'$\Delta$ VSCATTER')
        plots.plotp(ax[1,1],a['SNR'],a['VSCATTER']-b['VSCATTER'],yr=[-5,5],yt=r'$\Delta$ VSCATTER')
        if out is not None : 
            fig.savefig(out+'_2.png')
            plt.close()

    return a,b

def field(name,dr14=False,dir='./') :
    """ look at a single field
    """
    all=struct.concat([dir+'/apVisitSum*.fits'])
    alldr14=struct.concat([os.environ['APOGEE_REDUX']+'/r8/fields/apo25m/4162//apVisitSum*'])
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
    ax[0,0].hist(vscat,bins=np.arange(0.01,1,0.01),histtype='step',cumulative=True,normed=True,color='b')
    ax[0,0].hist(dr14vscat,bins=np.arange(0.01,1,0.01),histtype='step',cumulative=True,normed=True,color='r')
    gd=np.where(verr < 0.2)[0]
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

def all() :
    """ Do a series of RV comparisons for DR14 and proto DR16
    """
    grid=[]
    xtit=[]
    a=fits.open('allField.fits')[1].data
    load=apload.ApLoad(dr='dr14')
    b=load.allStar()[1].data
    comp(a,b,domatch=False,out='plots/dr14all')
    grid.append(['dr14all_1.png',''])
    xtit.append('all stars: DR14 (dotted) and test DR16 (solid)')
    comp(a,b,domatch=True,out='plots/dr14match')
    grid.append(['dr14match_1.png','dr14match_2.png'])
    xtit.append('same stars: DR14 (dotted) and test DR16 (solid)')
    apolco(a,out='plots/apolco')
    grid.append(['apolco.png',''])
    xtit.append('testDR16 : APO (solid) and LCO (dotted)')
    apolco(a,out='plots/apolco_nolowz',minfeh=-0.6)
    grid.append(['apolco_nolowz.png',''])
    xtit.append('testDR16, no low [Fe/H]: APO (solid) and LCO (dotted)')
    html.htmltab(grid,ytitle=xtit,file='plots/rv.html')

def visitspec(load,plate,mjd,fiber,gridfile='apg_rvsynthgrid',apstar=False) :

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

