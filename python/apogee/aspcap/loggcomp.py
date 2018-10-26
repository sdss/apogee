# routines for comparing gravities with asteroseismic sample

from apogee.utils import apload
from apogee.utils import apselect
from astropy.io import fits
#from holtz.gal import isochrones
#from holtz.gal import stars
from tools import match
from tools import plots
from tools import fit
from apogee.utils import bitmask
from apogee.aspcap import cal
import pdb
import matplotlib.pyplot as plt
import numpy as np
import os
import astropy

def rcrgb(allstar,apokasc='APOKASC_cat_v3.6.0.fits',logg='LOGG_SYD_SCALING',rclim=np.array([2.38,3.5]),out='rcrgbsep') :
    '''
    asteroseismic log g comparisons for input allStar structure
    '''

    gd=apselect.select(allstar,badval=['STAR_BAD'],logg=[0,3.8],teff=[3500,5500],raw=True)
    allstar=allstar[gd]

    # match ASPCAP with APOKASC, and get RC/RGB stars
    apokasc=fits.open(os.environ['APOGEE_DIR']+'/data/apokasc/'+apokasc)[1].data
    i1,i2=match.match(allstar['APOGEE_ID'],apokasc['2MASS_ID'])
    rgb=np.where(apokasc['CONS_EVSTATES'][i2] == 'RGB')[0]
    rc=np.where(apokasc['CONS_EVSTATES'][i2] == 'RC')[0]
    notrc=np.where(apokasc['CONS_EVSTATES'][i2] != 'RC')[0]
    rc2=np.where((apokasc['CONS_EVSTATES'][i2] == '2CL') & (apokasc[logg][i2] > -1))[0]
    rc2cl=np.where((apokasc['CONS_EVSTATES'][i2] == 'RC/2CL') & (apokasc[logg][i2] > -1))[0]
    rcall=np.append(rc,rc2)
    rcall=np.append(rcall,rc2cl)
    rc=i1[rc]
    rgb=i1[rgb]
    rc2=i1[rc2]
    rc2cl=i1[rc2cl]
    rcall=i1[rcall]

    # 2D fit for RGB Teff as f(log g, [M/H])
    fig,ax=plots.multi(3,2,figsize=(12,8))
    fit2d=fit.fit2d(allstar['FPARAM'][rgb,1]-2.5,allstar['FPARAM'][rgb,3],allstar['FPARAM'][rgb,0],plot=ax[0,0],xt='log g - 2.5',yt='[M/H]',zt='Teff')
    rgbfit=fit2d.parameters

    #histogram of RC logg values
    ax[1,0].hist(allstar['FPARAM'][rc,1],color='b',bins=np.arange(0,5,0.05),log=True)
    ax[1,0].set_xlim(1.5,3.5)
    ax[1,0].set_ylim(0.1,1000)
    ax[1,0].set_xlabel('log g')
    #ax[1,0].hist(allstar['FPARAM'][rgb,1],color='r',bins=np.arange(0,5,0.05))
    print('RC min log g: ',allstar['FPARAM'][rc,1].min())
    print('RC max log g: ',allstar['FPARAM'][rc,1].max())
    # limit log g range for RC
    plots.plotl(ax[1,0],[rclim[0],rclim[0]],[0,1000],color='k')
    plots.plotl(ax[1,0],[rclim[1],rclim[1]],[0,1000],color='k')
    rclogg = np.where((allstar['FPARAM'][rc,1] > rclim[0]) & (allstar['FPARAM'][rc,1]<rclim[1]))[0]
    rgblogg = np.where((allstar['FPARAM'][rgb,1] > rclim[0]) & (allstar['FPARAM'][rgb,1]<rclim[1]))[0]

    dt=allstar['FPARAM'][:,0]-fit2d(allstar['FPARAM'][:,1]-2.5,allstar['FPARAM'][:,3])  
    nbest=10000
    for dtcrit in range(-500,500) :
        rcbd = np.where((dt[rc[rclogg]] < dtcrit))[0]
        rgbbd = np.where(dt[rgb[rgblogg]] > dtcrit)[0]
        nbd=len(rcbd)+len(rgbbd)
        if nbd < nbest : 
            dtbest=dtcrit
            nbest=nbd
    dtcrit=dtbest
    rcbd = np.where((dt[rc[rclogg]] < dtcrit))[0]
    rgbbd = np.where(dt[rgb[rgblogg]] > dtcrit)[0]
    print('dtcrit: ',dtcrit)
    print('bad fractions (rc, rgb): ',float(len(rcbd))/len(rclogg),float(len(rgbbd))/len(rgblogg))

    dt=allstar['FPARAM'][:,0]-(rgbfit[0]+rgbfit[1]*(allstar['FPARAM'][:,1]-2.5)+rgbfit[2]*allstar['FPARAM'][:,3])  
    cn=allstar['FPARAM'][:,4]-allstar['FPARAM'][:,5]
    plots.plotc(ax[0,1],dt[rc],allstar['FPARAM'][rc,1],allstar['FPARAM'][rc,3],marker='s',xr=[-500,500],yr=[4,1],size=20,zr=[-2.0,0.5],xt='dt',yt='log g',zt='[M/H]',colorbar=True)
    plots.plotc(ax[0,1],dt[rgb],allstar['FPARAM'][rgb,1],allstar['FPARAM'][rgb,3],marker='o',xr=[-500,500],yr=[4,1],size=20,zr=[-2.0,0.5])
    plots.plotl(ax[0,1],[-500,500],[rclim[0],rclim[0]],color='k')
    plots.plotl(ax[0,1],[-500,500],[rclim[1],rclim[1]],color='k')
    ax[1,1].hist(dt[rc],color='b',bins=np.arange(-500,500,10))
    ax[1,1].hist(dt[rgb],color='r',bins=np.arange(-500,500,10))
    ax[1,1].hist(dt[rc2],color='g',bins=np.arange(-500,500,10))
    ax[1,1].hist(dt[rc2cl],color='m',bins=np.arange(-500,500,10))
    ax[1,1].set_xlabel('dt')

    # plot dt vs C/N 
    #plots.plotc(ax[0,2],dt[rc[rclogg]],cn[rc[rclogg]],allstar['FPARAM'][rc[rclogg],1],marker='s',xr=[-500,500],yr=[-1.0,0.5],zr=[2,4],size=20,xt='dt',yt='[C/N]',zt='log g',colorbar=True)
    #plots.plotc(ax[1,2],dt[rgb[rgblogg]],cn[rgb[rgblogg]],allstar['FPARAM'][rgb[rgblogg],1],marker='o',xr=[-500,500],yr=[-1.0,0.5],zr=[2,4],size=20,xt='dt',yt='[C/N]',zt='log g',colorbar=True)
    plots.plotc(ax[0,2],dt[rc[rclogg]],cn[rc[rclogg]],allstar['FPARAM'][rc[rclogg],3],marker='s',xr=[-500,500],yr=[-1.0,0.5],zr=[-1.5,0.5],size=20,xt='dt',yt='[C/N]',zt='[M/H]',colorbar=True)
    plots.plotc(ax[1,2],dt[rgb[rgblogg]],cn[rgb[rgblogg]],allstar['FPARAM'][rgb[rgblogg],3],marker='o',xr=[-500,500],yr=[-1.0,0.5],zr=[-1.5,0.5],size=20,xt='dt',yt='[C/N]',zt='[M/H]',colorbar=True)

    cnslopebest=-0.2/100.
    cnintbest=0.
    nbest=10000
    slopearray=np.arange(cnslopebest-20*0.0001,cnslopebest+20*0.0001,0.0001)
    intarray=np.arange(cnintbest-10*0.02,cnintbest+10*0.02,0.02)
    offarray=np.arange(-0.9,-0.3,0.01)
    x=np.array([-200,400])
    for cnslope in slopearray :
     for cnoff in offarray :
      for cnint in intarray :
        cnfit=np.array([cnint,cnoff,cnslope])
        rgbbd=np.where(cn[rgb[rgblogg]] > cnfit[0]+cnfit[1]*allstar['FPARAM'][rgb[rgblogg],3]+cnfit[2]*dt[rgb[rgblogg]])[0]
        rcbd= np.where(cn[rc[rclogg]] < cnfit[0]+cnfit[1]*allstar['FPARAM'][rc[rclogg],3]+cnfit[2]*dt[rc[rclogg]])[0]
        nbd=float(len(rcbd))/len(rclogg)+float(len(rgbbd))/len(rgblogg)
        if nbd < nbest : 
            cnfitbest=cnfit
            nbest=nbd
    print(nbest)
    cnfit=cnfitbest
    x=np.array([-200,400])
    for i in [0,1] :
      for j in [2] :
        ax[i,j].plot(x,cnfit[0]+cnfit[1]*0.+x*cnfit[2])
        ax[i,j].plot(x,cnfit[0]+cnfit[1]*(-0.5)+x*cnfit[2])
        ax[i,j].plot(x,cnfit[0]+cnfit[1]*0.5+x*cnfit[2])

    rcbd = np.where((cn[rc[rclogg]] < cnfit[0]+cnfit[1]*allstar['FPARAM'][rc[rclogg],3]+cnfit[2]*dt[rc[rclogg]]))[0]
    rgbbd = np.where((cn[rgb[rgblogg]] > cnfit[0]+cnfit[1]*allstar['FPARAM'][rgb[rgblogg],3]+cnfit[2]*dt[rgb[rgblogg]]))[0]
    ax[0,2].text(0.98,0.98,'RC bad: {:5.3f}'.format(float(len(rcbd))/len(rclogg)),transform=ax[0,2].transAxes,va='top',ha='right')
    ax[1,2].text(0.98,0.98,'RGB bad: {:5.3f}'.format(float(len(rgbbd))/len(rgblogg)),transform=ax[1,2].transAxes,va='top',ha='right')
    print('bad fractions (rc, rgb): ',float(len(rcbd))/len(rclogg),float(len(rgbbd))/len(rgblogg),len(rcbd),len(rclogg),len(rgbbd),len(rgblogg))

    plt.tight_layout()
    if out is not None :
        plt.savefig(out+'.jpg')

    fig,ax=plots.multi(2,1)
    plots.plotp(ax[0],allstar['FPARAM'][rgb,0],allstar['FPARAM'][rgb,1],color='r',xr=[5500,3500],yr=[4,1],xt='Teff',yt='log g')
    plots.plotp(ax[0],allstar['FPARAM'][rc,0],allstar['FPARAM'][rc,1],color='b')
    plots.plotp(ax[0],allstar['FPARAM'][rc2,0],allstar['FPARAM'][rc2,1],color='g')
    plots.plotp(ax[0],allstar['FPARAM'][rc2cl,0],allstar['FPARAM'][rc2cl,1],color='m')

    x = -0.08 - 0.5*allstar['FPARAM'][:,3] - 0.0039*dt
    plots.plotp(ax[1],x[rgb],cn[rgb],color='r',xt='-0.08-0.5[M/H]-0.0039 dt',yt='[C/N]',xr=[-2.5,1.5],yr=[-2,1],nxtick=5)
    plots.plotp(ax[1],x[rc],cn[rc],color='b')
    plots.plotp(ax[1],x[rc2],cn[rc2],color='g')
    plots.plotp(ax[1],x[rc2cl],cn[rc2cl],color='m')
    ax[1].plot([-2,1.5],[-2,1.5])
    fig.tight_layout()
    if out is not None :
        plt.savefig(out+'_hr.pdf')

    return {'rclim' : rclim, 'rgbsep' : rgbfit, 'cnsep' : cnfit}
    

def apokasc(allstar,logg='LOGG_SYD_SCALING',apokasc_cat='APOKASC_cat_v3.6.0.fits',raw=True,plotcal=False,out='loggcomp',calloggrange=[-1.,3.8],loggrange=[1.,3.2],mhrange=[-2.5,0.5],teffrange=[3500,5500],calib=False) :
    '''
    asteroseismic log g comparisons for input allStar structure
    '''

    if calib :
        param = 'PARAM'
    else :
        param = 'FPARAM'

    gd=apselect.select(allstar,badval=['STAR_BAD'],mh=mhrange,logg=loggrange,teff=teffrange,raw=True)
    allstar=allstar[gd]

    # match ASPCAP with APOKASC, and get RC/RGB stars
    apokasc=fits.open(os.environ['APOGEE_IDR']+'/data/apokasc/'+apokasc_cat)[1].data
    i1,i2=match.match(allstar['APOGEE_ID'],apokasc['2MASS_ID'])
    rgb=np.where((apokasc['CONS_EVSTATES'][i2] == 'RGB') & (apokasc[logg][i2] > -1))[0]
    rc=np.where((apokasc['CONS_EVSTATES'][i2] == 'RC') & (apokasc[logg][i2] > -1))[0]
    notrc=np.where(apokasc['CONS_EVSTATES'][i2] != 'RC')[0]
    rc2=np.where((apokasc['CONS_EVSTATES'][i2] == '2CL') & (apokasc[logg][i2] > -1))[0]
    rc2cl=np.where((apokasc['CONS_EVSTATES'][i2] == 'RC/2CL') & (apokasc[logg][i2] > -1))[0]
    rcall=np.append(rc,rc2)
    rcall=np.append(rcall,rc2cl)

    # Do a 2D fit for RGB stars
    fig,ax=plots.multi(2,1,figsize=(12,6))
    rgbfit = fit.fit2d(allstar['FPARAM'][i1[rgb],1],allstar['FPARAM'][i1[rgb],3],
        allstar[param][i1[rgb],1]-apokasc[logg][i2[rgb]],zr=[-1,0.5],gdrange=[-2,2],yr=[-3,1],xr=[1,4],degree=1,
        plot=ax[0],yt='[M/H]',xt='log g',zt='$\Delta log g$',reject=0.3)
    rgbrms=(allstar[param][i1[rgb],1]-rgbfit(allstar['FPARAM'][i1[rgb],1],allstar['FPARAM'][i1[rgb],3])-apokasc[logg][i2[rgb]]).std()
    ax[0].text(0.98,0.98,'rms: {:5.3f}'.format(rgbrms),transform=ax[0].transAxes,va='top',ha='right')
    rgberrpar = cal.errfit(allstar[param][i1[rgb],0],allstar['SNR'][i1[rgb]],allstar[param][i1[rgb],3],
                        allstar[param][i1[rgb],1]-rgbfit(allstar['FPARAM'][i1[rgb],1],allstar['FPARAM'][i1[rgb],3])-apokasc[logg][i2[rgb]],
                        out=out+'_rgb',title='log g',zr=[0,0.2])

    rcfit = fit.fit2d(allstar['FPARAM'][i1[rc],1],allstar['FPARAM'][i1[rc],3],
        allstar[param][i1[rc],1]-apokasc[logg][i2[rc]],zr=[-1,0.5],gdrange=[-2,2],yr=[-3,1],xr=[1,4],degree=1,
        plot=ax[1],yt='[M/H]',xt='log g',zt='$\Delta log g$',reject=0.3)
    rcfit2 = fit.fit1d(allstar['FPARAM'][i1[rcall],1], allstar[param][i1[rcall],1]-apokasc[logg][i2[rcall]],zr=[-1,0.5],yr=[-3,1],xr=[1,4],degree=2,reject=0.3)
    rcrms=(allstar[param][i1[rc],1]-rcfit(allstar['FPARAM'][i1[rc],1],allstar['FPARAM'][i1[rc],3])-apokasc[logg][i2[rc]]).std()
    rcerrpar = cal.errfit(allstar[param][i1[rc],0],allstar['SNR'][i1[rc]],allstar[param][i1[rc],3],
                        allstar[param][i1[rc],1]-rcfit(allstar['FPARAM'][i1[rc],1],allstar['FPARAM'][i1[rc],3])-apokasc[logg][i2[rc]],
                        out=out+'_rc',title='log g',zr=[0,0.2])
    ax[1].text(0.98,0.98,'rms: {:5.3f}'.format(rcrms),transform=ax[1].transAxes,va='top',ha='right')
    fig.tight_layout()
    if out is not None :
        fig.savefig(out+'.jpg')

    # set up plots
    if raw and plotcal :
        fig,ax=plots.multi(2,3,hspace=0.5,wspace=0.001,figsize=(12,12))
    else :
        fig,tmpax=plots.multi(1,3,hspace=0.5,wspace=0.001,figsize=(8,8))

    # diff color-coded by gravity as f([M/H])
    # diff color-coded by [M/H] as f(log g)
    # RGB and RC as f(log g)
    if raw :
        if plotcal: tmpax=ax[:,0]
        plots.plotc(tmpax[0],allstar['FPARAM'][i1,3],allstar[param][i1,1]-apokasc[logg][i2],
           allstar['FPARAM'][i1,1],zr=[0,5],xr=[-2.5,0.5],yr=[-1,1],xt='[M/H]',yt='ASPCAP-seismic log g',zt='log g',size=15,colorbar=True)
        plots.plotc(tmpax[1],allstar['FPARAM'][i1[rgb],1],allstar[param][i1[rgb],1]-apokasc[logg][i2[rgb]],
           allstar['FPARAM'][i1[rgb],3],xr=[0,5],zr=[-2.5,0.5],yr=[-1,1],zt='[M/H]',yt='ASPCAP-seismic log g',xt='log g',size=15,colorbar=True)
        loggfit=np.arange(1,3.5,0.01)
        mhfit=loggfit*0.
        plots.plotl(tmpax[1],loggfit,rgbfit(loggfit,mhfit),color='orange',linewidth=1.5)
        mhfit=loggfit*0-1.
        plots.plotl(tmpax[1],loggfit,rgbfit(loggfit,mhfit),color='c',linewidth=1.5)
        mhfit=loggfit*0+0.5
        plots.plotl(tmpax[1],loggfit,rgbfit(loggfit,mhfit),color='r',linewidth=1.5)

        plots.plotp(tmpax[2],allstar['FPARAM'][i1[rgb],1],allstar[param][i1[rgb],1]-apokasc[logg][i2[rgb]],
           xr=[0,5],yr=[-1,1],xt='seismic log g',yt='ASPCAP-seismic log g',label=[0.9,0.8,'RGB'],color='r',size=15)
        plots.plotp(tmpax[2],allstar['FPARAM'][i1[rc],1],allstar[param][i1[rc],1]-apokasc[logg][i2[rc]],
           xr=[0,5],yr=[-1,1],xt='seismic log g',yt='ASPCAP-seismic log g',label=[0.9,0.6,'RC'],color='b',size=15)
        plots.plotp(tmpax[2],allstar['FPARAM'][i1[rc2],1],allstar[param][i1[rc2],1]-apokasc[logg][i2[rc2]],
           xr=[0,5],yr=[-1,1],xt='seismic log g',yt='ASPCAP-seismic log g',label=[0.9,0.6,'RC'],color='g',size=15)
        plots.plotp(tmpax[2],allstar['FPARAM'][i1[rc2cl],1],allstar[param][i1[rc2cl],1]-apokasc[logg][i2[rc2cl]],
           xr=[0,5],yr=[-1,1],xt='seismic log g',yt='ASPCAP-seismic log g',label=[0.9,0.6,'RC'],color='m',size=15)

        loggfit=np.arange(2.5,3.5,0.01)
        mhfit=loggfit*0.
        plots.plotl(tmpax[2],loggfit,rcfit(loggfit,mhfit),color='g',linewidth=2)
        plots.plotl(tmpax[2],loggfit,rcfit2(loggfit),color='k',linewidth=2)

        #plots.plotc(tmpax[3],allstar['FPARAM'][i1[rgb],1],allstar['PARAM'][i1[rgb],1]-allstar['FPARAM'][i1[rgb],1],
        #   allstar['FPARAM'][i1[rgb],3],xr=[0,5],yr=[-1,1],xt='seismic log g',yt='corrected-raw log g',label=[0.1,0.9,'allstar (Kurucz)'],zr=[-2,0.5])

    if plotcal :
        if raw: tmpax=ax[:,1]
        param=allstar['FPARAM'][:,1]-rgbfit(allstar['FPARAM'][:,1],allstar['FPARAM'][:,3])
        param[i1[rc]]=allstar['FPARAM'][i1[rc],1]-rcfit(allstar['FPARAM'][i1[rc],1],allstar['FPARAM'][i1[rc],3])

        plots.plotc(tmpax[0],allstar['FPARAM'][i1,3],param[i1]-apokasc[logg][i2],
           allstar['FPARAM'][i1,1],zr=[0,5],xr=[-2.5,0.5],yr=[-1,1],xt='[M/H]',colorbar=True,zt='log g',size=15)
        plots.plotc(tmpax[1],allstar['FPARAM'][i1,1],param[i1]-apokasc[logg][i2],
           allstar['FPARAM'][i1,3],xr=[1,4],zr=[-2.5,0.5],yr=[-1,1],zt='[M/H]',colorbar=True,xt='log g',size=15)
        plots.plotp(tmpax[2],allstar['FPARAM'][i1[rgb],1],param[i1[rgb]]-apokasc[logg][i2[rgb]],
           xr=[1,4],yr=[-1,1],xt='log g',color='r',size=15)
        plots.plotp(tmpax[2],allstar['FPARAM'][i1[rc],1],param[i1[rc]]-apokasc[logg][i2[rc]],
           xr=[1,4],yr=[-1,1],xt='log g',color='b',size=15)
        plots.plotp(tmpax[2],allstar['FPARAM'][i1[rc2],1],param[i1[rc2]]-apokasc[logg][i2[rc2]],
           xr=[1,4],yr=[-1,1],xt='log g',color='g',size=15)
        plots.plotp(tmpax[2],allstar['FPARAM'][i1[rc2cl],1],param[i1[rc2cl]]-apokasc[logg][i2[rc2cl]],
           xr=[1,4],yr=[-1,1],xt='log g',color='m',size=15)
        #plots.plotc(tmpax[3],allstar['FPARAM'][i1[rc],1],allstar['PARAM'][i1[rc],1]-allstar['FPARAM'][i1[rc],1],
        #    allstar['FPARAM'][i1[rc],3],xr=[0,5],yr=[-1,1],xt='seismic log g',zr=[-2,0.5])
    fig.tight_layout()
    if out is not None :
        plt.savefig(out+'_b.jpg')
        plt.savefig(out+'_b.pdf')

    return {'calloggmin' : calloggrange[0], 'calloggmax' : calloggrange[1], 'loggmin' : loggrange[0], 'loggmax' : loggrange[1], 
            'mhmin' : mhrange[0], 'mhmax' : mhrange[1], 'calteffmin': teffrange[0], 'calteffmax' : teffrange[1],
            'rgbfit' : rgbfit.parameters, 'rcfit' : rcfit.parameters, 'rcfit2' : rcfit2.parameters, 'rgbrms' : rgbrms, 'rcrms' : rcrms ,
            'rgberrpar': rgberrpar, 'rcerrpar': rcerrpar}


def clusters(allstar,xr=[-2.75,0.5],yr=[-1.,1.],zr=[3500,5500],apokasc='APOKASC_cat_v3.6.0.fits',firstgen=False) :
    '''
    Compare ASPCAP gravities in clusters to physical gravities
    '''
    fig,ax=plots.multi(1,2,hspace=0.001)

    # put APOKASC underneath
    apokasc=fits.open(os.environ['APOGEE_DIR']+'/data/apokasc/'+apokasc)[1].data
    i1,i2=match.match(allstar['APOGEE_ID'],apokasc['2MASS_ID'])
    plots.plotc(ax[0],allstar['FPARAM'][i1,3],allstar['FPARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],allstar['FPARAM'][i1,0],zr=zr)
    plots.plotc(ax[1],allstar['PARAM'][i1,3],allstar['PARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],allstar['PARAM'][i1,0],zr=zr)

    # physical gravities
    clust=apselect.clustdata()
    itext=0
    out=open('clust.txt','w')
    for cluster in ['M92','M15','M53','M2','M13','M3','M5','N2420','M67','N6819','N6791'] :
        i=np.where(clust.name == cluster)
        dist=clust[i].dist*1000.
        mh=clust[i].mh
        mass=clust[i].giant_mass
        ejk=0.452*clust[i].ebv
        ah=1.55*clust[i].ebv
        age=np.log10(clust[i].age*1.e9)
        name=clust[i].name
        ytext=0.85-itext%3*0.15
        if mass > 0 :
            # get cluster members
            j=np.array(apselect.clustmember(allstar,cluster,raw=True,firstgen=firstgen))

            # calculate physical gravities
            lum=10.**(-0.4*(allstar['H'][j]-ah+isochrones.bc(allstar['FPARAM'][j,0],filt='h',agerange=[age-0.05,age+0.05])-(5*np.log10(dist)-5)-4.74))*astropy.constants.L_sun.cgs.value
            logg=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*astropy.constants.sigma_sb.cgs.value*allstar['FPARAM'][j,0]**4/lum)
            plots.plotc(ax[0],allstar['FPARAM'][j,3]*0+mh,allstar['FPARAM'][j,1]-logg,allstar['FPARAM'][j,0],xr=xr,yr=yr,zr=zr,yt='ASPCAP-physical log g')
            ax[0].text(0.9,0.1,'raw',transform=ax[0].transAxes,ha='right')

            plots.plotp(ax[0],allstar['FPARAM'][j,3]*0+mh,allstar['FPARAM'][j,1]-logg,color='k')
            plots.plotp(ax[0],mh[0],np.median(allstar['FPARAM'][j,1]-logg),size=40,color='r')
            ax[0].text(mh[0],ytext,name[0],ha='center')

            out.write('{:<20s}{:8.3f}{:8.3f}{:8.3f}\n'.format(clust[i].name[0],clust[i].dist[0],clust[i].ebv[0],mass[0]))

            gd=np.where((allstar['PARAM'][j,3]>-9)&(allstar['PARAM'][j,1]>-9))[0]
            axim=plots.plotc(ax[1],allstar['PARAM'][j[gd],3]*0+mh,allstar['PARAM'][j[gd],1]-logg[gd],allstar['PARAM'][j[gd],0],xr=xr,yr=yr,zr=zr,xt='[M/H]',yt='ASPCAP-physical log g')
            ax[1].text(0.9,0.1,'calibrated',transform=ax[1].transAxes,ha='right')

            plots.plotp(ax[1],mh[0],np.median(allstar['PARAM'][j[gd],1]-logg[gd]),size=40)
            # apply a temperature correction for the physical gravities
            logg_new=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*astropy.constants.sigma_sb.cgs.value*(allstar['FPARAM'][j,0]-100.*allstar['FPARAM'][j,3])**4/lum)
            plots.plotp(ax[1],mh[0],np.median(allstar['PARAM'][j,1]-logg_new),size=40,color='b')
            # use a photometric temperature
            logg_phot=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*astropy.constants.sigma_sb.cgs.value*stars.ghb(allstar['J'][j]-allstar['K'][j]-ejk,allstar['FPARAM'][j,3])[0]**4/lum)
            plots.plotp(ax[1],mh[0],np.median(allstar['PARAM'][j,1]-logg_phot),size=40,color='g')
            ax[1].text(mh[0],ytext,name[0],ha='center')
            itext+=1

    # Now adding the colorbar
    cbaxes = fig.add_axes([0.9, 0.1, 0.03, 0.8]) 
    cb = plt.colorbar(axim, cax = cbaxes)  
    out.close()


def dr13dr12() :
    '''
    ASPCAP compared with physical and asteroseismic log g, DR13/DR12/l30i
    '''
    apokasc=fits.open(os.environ['APOGEE_DIR']+'/data/apokasc/'+apokasc)[1].data
    j=np.where(apokasc['LOGG_SYD_SCALING'] > -1)[0]
    apokasc=apokasc[j]

    dr12load=apload.ApLoad(dr='dr12')
    dr12=dr12load.allStar()[1].data
    dr13load=apload.ApLoad(dr='dr13')
    dr13=dr13load.allStar()[1].data
    dr13load.aspcap='l30i'
    dr13load.results='l30i'
    l30i=dr13load.allStar()[1].data

    fig,ax =plots.multi(3,2,wspace=0.001,hspace=0.001)

    # physical gravities
    clust=apselect.clustdata()
    for cluster in ['M92','M15','M53','M2','M13','M3','M5'] :
        i=np.where(clust.name == cluster)
        dist=clust[i].dist*1000.
        mh=clust[i].mh
        mass=0.85

        #DR12
        j=apselect.clustmember(dr12,cluster,raw=True)
        lum=10.**(-0.4*(dr12['H'][j]+isochrones.bc(dr12['FPARAM'][j,0],filt='h',agerange=[10,10.1])-(5*np.log10(dist)-5)-4.74))*astropy.constants.L_sun.cgs.value
        logg=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*astropy.constants.sigma_sb.cgs.value*dr12['FPARAM'][j,0]**4/lum)
        plots.plotc(ax[0,0],dr12['FPARAM'][j,3]*0+mh,dr12['FPARAM'][j,1]-logg,dr12['FPARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],yt='ASPCAP-physical log g',label=[0.1,0.9,'DR2 raw'])
        plots.plotc(ax[1,0],dr12['PARAM'][j,3]*0+mh,dr12['PARAM'][j,1]-logg,dr12['PARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',yt='ASPCAP-physical log g',label=[0.1,0.9,'DR12 cal'])

        #DR13
        j=apselect.clustmember(dr13,cluster,raw=True)
        lum=10.**(-0.4*(dr13['H'][j]+isochrones.bc(dr13['FPARAM'][j,0],filt='h',agerange=[10,10.1])-(5*np.log10(dist)-5)-4.74))*astropy.constants.L_sun.cgs.value
        logg=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*astropy.constants.sigma_sb.cgs.value*dr13['FPARAM'][j,0]**4/lum)
        plots.plotc(ax[0,1],dr13['FPARAM'][j,3]*0+mh,dr13['FPARAM'][j,1]-logg,dr13['FPARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],label=[0.1,0.9,'DR13 raw'])
        plots.plotc(ax[1,1],dr13['PARAM'][j,3]*0+mh,dr13['PARAM'][j,1]-logg,dr13['PARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],label=[0.1,0.9,'DR13 cal'],xt='[M/H]')

        #l30i
        j=apselect.clustmember(l30i,cluster,raw=True)
        lum=10.**(-0.4*(l30i['H'][j]+isochrones.bc(l30i['FPARAM'][j,0],filt='h',agerange=[10,10.1])-(5*np.log10(dist)-5)-4.74))*astropy.constants.L_sun.cgs.value
        logg=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*astropy.constants.sigma_sb.cgs.value*l30i['FPARAM'][j,0]**4/lum)
        plots.plotc(ax[0,2],l30i['FPARAM'][j,3]*0+mh,l30i['FPARAM'][j,1]-logg,l30i['FPARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],label=[0.1,0.9,'l30i raw'])
        plots.plotc(ax[1,2],l30i['PARAM'][j,3]*0+mh,l30i['PARAM'][j,1]-logg,l30i['PARAM'][j,0],xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],label=[0.1,0.9,'l30i cal'],xt='[M/H]')
    
    plt.show()
    pdb.set_trace()

    # plots vs asterseismic
    fig,ax =plots.multi(3,2,wspace=0.001,hspace=0.001)
    
    i1,i2=match.match(dr13['APOGEE_ID'],apokasc['2MASS_ID'])
    plots.plotc(ax[0,0],dr13['FPARAM'][i1,3],dr13['FPARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],dr13['FPARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',label=[0.1,0.9,'DR13 raw'])
    plots.plotc(ax[1,0],dr13['PARAM'][i1,3],dr13['PARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],dr13['PARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',yt='ASPCAP-seismic log g',label=[0.1,0.9,'DR13 cal'])

    i1,i2=match.match(dr12['APOGEE_ID'],apokasc['2MASS_ID'])
    plots.plotc(ax[0,1],dr12['FPARAM'][i1,3],dr12['FPARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],dr12['FPARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],label=[0.1,0.9,'DR12 raw'])
    plots.plotc(ax[1,1],dr12['PARAM'][i1,3],dr12['PARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],dr12['PARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',label=[0.1,0.9,'DR12 cal'])

    i1,i2=match.match(l30i['APOGEE_ID'],apokasc['2MASS_ID'])
    plots.plotc(ax[0,2],l30i['FPARAM'][i1,3],l30i['FPARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],l30i['FPARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],label=[0.1,0.9,'l30i raw'])
    plots.plotc(ax[1,2],l30i['PARAM'][i1,3],l30i['PARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],l30i['PARAM'][i1,0],xr=[-2.5,1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',label=[0.1,0.9,'l30i cal'])

    plt.show()

def kurucz_marcs(logg='LOGG_SYD_SCALING',apokasc='APOKASC_cat_v3.6.0.fits') :
    '''
    asteroseismic log g comparisons for Kurucz and MARCS results
    '''
    # read APOKASC
    apokasc=fits.open('APOKASC_cat_v3.6.0.fits')[1].data
    #j=np.where((apokasc['TEFF_FIT'] < 4000) & (apokasc[logg] > -500))[0]
    #j=np.where((apokasc['CONS_EVSTATES'] == 'RGB'))[0]
    #apokasc=apokasc[j]

    # read DR13 and l30i
    dr13load=apload.ApLoad(dr='dr13')
    dr13=dr13load.allStar()[1].data
    dr13load.aspcap='l30i'
    dr13load.results='l30i'
    l30i=dr13load.allStar()[1].data

    # define axes
    fig=plt.figure()
    ax1=fig.add_subplot(211)
    ax2=fig.add_subplot(212)
    #ax3=fig.add_subplot(223)
    #ax4=fig.add_subplot(224)

    # match l30i with APOKASC
    i1,i2=match.match(l30i['APOGEE_ID'],apokasc['2MASS_ID'])
    warn=np.where(l30i['ASPCAPFLAG'][i1] & bitmask.aspcapflagval('ATMOS_HOLE_WARN'))[0]
    bad=np.where(l30i['ASPCAPFLAG'][i1] & bitmask.aspcapflagval('ATMOS_HOLE_BAD'))[0]
    rgb=np.where(apokasc['CONS_EVSTATES'][i2] == 'RGB')[0]
    rc=np.where(apokasc['CONS_EVSTATES'][i2] == 'RC')[0]
    #plt.plot(apokasc[logg][i2],l30i['FPARAM'][i1,1],'ro')
    # color code by [M/H]
    #plots.plotc(ax1,l30i['FPARAM'][i1,0],l30i['FPARAM'][i1,1]-apokasc[logg][i2],l30i['FPARAM'][i1,3],zr=[-1,0.5],xr=[3500,5000],yr=[-1,1],xt='Teff',yt='ASPCAP-seismic log g',label=[0.1,0.9,'l30i (MARCS)'],colorbar=True,zt='[M/H]')
    # color code by [alpha/M]
    plots.plotc(ax1,l30i['FPARAM'][i1,0],l30i['FPARAM'][i1,1]-apokasc[logg][i2],l30i['FPARAM'][i1,6],zr=[-0.1,0.4],xr=[3500,5000],yr=[-1,1],xt='seismic log g',yt='ASPCAP-seismic log g',label=[0.1,0.9,'l30i (MARCS)'],colorbar=True,zt='[alpha/M]')
    # plot ATMOS_HOLE_WARN and BAD points 
    plots.plotp(ax1,l30i['FPARAM'][i1[warn],0],l30i['FPARAM'][i1[warn],1]-apokasc[logg][i2[warn]],color='y')
    plots.plotp(ax1,l30i['FPARAM'][i1[bad],0],l30i['FPARAM'][i1[bad],1]-apokasc[logg][i2[bad]],color='r')
    # colod code by Teff
    #plots.plotc(ax2,l30i['FPARAM'][i1,1],l30i['FPARAM'][i1,1]-apokasc[logg][i2],l30i['FPARAM'][i1,0],zr=[3500,5000],xr=[0,5],yr=[-1,1],xt='seismic log g',yt='ASPCAP-seismic log g',label=[0.1,0.9,'l30i'],colorbar=True,zt='Teff')

    # match dr13 with APOKASC
    i1,i2=match.match(dr13['APOGEE_ID'],apokasc['2MASS_ID'])
    #plt.plot(apokasc[logg][i2],dr13['FPARAM'][i1,1],'bo')
    #plt.xlim(1.60,1.15)
    #plt.ylim(1.15,1.85)
    #plots.plotc(ax2,dr13['FPARAM'][i1,0],dr13['FPARAM'][i1,1]-apokasc[logg][i2],dr13['FPARAM'][i1,3],zr=[-1,0.5],xr=[3500,5000],yr=[-1,1],xt='Teff',yt='ASPCAP-seismic log g',label=[0.1,0.9,'dr13 (Kurucz)'],colorbar=True,zt='[M/H]')
    plots.plotc(ax2,dr13['FPARAM'][i1,0],dr13['FPARAM'][i1,1]-apokasc[logg][i2],dr13['FPARAM'][i1,6],zr=[-0.1,0.4],xr=[3500,5000],yr=[-1,1],xt='seismic log g',yt='ASPCAP-seismic log g',label=[0.1,0.9,'dr13 (Kurucz)'],colorbar=True,zt='[alpha/M]')
    ##plots.plotc(ax4,dr13['FPARAM'][i1,1],dr13['FPARAM'][i1,1]-apokasc[logg][i2],dr13['FPARAM'][i1,0],zr=[3500,5000],xr=[0,5],yr=[-1,1],xt='seismic log g',yt='ASPCAP-seismic log g',label=[0.1,0.9,'dr13'],colorbar=True,zt='Teff')

    plt.tight_layout()
    plt.show()

