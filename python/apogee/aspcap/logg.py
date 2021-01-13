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
try: from apogee.aspcap import cal
except: pass
from apogee.aspcap import err, teff
from apogee.speclib import isochrones
import pdb
import matplotlib.pyplot as plt
import numpy as np
import os
import astropy
from scipy import interpolate

def rcrgb(allstar,apokasc='APOKASC_cat_v3.6.0.fits',logg='LOGG_SYD_SCALING',rclim=np.array([2.38,3.5]),out='rcrgbsep') :
    '''
    asteroseismic log g comparisons for input allStar structure
    '''

    gd=apselect.select(allstar,badval=['STAR_BAD'],logg=[0,3.8],teff=[3500,5500],raw=True)
    allstar=allstar[gd]

    # match ASPCAP with APOKASC, and get RC/RGB stars
    apokasc=fits.open(os.environ['APOGEE_DIR']+'/data/apokasc/'+apokasc)[1].data
    # strip off .XXXX if we have it, e.g. from calibration fields where we have added .FIELD
    apogee_id = np.array(np.core.defchararray.split(allstar['APOGEE_ID'],'.').tolist())[:,0]
    i1,i2=match.match(apogee_id,apokasc['2MASS_ID'])
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
        plt.savefig(out+'.png')
        plt.close(fig)

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
        plt.close(fig)

    return {'rclim' : rclim, 'rgbsep' : rgbfit, 'cnsep' : cnfit}
    
def rcrgb_plot(a,out=None) :
    """ Plot logg classification from bitmask
    """
    b=bitmask.ParamBitMask()
    rgb=np.where((a['PARAMFLAG'][:,1] & b.getval('LOGG_CAL_RGB')) > 0)[0]
    rc=np.where((a['PARAMFLAG'][:,1] & b.getval('LOGG_CAL_RC')) > 0)[0]
    ms=np.where((a['PARAMFLAG'][:,1] & b.getval('LOGG_CAL_MS')) > 0)[0]
    rgb_ms=np.where((a['PARAMFLAG'][:,1] & b.getval('LOGG_CAL_RGB_MS')) > 0)[0]

    fig,ax = plots.multi(1,1)
    plots.plotp(ax,a['FPARAM'][rgb,0],a['FPARAM'][rgb,1],color='r',size=1,
                xr=[8000,3000],yr=[6,-1],xt='$T_{eff}$',yt='log g')
    plots.plotp(ax,a['FPARAM'][rc,0],a['FPARAM'][rc,1],color='b',size=1)
    plots.plotp(ax,a['FPARAM'][ms,0],a['FPARAM'][ms,1],color='g',size=1)
    plots.plotp(ax,a['FPARAM'][rgb_ms,0],a['FPARAM'][rgb_ms,1],color='m',size=1)
    if out is not None :
        fig.savefig(out+'.png')
        plt.close()


def dwarf(allstar,mhrange=[-2.5,1.0],loggrange=[3.8,5.5],teffrange=[3000,7500],apokasc_cat='APOKASC_cat_v4.4.2.fits',out='logg',calib=False) :
    """ logg calibration for dwarfs, from asteroseismic and isochrones
    """
    if calib :
        param = 'PARAM'
    else :
        param = 'FPARAM'

    gd=apselect.select(allstar,badval=['STAR_BAD'],mh=mhrange,logg=loggrange,teff=teffrange,raw=True)
    allstar=allstar[gd]
    try:
        gd=np.where(allstar['VISIT'] == 0)[0]
        allstar=allstar[gd]
    except: pass
 
    # match ASPCAP with APOKASC, and get RC/RGB stars
    apokasc=fits.open(os.environ['APOGEE_DIR']+'/data/apokasc/'+apokasc_cat)[1].data
    # strip off .XXXX if we have it, e.g. from calibration fields where we have added .FIELD
    apogee_id = np.array(np.core.defchararray.split(allstar['APOGEE_ID'],'.').tolist())[:,0]
    i1,i2=match.match(apogee_id,apokasc['2MASS_ID'])

    # now get isochrone logg from lower main sequence
    isologg=isochrone(allstar,snrbd=50)
    isochrone_id = np.array(np.core.defchararray.split(isologg['APOGEE_ID'],'.').tolist())[:,0]
    j1,j2=match.match(apogee_id,isochrone_id)

    # plots of gravity differences
    fig,ax=plots.multi(2,2)
    plots.plotc(ax[0,0],allstar[param][i1,1],allstar[param][i1,1]-apokasc['LOGG_DW'][i2],allstar[param][i1,0],yr=[-1,1],
                xt='log g',yt=r'$\Delta$logg',zt='Teff',colorbar=True,xr=[3,6],zr=[4000,7000])
    plots.plotc(ax[0,1],allstar[param][i1,3],allstar[param][i1,1]-apokasc['LOGG_DW'][i2],allstar[param][i1,0],yr=[-1,1],
                xt='[M/H]',yt=r'$\Delta$logg',zt='Teff',colorbar=True,xr=[-2,0.5],zr=[4000,7000])
    plots.plotc(ax[1,0],allstar[param][i1,0],allstar[param][i1,1]-apokasc['LOGG_DW'][i2],10.**allstar[param][i1,2],yr=[-1,1],
                xt='Teff',yt=r'$\Delta$logg',zt='vmicro',colorbar=True,xr=[3000,8000],zr=[0.5,2.5])
    plots.plotc(ax[1,1],allstar[param][i1,0],allstar[param][i1,1]-apokasc['LOGG_DW'][i2],allstar[param][i1,3],yr=[-1,1],
                xt='Teff',yt=r'$\Delta$logg',zt='[M/H]',colorbar=True,xr=[3000,8000],zr=[-2,0.5])
    # only add main sequence in Teff plot
    plots.plotc(ax[1,1],allstar[param][j1,0],allstar[param][j1,1]-isologg['ISOLOGG'][j2],allstar[param][j1,3],zr=[-2,0.5])
    plt.tight_layout()

    # 2D fit as f(Teff,[M/H]), using both APOKASC and isochrone log g
    gd=np.where(apokasc['LOGG_DW'][i2] > -99)[0]
    tfit=allstar[param][i1[gd],0]
    mhfit=allstar[param][i1[gd],3]
    diff=allstar[param][i1[gd],1]-apokasc['LOGG_DW'][i2[gd]]
    snrfit=allstar['SNR'][i1[gd]]

    # do fit from high S/N, but get uncertainties from all
    gd=np.where(allstar['SNR'][j1] > 300)[0]
    msfit = fit.fit2d(np.append(tfit,allstar[param][j1[gd],0]),
                      np.append(mhfit,allstar[param][j1[gd],3]),
                      np.append(diff,allstar[param][j1[gd],1]-isologg['ISOLOGG'][j2[gd]]),degree=1,reject=0.3)

    # for uncertainties, all all S/N
    tfit=np.append(tfit,allstar[param][j1,0])
    mhfit=np.append(mhfit,allstar[param][j1,3])
    diff=np.append(diff,allstar[param][j1,1]-isologg['ISOLOGG'][j2])
    snrfit=np.append(snrfit,allstar['SNR'][j1])
    mserrpar = err.errfit(tfit,np.clip(snrfit,0.,249.),mhfit,diff-msfit(tfit,mhfit),
                          out=out+'_ms',title='log g',zr=[0,0.2])
    #mserrpar=np.zeros([4])

    # plot the relation
    tfit=np.arange(teffrange[0],teffrange[1],10)
    mhfit=tfit*0.
    plots.plotl(ax[1,1],tfit,msfit(tfit,mhfit),color='orange',linewidth=1.5)
    mhfit=tfit*0-1.
    plots.plotl(ax[1,1],tfit,msfit(tfit,mhfit),color='c',linewidth=1.5)
    mhfit=tfit*0+0.5
    plots.plotl(ax[1,1],tfit,msfit(tfit,mhfit),color='r',linewidth=1.5)
    ax[1,1].grid()
    if out is not None:
        fig.savefig(out+'_dwarfs.png')
        plt.close()
    
    # HR diagram plot color coded by asteroseismic gravity differences
    hrfig,hrax=plots.multi(1,2,hspace=0.001)
    gd=np.where(apokasc['APOKASC2_LOGG'][i2] > -99)[0]
    plots.plotc(hrax[0],allstar[param][i1[gd],0],allstar[param][i1[gd],1],allstar[param][i1[gd],1]-apokasc['APOKASC2_LOGG'][i2[gd]],
                xr=[8000,3000],yr=[6,0],zr=[-0.5,0.5],colorbar=True,zt=r'$\Delta$ logg',xt='Teff',yt='logg')
    plots.plotc(hrax[0],allstar[param][j1,0],allstar[param][j1,1],allstar[param][j1,1]-isologg['ISOLOGG'][j2],zr=[-0.5,0.5])
    plots.plotc(hrax[1],allstar[param][i1[gd],0],apokasc['APOKASC2_LOGG'][i2[gd]],allstar[param][i1[gd],1]-apokasc['APOKASC2_LOGG'][i2[gd]],
                xr=[8000,3000],yr=[6,0],zr=[-0.5,0.5],colorbar=True,zt=r'$\Delta$ logg',xt='Teff',yt='APOKASC logg')
    # use asteroseismic logg on y axis
    gd=np.where(apokasc['LOGG_DW'][i2] > -99)[0]
    plots.plotc(hrax[0],allstar[param][i1[gd],0],allstar[param][i1[gd],1],allstar[param][i1[gd],1]-apokasc['LOGG_DW'][i2[gd]],
                xr=[8000,3000],yr=[6,0],zr=[-0.5,0.5])
    plots.plotc(hrax[1],allstar[param][i1[gd],0],apokasc['LOGG_DW'][i2[gd]],allstar[param][i1[gd],1]-apokasc['LOGG_DW'][i2[gd]],
                xr=[8000,3000],yr=[6,0],zr=[-0.5,0.5])

    if out is not None:
        hrfig.savefig(out+'_all.png')
        plt.close()

    return {'calloggmin' : loggrange[0], 'calloggmax' : loggrange[1], 'loggmin' : loggrange[0], 'loggmax' : loggrange[1], 
            'mhmin' : mhrange[0], 'mhmax' : mhrange[1], 'temin': teffrange[0], 'temax' : teffrange[1],
            'msfit' : msfit.parameters, 'errpar' : mserrpar }

 
def apokasc(allstar,apokasc_cat='APOKASC_cat_v4.4.2.fits',raw=True,plotcal=False,out='loggcomp',calloggrange=[-1.,3.8],loggrange=[-1.,3.8],mhrange=[-2.5,0.5],teffrange=[3500,5500],calteffrange=[3000,6000],calib=False) :
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
    apokasc=fits.open(os.environ['APOGEE_DIR']+'/data/apokasc/'+apokasc_cat)[1].data
    # strip off .XXXX if we have it, e.g. from calibration fields where we have added .FIELD
    apogee_id = np.array(np.core.defchararray.split(allstar['APOGEE_ID'],'.').tolist())[:,0]
    i1,i2=match.match(apogee_id,apokasc['2MASS_ID'])
    try:
        print('trying APOKASC2 catalog tags...')
        logg='APOKASC2_LOGG'
        rgb=np.where((apokasc['CONS_EVSTATES'][i2] == 'RGB') & (apokasc[logg][i2] > -1))[0]
        rc=np.where((apokasc['CONS_EVSTATES'][i2] == 'RC') & (apokasc[logg][i2] > -1))[0]
        notrc=np.where(apokasc['CONS_EVSTATES'][i2] != 'RC')[0]
        rc2=np.where((apokasc['CONS_EVSTATES'][i2] == '2CL') & (apokasc[logg][i2] > -1))[0]
        rc2cl=np.where((apokasc['CONS_EVSTATES'][i2] == 'RC/2CL') & (apokasc[logg][i2] > -1))[0]
        # use LOGG_DW if we have it
        dw=np.where((apokasc[logg][i2] < -99) & (apokasc['LOGG_DW'][i2] >-99) )[0]
        apokasc[logg][i2[dw]] = apokasc['LOGG_DW'][i2[dw]]
        rgb=np.append(rgb,dw)
    except :
        # DR14 used APOKASC_cat_v3.6.0
        print('trying older APOKASC catalog tags...')
        logg='LOGG_SYD_SCALING'
        rgb=np.where((apokasc['CONS_EVSTATES'][i2] == 'RGB') & (apokasc[logg][i2] > -1))[0]
        rc=np.where((apokasc['CONS_EVSTATES'][i2] == 'RC') & (apokasc[logg][i2] > -1))[0]
        notrc=np.where(apokasc['CONS_EVSTATES'][i2] != 'RC')[0]
        rc2=np.where((apokasc['CONS_EVSTATES'][i2] == '2CL') & (apokasc[logg][i2] > -1))[0]
        rc2cl=np.where((apokasc['CONS_EVSTATES'][i2] == 'RC/2CL') & (apokasc[logg][i2] > -1))[0]
    rcall=np.append(rc,rc2)
    rcall=np.append(rcall,rc2cl)

    # Do some 2D fits for RGB stars
    fig,ax=plots.multi(2,1,figsize=(12,6))
    # linear in logg and [M/H]
    print('rgbfit', len(rgb))
    rgbfit = fit.fit2d(allstar['FPARAM'][i1[rgb],1],allstar['FPARAM'][i1[rgb],3],
        allstar[param][i1[rgb],1]-apokasc[logg][i2[rgb]],zr=[-1,0.5],gdrange=[-2,2],yr=[-3,1],xr=[1,4],degree=1,
        plot=ax[0],yt='[M/H]',xt='log g',zt='$\Delta log g$',reject=0.3)
    # cubic in logg, linear in [M/H]
    data=allstar[param][i1[rgb],1]-apokasc[logg][i2[rgb]]
    design=np.ones([5,len(rgb)])
    design[1,:]=allstar['FPARAM'][i1[rgb],1]
    design[2,:]=allstar['FPARAM'][i1[rgb],1]**2
    design[3,:]=allstar['FPARAM'][i1[rgb],1]**3
    design[4,:]=allstar['FPARAM'][i1[rgb],3]
    params=fit.linear(data,design)[0]
    rgbrms=(allstar[param][i1[rgb],1]-rgbfit(allstar['FPARAM'][i1[rgb],1],allstar['FPARAM'][i1[rgb],3])-apokasc[logg][i2[rgb]]).std()
    ax[0].text(0.98,0.98,'rms: {:5.3f}'.format(rgbrms),transform=ax[0].transAxes,va='top',ha='right')
    print('errfit')
    rgberrpar = err.errfit(allstar[param][i1[rgb],0],allstar['SNR'][i1[rgb]],allstar[param][i1[rgb],3],
                        allstar[param][i1[rgb],1]-rgbfit(allstar['FPARAM'][i1[rgb],1],allstar['FPARAM'][i1[rgb],3])-apokasc[logg][i2[rgb]],
                        out=out+'_rgb',title='log g',zr=[0,0.2])
    loggmin=allstar['FPARAM'][i1[rgb],1].min()
    loggmax=allstar['FPARAM'][i1[rgb],1].max()

    # RC fits
    # linear in logg and [M/H]
    print('rcfit',len(rc))
    rcfit = fit.fit2d(allstar['FPARAM'][i1[rc],1],allstar['FPARAM'][i1[rc],3],
        allstar[param][i1[rc],1]-apokasc[logg][i2[rc]],zr=[-1,0.5],gdrange=[-2,2],yr=[-3,1],xr=[1,4],degree=1,
        plot=ax[1],yt='[M/H]',xt='log g',zt='$\Delta log g$',reject=0.3)
    # quadratic in logg
    print('rcfit 2')
    rcfit2 = fit.fit1d(allstar['FPARAM'][i1[rcall],1], allstar[param][i1[rcall],1]-apokasc[logg][i2[rcall]],zr=[-1,0.5],yr=[-3,1],xr=[1,4],degree=2,reject=0.3)
    rcrms=(allstar[param][i1[rc],1]-rcfit(allstar['FPARAM'][i1[rc],1],allstar['FPARAM'][i1[rc],3])-apokasc[logg][i2[rc]]).std()
    rcerrpar = err.errfit(allstar[param][i1[rc],0],allstar['SNR'][i1[rc]],allstar[param][i1[rc],3],
                        allstar[param][i1[rc],1]-rcfit(allstar['FPARAM'][i1[rc],1],allstar['FPARAM'][i1[rc],3])-apokasc[logg][i2[rc]],
                        out=out+'_rc',title='log g',zr=[0,0.2])
    ax[1].text(0.98,0.98,'rms: {:5.3f}'.format(rcrms),transform=ax[1].transAxes,va='top',ha='right')
    fig.tight_layout()
    if out is not None :
        fig.savefig(out+'.png')
        plt.close()

    # set up plots
    if raw and plotcal :
        fig,ax=plots.multi(2,3,hspace=0.5,wspace=0.001,figsize=(12,12))
    else :
        fig,tmpax=plots.multi(1,4,hspace=0.5,wspace=0.001,figsize=(8,10))
        fig2,ax2=plots.multi(1,1)

    # diff color-coded by gravity as f([M/H])
    # diff color-coded by [M/H] as f(log g)
    # RGB and RC as f(log g)
    if raw :
        if plotcal: tmpax=ax[:,0]
        plots.plotc(tmpax[0],allstar['FPARAM'][i1,3],allstar[param][i1,1]-apokasc[logg][i2],
           allstar['FPARAM'][i1,1],zr=[0,5],xr=[-2.5,0.5],yr=[-0.75,0.75],xt='[M/H]',yt='ASPCAP-seismic log g',zt='log g',size=15,colorbar=True)
        plots.plotc(tmpax[1],allstar['FPARAM'][i1[rgb],1],allstar[param][i1[rgb],1]-apokasc[logg][i2[rgb]],
           allstar['FPARAM'][i1[rgb],3],xr=[0,5],zr=[-2.5,0.5],yr=[-0.75,0.75],zt='[M/H]',yt='ASPCAP-seismic log g',xt='log g',size=15,colorbar=True)
        loggfit=np.arange(1,3.5,0.01)
        mhfit=loggfit*0.
        plots.plotl(tmpax[1],loggfit,rgbfit(loggfit,mhfit),color='orange',linewidth=1.5)
        plots.plotl(tmpax[1],loggfit,params[0]+params[1]*loggfit+params[2]*loggfit**2+params[3]*loggfit**3+params[4]*mhfit,color='orange',linewidth=1.5)
        mhfit=loggfit*0-2.
        plots.plotl(tmpax[1],loggfit,rgbfit(loggfit,mhfit),color='b',linewidth=1.5)
        plots.plotl(tmpax[1],loggfit,params[0]+params[1]*loggfit+params[2]*loggfit**2+params[3]*loggfit**3+params[4]*mhfit,color='b',linewidth=1.5)
        mhfit=loggfit*0-1.
        plots.plotl(tmpax[1],loggfit,rgbfit(loggfit,mhfit),color='c',linewidth=1.5)
        plots.plotl(tmpax[1],loggfit,params[0]+params[1]*loggfit+params[2]*loggfit**2+params[3]*loggfit**3+params[4]*mhfit,color='c',linewidth=1.5)
        mhfit=loggfit*0+0.5
        plots.plotl(tmpax[1],loggfit,rgbfit(loggfit,mhfit),color='r',linewidth=1.5)
        plots.plotl(tmpax[1],loggfit,params[0]+params[1]*loggfit+params[2]*loggfit**2+params[3]*loggfit**3+params[4]*mhfit,color='r',linewidth=1.5)
        tmpax[0].grid()
        tmpax[1].grid()

        iax=tmpax[2]
        plots.plotp(iax,allstar['FPARAM'][i1[rgb],1],allstar[param][i1[rgb],1]-apokasc[logg][i2[rgb]],
           xr=[0,5],yr=[-0.5,0.5],xt='seismic log g',yt='ASPCAP-seismic log g',label=[0.9,0.8,'RGB'],color='r',size=5)
        plots.plotp(iax,allstar['FPARAM'][i1[rc],1],allstar[param][i1[rc],1]-apokasc[logg][i2[rc]],
           xr=[0,5],yr=[-0.5,0.5],xt='seismic log g',yt='ASPCAP-seismic log g',label=[0.9,0.6,'RC'],color='b',size=5)
        plots.plotp(iax,allstar['FPARAM'][i1[rc2],1],allstar[param][i1[rc2],1]-apokasc[logg][i2[rc2]],
           xr=[0,5],yr=[-0.5,0.5],xt='seismic log g',yt='ASPCAP-seismic log g',label=[0.9,0.6,'RC'],color='g',size=5)
        plots.plotp(iax,allstar['FPARAM'][i1[rc2cl],1],allstar[param][i1[rc2cl],1]-apokasc[logg][i2[rc2cl]],
           xr=[0,5],yr=[-0.5,0.5],xt='seismic log g',yt='ASPCAP-seismic log g',label=[0.9,0.6,'RC'],color='m',size=5)
        # single plot as f(Teff)
        iax=ax2
        plots.plotp(iax,allstar['FPARAM'][i1[rgb],0],allstar[param][i1[rgb],1]-apokasc[logg][i2[rgb]],
           xr=[4500,5200],yr=[-0.5,0.5],xt='Teff',yt='ASPCAP-seismic log g',label=[0.9,0.8,'RGB'],color='r',size=5)
        plots.plotp(iax,allstar['FPARAM'][i1[rc],0],allstar[param][i1[rc],1]-apokasc[logg][i2[rc]],
           xr=[4500,5200],yr=[-0.5,0.5],xt='Teff',yt='ASPCAP-seismic log g',label=[0.9,0.6,'RC'],color='b',size=5)
        plots.plotp(iax,allstar['FPARAM'][i1[rc2],0],allstar[param][i1[rc2],1]-apokasc[logg][i2[rc2]],
           xr=[4500,5200],yr=[-0.5,0.5],xt='Teff',yt='ASPCAP-seismic log g',label=[0.9,0.6,'RC'],color='g',size=15)
        plots.plotp(iax,allstar['FPARAM'][i1[rc2cl],0],allstar[param][i1[rc2cl],1]-apokasc[logg][i2[rc2cl]],
           xr=[4500,5200],yr=[-0.5,0.5],xt='Teff',yt='ASPCAP-seismic log g',label=[0.9,0.6,'RC'],color='m',size=5)
          

        loggfit=np.arange(2.5,3.5,0.01)
        mhfit=loggfit*0.
        plots.plotl(tmpax[2],loggfit,rcfit(loggfit,mhfit),color='g',linewidth=2)
        plots.plotl(tmpax[2],loggfit,rcfit2(loggfit),color='k',linewidth=2)
        tmpax[2].grid()

        #plots.plotp(tmpax[3],allstar['FPARAM'][i1[rgb],0],allstar[param][i1[rgb],1]-apokasc[logg][i2[rgb]],
        #   xr=[3500,7000],yr=[-0.75,0.75],xt='Teff',yt='ASPCAP-seismic log g',label=[0.9,0.8,'RGB'],color='r',size=15)
        plots.plotc(tmpax[3],allstar['FPARAM'][i1[rgb],0],allstar[param][i1[rgb],1]-apokasc[logg][i2[rgb]],allstar[param][i1[rgb],3],
           xr=[3500,7000],yr=[-0.75,0.75],xt='Teff',yt='ASPCAP-seismic log g',label=[0.9,0.8,'RGB'],zr=[-2,0.5],size=15,colorbar=True)
        plots.plotp(tmpax[3],allstar['FPARAM'][i1[rc],0],allstar[param][i1[rc],1]-apokasc[logg][i2[rc]],
           xr=[3500,7000],yr=[-0.75,0.75],xt='Teff',yt='ASPCAP-seismic log g',label=[0.9,0.6,'RC'],color='b',size=15)
        plots.plotp(tmpax[3],allstar['FPARAM'][i1[rc2],0],allstar[param][i1[rc2],1]-apokasc[logg][i2[rc2]],
           xr=[3500,7000],yr=[-0.75,0.75],xt='Teff',yt='ASPCAP-seismic log g',label=[0.9,0.6,'RC'],color='g',size=15)
        plots.plotp(tmpax[3],allstar['FPARAM'][i1[rc2cl],0],allstar[param][i1[rc2cl],1]-apokasc[logg][i2[rc2cl]],
           xr=[3500,7000],yr=[-0.75,0.75],xt='Teff',yt='ASPCAP-seismic log g',label=[0.9,0.6,'RC'],color='m',size=15)
        tmpax[3].grid()
        #plots.plotc(tmpax[3],allstar['FPARAM'][i1[rgb],1],allstar['PARAM'][i1[rgb],1]-allstar['FPARAM'][i1[rgb],1],
        #   allstar['FPARAM'][i1[rgb],3],xr=[0,5],yr=[-1,1],xt='seismic log g',yt='corrected-raw log g',label=[0.1,0.9,'allstar (Kurucz)'],zr=[-2,0.5])

    if plotcal :
        if raw: tmpax=ax[:,1]
        param=allstar['FPARAM'][:,1]-rgbfit(allstar['FPARAM'][:,1],allstar['FPARAM'][:,3])
        param[i1[rc]]=allstar['FPARAM'][i1[rc],1]-rcfit(allstar['FPARAM'][i1[rc],1],allstar['FPARAM'][i1[rc],3])

        plots.plotc(tmpax[0],allstar['FPARAM'][i1,3],param[i1]-apokasc[logg][i2],
           allstar['FPARAM'][i1,1],zr=[0,5],xr=[-2.5,0.5],yr=[-0.75,0.75],xt='[M/H]',colorbar=True,zt='log g',size=15)
        plots.plotc(tmpax[1],allstar['FPARAM'][i1,1],param[i1]-apokasc[logg][i2],
           allstar['FPARAM'][i1,3],xr=[1,4],zr=[-2.5,0.5],yr=[-0.75,0.75],zt='[M/H]',colorbar=True,xt='log g',size=15)
        plots.plotp(tmpax[2],allstar['FPARAM'][i1[rgb],1],param[i1[rgb]]-apokasc[logg][i2[rgb]],
           xr=[1,4],yr=[-0.75,0.75],xt='log g',color='r',size=15)
        plots.plotp(tmpax[2],allstar['FPARAM'][i1[rc],1],param[i1[rc]]-apokasc[logg][i2[rc]],
           xr=[1,4],yr=[-0.75,0.75],xt='log g',color='b',size=15)
        plots.plotp(tmpax[2],allstar['FPARAM'][i1[rc2],1],param[i1[rc2]]-apokasc[logg][i2[rc2]],
           xr=[1,4],yr=[-0.75,0.75],xt='log g',color='g',size=15)
        plots.plotp(tmpax[2],allstar['FPARAM'][i1[rc2cl],1],param[i1[rc2cl]]-apokasc[logg][i2[rc2cl]],
           xr=[1,4],yr=[-0.75,0.75],xt='log g',color='m',size=15)
        #plots.plotc(tmpax[3],allstar['FPARAM'][i1[rc],1],allstar['PARAM'][i1[rc],1]-allstar['FPARAM'][i1[rc],1],
        #    allstar['FPARAM'][i1[rc],3],xr=[0,5],yr=[-1,1],xt='seismic log g',zr=[-2,0.5])
    fig.tight_layout()
    if out is not None :
        fig.savefig(out+'_b.png')
        plt.close(fig)
        fig2.savefig(out+'_c.png')
        plt.close(fig2)

    return {'calloggmin' : calloggrange[0], 'calloggmax' : calloggrange[1], 'loggmin' : loggmin, 'loggmax' : loggmax, 
            'mhmin' : mhrange[0], 'mhmax' : mhrange[1], 'calteffmin': calteffrange[0], 'calteffmax' : calteffrange[1],
            'rgbfit' : rgbfit.parameters, 'rgbfit2' : params, 'rcfit' : rcfit.parameters, 'rcfit2' : rcfit2.parameters, 'rgbrms' : rgbrms, 'rcrms' : rcrms ,
            'rgberrpar': rgberrpar, 'rcerrpar': rcerrpar}

def isochrone(allstar,snrbd=300) :
    """ logg correction for cool dwarfs based on isochrones
        returns structured array with APOGEE_ID, ISOLOGG
    """

    print('getting isochrone log g')
    # restrict the sample to good high S/N stars
    aspcapmask=bitmask.AspcapBitMask()
    starmask=bitmask.StarBitMask()
    gd=np.where( ((allstar['ASPCAPFLAG']&aspcapmask.badval()) == 0) &
                 ((allstar['STARFLAG']&starmask.badval()) == 0) &
                  (allstar['SNR']>=snrbd) ) [0]
    allstar=allstar[gd]
    if 'TARGFLAGS' in allstar.columns.names : badtarg=['YOUNG','EMBEDDED','EXTENDED','M31','M33','EMISSION','RRLYR','DSPH','MAGCLOUD']
    else : badtarg = None

    gd=apselect.select(allstar,raw=True,teff=[3000,5000],logg=[4.0,5.5],badtarg=badtarg)
    allstar=allstar[gd]
    print(len(allstar))

    # loop through isochrones, reading, finding matches, and calculating expected isochrone logg given Teff
    first=True
    for z in np.arange(-1.0,0.3,0.1) :
        if z<-0.01 : name='zm{:02d}'.format(int(abs(z)*10))
        else :name='zp{:02d}'.format(int(abs(z)*10))
        j=np.where(abs(allstar['FPARAM'][:,3]-z) <0.05)[0]
        if len(j) > 0: 
            print(z,len(j),name)
            isodata=isochrones.read(os.environ['ISOCHRONE_DIR']+'/'+name+'.dat',agerange=[9.29,9.31])
            mdiff = isodata['mini'][0:-1]-isodata['mini'][1:]
            use=np.where(abs(mdiff) < 1.e-8)[0]
            if len(use) > 0 : use=use[0]
            else : use=len(isodata)
            if use < 10 : pdb.set_trace()
            gd=np.where(isodata['logg'][0:use]>4)[0]
            f = interpolate.interp1d(isodata['teff'][gd], isodata['logg'][gd],bounds_error=False)
            isologg = f(allstar['FPARAM'][j,0])
            if first :
                out_id=allstar['APOGEE_ID'][j]
                out_isologg=isologg
                first= False
            else :
                out_id=np.append(out_id,allstar['APOGEE_ID'][j])
                out_isologg=np.append(out_isologg,isologg)

    # output structured array
    outtype=np.dtype([('APOGEE_ID',out_id.dtype),('ISOLOGG',isologg.dtype)])
    outdata=np.empty(len(out_id),dtype=outtype)
    outdata['APOGEE_ID']=out_id
    outdata['ISOLOGG']=out_isologg
    return outdata
            

def clusters(allstar,xr=[-2.75,0.5],yr=[-1.,1.],zr=[3500,5500],apokasc='APOKASC_cat_v3.6.0.fits',firstgen=False,
             clusters=['M92','M15','M53','M2','M13','M3','M5','N2420','M67','N6819','N6791'],calib=True,title=None,
             out='cluster_logg') :
    '''
    Compare ASPCAP gravities in clusters to physical gravities
    '''
    if calib : fig,ax=plots.multi(1,2,hspace=0.001)
    else : fig,ax=plots.multi(1,1,hspace=0.001)
    ax=np.atleast_1d(ax)

    """
    # put APOKASC underneath
    apokasc=fits.open(os.environ['APOGEE_DIR']+'/data/apokasc/'+apokasc)[1].data
    # strip off .XXXX if we have it, e.g. from calibration fields where we have added .FIELD
    apogee_id = np.array(np.core.defchararray.split(allstar['APOGEE_ID'],'.').tolist())[:,0]
    i1,i2=match.match(apogee_id,apokasc['2MASS_ID'])
    plots.plotc(ax[0],allstar['FPARAM'][i1,3],allstar['FPARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],allstar['FPARAM'][i1,0],zr=zr)
    plots.plotc(ax[1],allstar['PARAM'][i1,3],allstar['PARAM'][i1,1]-apokasc['LOGG_SYD_SCALING'][i2],allstar['PARAM'][i1,0],zr=zr)
    """

    # physical gravities
    clust=apselect.clustdata()
    itext=0
    txt=open('clust.txt','w')
    for cluster in clusters :
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
            axim=plots.plotc(ax[0],allstar['FPARAM'][j,3]*0+mh,allstar['FPARAM'][j,1]-logg,allstar['FPARAM'][j,0],xr=xr,yr=yr,zr=zr,yt='ASPCAP-physical log g')
            ax[0].text(0.9,0.1,'raw',transform=ax[0].transAxes,ha='right')

            #plots.plotp(ax[0],allstar['FPARAM'][j,3]*0+mh,allstar['FPARAM'][j,1]-logg,color='k')
            plots.plotp(ax[0],mh[0],np.median(allstar['FPARAM'][j,1]-logg),size=50,color='k')
            ax[0].text(mh[0],ytext,name[0],ha='center')
            ax[0].grid()

            txt.write('{:<20s}{:8.3f}{:8.3f}{:8.3f}\n'.format(clust[i].name[0],clust[i].dist[0],clust[i].ebv[0],mass[0]))

            if calib :
                gd=np.where((allstar['PARAM'][j,3]>-9)&(allstar['PARAM'][j,1]>-9))[0]
                axim=plots.plotc(ax[1],allstar['PARAM'][j[gd],3]*0+mh,allstar['PARAM'][j[gd],1]-logg[gd],allstar['PARAM'][j[gd],0],xr=xr,yr=yr,zr=zr,xt='[M/H]',yt='ASPCAP-physical log g')
                ax[1].text(0.9,0.1,'calibrated',transform=ax[1].transAxes,ha='right')

                plots.plotp(ax[1],mh[0],np.median(allstar['PARAM'][j[gd],1]-logg[gd]),size=40)
                # apply a temperature correction for the physical gravities
                logg_new=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*astropy.constants.sigma_sb.cgs.value*(allstar['FPARAM'][j,0]-100.*allstar['FPARAM'][j,3])**4/lum)
                plots.plotp(ax[1],mh[0],np.median(allstar['PARAM'][j,1]-logg_new),size=40,color='b')
                # use a photometric temperature
                logg_phot=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*astropy.constants.sigma_sb.cgs.value*teff.cte_ghb(allstar['J'][j]-allstar['K'][j]-ejk,allstar['FPARAM'][j,3])[0]**4/lum)
                plots.plotp(ax[1],mh[0],np.median(allstar['PARAM'][j,1]-logg_phot),size=40,color='g')
                ax[1].text(mh[0],ytext,name[0],ha='center')
            itext+=1

    if title is not None : fig.suptitle(title)

    # Now adding the colorbar
    cbaxes = fig.add_axes([0.9, 0.1, 0.03, 0.8]) 
    cb = plt.colorbar(axim, cax = cbaxes)  
    txt.close()
    if out is not None : fig.savefig(out+'.png')


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
        lum=10.**(-0.4*(dr12['H'][j]+isochrones.bc(dr12['FPARAM'][j,0],filt='h',
            agerange=[10,10.1])-(5*np.log10(dist)-5)-4.74))*astropy.constants.L_sun.cgs.value
        logg=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*
            astropy.constants.sigma_sb.cgs.value*dr12['FPARAM'][j,0]**4/lum)
        plots.plotc(ax[0,0],dr12['FPARAM'][j,3]*0+mh,dr12['FPARAM'][j,1]-logg,dr12['FPARAM'][j,0],
                    xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],yt='ASPCAP-physical log g',label=[0.1,0.9,'DR2 raw'])
        plots.plotc(ax[1,0],dr12['PARAM'][j,3]*0+mh,dr12['PARAM'][j,1]-logg,dr12['PARAM'][j,0],
                    xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],xt='[M/H]',yt='ASPCAP-physical log g',label=[0.1,0.9,'DR12 cal'])

        #DR13
        j=apselect.clustmember(dr13,cluster,raw=True)
        lum=10.**(-0.4*(dr13['H'][j]+isochrones.bc(dr13['FPARAM'][j,0],filt='h',
            agerange=[10,10.1])-(5*np.log10(dist)-5)-4.74))*astropy.constants.L_sun.cgs.value
        logg=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*
            astropy.constants.sigma_sb.cgs.value*dr13['FPARAM'][j,0]**4/lum)
        plots.plotc(ax[0,1],dr13['FPARAM'][j,3]*0+mh,dr13['FPARAM'][j,1]-logg,dr13['FPARAM'][j,0],
            xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],label=[0.1,0.9,'DR13 raw'])
        plots.plotc(ax[1,1],dr13['PARAM'][j,3]*0+mh,dr13['PARAM'][j,1]-logg,dr13['PARAM'][j,0],
            xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],label=[0.1,0.9,'DR13 cal'],xt='[M/H]')

        #l30i
        j=apselect.clustmember(l30i,cluster,raw=True)
        lum=10.**(-0.4*(l30i['H'][j]+isochrones.bc(l30i['FPARAM'][j,0],filt='h',
             agerange=[10,10.1])-(5*np.log10(dist)-5)-4.74))*astropy.constants.L_sun.cgs.value
        logg=np.log10(4*np.pi*astropy.constants.G.cgs.value*mass*astropy.constants.M_sun.cgs.value*
             astropy.constants.sigma_sb.cgs.value*l30i['FPARAM'][j,0]**4/lum)
        plots.plotc(ax[0,2],l30i['FPARAM'][j,3]*0+mh,l30i['FPARAM'][j,1]-logg,l30i['FPARAM'][j,0],
             xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],label=[0.1,0.9,'l30i raw'])
        plots.plotc(ax[1,2],l30i['PARAM'][j,3]*0+mh,l30i['PARAM'][j,1]-logg,l30i['PARAM'][j,0],
             xr=[-2.75,-1.],yr=[-2,2],zr=[3500,5500],label=[0.1,0.9,'l30i cal'],xt='[M/H]')
    
    plt.show()
    pdb.set_trace()

    # plots vs asteroseismic
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

def cal(a,caldir='cal/') :
    """ apply log g calibration
    """

    # get bitmask, and select stars without STAR_BAD
    aspcapmask=bitmask.AspcapBitMask()
    parammask=bitmask.ParamBitMask()
    starmask=bitmask.StarBitMask()
    gd=np.where( ((a['ASPCAPFLAG']&aspcapmask.badval()) == 0) )[0]

    # start with CALRANGE_BAD
    a['PARAM'][:,1] = np.nan
    a['PARAMFLAG'][gd,1] |= parammask.getval('CALRANGE_BAD')

    # get calibration data
    calpars=fits.open(caldir+'/giant_loggcal.fits')[1].data
    rgbsep=calpars['rgbsep'][0]
    cnsep=calpars['cnsep'][0]
    rclim=calpars['rclim'][0]
    rcfit2=calpars['rcfit2'][0]
    rgbfit2=calpars['rgbfit2'][0]
    calloggmin=calpars['calloggmin']
    calloggmax=calpars['calloggmax']
    calteffmin=calpars['calteffmin']
    calteffmax=calpars['calteffmax']


    # for stars that aren't bad, get cn and dt
    cn=a['FPARAM'][gd,4]-a['FPARAM'][gd,5]
    dt=a['FPARAM'][gd,0] - (rgbsep[0] + rgbsep[1]*(a['FPARAM'][gd,1]-2.5) +rgbsep[2]*a['FPARAM'][gd,3])
    try: snr=np.clip(a['SNREV'][gd],0,200.)
    except:
        print('No SNREV, continnue with SNR?')
        pdb.set_trace()
        snr=np.clip(a['SNR'][gd],0,200.)

    # select RC
    rc=np.where((a['FPARAM'][gd,1]<rclim[1])&(a['FPARAM'][gd,1]>rclim[0])&
                (cn>cnsep[0]+cnsep[1]*a['FPARAM'][gd,3] + cnsep[2]*dt)&
                (a['FPARAM'][gd,1]<calloggmax)&(a['FPARAM'][gd,1]>calloggmin) &
                (a['FPARAM'][gd,0]<calteffmax)&(a['FPARAM'][gd,0]>calteffmin))[0]
    rccorr=rcfit2[0] + rcfit2[1]*a['FPARAM'][gd,1] + rcfit2[2]*a['FPARAM'][gd,1]**2
    a['PARAM'][gd[rc],1]=a['FPARAM'][gd[rc],1]-rccorr[rc]
    a['PARAM_COV'][gd[rc],1,1]=err.elemerr(calpars['rcerrpar'][0],a['FPARAM'][gd[rc],0]-4500,snr[rc]-100,a['FPARAM'][gd[rc],3])**2
    a['PARAMFLAG'][gd[rc],1] &= ~parammask.getval('CALRANGE_BAD')
    a['PARAMFLAG'][gd[rc],1] |= parammask.getval('LOGG_CAL_RC')

    # select RGB
    rgb=np.where(((a['FPARAM'][gd,1]>rclim[1])|(a['FPARAM'][gd,1]<rclim[0])|
                (cn<cnsep[0]+cnsep[1]*a['FPARAM'][gd,3] + cnsep[2]*dt)) &
                (a['FPARAM'][gd,1]<calloggmax)&(a['FPARAM'][gd,1]>calloggmin) &
                (a['FPARAM'][gd,0]<calteffmax)&(a['FPARAM'][gd,0]>calteffmin))[0]
    #clip logg at loggmin and loggmax
    logg=np.clip(a['FPARAM'][gd,1],calpars['loggmin'],calpars['loggmax'])
    mh=np.clip(a['FPARAM'][gd,3],calpars['mhmin'],calpars['mhmax'])
    # get correction
    rgbcorr=(rgbfit2[0] + rgbfit2[1]*logg + rgbfit2[2]*logg**2 +
                       rgbfit2[3]*logg**3 + rgbfit2[4]*mh )
    a['PARAM'][gd[rgb],1]=a['FPARAM'][gd[rgb],1]-rgbcorr[rgb]
    a['PARAM_COV'][gd[rgb],1,1]=err.elemerr(calpars['rgberrpar'][0],a['FPARAM'][gd[rgb],0]-4500,snr[rgb]-100,a['FPARAM'][gd[rgb],3])**2
    a['PARAMFLAG'][gd[rgb],1] &= ~parammask.getval('CALRANGE_BAD')
    a['PARAMFLAG'][gd[rc],1] |= parammask.getval('LOGG_CAL_RGB')

    # dwarfs
    calpars=fits.open(caldir+'/dwarf_loggcal.fits')[1].data
    teff=np.clip(a['FPARAM'][gd,0],calpars['temin'],calpars['temax'])
    logg=np.clip(a['FPARAM'][gd,1],calpars['loggmin'],calpars['loggmax'])
    mh=np.clip(a['FPARAM'][gd,3],calpars['mhmin'],calpars['mhmax'])
    msfit=calpars['msfit'][0]
    mscorr=msfit[0]+msfit[1]*teff+msfit[2]*mh
    ms=np.where(a['FPARAM'][gd,1] > calpars['calloggmin'])[0]
    a['PARAM'][gd[ms],1]=a['FPARAM'][gd[ms],1]-mscorr[ms]
    a['PARAM_COV'][gd[ms],1,1]=err.elemerr(calpars['errpar'][0],a['FPARAM'][gd[ms],0]-4500,snr[ms]-100,a['FPARAM'][gd[ms],3])**2
    a['PARAMFLAG'][gd[ms],1] &= ~parammask.getval('CALRANGE_BAD')
    a['PARAMFLAG'][gd[rc],1] |= parammask.getval('LOGG_CAL_MS')

    # dwrarf-RGB transition
    trans=np.where((a['FPARAM'][gd,1] < 4) & (a['FPARAM'][gd,1] > 3.5) &
                (a['FPARAM'][gd,0] < calteffmax) )[0]
    ms_weight=(a['FPARAM'][gd[trans],1]-3.5)/0.5
    a['PARAM'][gd[trans],1] = a['FPARAM'][gd[trans],1]-(mscorr[trans]*ms_weight+rgbcorr[trans]*(1-ms_weight))
    a['PARAMFLAG'][gd[trans],1] &= ~parammask.getval('CALRANGE_BAD')
    a['PARAMFLAG'][gd[rc],1] |= parammask.getval('LOGG_CAL_RGB_MS')
    a['PARAMFLAG'][gd[trans],1] &= ~parammask.getval('LOGG_CAL_RC')
    a['PARAMFLAG'][gd[trans],1] &= ~parammask.getval('LOGG_CAL_RGB')
    a['PARAMFLAG'][gd[trans],1] &= ~parammask.getval('LOGG_CAL_MS')

    return 


