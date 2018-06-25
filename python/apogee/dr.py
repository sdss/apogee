#
# Routine to make plots comparing parameters and abundances for dr papers
#
import pdb
import numpy as np
import matplotlib.pyplot as plt
from apogee.utils import apload
from tools import match
from tools import plots

def dr_compare() :
    # load the DRs, select stars with SN>150
    apload.dr12()
    dr12=apload.allStar()[1].data
    gd=np.where(dr12['SNR'] > 150)[0]
    dr12=dr12[gd]
    apload.dr13()
    dr13=apload.allStar()[1].data
    gd=np.where(dr13['SNR'] > 150)[0]
    dr13=dr13[gd]
    apload.dr14()
    dr14=apload.allStar()[1].data
    gd=np.where(dr14['SNR'] > 150)[0]
    dr14=dr14[gd]
    c=apload.allStar()[3].data

    # match them
    m1a,m2a=match.match(dr12['APOGEE_ID'],dr13['APOGEE_ID'])
    m1b,m2b=match.match(dr12['APOGEE_ID'],dr14['APOGEE_ID'])
    m1c,m2c=match.match(dr13['APOGEE_ID'],dr14['APOGEE_ID'])

    # parameter figures
    figu,axu=plots.multi(3,7,hspace=0.001,wspace=0.001)
    figc,axc=plots.multi(3,7,hspace=0.001,wspace=0.001)

    tit=[r'T$_{\rm eff}$','log g',r'V$_{\rm micro}$','[M/H]','[C/M]','[N/M]',r'[$\alpha$/M]']
    for iparam in range(7) :
    
      print(iparam)
      for iy,param in enumerate(['FPARAM','PARAM']) :
        if iy == 0 :
          ax=axu
        else :
          ax=axc
        yt=r'$\Delta$'+tit[iparam]
        if iparam == 6 : xt=r'T$_{\rm eff}$'
        else : xt=None
        if iparam == 0 :
          ax[iparam,0].text(0.5,1.0,'DR13-DR12',transform=ax[iparam,0].transAxes,ha='center',va='bottom')
          ax[iparam,1].text(0.5,1.0,'DR14-DR12',transform=ax[iparam,1].transAxes,ha='center',va='bottom')
          ax[iparam,2].text(0.5,1.0,'DR14-DR13',transform=ax[iparam,2].transAxes,ha='center',va='bottom')

        if iparam == 0 :
            yr=[-300,300]
        elif iparam == 1 : 
            yr=[-0.5,0.5]
        else :
            yr=[-0.3,0.3]

        xr=[3500,6000]
 
        axim = plots.plotc(ax[iparam,0],dr12['TEFF'][m1a],dr13[param][m2a,iparam]-dr12[param][m1a,iparam],dr12[param][m1a,3],size=1,xr=xr,yr=yr,zr=[-1,0.5],yt=yt,xt=xt,rasterized=True)
        plots.plotl(ax[iparam,0],xr,[0.,0.],ls=':')
        plots.plotc(ax[iparam,1],dr12['TEFF'][m1b],dr14[param][m2b,iparam]-dr12[param][m1b,iparam],dr12[param][m1b,3],size=1,xr=xr,yr=yr,zr=[-1,0.5],xt=xt,rasterized=True)
        plots.plotl(ax[iparam,1],xr,[0.,0.],ls=':')
        plots.plotc(ax[iparam,2],dr13['TEFF'][m1c],dr14[param][m2c,iparam]-dr13[param][m1c,iparam],dr13[param][m1c,3],size=1,xr=xr,yr=yr,zr=[-1,0.5],xt=xt,rasterized=True)
        plots.plotl(ax[iparam,2],xr,[0.,0.],ls=':')
        for iax in range(3) :
          ax[iparam,iax].tick_params(axis='both',labelsize=8)

    # add colorbar
    for fig in [figu, figc] :
        cbaxes = fig.add_axes([0.91, 0.1, 0.01, 0.8])
        cb = plt.colorbar(axim, cax = cbaxes)
        cb.set_label('[M/H]')
        cbaxes.tick_params(axis='both',labelsize=8)

    figu.savefig('drcomp_uncal.pdf')
    figc.savefig('drcomp_cal.pdf')
    plots.close()

    # abundance figure
    fig,ax=plots.multi(3,14,hspace=0.001,wspace=0.001,figsize=(8,32))

    for ielem,elem in enumerate(['C','N','O','Na','Mg','Al','Si','S','K','Ca','Ti','V','Mn','Ni']) :
        print(elem)
        yt=r'$\Delta$'+elem
        if ielem == 13 : xt=r'T$_{\rm eff}$'
        else : xt=None
        if ielem == 0 :
          ax[ielem,0].text(0.5,1.0,'DR13-DR12',transform=ax[ielem,0].transAxes,ha='center',va='bottom')
          ax[ielem,1].text(0.5,1.0,'DR14-DR12',transform=ax[ielem,1].transAxes,ha='center',va='bottom')
          ax[ielem,2].text(0.5,1.0,'DR14-DR13',transform=ax[ielem,2].transAxes,ha='center',va='bottom')

        yr=[-0.5,0.5]

        dr12elem=dr12[elem.upper()+'_H'][m1a]-dr12['FE_H'][m1a]
        dr13elem=dr13[elem.upper()+'_FE'][m2a]
        gd=np.where((dr12elem > -99) & (dr13elem>-99))[0]
        plots.plotc(ax[ielem,0],dr12['TEFF'][m1a[gd]],dr13elem[gd]-dr12elem[gd],dr12['PARAM'][m1a[gd],3],size=1,xr=[3500,6000],yr=yr,zr=[-1,0.5],yt=yt,xt=xt,nytick=5,rasterized=True)
        plots.plotl(ax[ielem,0],xr,[0.,0.],ls=':')
        ax[ielem,0].tick_params(axis='both',labelsize=8)
  
        dr12elem=dr12[elem.upper()+'_H'][m1b]-dr12['FE_H'][m1b]
        dr14elem=dr14[elem.upper()+'_FE'][m2b]
        gd=np.where((dr12elem > -99) & (dr14elem>-99))[0]
        plots.plotc(ax[ielem,1],dr12['TEFF'][m1b[gd]],dr14elem[gd]-dr12elem[gd],dr12['PARAM'][m1b[gd],3],size=1,xr=[3500,6000],yr=yr,zr=[-1,0.5],xt=xt,nytick=5,rasterized=True)
        plots.plotl(ax[ielem,1],xr,[0.,0.],ls=':')
        ax[ielem,1].tick_params(axis='both',labelsize=8)

        dr13elem=dr13[elem.upper()+'_FE'][m1c]
        dr14elem=dr14[elem.upper()+'_FE'][m2c]
        gd=np.where((dr13elem > -99) & (dr14elem>-99))[0]
        plots.plotc(ax[ielem,2],dr13['TEFF'][m1c[gd]],dr14elem[gd]-dr13elem[gd],dr13['PARAM'][m1c[gd],3],size=1,xr=[3500,6000],yr=yr,zr=[-1,0.5],xt=xt,nytick=5,rasterized=True)
        plots.plotl(ax[ielem,2],xr,[0.,0.],ls=':')
        ax[ielem,2].tick_params(axis='both',labelsize=8)

    cbaxes = fig.add_axes([0.91, 0.1, 0.01, 0.8])
    cb = plt.colorbar(axim, cax = cbaxes)
    cb.set_label('[M/H]')
    cbaxes.tick_params(axis='both',labelsize=8)

    for item in (cbaxes.get_xticklabels() + cbaxes.get_yticklabels()) : item.set_fontsize(8)
    fig.savefig('drcomp_elem.pdf')

def kurucz_marcs() :

    apload.dr13()
    dr13=apload.allStar()[1].data
    gd=np.where(dr13['SNR'] > 150)[0]
    dr13=dr13[gd]

    apload.aspcap = 'l30g'
    dr13_marcs=apload.allStar()[1].data
    gd=np.where(dr13_marcs['SNR'] > 150)[0]
    dr13_marcs=dr13_marcs[gd]

    fig,ax=plots.multi(2,1,wspace=0.001)
    plots.plotc(ax[0],dr13['FPARAM'][:,0],dr13['FPARAM'][:,1],dr13['FPARAM'] [:,3],
                xr=[4200,3000],yr=[5,-1],zr=[-2,0.5],xt=r'T$_{\rm eff}$',yt='log g',rasterized=True)
    plots.plotc(ax[1],dr13_marcs['FPARAM'][:,0],dr13_marcs['FPARAM'][:,1],dr13_marcs['FPARAM'] [:,3],
                xr=[4200,3000],yr=[5,-1],zr=[-2,0.5],xt='Teff',rasterized=True)
    fig.savefig('kurucz_marcs.pdf')

