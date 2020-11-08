import matplotlib.pyplot as plt
import numpy as np
from tools import plots
import pdb
import glob
import os
from astropy.io import fits
from tools import html

def fwhm(mjd=[0,99999],apred='current',telescope='apo25m'):
    """ FWHM as a function of MJD from Sci.fits summary files
    """
    if telescope == 'lco25m' :
        a=fits.open(os.environ['APOGEE_REDUX']+'/'+apred+'/apogee-nSci'.fits)
    else :
        a=fits.open(os.environ['APOGEE_REDUX']+'/'+apred+'/apogee-sSci'.fits)

    thar = np.where((a['THAR'] == 1) )[0]
    gd = np.where((a['THAR'] == 1) & (a['JD']>mjd[0]+2400000) & (a['JD']<=mjd[1]+2400000) )[0]
    xr= [a['JD'][gd].min()-2400000-10,a['JD'][gd].max()-2400000+10]
    print(len(gd))
    colors = ['r','g','b','m','c']
    fibers=[290,220,150,80,10]
    fig,ax = plots.multi(2,3,hspace=0.001,wspace=0.001,sharex=True,figsize=(14,8))
    for iline in range(2) :
        if iline == 0 : yt='FWHM'
        else : yt=''
        for ichip in range(3) :
            # use same yr for both lines
            med = np.median(a['GAUSS'][gd,0,ichip,:,2]*2.354)
            yr= [med-1,med+1]
            for ifiber in range(5) :
                plots.plotp(ax[ichip,iline],a['JD'][gd]-2400000,a['GAUSS'][gd,iline,ichip,ifiber,2]*2.354,color=colors[ifiber],yr=yr,xr=xr,
                            label='Fiber {:d}'.format(fibers[ifiber]),xt='MJD',yt=yt)
                allgd = np.where((a['GAUSS'][thar,iline,ichip,ifiber,2]*2.354 > 0.5) & (a['GAUSS'][thar,iline,ichip,ifiber,2]*2.354 < 5) )[0]
                allmed = np.median(a['GAUSS'][thar[allgd],iline,ichip,ifiber,2]*2.354)
                print(iline,ichip,ifiber,allmed)
                plots.plotl(ax[ichip,iline],ax[ichip,iline].get_xlim(),[allmed,allmed],color=colors[ifiber],xr=xr,yr=yr)
            ax[ichip,iline].legend(fontsize='xx-small',loc='upper left')
        plt.draw()
    fig.suptitle('FWHM of two (left/right) ThArNe lines ')
    pdb.set_trace()

def flats(mjd=[58000,59155],apred='current',inst='apogee-n',clobber=False) :
    """ Look at flux frames (from petal and internal flats) levels as a function  of fiber
    """
    grid=[]
    if inst == 'apogee-s' : prefix = 'as'
    else : prefix = 'ap'
    rms=[]
    jd=[]
    cartid=[]
    for m in range(mjd[0],mjd[1]) :
        print(m)
        dayno = m-55562
        files=np.sort(glob.glob(os.environ['APOGEE_REDUX']+'/'+apred+'/cal/flux/'+prefix+'Flux-b-'+'{:04d}'.format(dayno)+'*.fits'))
        if len(files) > 0 :
            for i,file in enumerate(files) :
                a=fits.open(file)
                ngd=np.where(np.isfinite(a[1].data[:,1000]))[0]
                if len(ngd) > 150 :
                    jd.append(a[0].header['JD-MID'])
                    rms.append(np.nanstd(a[1].data[:,1000]))
                    try:
                        if a[0].header['EXPTYPE'] == 'QUARTZFLAT' : cart = 0
                        else : cart = a[0].header['CARTID']
                    except: cart=0
                    cartid.append(cart)
           
            if clobber or not os.path.exists(os.environ['APOGEE_REDUX']+'/'+apred+'/cal/flux/plots/'+prefix+'Flux-'+str(m)+'.png') :
                fig,ax=plots.multi(1,len(files),figsize=(6,2+len(files)),wspace=0.001,hspace=0.001)
                ax=np.atleast_1d(ax)
                for i,file in enumerate(files) :
                    a=fits.open(file)
                    print(file)
                    plots.plotp(ax[i],300-np.arange(300),a[1].data[:,1000],yr=[0.25,1.75],xt='Fiber')
                    try:
                        if a[0].header['EXPTYPE'] == 'QUARTZFLAT' :
                            plate = 'QUARTZFLAT'
                            cart = 0
                        else :
                            plate = 'PLATE: '+str(a[0].header['PLATEID'])
                            cart =  'CART: '+str(a[0].header['CARTID'])
                        ax[i].text(0.1,0.1,plate,transform=ax[i].transAxes)
                        ax[i].text(0.8,0.1,cart,transform=ax[i].transAxes)
                    except: pass
                fig.suptitle(m)
                fig.savefig(os.environ['APOGEE_REDUX']+'/'+apred+'/cal/flux/plots/'+prefix+'Flux-'+str(m)+'.png')
            grid.append([prefix+'Flux-'+str(m)+'.png'])
            plt.close()
    fig,ax=plots.multi(1,1,figsize=(16,8))
    plots.plotc(ax,np.array(jd),np.array(rms),np.array(cartid),xt='JD',yt='flat rms',yr=[0,1],colorbar=True,zt='Cart')
    fig.savefig(os.environ['APOGEE_REDUX']+'/'+apred+'/monitor/'+inst+'/'+prefix+'Flux.png')
    html.htmltab(grid,file=os.environ['APOGEE_REDUX']+'/'+apred+'/cal/flux/plots/'+prefix+'Flux.html')
    return jd,rms,cartid

