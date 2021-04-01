import os
import numpy as np
import matplotlib.pyplot as plt
import pdb
from numpy.random import normal

from tools import plots
from tools import html
from tools import fit
from apogee.utils import apselect, bitmask
from apogee.aspcap import aspcap
from astropy.io import fits

def apply(data,caldir=None) :
    """ Add empirical uncertainties to data
    """

    if caldir is None or caldir == 'none' :
        print('no uncertainty calibration specified, not populating')
        return

    aspcapmask=bitmask.AspcapBitMask()
    parammask=bitmask.ParamBitMask()
    gd=np.where( ((data['ASPCAPFLAG']&aspcapmask.badval()) == 0) )[0]

    # separate giants and dwarfs
    giant = np.where( (data['FPARAM'][gd,1] < 2./1300.*(data['FPARAM'][gd,0]-3500)+2.) &
                      (data['FPARAM'][gd,1] < 3.8) & (data['FPARAM'][gd,0] < 5500) )[0]
    tmp = np.zeros(len(gd),dtype=bool)
    tmp[giant] = True
    dwarf = np.where(~tmp)[0]
    giant = gd[giant]
    dwarf = gd[dwarf]

    # cap S/N at 200.
    try: snr=np.clip(data['SNREV'],0,200.)
    except:
        print('No SNREV, continnue with SNR?')
        pdb.set_trace()
        snr=np.clip(a['SNR'],0,200.)

    # read uncertainty fit parameters
    giant_errfit=fits.open(caldir+'/giant_errfit.fits')
    dwarf_errfit=fits.open(caldir+'/dwarf_errfit.fits')

    # add uncertainties for parameters
    fig,ax=plots.multi(1,1)
    for iparam,param in enumerate(aspcap.params()[1]) :
        data['PARAM_COV'][giant,iparam,iparam] = elemerr(giant_errfit[1].data['errfit'][iparam],data['FPARAM'][giant,0]-4500.,snr[giant]-100.,data['FPARAM'][giant,3])**2
        data['PARAM_COV'][dwarf,iparam,iparam] = elemerr(dwarf_errfit[1].data['errfit'][iparam],data['FPARAM'][dwarf,0]-4500.,snr[dwarf]-100.,data['FPARAM'][dwarf,3])**2
        j=np.where(data['FPARAM_COV'][:,iparam,iparam] > data['PARAM_COV'][:,iparam,iparam])[0]
        data['PARAM_COV'][j,iparam,iparam] = data['FPARAM_COV'][j,iparam,iparam]
        data['PARAMFLAG'][j,iparam] |= parammask.getval('FERRE_ERR_USED')
        print(param,'FERRE ERR used: ',len(j))

    # add uncertainties for abundances
    for ielem,elem in enumerate(aspcap.elems()[0]) :
        data['X_H_ERR'][giant,ielem] = elemerr(giant_errfit[2].data['errfit'][ielem],data['FPARAM'][giant,0]-4500.,snr[giant]-100.,data['FPARAM'][giant,3])
        data['X_H_ERR'][dwarf,ielem] = elemerr(dwarf_errfit[2].data['errfit'][ielem],data['FPARAM'][dwarf,0]-4500.,snr[dwarf]-100.,data['FPARAM'][dwarf,3])
        # use FERRE uncertainty if larger
        j=np.where(data['FELEM_ERR'][:,ielem] > data['X_H_ERR'][:,ielem])[0]
        data['X_H_ERR'][j,ielem] = data['FELEM_ERR'][j,ielem]
        data['ELEMFLAG'][j,ielem] |= parammask.getval('FERRE_ERR_USED')
        print(elem,'FERRE ERR used: ',len(j))
        data['X_M_ERR'][:,ielem] = data['X_H_ERR'][:,ielem]


def elemerr(soln,te,sn,fe, quad=False, log=True, fact=1.0) :
    ''' 
    Function to evaluate fit for uncertainty
    '''
    out=soln[0]+soln[1]*te+soln[2]*sn
    if len(soln) > 3: out+= soln[3]*fe
    if quad : out +=soln[4]*te**2
    if log : return np.exp(out)*fact
    else : return out

def errfit(te, snr, mh, val, snbins=np.arange(50,250,50), tebins=np.arange(3500,6000,250), mhbins=np.arange(-2.25,0.75,0.5),verbose=False,
           out=None,title='', zr=[0,0.1], snplot=True, meanerr=None,quad=False,mkhtml=True ) :
    '''
    Fits for empirical uncertainty as function of Teff, S/N, and [M/H] given measured values of a quantity
    '''
    if out is not None :
        fig,ax=plots.multi(len(snbins),2,wspace=0.001,figsize=(3*len(snbins),5))

    # bin sizes and initialize data arrays
    dte = tebins[1]-tebins[0]
    dmh = mhbins[1]-mhbins[0]
    dsn = snbins[1]-snbins[0]
    rmsdata=[]
    rmsderiv=[]
    nbin=[]
    # accumulate the data: rms in bins of Teff, S/N, and [M/H]
    snmin=snbins[-1]
    snmax=snbins[0]
    temin=tebins[-1]
    temax=tebins[0]
    mhmin=mhbins[-1]
    mhmax=mhbins[0]
    npts=np.zeros([len(tebins),len(mhbins),len(snbins)])
    for imhbin,mhbin in enumerate(mhbins) :
        for itebin,tebin in enumerate(tebins) :
            for isnbin,snbin in enumerate(snbins) :
                ibin = np.where(( te > tebin) & (te <= tebin+dte) &
                                ( mh > mhbin ) & (mh <= mhbin+dmh) &
                                ( snr > snbin) & (snr <= snbin+dsn) & (val > -9990.) )[0]
                if len(ibin) > 3 :
                    npts[itebin,imhbin,isnbin] = len(ibin)
                    if meanerr is not None :
                        err = np.sqrt(np.clip(val[ibin].std()**2 - np.median(meanerr[ibin])**2,0.001,10000000.))
                    else :
                        err = val[ibin].std()
                    rmsdata.append(err)
                    if quad :
                        rmsderiv.append([1.,tebin+dte/2.-4500.,snbin+dsn/2.-100.,mhbin+dmh/2.,(tebin+dte/2.-4500.)**2])
                    else :
                        rmsderiv.append([1.,tebin+dte/2.-4500.,snbin+dsn/2.-100.,mhbin+dmh/2.])
                    if verbose :
                        print(tebin+dte/2.,snbin+dsn/2.,mhbin+dmh/2.,err,len(ibin))
                    snmin=np.array([snmin,snbin]).min()
                    snmax=np.array([snmax,snbin]).max()
                    temin=np.array([temin,tebin]).min()
                    temax=np.array([temax,tebin]).max()
                    mhmin=np.array([mhmin,mhbin]).min()
                    mhmax=np.array([mhmax,mhbin]).max()
                    if out is not None :
                        iplt=np.where(snbins == snbin)[0][0]
                        plots.plotc(ax[0,iplt],np.array([mhbin+dmh/2.]),np.array([tebin+dte/2.]),np.array([err]),
                                    xr=[mhbins[0],mhbins[-1]],yr=[tebins[0],tebins[-1]],zr=zr,size=30,linewidth=1)

    # do the fit in log(rms) so that empirical uncertainty is positive-definite
    rmsdata=np.log(np.array(rmsdata))
    rmsderiv=np.array(rmsderiv)
    if len(rmsdata) > 5 :
        soln,inv = fit.linear(rmsdata,rmsderiv.transpose())

    figs=[]
    if out is not None :
        y, x = np.mgrid[tebins[0]:tebins[-1]:200j,mhbins[0]:mhbins[-1]:200j]
        for iplt in range(len(snbins)) :
            sn = snbins[iplt]+dsn/2.
            try :
                ax[0,iplt].imshow(elemerr(soln,y-4500.,sn-100.,x, quad=quad),extent=[mhbins[0],mhbins[-1],tebins[0],tebins[-1]], 
                                  aspect='auto',vmin=zr[0],vmax=zr[1], origin='lower',cmap='rainbow')
                ax[1,iplt].imshow(npts[:,:,iplt],extent=[mhbins[0],mhbins[-1],tebins[0],tebins[-1]], 
                                  aspect='auto',vmin=0,vmax=50, origin='lower',cmap='rainbow')
            except: pass

            ax[0,iplt].text(0.98,0.98,title+' S/N={:4.0f}'.format(sn),va='top',ha='right',transform=ax[0,iplt].transAxes)

        fig.savefig(out+'_err.png')
        plt.close()
        figs.append([os.path.basename(out+'_err.png')])

        if snplot :
            fig,ax=plots.multi(len(tebins),len(mhbins),wspace=0.001,hspace=0.001,figsize=(2*len(tebins),2*len(mhbins)))
            x=np.arange(0,250)
            for ix in range(len(tebins)) :
              if ix == 0 : yt=r'$\sigma$'
              else : yt=''
              for iy in range(len(mhbins)) :
                try :
                    gdplt=np.where((np.isclose(rmsderiv[:,1]+4500,tebins[ix]+dte/2.)) & (np.isclose(rmsderiv[:,3],mhbins[iy]+dmh/2.)))[0]
                    plots.plotc(ax[iy,ix],rmsderiv[gdplt,2]+100,np.exp(rmsdata[gdplt]),rmsderiv[gdplt,3],size=30,zr=[-2,0.5],
                                yr=zr,xr=[snbins[0],snbins[-1]],xt='S/N',yt=yt)
                    ax[iy,ix].plot(x,elemerr(soln,tebins[ix]+dte/2.-4500,x-100,mhbins[iy]+dmh/2., quad=quad))
                except: pass
                ax[iy,ix].text(0.98,0.98,'{:8.0f} {:8.2f}'.format(tebins[ix]+dte/2.,mhbins[iy]+dmh/2.),ha='right',va='top',transform=ax[iy,ix].transAxes)
            fig.savefig(out+'_err_sn.png')
            plt.close()
            figs.append([os.path.basename(out+'_err_sn.png')])

        if mkhtml : html.htmltab(figs,file=out+'_err.html',header=title+' empirical uncertainties')

    if verbose : 
        print(soln)
        print(snmin,snmax,temin,temax,mhmin,mhmax)
        pdb.set_trace()

    try : return soln
    except : return 0.

def repeat(data=None,params=None,elems=None,inds=None,stars=None, out='./',elem=True,logg=[-1,6], teff=[3000,8000], log=True, fact=1.0) :
    """ Fits for empirical uncertainty given repeat observations of objects
    """

    # duplicates
    # get indices in file for each degree of RA, to shorten search
    gd=apselect.select(data,badval='STAR_BAD',raw=True,logg=logg,teff=teff)
    a=data[gd]
    ind = np.zeros(360,dtype=int)
    for i in range(1,360) :
        ind[i] = np.where(a['RA']>i)[0][0]
        print(i,ind[i])

    # bins for plots only
    snbins=np.arange(50,300,50)
    tebins=np.arange(3500,7500,250)
    tebins=np.arange(3500,6000,250)
    mhbins=np.arange(-2.25,1.00,0.5)
    dte = tebins[1]-tebins[0]
    dmh = mhbins[1]-mhbins[0]
    dsn = snbins[1]-snbins[0]

    # output data files
    if params is None : params=aspcap.params()[1]
    if elems is None : elems=aspcap.elems()[0]
    of=[]
    for el in elems: 
        of.append(open(out+el+'.txt','w'))

    # loop over stars looking for duplications
    rmsderiv=[]
    rmsparam=[]
    rmselem=[]
    quad= True
    fmt='{:<20s}'+6*' {:9.2f}'+3*' {:10.3f}'+'\n'
    for istar,star in enumerate(a) :
        ira=int(star['RA'])
        i1 = ind[ira]
        if ira < 359 : i2 = ind[ira+1 ]
        else : i2 = len(a)
        jj = i1 + np.where(a['APOGEE_ID'][i1:i2] == star['APOGEE_ID'])[0]
        n = len(jj)
        if n > 1 :
            print(istar,n)
            isort = np.argsort(a['SNR'][jj])
            j=jj[isort]
            i=0
            while i < n-1 :
                if a['SNR'][j[i+1]]/a['SNR'][j[i]] < 1.2 :
                   rmsderiv.append([1.,
                                   a['FPARAM'][j[i:i+2],0].mean()-4500.,
                                   np.min([250.,a['SNR'][j[i:i+2]].mean()])-100.,  # cap S/N at 250
                                   a['FPARAM'][j[i:i+2],3].mean(),
                                   (a['FPARAM'][j[i:i+2],0].mean()-4500.)**2 ])
                   rmsparam.append(np.abs(a['FPARAM'][j[i],:]-a['FPARAM'][j[i+1],:])*np.sqrt(np.pi)/2.)
                   try: rmselem.append(np.abs(a['FELEM'][j[i],0,:]-a['FELEM'][j[i+1],0,:])*np.sqrt(np.pi)/2.)
                   except: rmselem.append(np.abs(a['FELEM'][j[i],:]-a['FELEM'][j[i+1],:])*np.sqrt(np.pi)/2.)
                   for iel in range(len(elems)) :
                       try : 
                           felem=a['FELEM'][j[i],0,iel]
                           felem1=a['FELEM'][j[i+1],0,iel]
                       except : 
                           felem=a['FELEM'][j[i],iel]
                           felem1=a['FELEM'][j[i+1],iel]
                       of[iel].write(fmt.format(star['APOGEE_ID'],a['FPARAM'][j[i],0],a['FPARAM'][j[i+1],0],
                                                a['SNR'][j[i]],a['SNR'][j[i+1]],
                                                a['FPARAM'][j[i],3],a['FPARAM'][j[i+1],3],
                                                felem,felem1,
                                                np.abs(felem-felem1)))
                   i+=2
                else :
                   i+=1
    for iel in range(len(elems)) : of[iel].close()
    rmsderiv=np.array(rmsderiv)
    rmsparam=np.array(rmsparam)
    rmselem=np.array(rmselem)

    # default grids for showing 2D fits
    y, x = np.mgrid[tebins[0]:tebins[-1]:200j,mhbins[0]:mhbins[-1]:200j]

    # parameters
    grid=[]
    ytit=[]
    outtype=np.dtype([('PARAM',params.dtype),('ERRFIT','5f4')])
    outparam=np.empty(len(params),dtype=outtype)

    for i,param in enumerate(params) :
        print(param)
        outparam[i]['PARAM']=param
        if param == 'TEFF' : 
            gd=np.where((rmsparam[:,i] < 500.) & (rmsparam[:,i] > 0.) )[0]
            zr=[0,100] 
        else : 
            gd=np.where((rmsparam[:,i] < 1.) & (rmsparam[:,i] > 0.) )[0]
            zr=[0,0.2]
        if len(gd)<5 : continue

        if log : soln,inv = fit.linear(np.log(rmsparam[gd,i]),rmsderiv[gd,:].transpose())
        else : soln,inv = fit.linear(rmsparam[gd,i],rmsderiv[gd,:].transpose())
        outparam[i]['ERRFIT']=soln

        # plots of 2D fits in bins of S/N
        fig,ax=plots.multi(len(snbins),1,wspace=0.001,figsize=(3*len(snbins),4))
        for iplt in range(len(snbins)) :
            sn = snbins[iplt]+dsn/2.
            ax[iplt].imshow(elemerr(soln,y-4500.,sn-100.,x, quad=quad, log=log, fact=fact),extent=[mhbins[0],mhbins[-1],tebins[0],tebins[-1]], 
                              aspect='auto',vmin=zr[0],vmax=zr[1], origin='lower',cmap='rainbow')
            cs=ax[iplt].contour(np.linspace(mhbins[0],mhbins[-1],200),np.linspace(tebins[0],tebins[-1],200),
                             elemerr(soln,y-4500.,sn-100.,x, quad=quad, log=log, fact=fact),extent=[mhbins[0],mhbins[-1],tebins[0],tebins[-1]], 
                             levels=[0.05,0.1,0.15,0.2],colors='k')
            ax[iplt].clabel(cs,fontsize=8)
            ax[iplt].set_xlabel('[M/H]')
            ax[iplt].set_ylabel('Teff')

            ax[iplt].text(0.98,0.98,param+' S/N={:4.0f}'.format(sn),va='top',ha='right',transform=ax[iplt].transAxes)

        # plots of rms
        snfig,snax=plots.multi(len(tebins)-1,len(mhbins)-1,wspace=0.001,hspace=0.001,figsize=(3*len(tebins),2*len(mhbins)),xtickrot=60)
        xx=np.arange(0,250)
        for ix in range(len(tebins)-1) :
            if ix == 0 : yt=r'$\sigma$'
            else : yt=''
            for iy in range(len(mhbins)-1) :
                gdplt = np.where((rmsderiv[:,1]+4500.>tebins[ix]) & (rmsderiv[:,1]+4500.<tebins[ix+1]) &
                                  (rmsderiv[:,3]>mhbins[iy]) & (rmsderiv[:,3]<mhbins[iy+1]) )[0]
                if len(gdplt) > 1 :
                    plots.plotc(snax[iy,ix],rmsderiv[gdplt,2]+100,rmsparam[gdplt,i],rmsderiv[gdplt,3],size=10,zr=[-2,0.5],
                                yr=zr,xr=[snbins[0],snbins[-1]+10],xt='S/N',yt=yt)
                snax[iy,ix].set_ylim(zr)
                snax[iy,ix].plot(xx,elemerr(soln,tebins[ix]+dte/2.-4500,xx-100,mhbins[iy]+dmh/2., quad=quad, log=log, fact=fact))
                snax[iy,ix].text(0.98,0.98,'{:8.0f}{:6.2f}'.format(tebins[ix]+dte/2.,mhbins[iy]+dmh/2.),ha='right',va='top',transform=snax[iy,ix].transAxes,fontsize='x-small')
                for iz in range(len(snbins)-1) :
                    gd= np.where((rmsderiv[gdplt,2]+100.>snbins[iz]) & (rmsderiv[gdplt,2]+100.<snbins[iz+1]) ) [0]
                    if len(gd) > 0 :
                        snax[iy,ix].text(0.98,0.88-iz*0.08,'{:8.2f}'.format(
                                    rmsparam[gdplt[gd],i].mean()),transform=snax[iy,ix].transAxes,fontsize='x-small',ha='right',va='top')
        fig.savefig(out+param+'.png')
        plt.close(fig)
        snfig.savefig(out+param+'_sn.png')
        plt.close(snfig)
        grid.append([os.path.basename(out+param+'.png'),os.path.basename(out+param+'_sn.png')])
        ytit.append(param)
    html.htmltab(grid,file=out+'repeat_param.html',ytitle=ytit)
  
    # elements 
    grid=[]
    ytit=[]
    outtype=np.dtype([('ELEM',elems.dtype),('ERRFIT','5f4')])
    outelem=np.empty(len(elems),dtype=outtype)
    allfig,allax=plots.multi(5,5,hspace=0.001,wspace=0.001,xtickrot=60)
    for i,el in enumerate(elems) :
        print(el)
        outelem[i]['ELEM']=el
        gd=np.where((rmselem[:,i] < 1.) & (rmselem[:,i] > 0.) )[0]
        zr=[0,0.19]
        if len(gd)<5 : continue

        if log :soln,inv = fit.linear(np.log(rmselem[gd,i]),rmsderiv[gd,:].transpose())
        else :soln,inv = fit.linear(rmselem[gd,i],rmsderiv[gd,:].transpose())
        outelem[i]['ERRFIT']=soln
        fig,ax=plots.multi(len(snbins),1,wspace=0.001,figsize=(3*len(snbins),4))
        for iplt in range(len(snbins)) :
            sn = snbins[iplt]+dsn/2.
            cs=ax[iplt].contour(np.linspace(mhbins[0],mhbins[-1],200),np.linspace(tebins[0],tebins[-1],200),
                             elemerr(soln,y-4500.,sn-100.,x, quad=quad, log=log, fact=fact),
                             levels=[0.05,0.1,0.15,0.2],colors='k') 
            ax[iplt].clabel(cs,fontsize=8)
            ax[iplt].imshow(elemerr(soln,y-4500.,sn-100.,x, quad=quad, log=log, fact=fact),
                              extent=[mhbins[0],mhbins[-1],tebins[0],tebins[-1]], 
                              aspect='auto',vmin=zr[0],vmax=zr[1], origin='lower',cmap='rainbow')
            ax[iplt].text(0.98,0.98,el+' S/N={:4.0f}'.format(sn),va='top',ha='right',transform=ax[iplt].transAxes)
            if i<25 and sn == 125 :
                allax[i//5,i%5].text(0.05,0.95,el,va='top',ha='left',transform=allax[i//5,i%5].transAxes,color='w')
                allax[i//5,i%5].set_xlabel('[M/H]')
                if i%5 == 0 : allax[i//5,i%5].set_ylabel(r'$T_{eff}$')
                cs=allax[i//5,i%5].contour(np.linspace(mhbins[0],mhbins[-1],200),np.linspace(tebins[0],tebins[-1],200),
                                        elemerr(soln,y-4500.,sn-100.,x, quad=quad, log=log, fact=fact),
                                        levels=[0.05,0.1,0.15,0.2],colors='k') 
                allax[i//5,i%5].clabel(cs,fontsize=6)
                cm= allax[i//5,i%5].imshow(elemerr(soln,y-4500.,sn-100.,x, quad=quad, log=log, fact=fact),
                                       extent=[mhbins[0],mhbins[-1],tebins[0],tebins[-1]], 
                                       aspect='auto',vmin=zr[0],vmax=0.1, origin='lower',cmap='rainbow')

        snfig,snax=plots.multi(len(tebins)-1,len(mhbins)-1,wspace=0.001,hspace=0.001,figsize=(2.2*len(tebins),3*len(mhbins)),xtickrot=60)
        fig2,ax2=plots.multi(len(tebins)-1,len(mhbins)-1,wspace=0.001,hspace=0.001,figsize=(2*len(tebins),3*len(mhbins)),xtickrot=60)
        xx=np.arange(0,250)
        for ix in range(len(tebins)-1) :
            for iy in range(len(mhbins)-1) :
                if ix == 0 : snax[iy,ix].set_ylabel(r'$\sigma$('+el+')')
                gdplt = np.where((rmsderiv[:,1]+4500.>tebins[ix]) & (rmsderiv[:,1]+4500.<tebins[ix+1]) &
                                  (rmsderiv[:,3]>mhbins[iy]) & (rmsderiv[:,3]<mhbins[iy+1]) )[0]
                if len(gdplt) > 1 :
                    #plots.plotc(snax[iy,ix],rmsderiv[gdplt,2]+100,rmselem[gdplt,i],rmsderiv[gdplt,3],size=5,zr=[-2,0.5],
                    #            yr=zr,xr=[snbins[0],snbins[-1]],xt='S/N')
                    plots.plotp(snax[iy,ix],rmsderiv[gdplt,2]+100,rmselem[gdplt,i],color='k',size=10,zr=[-2,0.5],
                                yr=zr,xr=[snbins[0],snbins[-1]+10],xt='S/N')
                snax[iy,ix].set_xlim(snbins[0]+1,snbins[-1]-1)
                snax[iy,ix].set_xlabel('S/N')
                snax[iy,ix].set_ylim(zr)
                snax[iy,ix].plot(xx,elemerr(soln,tebins[ix]+dte/2.-4500,xx-100,mhbins[iy]+dmh/2., quad=quad, log=log, fact=fact))
                snax[iy,ix].text(0.02,0.95,r'$T_{eff}:$'+'{:6.0f}'.format(tebins[ix]+dte/2.),
                                 ha='left',va='top',transform=snax[iy,ix].transAxes,fontsize='large')
                snax[iy,ix].text(0.02,0.80,'[M/H]: {:6.2f}'.format(mhbins[iy]+dmh/2.),
                                 ha='left',va='top',transform=snax[iy,ix].transAxes,fontsize='large')
                for iz in range(len(snbins)-1) :
                    gd= np.where((rmsderiv[gdplt,2]+100.>snbins[iz]) & (rmsderiv[gdplt,2]+100.<snbins[iz+1]) ) [0]
                    #snax[iy,ix].text(0.98,0.88-iz*0.08,'{:8.3f}'.format(rmselem[gdplt[gd],i].mean()),
                    #                 transform=snax[iy,ix].transAxes,fontsize='x-small',ha='right',va='top')
                    ax2[iy,ix].text(0.98,0.98,'{:8.0f}{:6.2f}'.format(tebins[ix]+dte/2.,mhbins[iy]+dmh/2.),
                                    ha='right',va='top',transform=ax2[iy,ix].transAxes,fontsize='x-small')
                    ax2[iy,ix].set_xlim(0,0.16)
                    bins=np.arange(0,0.161,0.02)
                    ax2[iy,ix].set_ylim(0,10)
                    if len(gd) > 0 :
                        if iz == 1 and len(gd) > 3 :
                            ax2[iy,ix].hist(np.clip(rmselem[gdplt[gd],i],bins[0],bins[-1]-0.001),bins=bins,histtype='step',color='k')
                            ylim=ax2[iy,ix].get_ylim()
                            m=rmselem[gdplt[gd],i].mean()
                            ax2[iy,ix].plot([m,m],[ylim[0],ylim[1]],color='r')
                            m=np.median(rmselem[gdplt[gd],i])
                            ax2[iy,ix].plot([m,m],[ylim[0],ylim[1]],color='g')
                            m=elemerr(soln,tebins[ix]+dte/2.-4500,snbins[iz]+dsn/2.-100,mhbins[iy]+dmh/2.,
                                      quad=quad,log=log,fact=fact)
                            ax2[iy,ix].plot([m,m],[ylim[0],ylim[1]],color='b')
        fig.savefig(out+el+'.png')
        plt.close(fig)
        snfig.savefig(out+el+'_sn.png')
        plt.close(snfig)
        fig2.savefig(out+el+'_hist2.png')
        plt.close(fig2)
        fig,ax=plots.multi(1,1)
        ax.hist(rmselem[:,i]-elemerr(soln,rmsderiv[:,1],rmsderiv[:,2],rmsderiv[:,3],quad=quad,log=log,fact=fact),bins=np.arange(-0.4,0.4,0.005))
        ylim=ax.get_ylim()
        ax.plot([0.,0.],[ylim[0],ylim[1]])
        ax.set_xlabel('diff*sqrt(pi)/2.-fit')
        fig.savefig(out+el+'_hist.png')
        plt.close(fig)
        grid.append([os.path.basename(out+el+'.png'),os.path.basename(out+el+'_sn.png'),os.path.basename(out+el+'_hist.png'),
                     os.path.basename(out+el+'_hist2.png')])
        ytit.append('<A HREF='+os.path.basename(out)+el+'.txt>'+el+'</A>')

    cb=allfig.colorbar(cm, ax=list(allax[:,-1]))
    cb.ax.set_ylabel('$\sigma$')

    allfig.savefig(out+'allel.png')
    html.htmltab(grid,file=out+'repeat_elem.html',ytitle=ytit)
    hdulist=fits.HDUList()
    hdulist.append(fits.BinTableHDU(outparam))
    hdulist.append(fits.BinTableHDU(outelem))
    hdulist.append(fits.ImageHDU(rmsderiv))
    hdulist.append(fits.ImageHDU(rmselem))
    hdulist.writeto(out+'errfit.fits',overwrite=True)

    return outparam,outelem

def gauss() :
    """ simple simulator to assess statistics of drawing repeats from a Gaussian distribution
    """
    # draw pairs from a unit-Gaussian and accumulate differences
    diff=[]
    for i in range(1000000) :
        a=normal(size=2)
        diff.append(a[0]-a[1])

    fig,ax=plots.multi(1,1,figsize=(10,5))

    diff=np.array(diff)
    ax.hist(diff,bins=51)
    ax.set_xlabel('Difference between pair')
    ax.text(0.02,0.95,'<diff>: {:7.3f}'.format(diff.mean()), transform=ax.transAxes,size='small')
    # Vijith says std of differences is sqrt(2), which is verified
    ax.text(0.02,0.90,'std: {:7.3f}'.format(diff.std()), transform=ax.transAxes,size='small')
    # biased estimator of standard deviation
    ax.text(0.02,0.85,'<abs(diff/sqrt(2))>: {:7.3f}'.format(abs(diff/np.sqrt(2.)).mean()), transform=ax.transAxes,size='small')
    # is unbiased estimator of variance!
    ax.text(0.02,0.80,'<diff**2/2>: {:7.3f}'.format(((diff**2/2)).mean()), transform=ax.transAxes,size='small')
    # corrected biased estimator by sqrt(pi/2.)
    ax.text(0.02,0.75,'<abs(diff/sqrt(2)*sqrt(pi/2)>: {:7.3f}'.format(abs(diff/np.sqrt(2)*np.sqrt(np.pi/2.)).mean()), transform=ax.transAxes,size='small')
    # but what if we take the log, compare to exp(<log>)
    ax.text(0.02,0.70,'exp(<log(abs(diff/sqrt(2)*sqrt(pi/2))>): {:7.3f}'.format(np.exp(np.log(abs(diff/np.sqrt(2)*np.sqrt(np.pi/2.))).mean())), transform=ax.transAxes,size='small')
    ax.set_xlim(-10,8)
    print(diff.mean(),diff.std(),abs(diff).mean(),(diff**2).mean(),(diff/np.sqrt(2.)).std(),abs(diff*np.sqrt(np.pi/2.)).mean())
    print(abs(diff*np.sqrt(np.pi/2.)).mean(),np.exp(np.log(abs(diff*np.sqrt(np.pi/2.))).mean()))

    """
    two observations x_1 and x_2, separated by dx
    x2=x1+dx
    x1-xmean= dx/2
    x2-xmean=dx/2

    for 2 points (n=2) :
    std = sqrt ( (x1 - x_mean)**2 + (x2-x_mean)**2 )
        = sqrt ( (dx/2)**2 + (dx/2)**2 )
        = sqrt ( dx**2 / 2. )
        = dx / sqrt(2)
 
    but this is a biased estimator
    """

