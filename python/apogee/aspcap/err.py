import os
import numpy as np
import matplotlib.pyplot as plt

from tools import plots
from tools import html
from tools import fit

def errfit(te, snr, mh, val, snbins=np.arange(50,250,50), tebins=np.arange(3500,6000,250), mhbins=np.arange(-2.25,0.75,0.5),verbose=False,out=None,title='', zr=[0,0.1], snplot=True, meanerr=None ) :
    '''
    Fits for empirical uncertainty as function of Teff, S/N, and [M/H]
    '''
    if out is not None :
        fig,ax=plots.multi(len(snbins),1,wspace=0.001,figsize=(3*len(snbins),2))

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
    for mhbin in mhbins :
        for tebin in tebins :
            for snbin in snbins :
                ibin = np.where(( te > tebin) & (te <= tebin+dte) &
                                ( mh > mhbin ) & (mh <= mhbin+dmh) &
                                ( snr > snbin) & (snr <= snbin+dsn) & (val > -9990.) )[0]
                if len(ibin) > 3 :
                    if meanerr is not None :
                        err = np.sqrt(np.clip(val[ibin].std()**2 - np.median(meanerr[ibin])**2,0.001,10000000.))
                    else :
                        err = val[ibin].std()
                    rmsdata.append(err)
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
                        plots.plotc(ax[iplt],mhbin+dmh/2.,tebin+dte/2.,err,xr=[mhbins[0],mhbins[-1]],yr=[tebins[0],tebins[-1]],zr=zr,size=30,linewidth=1)

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
                ax[iplt].imshow(elemerr(soln,y-4500.,sn-100.,x),extent=[mhbins[0],mhbins[-1],tebins[0],tebins[-1]], 
                                aspect='auto',vmin=zr[0],vmax=zr[1], origin='lower',cmap='rainbow')
            except: pass

            ax[iplt].text(0.98,0.98,title+' S/N={:4.0f}'.format(sn),va='top',ha='right',transform=ax[iplt].transAxes)

        fig.savefig(out+'_err.jpg')
        plt.close()
        figs.append([os.path.basename(out+'_err.jpg')])

        if snplot :
            fig,ax=plots.multi(len(tebins),len(mhbins),wspace=0.001,hspace=0.001,figsize=(2*len(tebins),2*len(mhbins)))
            for ix in range(len(tebins)) :
              if ix == 0 : yt=r'$\sigma$'
              else : yt=''
              for iy in range(len(mhbins)) :
                try :
                    gdplt=np.where((np.isclose(rmsderiv[:,1]+4500,tebins[ix]+dte/2.)) & (np.isclose(rmsderiv[:,3],mhbins[iy]+dmh/2.)))[0]
                    plots.plotc(ax[iy,ix],rmsderiv[gdplt,2]+100,np.exp(rmsdata[gdplt]),rmsderiv[gdplt,3],size=30,zr=[-2,0.5],
                                yr=zr,xr=[snbins[0],snbins[-1]],xt='S/N',yt=yt)
                except: pass
                ax[iy,ix].text(0.98,0.98,'{:8.0f} {:8.2f}'.format(tebins[ix]+dte/2.,mhbins[iy]+dmh/2.),ha='right',va='top',transform=ax[iy,ix].transAxes)
            fig.savefig(out+'_err_sn.jpg')
            plt.close()
            figs.append([os.path.basename(out+'_err_sn.jpg')])

        html.htmltab(figs,file=out+'_err.html',header=title+' empirical uncertainties')

    if verbose : 
        print(soln)
        print(snmin,snmax,temin,temax,mhmin,mhmax)
        pdb.set_trace()

    try : return soln
    except : return 0.


def elemerr(soln,te,sn,fe) :
    ''' 
    Function to evaluate fit for uncertainty
    '''
    out=soln[0]+soln[1]*te+soln[2]*sn
    if len(soln) > 3: out+= soln[3]*fe
    return np.exp(out)

