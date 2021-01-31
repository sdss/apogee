import numpy as np
from tools import plots, html
from apogee.utils import apload
from apogee.aspcap import aspcap
from astropy.io import fits
import pdb
import yaml
import matplotlib.pyplot as plt

def check(allstar, allvisit ) :
    """ Calculate for each star fraction of visits and flux with MTPFLUX less than several thresholds
    """

    # get indices in allVisit file for each degree of RA, to shorten search
    ind = np.zeros(360,dtype=int)
    for i in range(1,360) :
        ind[i] = np.where(allvisit['RA']>i)[0][0]
        print(i,ind[i])

    # loop over all stars and get corresponding visits
    nvisits=[]
    nmax=0
    frac=[]
    fluxfrac=[]
    threshes=[0.3,0.5,0.7]
    fig,ax=plots.multi(len(threshes),2,hspace=0.001,sharey=True)
    apo25m = np.where(allstar['TELESCOPE'] == 'apo25m')[0]
    for i,star in enumerate(allstar[apo25m]) :
        ira=int(star['RA'])
        i1 = ind[ira]
        if ira < 359 : i2 = ind[ira+1 ]
        else : i2 = len(allvisit)
        j = i1+np.where((allvisit['APOGEE_ID'][i1:i2] == star['APOGEE_ID']) & (allvisit['FIELD'][i1:i2] == star['FIELD']))[0]

        if len(j) > 0 :
            tmp=[]
            fluxtmp=[]
            for thresh in threshes :
                bd=np.where(allvisit['MTPFLUX'][j] < thresh)[0]
                tmp.append( len(bd)/len(j) )
                fluxtmp.append( allvisit['MTPFLUX'][j[bd]].sum()/allvisit['MTPFLUX'][j].sum() )
            frac.append(tmp)
            fluxfrac.append(fluxtmp)
        else :
            print('{:s}: no visits! {:s} '.format(star['APOGEE_ID'],star['FIELD']))

    frac=np.array(frac)
    fluxfrac=np.array(fluxfrac)
    for i,thresh in enumerate(threshes) :
        ax[0,i].hist(frac[:,i],bins=np.arange(0,1,0.05))    
        ax[0,i].hist(frac[:,i],bins=np.arange(0,1,0.05),cumulative=True,histtype='step')    
        ax[0,i].set_xlabel('Frac of bad visits')
        ax[0,i].text(0.4,0.9,'Thresh: {:5.1f}'.format(thresh),transform=ax[0,i].transAxes)

        ax[1,i].hist(fluxfrac[:,i],bins=np.arange(0,1,0.05))    
        ax[1,i].hist(fluxfrac[:,i],bins=np.arange(0,1,0.05),cumulative=True,histtype='step')    
        ax[1,i].set_xlabel('Flux frac of bad visits')
        #ax[0,i].set_yscale('log')
        #ax[1,i].set_yscale('log')
    fig.tight_layout()
    stats(allstar,(frac,fluxfrac))
    return frac, fluxfrac


def stats(allstar,out,fluxfrac=True) :
    """ Statistics of low flux
    """
  
    apo25m = np.where(allstar['TELESCOPE'] == 'apo25m')[0]
    print('exceeding thresh')
    if fluxfrac : iout=1
    else : iout=0
    for frac in [0.5, 0.99] :
        print(frac)
        for i,thresh in enumerate([0.3,0.5,0.7]) :
            n=np.where(out[iout][:,i]>frac)[0]
            print('{:8.1f}{:8d}{:8.2f}'.format(thresh,len(n),len(n)/len(apo25m)*100))
    print('below thresh')
    for frac in [0.99, 0.5, 0.25, 0.1, 0.01] :
        print(frac)
        for i,thresh in enumerate([0.3,0.5,0.7]) :
            n=np.where(out[1][:,i]<frac)[0]
            print('{:8.1f}{:8d}{:8.2f}'.format(thresh,len(n),len(n)/len(apo25m)*100))

def spec(a,v,out) :
    """ Plot spectra of different visits labeleed by MTPFLUX
    """

    # select bright stars to get good S/N per visit
    apo25m = np.where(a['TELESCOPE'] == 'apo25m')[0]
    jj=np.where((out[1][:,1]>0) & (a['H'][apo25m]<10) & 
                (a['RV_TEFF'][apo25m]<4500) & (a['RV_TEFF'][apo25m]>4000) &
                (a['RV_LOGG'][apo25m]<3.5)) [0]

    # see how many stars per field, so we can choose a good field for ASPCAP test
    fields = set(a['FIELD'][apo25m[jj]])
    for field in fields :
        j=np.where(a['FIELD'][apo25m[jj]] == field)[0]
        print(field,len(j))

    # choose 120-08-RV
    jj=np.where((out[1][:,1]>0) & (a['H'][apo25m]<10) & 
                (a['RV_TEFF'][apo25m]<4500) & (a['RV_TEFF'][apo25m]>4000) &
                (a['RV_LOGG'][apo25m]<3.5) & (a['FIELD'][apo25m] == '120-08-RV')) [0]

    # loop over stars and plot spectra
    load=apload.ApLoad(apred='dr17')
    wave=aspcap.apStarWave()

    grid=[]
    print(a['APOGEE_ID'][apo25m[jj]])
    for iy,j in enumerate(apo25m[jj]) :
        jv = np.where(v['APOGEE_ID'] == a['APOGEE_ID'][j])[0]
        spec=load.apStar(a['FIELD'][j],a['APOGEE_ID'][j],load=True)
        n=len(spec.flux)-2
        fig,ax=plots.multi(1,1,hspace=0.001,figsize=(12,6))

        # get the MTPFLUX for individual visit, matching spectra in apStar by JD
        mtpflux=[]
        for i in range(n) :
            jjj=np.where(abs(v['JD'][jv]-spec.header['JD{:d}'.format(i+1)]) < 1)[0][0]
            mtpflux.append(v['MTPFLUX'][jv[jjj]])
        ind=np.argsort(np.array(mtpflux))
        for i in range(n) :
            jjj=np.where(abs(v['JD'][jv]-spec.header['JD{:d}'.format(ind[i]+1)]) < 1)[0]
            plots.plotl(ax,wave,spec.flux[ind[i]+2]/spec.flux[0]+i*0.1)
            ax.text(17000,1+i*0.1,'{:7.2f}'.format(v['MTPFLUX'][jv[jjj[0]]]),va='center')
        ax.set_ylim(0.9, 1.1+n*0.1)
        ax.set_title(a['APOGEE_ID'][j])
        fig.savefig(a['APOGEE_ID'][j]+'.png')
        plt.close()
        grid.append([a['APOGEE_ID'][j]+'.png'])
    html.htmltab(grid,file='spec.html')

def param(field='120-08-RV') :
    """ Plot parameters of individual visits vs combined
    """
    a=fits.open('apo25m/{:s}/aspcapField-{:s}.fits'.format(field,field))
    config=yaml.safe_load(open('test.yml','r'))
    stars = config['obj']
    grid=[]
    for star in stars :
        j=np.where(np.core.defchararray.find(a[1].data['APOGEE_ID'],star)>=0)[0]

        fig,ax=plots.multi(1,8,hspace=0.001,figsize=(8,6))
        for i,par in enumerate(aspcap.params()[2]) :
            if i < 7 :
                ref=a[1].data['FPARAM'][j[0],i]
                plots.plotp(ax[i],a[1].data['SNR'][j[1:]],a[1].data['FPARAM'][j[1:],i]-ref,xt='SNR',yt=par,size=10)
                ax[i].plot(ax[i].get_xlim(),[0,0],ls=':')
            elif i==7 :
                plots.plotp(ax[i],a[1].data['SNR'][j[1:]],a[1].data['ASPCAP_CHI2'][j[1:]],xt='SNR',yt='CHI2',size=10)
        fig.savefig(star+'_param.png')
        plt.close()

        fig,ax=plots.multi(1,26,hspace=0.001,figsize=(8,18))
        for i,par in enumerate(aspcap.elems()[0]) :
            ref=a[1].data['FELEM'][j[0],i]
            plots.plotp(ax[i],a[1].data['SNR'][j[1:]],a[1].data['FELEM'][j[1:],i]-ref,xt='SNR',yt=par,size=10)
            ax[i].plot(ax[i].get_xlim(),[0,0],ls=':')
        fig.savefig(star+'_elem.png')
        plt.close()

        grid.append([star+'.png',star+'_param.png',star+'_elem.png'])

    html.htmltab(grid,file='mtp.html')
