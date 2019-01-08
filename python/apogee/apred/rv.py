import os
import pdb
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from apogee.utils import apload
from tools import plots
from tools import html
from tools import match
from tools import struct

colors=['r','g','b','c','m','y','k']


def allField(files=['apo*/*/apField-*.fits','apo*/*/apFieldC-*.fits','lco*/*/apField-*.fits'],out='allField.fits',verbose=False) :
    '''
    Concatenate set of apField files
    '''
    # concatenate the structures
    all=struct.concat(files,verbose=verbose)

    # write out the file
    if out is not None:
        print('writing',out)
        struct.wrfits(all,out)

    return all

def allFieldVisits(files=['apo*/*/apFieldVisits-*.fits','apo*/*/apFieldC-*.fits','lco*/*/apFieldVisits-*.fits'],out='allFieldVisits.fits',verbose=False) :
    '''
    Concatenate set of apField files
    '''
    # concatenate the structures
    all=struct.concat(files,verbose=verbose)

    # write out the file
    if out is not None:
        print('writing',out)
        struct.wrfits(all,out)

    return all


def vscat(a,fig=None,ls=None,marker='o') :
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
    gd=np.where((a['TELESCOPE'] == 'apo25m') & (a['RV_FEH']>minfeh) )[0]
    fig=vscat(a[gd],marker='o')
    gd=np.where((a['TELESCOPE'] == 'lco25m') & (a['RV_FEH']>minfeh) )[0]
    vscat(a[gd],fig=fig,ls=':',marker='+')
    if out is not None : 
        fig[0].savefig(out+'.png')
        plt.close()

def comp(a,b,domatch=True,out=None) :

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

def all() :
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
