import numpy as np
from tools import plots
import os
import pdb
from brokenaxes import brokenaxes
from matplotlib.pyplot import GridSpec
from apogee.aspcap import aspcap
from apogee.utils import apload
import matplotlib.pyplot as plt
from tools import html

def chip(a,row=150,ax=None,pixel=False,visit=False,color=None) :
    """ Routine to plot 3 chips in 3 panels
    """
    if ax is None : fig,ax=plots.multi(1,3,hspace=0.3)
    if visit : rows=range(3)
    else : rows=[row,row,row]
    chips=['a','b','c']
    x=np.arange(a['a'][1].header['NAXIS1'])
    for ichip,chip in enumerate(chips) :
        if pixel :
            plots.plotl(ax[ichip],x,a[chip][1].data[rows[ichip],:],color=color)
        else :
            plots.plotl(ax[ichip],a[chip][4].data[rows[ichip],:],a[chip][1].data[rows[ichip],:],color=color)

    try : return fig,ax
    except : return

def elem(data,inds,els=['all'],out=None,apred_vers='dr17',aspcap_vers='l33') :
    """ Plots spectra in windows for requested elements of requested stars
            data: allStar table
            inds: list of indices of stars to plot
            els: list of elements to create plots for
            out: prefix for hardcopy if desired
            apred_vers, aspcap_vers : version identifier for aspcapStar files
    """

    # lists of elements and indices
    if not isinstance(els,list) : els=[els]
    nel=len(els)
    if not isinstance(inds,list) : inds=[inds]
    nspec=len(inds)

    # setup reader
    load=apload.ApLoad(apred=apred_vers,aspcap=aspcap_vers)

    # get wavelengths
    wgrid=np.hstack(aspcap.gridWave())
    wave=aspcap.apStarWave()

    # loop over elements
    elems=aspcap.elems()
    grid=[]
    for iel,el in enumerate(els) :
        if el == 'all' :
            wplot=[(15150,15800),(15865,16425),(16480,16950)]
            filt=np.loadtxt(os.environ['APOGEE_DIR']+'/data/windows/dr17/global_mask_v02.txt')
            index=-1
        else :
            index=np.where(elems[0] == el)[0][0]
            # read windows and filt file, set up wavelength ranges
            wind=np.loadtxt(os.environ['APOGEE_DIR']+'/data/windows/dr17/'+el+'.wind')
            wplot=[]
            for w in wind :
                if w[2]>0 : wplot.append((w[0],w[1]))
            filt=np.loadtxt(os.environ['APOGEE_DIR']+'/data/windows/dr17/'+el+'.filt')

        # do the plots
        nw=len(wplot)
        if el=='all' :
            fig,ax=plots.multi(nw,nspec,wspace=0.1,figsize=(16*nw,4*nspec),hspace=0.2,brokenx=True)
            ax=ax.reshape([nspec,nw])
        else :
            fig,ax=plots.multi(nw,nspec,wspace=0.1,figsize=(3*nw,4*nspec),hspace=0.2,brokenx=True)
            fig.suptitle(el) 
            ax=np.atleast_2d(ax)
        for iy,ind in enumerate(inds) :
            load.telescope=data['TELESCOPE'][ind]
            spec=load.aspcapStar(data['FIELD'][ind],data['APOGEE_ID'][ind])
            for i,w in enumerate(wplot) :
                plots.plotl(ax[iy,i],wave,spec[1].data,xr=w,yr=[-0.1,1.2],linewidth=2)
                ax[iy,i].plot(wave,spec[2].data,linewidth=1)
                ax[iy,i].plot(wave,spec[3].data,linewidth=1)
                #plots.plotl(ax[iy,i],wave,spec['SPEC_BESTFIT'],xr=w,yr=[0,1.2],linewidth=2)
                ax[iy,i].plot(wgrid,filt*0.5,linewidth=1)
                ax[iy,i].ticklabel_format(style='plain',useOffset=False)
                if index >= 0:
                    if elems[1][index] == 1 : xm = data['FELEM'][ind,index]
                    else : xm = data['FELEM'][ind,index] +data['FPARAM'][ind,3]
                else : xm=-99
                lab='{:s} Teff: {:6.0f} log g: {:6.2f} [M/H] :{:6.2f} [{:s}/H] :{:6.2f} CHI2:{:6.2f}'.format(
                     data['APOGEE_ID'][ind],data['FPARAM'][ind,0],data['FPARAM'][ind,1],data['FPARAM'][ind,3],el,xm,data['ASPCAP_CHI2'][ind]) 
                if i == 0 : ax[iy,i].text(-0.05,1.02,lab,transform=ax[iy,i].transAxes)
                if iy < nspec-1 : ax[iy,i].tick_params(labelbottom=False) 
                if i == nw-1 : ax[iy,i].text(0.9,0.5,el,transform=ax[iy,i].transAxes,fontsize='large')

        #hardcopy
        if out is not None :
            fig.savefig(out+'_'+el+'.png')
            plt.close()
            grid.append(os.path.basename(out+'_'+el+'.png'))

    # html
    if out is not None : html.htmltab([grid],file=out+'.html')
