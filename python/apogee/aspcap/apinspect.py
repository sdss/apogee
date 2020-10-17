from apogee.aspcap import aspcap
from apogee.aspcap import qa
from apogee.utils import apload
from tools import plots
import glob
import numpy as np
import os
import pdb
import matplotlib.pyplot as plt

def tostr(dat) :
    if type(dat) is str : return dat
    else : return dat.decode()

def getaxis(a,axis,raw=True) :
    if axis == 'Teff' :
        if raw : x=a['FPARAM'][:,0]
        else : x=a['PARAM'][:,0]
        xlim=[8000,3000]
    elif axis == 'log g' :
        if raw : x=a['FPARAM'][:,1]
        else : x=a['PARAM'][:,1]
        xlim=[6,-1]
    elif axis == 'M_H' :
        if raw : x=a['FPARAM'][:,3]
        else : x=a['PARAM'][:,3]
        xlim=[-2.5,0.5]
    elif axis == 'C_M' :
        if raw : x=a['FPARAM'][:,4]
        else : x=a['PARAM'][:,4]
        xlim=[-0.5,0.5]
    elif axis == 'N_M' :
        if raw : x=a['FPARAM'][:,5]
        else : x=a['PARAM'][:,5]
        xlim=[-0.25,1.5]
    elif axis == 'A_M' :
        if raw : x=a['FPARAM'][:,6]
        else : x=a['PARAM'][:,6]
        xlim=[-0.5,0.5]
    else :
        j=np.where(aspcap.elems()[0] == axis)[0]
        if len(j) > 0 :
            if raw : 
                x=a['FELEM'][:,j]
                xlim=[np.median(x)-0.5,np.median(x)+0.5]
            else : 
                x=a['X_M'][:,j]
                xlim=[-0.5,0.5]
        else :
            print('Unknown axis!')
    return x,xlim

def inspect(a=None,xaxis='Teff',yaxis='log g',zaxis='M_H',raw=True,param='FPARAM',indir='cal',apred='r14',aspcap_vers='l33',verbose=False,
            xr=None, yr=None, zr=None) :
    """ Given input structure, plot HR diagram, and enter event loop to mark stars to plot spectra
    """

    x,xlim = getaxis(a,xaxis,raw=raw)
    y,ylim = getaxis(a,yaxis,raw=raw)
    z,zlim = getaxis(a,zaxis,raw=raw)
    if xr == None : xr=xlim
    if yr == None : yr=ylim
    if zr == None : zr=zlim

    load=apload.ApLoad(apred=apred,aspcap=aspcap_vers,verbose=verbose)
    if a is None : a=load.allCal()[1].data

    fig,ax = plots.multi(1,1)
    plots.plotc(ax,x,y,z,xr=xr,yr=yr,zr=zr)
    plots.event(fig)
    sf,sa=plots.multi(1,1)
    nplot = 11
    hf,ha=plots.multi(1,nplot,figsize=(8.5,11),hspace=0.2)
    ha2=[]
    for i in range(nplot) : ha2.append(ha[i].twinx())
    print('hit any key near object to plot in HR diagram, q to quit')
    while (1) :
        ret=plots.mark(fig)
        if ret[2] == 'q' : break
        ind=plots._index[0]
        load.settelescope('apo25m')
        try : f=load.aspcapField(tostr(a['ALTFIELD'][ind]))
        except : f=load.aspcapField(tostr(a['FIELD'][ind]))
        if f is None :
            load.settelescope('lco25m')
            try : f=load.aspcapField(tostr(a['ALTFIELD'][ind]))
            except : f=load.aspcapField(tostr(a['FIELD'][ind]))
        if f is None :
            f=glob.glob(indir+'/*'+tostr(a['APOGEE_ID'][plots._index[0]])+'*')
            dir=os.path.dirname(f[0])
            f=glob.glob(dir+'/*aspcapField*.fits')
            f=fits.open(f[0])
        f[3].data['WAVE']=np.hstack(aspcap.gridWave())
        print('f: ',f)
        data=f[1].data
        j=np.where(data['APOGEE_ID'] == tostr(a['APOGEE_ID'][plots._index[0]]))[0][0]
        sa.cla()
        for i in range(11) :
            ha[i].cla()
            ha2[i].cla()
            ha[i].set_ylabel('Flux')
            ha2[i].set_ylabel(r'$\chi^2$')
            ha[i].set_ylim(0.5,1.3)
            ha2[i].set_ylim(0.,20.)
        qa.plot(f[3].data['WAVE'][0],f[2].data['SPEC'][j,:],ax=ha,sum=True,color='k')
        qa.plot(f[3].data['WAVE'][0],f[2].data['SPEC_BESTFIT'][j,:],ax=ha,sum=True,color='b')
        chi2 = (f[2].data['SPEC'][j,:]-f[2].data['SPEC_BESTFIT'][j,:])**2/f[2].data['ERR'][j,:]**2
        qa.plot(f[3].data['WAVE'][0],chi2,ax=ha2,sum=True,alpha=0.4)
   
        plots.plotl(sa,f[3].data['WAVE'][0],f[2].data['SPEC'][j,:])
        plots.plotl(sa,f[3].data['WAVE'][0],f[2].data['SPEC_BESTFIT'][j,:])
        try: chi2 = data['ASPCAP_CHI2']
        except: chi2 = data['PARAM_CHI2']
        text1=r'ID: {:s} FIELD: {:s} SNR: {:6.1f} $\chi^2$: {:6.1f}'.format(
             data['APOGEE_ID'][j],data['FIELD'][j],data['SNR'][j],chi2)
        text2=r'Teff: {:5.0f} logg: {:5.1f} [M/H]: {:5.2f} [$\alpha$/M]: {:5.2f} [C/M]: {:5.2f} [N/M]: {:5.2f}'.format(
             data[param][j,0],data[param][j,1],data[param][j,3],data[param][j,6],data[param][j,4],data[param][j,5])
        sf.suptitle(text1+'\n'+text2)
        hf.suptitle(text1+'\n'+text2)
        plt.draw()
        plt.show()
    plt.close(hf)
    plt.close(sf)
    plt.close(fig)

