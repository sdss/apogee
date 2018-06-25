from astropy.io import fits
from apogee.utils import bitmask
from apogee.utils import apselect
from tools import plots
import os
import numpy as np
import matplotlib.pyplot as plt

a=fits.open(os.environ['APOGEE_REDUX']+'/r8/stars/l31c/l31c.1/allStar-l31c.1.fits')[1].data
c=fits.open(os.environ['APOGEE_REDUX']+'/r8/stars/l31c/l31c.1/allStar-l31c.1.fits')[1].data

gd=apselect.select(a,badval=['STAR_BAD'],logg=[0,3.8])
a=a[gd]

hmin=12
for mag in [0] :
    if mag == 0 :
        pand = np.where( (a['andflag'] & bitmask.starflagval('PERSIST_HIGH') > 0 ) &
                         (a['h'] > hmin) & (a['commiss'] == 0) )[0]
        por = np.where( (a['andflag'] & bitmask.starflagval('PERSIST_HIGH') == 0 ) &
                        (a['starflag'] & bitmask.starflagval('PERSIST_HIGH') > 0 ) &
                        (a['h'] > hmin) & (a['commiss'] == 0) )[0]
        pno = np.where( (a['starflag'] & bitmask.starflagval('PERSIST_HIGH') == 0 ) &
                       (a['starflag'] & bitmask.starflagval('PERSIST_MED') == 0 ) &
                       (a['starflag'] & bitmask.starflagval('PERSIST_LOW') == 0 ) &
                       (a['h'] > hmin) & (a['commiss'] == 0) )[0]
    else :
        pand = np.where( (a['andflag'] & bitmask.starflagval('PERSIST_HIGH') > 0 ) &
                         (a['h'] < hmin) & (a['commiss'] == 0) )[0]
        por = np.where( (a['andflag'] & bitmask.starflagval('PERSIST_HIGH') == 0 ) &
                        (a['starflag'] & bitmask.starflagval('PERSIST_HIGH') > 0 ) &
                        (a['h'] < hmin) & (a['commiss'] == 0) )[0]
        pno = np.where( (a['starflag'] & bitmask.starflagval('PERSIST_HIGH') == 0 ) &
                       (a['starflag'] & bitmask.starflagval('PERSIST_MED') == 0 ) &
                       (a['starflag'] & bitmask.starflagval('PERSIST_LOW') == 0 ) &
                       (a['h'] < hmin) & (a['commiss'] == 0) )[0]

    fig,ax = plots.multi(3,1,wspace=0.001,hspace=0.001)
    axim=plots.plotc(ax[0],a['TEFF'][pno],a['LOGG'][pno],a['M_H'][pno],size=2,xr=[5500,3500],yr=[5,0],zr=[-2,0.5],xt='Teff',yt='logg',rasterized=True)
    plots.plotc(ax[1],a['TEFF'][por],a['LOGG'][por],a['M_H'][por],size=2,xr=[5500,3500],yr=[5,0],zr=[-2,0.5],xt='Teff',rasterized=True)
    plots.plotc(ax[2],a['TEFF'][pand],a['LOGG'][pand],a['M_H'][pand],size=2,xr=[5500,3500],yr=[5,0],zr=[-2,0.5],xt='Teff',rasterized=True)
    cbaxes = fig.add_axes([0.91, 0.1, 0.01, 0.8])
    cb = plt.colorbar(axim, cax = cbaxes)
    cb.set_label('[M/H]')
    cbaxes.tick_params(axis='both',labelsize=8)
    fig.savefig('persist_f_hr.pdf')
   
    tags=['C_FE','N_FE','O_FE']
    fig,ax = plots.multi(3,len(tags),wspace=0.001,hspace=0.001)
    for i,tag in enumerate(['C_FE','N_FE','O_FE']) :
        el=tag.split('_')[0]
        axim=plots.plotc(ax[i,0],a['M_H'][pno],a[tag][pno],a['TEFF'][pno],size=2,xr=[-2.5,1],yr=[-0.9,0.9],zr=[3500,5500],xt='[M/H]',yt='['+el.capitalize()+'/Fe]',rasterized=True)
        plots.plotc(ax[i,1],a['M_H'][por],a[tag][por],a['TEFF'][por],size=2,xr=[-2.5,1],yr=[-0.9,0.9],zr=[3500,5500],xt='[M/H]',rasterized=True)
        plots.plotc(ax[i,2],a['M_H'][pand],a[tag][pand],a['TEFF'][pand],size=2,xr=[-2.5,1],yr=[-0.9,0.9],zr=[3500,5500],xt='[M/H]',rasterized=True)
    ax[2,0].tick_params(axis='x',labelsize=8)
    ax[2,1].tick_params(axis='x',labelsize=8)
    ax[2,2].tick_params(axis='x',labelsize=8)
    cbaxes = fig.add_axes([0.91, 0.1, 0.01, 0.8])
    cb = plt.colorbar(axim, cax = cbaxes)
    cb.set_label(r'T$_{\rm eff}$')
    cbaxes.tick_params(axis='both',labelsize=8)
    fig.savefig('persist_f_cno.pdf')

    tags=['MG_FE','SI_FE','S_FE','CA_FE','TI_FE']
    fig,ax = plots.multi(3,len(tags),wspace=0.001,hspace=0.001)
    for i,tag in enumerate(tags) :
        el=tag.split('_')[0]
        axim=plots.plotc(ax[i,0],a['M_H'][pno],a[tag][pno],a['TEFF'][pno],size=2,xr=[-2.5,1],yr=[-0.9,0.9],zr=[3500,5500],xt='[M/H]',yt='['+el.capitalize()+'/Fe]',rasterized=True)
        plots.plotc(ax[i,1],a['M_H'][por],a[tag][por],a['TEFF'][por],size=2,xr=[-2.5,1],yr=[-0.9,0.9],zr=[3500,5500],xt='[M/H]',rasterized=True)
        plots.plotc(ax[i,2],a['M_H'][pand],a[tag][pand],a['TEFF'][pand],size=2,xr=[-2.5,1],yr=[-0.9,0.9],zr=[3500,5500],xt='[M/H]',rasterized=True)
    ax[4,0].tick_params(axis='x',labelsize=8)
    ax[4,1].tick_params(axis='x',labelsize=8)
    ax[4,2].tick_params(axis='x',labelsize=8)
    cbaxes = fig.add_axes([0.91, 0.1, 0.01, 0.8])
    cb = plt.colorbar(axim, cax = cbaxes)
    cb.set_label(r'T$_{\rm eff}$')
    cbaxes.tick_params(axis='both',labelsize=8)
    fig.savefig('persist_f_alpha.pdf')

    tags=['AL_FE','K_FE','MN_FE','NI_FE']
    fig,ax = plots.multi(3,len(tags),wspace=0.001,hspace=0.001)
    for i,tag in enumerate(tags) :
        el=tag.split('_')[0]
        axim=plots.plotc(ax[i,0],a['M_H'][pno],a[tag][pno],a['TEFF'][pno],size=2,xr=[-2.5,1],yr=[-0.9,0.9],zr=[3500,5500],xt='[M/H]',yt='['+el.capitalize()+'/Fe]',rasterized=True)
        plots.plotc(ax[i,1],a['M_H'][por],a[tag][por],a['TEFF'][por],size=2,xr=[-2.5,1],yr=[-0.9,0.9],zr=[3500,5500],xt='[M/H]',rasterized=True)
        plots.plotc(ax[i,2],a['M_H'][pand],a[tag][pand],a['TEFF'][pand],size=2,xr=[-2.5,1],yr=[-0.9,0.9],zr=[3500,5500],xt='[M/H]',rasterized=True)
    ax[3,0].tick_params(axis='x',labelsize=8)
    ax[3,1].tick_params(axis='x',labelsize=8)
    ax[3,2].tick_params(axis='x',labelsize=8)
    cbaxes = fig.add_axes([0.91, 0.1, 0.01, 0.8])
    cb = plt.colorbar(axim, cax = cbaxes)
    cb.set_label(r'T$_{\rm eff}$')
    cbaxes.tick_params(axis='both',labelsize=8)
    fig.savefig('persist_f_fe.pdf')
