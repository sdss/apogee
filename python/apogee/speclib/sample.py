# encoding: utf-8
#
# @Author: Jon Holtzman
# @Date: March 2018
# @Filename: synth.py
# @License: BSD 3-Clause
# @Copyright: Jon Holtzman

# Routines for making a "representative" sample of stellar parameter/abundance combinations
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import numpy as np
import os
import glob
from tools import plots
from apogee.speclib import isochrones
from astropy.io import ascii

def sample(name='test',gridclass=None,eps=0.01,tefflim=[3000,8000],dtlo=100.,logglim=[-0.5,5.5],mhlim=[-2.5,0.75],nmlim=[-0.5,2.],cmlim=[-1.5,1.],emlim=[-0.5,1.],vmicrolim=[0.5,8.],amlim=[-0.5,1.],rot=True,nsamp=1) :
    """ Generate a test sample of parameters and abundances from isochrones
    """

    # set output limits
    if gridclass == 'GK' :
        tefflim=[3500,6000]
        logglim=[0,4.5]
        dtlo=250.
    elif gridclass == 'M' :
        tefflim=[3000,4000]
        logglim=[-0.5,3.0]
        dtlo=100.
    elif gridclass == 'F' :
        tefflim=[5500,8000]
        logglim=[2.0,5.5]
        dtlo=250.
    grid=[]

    # loop through isochrone data and take grid points nearest and +/- 1
    # accumulate unique set of these
    files = glob.glob(os.environ['ISOCHRONE_DIR']+'/z*.dat')
    for file in files :
        a = isochrones.read(file,agerange=[7,20])
        print(file)
        for i in range(len(a)) :
            if a['teff'][i] < 4000 : dt=dtlo
            else : dt = 250.
            for j in range(-1,2) :
              teff = (int(round(a['teff'][i]/dt))+j)*int(dt)
              logg = (int(round(a['logg'][i]/0.5))+j)*0.5
              mh = (int(round(a['feh'][i]/0.25))+j)*0.25
              # clip to stay within requested grid limits
              teff=clip(teff,tefflim)
              logg=clip(logg,logglim)
              mh=clip(mh,mhlim)
              grid.append(tuple([teff,logg,mh]))
        grid = list(set(grid))
        print(len(grid))

    # output file
    f=open(name,'w')
    finp=open(name+'.inp','w')
    f.write("#   Teff   logg    [M/H] [alpha/M] [C/M]   [N/M]  vmicro  vrot")
    allteff=[]
    alllogg=[]
    allmh=[]
    allvmic=[]
    allvrot=[]
    allam=[]
    allcm=[]
    allnm=[]
    els = np.array(['O','Na','Mg','Al','Si','P','S','K','Ca','Ti','V','Cr','Mn','Co','Ni','Cu','Ge','Rb','Ce','Nd'])
    els_alpha = np.where((els == 'O') | (els == 'Mg') | (els == 'Si') | (els == 'S') | (els == 'Ca') | (els == 'Ti'))[0]
    for el in els: f.write('{:>7s}'.format(el))
    f.write('\n')
    nel=len(els)
    for i,x in enumerate(grid) :
      for j in range(nsamp) :
        teff=x[0]
        logg=x[1]
        mh=x[2]
        vmicro=10.**(0.226-0.0228*logg+0.0297*logg**2-0.0113*logg**3)+np.random.normal(0.,0.3)
        vmicro=clip(vmicro,vmicrolim)
        vrot=0.
        if (logg < 3) & (teff<6000) :
            # for giants, use vmacro relation + small rotation
            if rot : vrot = np.max([0.,np.random.normal(1.5,0.5)])
            # carbon and nitrogen with significant range
            cm=np.random.normal(0.,0.5)
            cm = (int(round(cm/0.25)))*0.25
            nm=np.random.normal(0.3,1.0)
            # no need to pin [N/M] to grid since it is varied in synthesis!
            #nm = (int(round(nm/0.5)))*0.5
        else :
            # for dwarfs, use significant rotation
            if rot : vrot=abs(np.random.normal(0.,30))
            # carbon and nitrogen with small range
            cm=np.random.normal(0.,0.3)
            cm = (int(round(cm/0.25)))*0.25
            nm=np.random.normal(0.,0.3)
        cm=clip(cm,cmlim)
        nm=clip(nm,nmlim)
        am=np.random.uniform(-0.25,0.5)
        am = (round(am/0.25))*0.25
        am=clip(am,amlim)
        allteff.append(teff)
        alllogg.append(logg)
        allmh.append(mh)
        allvmic.append(vmicro)
        allvrot.append(vrot)
        allam.append(am)
        allcm.append(cm)
        allnm.append(nm)

        out = '{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}'.format(teff,logg,mh,am,cm,nm,vmicro,vrot)      
        # individual elemental abundances
        el=np.random.normal(0.,0.2,size=nel)
        for i,e in enumerate(el): el[i]=clip(e,emlim)
        el[els_alpha] += am
        for e in el :
          # add element abundances
          out = out + '{:7.2f}'.format(e)      # other elements
        print(out)
        f.write(out+'\n')

        # clip to adjust slightly off grid edges for FERRE input file
        teff=clip(teff,tefflim,eps=eps)
        logg=clip(logg,logglim,eps=eps)
        mh=clip(mh,mhlim,eps=eps)
        cm=clip(cm,cmlim,eps=eps)
        nm=clip(nm,nmlim,eps=eps)
        am=clip(am,amlim,eps=eps)

        inp = '{:s}{:d} {:7.2f} {:7.2f} {:7.2f} {:7.2f} {:7.2f} {:7.2f} {:8.2f}'.format(
               name,i+1,np.log10(vmicro),cm,nm,am,mh,logg,teff)
        finp.write(inp+'\n')

    f.close()
    finp.close()

    # plots of sample
    t=[x[0] for x in grid]
    g=[x[1] for x in grid]
    m=[x[2] for x in grid]
    fig,ax=plots.multi(2,3)
    plots.plotc(ax[0,0],allteff+np.random.normal(0.,25.,size=len(allteff)),alllogg+np.random.normal(0.,0.05,size=len(alllogg)),allmh,
                xr=[8000,2500],yr=[6.,-1],zr=mhlim,zt='[M/H]',colorbar=True)
    plots.plotc(ax[0,1],allteff+np.random.normal(0.,25.,size=len(allteff)),alllogg+np.random.normal(0.,0.05,size=len(alllogg)),allam,
                xr=[8000,2500],yr=[6.,-1],zr=amlim,zt='[alpha/M]',colorbar=True)
    plots.plotc(ax[1,0],allteff+np.random.normal(0.,25.,size=len(allteff)),alllogg+np.random.normal(0.,0.05,size=len(alllogg)),allvmic,
                xr=[8000,2500],yr=[6.,-1],zr=vmicrolim,zt='vmicro',colorbar=True)
    plots.plotc(ax[1,1],allteff+np.random.normal(0.,25.,size=len(allteff)),alllogg+np.random.normal(0.,0.05,size=len(alllogg)),allvrot,
                xr=[8000,2500],yr=[6.,-1],zr=[0,30],zt='vrot',colorbar=True)
    plots.plotc(ax[2,0],allteff+np.random.normal(0.,25.,size=len(allteff)),alllogg+np.random.normal(0.,0.05,size=len(alllogg)),allcm,
                xr=[8000,2500],yr=[6.,-1],zr=cmlim,zt='[C/M]',colorbar=True)
    plots.plotc(ax[2,1],allteff+np.random.normal(0.,25.,size=len(allteff)),alllogg+np.random.normal(0.,0.05,size=len(alllogg)),allnm,
                xr=[8000,2500],yr=[6.,-1],zr=nmlim,zt='[N/M]',colorbar=True)
    fig.tight_layout()
    fig.savefig(name+'.png')

def comp(file,true='test.inp',hard=False) :
    """ Compare input parameters with output results
    """
    true=ascii.read(true,names=['id','vmicro','cm','nm','am','mh','logg','teff'])
    ##spec=np.loadtxt('test.dat')

    obs=ascii.read(file+'.spm',names=['id','vmicro','cm','nm','am','mh','logg','teff','evm','ecm','enm','eam','emh','elogg','eteff','a','b','c'])
    #mdl=np.loadtxt(file+'.out')
    i1,i2=match.match(true['id'],obs['id'])

    fig,ax=plots.multi(2,7,hspace=0.001,wspace=0.5)
    plots.plotc(ax[0,0],obs['teff'][i2],obs['teff'][i2]-true['teff'][i1],true['mh'][i1],xt='Teff',yt=r'$\Delta$Teff',yr=[-200,200])
    plots.plotc(ax[1,0],obs['teff'][i2],obs['logg'][i2]-true['logg'][i1],true['mh'][i1],xt='Teff',yt=r'$\Delta$logg',yr=[-0.5,0.5])
    plots.plotc(ax[2,0],obs['teff'][i2],obs['mh'][i2]-true['mh'][i1],true['mh'][i1],xt='Teff',yt=r'$\Delta$[M/H]',yr=[-0.5,0.5])
    plots.plotc(ax[3,0],obs['teff'][i2],obs['am'][i2]-true['am'][i1],true['mh'][i1],xt='Teff',yt=r'$\Delta$[a/M]',yr=[-0.5,0.5])
    plots.plotc(ax[4,0],obs['teff'][i2],obs['cm'][i2]-true['cm'][i1],true['mh'][i1],xt='Teff',yt=r'$\Delta$[C/M]',yr=[-0.5,0.5])
    plots.plotc(ax[5,0],obs['teff'][i2],obs['nm'][i2]-true['nm'][i1],true['mh'][i1],xt='Teff',yt=r'$\Delta$[N/M]',yr=[-0.5,0.5])
    plots.plotc(ax[6,0],obs['teff'][i2],10.**obs['vmicro'][i2]-10.**true['vmicro'][i1],true['mh'][i1],xt='Teff',yt=r'$\Delta$vmicro',yr=[-0.5,0.5])
    ax[0,1].hist(obs['teff'][i2]-true['teff'][i1],bins=np.arange(-200,200,10),histtype='step')
    ax[1,1].hist(obs['logg'][i2]-true['logg'][i1],bins=np.arange(-0.5,0.5,0.01),histtype='step')
    ax[2,1].hist(obs['mh'][i2]-true['mh'][i1],bins=np.arange(-0.5,0.5,0.01),histtype='step')
    ax[3,1].hist(obs['am'][i2]-true['am'][i1],bins=np.arange(-0.5,0.5,0.01),histtype='step')
    ax[4,1].hist(obs['cm'][i2]-true['cm'][i1],bins=np.arange(-0.5,0.5,0.01),histtype='step')
    ax[5,1].hist(obs['nm'][i2]-true['nm'][i1],bins=np.arange(-0.5,0.5,0.01),histtype='step')
    ax[6,1].hist(obs['vmicro'][i2]-true['vmicro'][i1],bins=np.arange(-0.5,0.5,0.01),histtype='step')
    fig.suptitle(file)
    plt.show()
    if hard :
        fig.savefig(file+'.png')

def clip(x,lim,eps=None) :
    """ Utility routine to clip values within limits, and move slightly off edges if requested
    """
    # set negative zero to zero
    if np.isclose(x,0.) : x=0.
    # clip to limits
    tmp=np.max([lim[0],np.min([lim[1],x])])
    # move off limit if requested
    if eps is not None :
        if np.isclose(x,lim[0]) : tmp+=eps
        if np.isclose(x,lim[1]) : tmp-=eps
    return tmp

