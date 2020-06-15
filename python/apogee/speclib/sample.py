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

import matplotlib
matplotlib.use('Agg')
import numpy as np
import math
import os
import glob
import pdb
import matplotlib.pyplot as plt
#import thread
from tools import plots
from tools import match
from apogee.speclib import isochrones
from astropy.io import ascii

def elemsens(teffs=[3500,4500,5500],loggs=[1.0,3.0,5.0],mhs=[0.0]) :
    """ create sample with small delta of each element at grid of [teff,logg,mh] to see sensitivities
    """
    els = np.array(['O','Na','Mg','Al','Si','P','S','K','Ca','Ti','V','Cr','Mn','Co','Fe','Ni','Cu','Ge','Rb','Ce','Nd'])
    els_alpha = np.where((els == 'O') | (els == 'Mg') | (els == 'Si') | (els == 'S') | (els == 'Ca') | (els == 'Ti'))[0]

    cm=0.
    nm=0.
    am=0.
    vrot=0.
    files=[]
    for el in np.append(['','C','N'],els):
        if el == '' : name='ref.dat'
        else : name=el+'.dat'
        files.append(name)
        f=open(name,'w')
        f.write("#   Teff   logg    [M/H] [alpha/M] [C/M]   [N/M]  vmicro  vrot")
        for e in els: f.write('{:>7s}'.format(e))
        f.write('\n')
        for teff in teffs :
            for logg in loggs :
                vmicro=10.**(0.226-0.0228*logg+0.0297*logg**2-0.0113*logg**3)
                for mh in mhs :
                    # model with enhanced abundance in desired element
                    if el == 'C' : dcm=0.1
                    else : dcm=0.
                    if el == 'N' : dnm=0.1
                    else : dnm=0.
                    out = '{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}'.format(teff,logg,mh,am,cm+dcm,nm+dnm,vmicro,vrot)      
                    for e in els : 
                      if e == el : ab=0.1 
                      else: ab=0.0
                      out = out + '{:7.2f}'.format(ab)      # other elements
                    f.write(out+'\n')
    return files

def sample(name='test',gridclass=None,eps=0.01,tefflim=[3000,8000],dtlo=100.,logglim=[-0.5,5.5],mhlim=[-2.5,0.75],nmlim=[-0.5,2.],cmlim=[-1.5,1.],emlim=[-0.5,1.],vmicrolim=[0.3,4.8],amlim=[-0.5,1.],vrotlim=[1.5,96.],rot=True,nsamp=1,niso=None,elems='all',fact=1.0,offgrid=False,extracool=1.) :
    """ Generate a test sample of parameters and abundances from isochrones
    """

    # set range around isochrones
    isorange = [-1,0,1]
    # set output limits
    dthi=250
    dthot=250
    dtvhot=250
    dlogg=0.5
    dmh=0.25
    dam=0.25
    if gridclass == 'GKg' :
        tefflim=[3500,6000]
        logglim=[0,4.5]
        dtlo=250.
        rot=False
    elif gridclass == 'Mg' :
        tefflim=[3000,4000]
        logglim=[-0.5,3.0]
        dtlo=100.
        rot=False
    elif gridclass == 'GKd' :
        tefflim=[3500,6000]
        logglim=[2.5,5.5]
        cmlim=[-0.5,0.5]
        nmlim=[-0.5,1.5]
        dtlo=250.
        rot=True
    elif gridclass == 'Md' :
        tefflim=[3000,4000]
        logglim=[2.5,5.5]
        cmlim=[-0.5,0.5]
        nmlim=[-0.5,1.5]
        dtlo=100.
        rot=True
    elif gridclass == 'Fd' :
        tefflim=[5500,8000]
        logglim=[2.5,5.5]
        cmlim=[-0.5,0.5]
        nmlim=[-0.5,1.5]
        dtlo=250.
        rot=True
    elif gridclass == 'coarse' :
        tefflim=[3500,8000]
        logglim=[0.5,5. ]
        cmlim=[-0.5,0.5]
        nmlim=[-0.5,1.5]
        dtlo=500.
        dthi=500.
        dlogg=1.0
        dmh=0.5
        dam=0.5
        rot=True
    elif gridclass == 'rv' :
        tefflim=[3000,20000]
        logglim=[0.5,5. ]
        cmlim=[0.,0.]
        nmlim=[0.,0.]
        emlim=[0.,0.]
        dtlo=200.
        dthi=250.
        dthot=500.
        dtvhot=1000.
        dlogg=1.0
        dmh=0.5
        dam=0.5
        rot=False
        isorange=[0]
    if gridclass is not None : name = name+'_'+gridclass
    grid=[]

    # loop through isochrone data and take grid points nearest and +/- 1 (isorange)
    # accumulate unique set of these
    files = glob.glob(os.environ['ISOCHRONE_DIR']+'/z*.dat')
    if niso is None : niso = len(files)
    for file in files[0:niso] :
        a = isochrones.read(file,agerange=[7,20])
        print(file)
        for i in range(len(a)) :
            if a['teff'][i] < 4000 : dt=dtlo
            elif a['teff'][i] > 8000 : dt=dtvhot
            elif a['teff'][i] > 5500 : dt=dthot
            else : dt = dthi
            for j in isorange :
              teff = (int(round(a['teff'][i]/dt))+j)*int(dt)
              logg = (int(round(a['logg'][i]/dlogg))+j)*dlogg
              mh = (int(round(a['feh'][i]/dmh))+j)*dmh
              # clip to stay within requested grid limits
              teff=clip(teff,tefflim)
              logg=clip(logg,logglim)
              mh=clip(mh,mhlim)
              grid.append(tuple([teff,logg,mh]))
        grid = list(set(grid))
        print(len(grid))

    # output file
    f=open(name,'w')
    fipf=open(name+'.ipf','w')
    pars = ['Teff','logg','[M/H]','[alpha/M]','[C/M]','[N/M]','vmicro','vrot']
    f.write("# ")
    for par in pars: f.write(' {:s}'.format(par))
    allteff=[]
    alllogg=[]
    allmh=[]
    allvmic=[]
    allvrot=[]
    allam=[]
    allcm=[]
    allnm=[]
    els = np.array(['O','Na','Mg','Al','Si','P','S','K','Ca','Ti','V','Cr','Mn','Co','Fe','Ni','Cu','Ge','Rb','Ce','Nd'])

    els_alpha = np.where((els == 'O') | (els == 'Mg') | (els == 'Si') | (els == 'S') | (els == 'Ca') | (els == 'Ti'))[0]
    for el in els: f.write('{:>7s}'.format(el))
    f.write('\n')
    nel=len(els)
    for i,x in enumerate(grid) :
      if a['teff'][i] <= 4000 : nf=extracool
      else : nf=1
      for j in range(nsamp*nf) :
        if offgrid :
            if a['teff'][i] < 4000 : dt=dtlo
            elif a['teff'][i] > 8000 : dt=dtvhot
            elif a['teff'][i] > 5500 : dt=dthot
            teff=x[0]+np.random.uniform(-dt/2.,dt/2.)
            logg=x[1]+np.random.uniform(-dlogg/2.,dlogg/2.)
            mh=x[2]+np.random.uniform(-dmh/2.,dmh/2.)
            teff=clip(teff,tefflim)
            logg=clip(logg,logglim)
            mh=clip(mh,mhlim)
        else :
            teff=x[0]
            logg=x[1]
            mh=x[2]
        vmicro=10.**(0.226-0.0228*logg+0.0297*logg**2-0.0113*logg**3)+np.random.normal(0.,0.3)
        vmicro=clip(vmicro,vmicrolim)
        if (gridclass == 'rv') : 
            vmicro=0.226-0.0228*logg+0.0297*logg**2-0.0113*logg**3
            vmicro=10.**(int(round(vmicro/0.30103))*0.30103 - 0.522878)

        vrot=0.
        if (logg < 3.5) & (teff<6000) :
            # for giants, use vmacro relation + small rotation
            if rot : vrot = np.max([0.,np.random.normal(1.5,0.5*fact)])
            # carbon and nitrogen with significant range
            cm=np.random.normal(-0.20,.5*fact)
            if not offgrid: cm = (int(round(cm/0.25)))*0.25
            nm=np.random.normal(0.3,0.7*fact)
            # no need to pin [N/M] to grid since it is varied in synthesis!
            #nm = (int(round(nm/0.5)))*0.5
        else :
            # for dwarfs, use significant rotation
            if rot : vrot=abs(np.random.normal(0.,30*fact))
            # carbon and nitrogen with small range
            cm=np.random.normal(0.,0.25*fact)
            if not offgrid: cm = (int(round(cm/0.25)))*0.25
            nm=np.random.normal(0.,0.25*fact)
        cm=clip(cm,cmlim)
        nm=clip(nm,nmlim)
        # for RV grid, enhance N for giants
        if (gridclass == 'rv') & (logg < 3) & (teff<6000) : nm=0.25
        # draw a random alpha/M
        am=np.random.uniform(-0.25,0.5)
        if not offgrid: am = (round(am/dam))*dam
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
        el=np.random.normal(0.,0.2*fact,size=nel)
        for ie,e in enumerate(el): 
            if elems == 'all' or els[ie] in elems : el[ie]=clip(e,emlim)
            else : el[ie] = 0.
        el[els_alpha] += am
        for e in el :
          # add element abundances
          out = out + '{:7.2f}'.format(e)      # other elements
        print(out)
        f.write(out+'\n')
        if gridclass == 'rv' :
            # add some alpha-enhanced and carbon-enhanced models
            out = None
            if (mh < -0.5) & (teff < 5000) :
                out = '{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}'.format(teff,logg,mh,am+0.25,cm,nm,vmicro,vrot)      
            elif (mh > -0.5) & (teff < 3500) :
                out = '{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}'.format(teff,logg,mh,am+0.25,cm+0.5,nm,vmicro,vrot)      
            if out is not None :
                for e in el : out = out + '{:7.2f}'.format(e)      # other elements
                f.write(out+'\n')

        # clip to adjust slightly off grid edges for FERRE input file
        teff=clip(teff,tefflim,eps=eps)
        logg=clip(logg,logglim,eps=eps)
        mh=clip(mh,mhlim,eps=eps)
        cm=clip(cm,cmlim,eps=eps)
        if not offgrid : nm=(round(nm/0.5))*0.5
        nm=clip(nm,nmlim,eps=eps)
        am=clip(am,amlim,eps=eps)
        vmicro=round((np.log10(vmicro)-math.log10(vmicrolim[0]))/math.log10(2.))*math.log10(2.)+math.log10(vmicrolim[0])
        vmicro=clip(vmicro,np.log10(np.array(vmicrolim)),eps=eps)
        vrot=round((np.log10(vrot)-math.log10(vrotlim[0]))/math.log10(2.))*math.log10(2.)+math.log10(vrotlim[0])
        vrot=clip(vrot,np.log10(np.array(vrotlim)),eps=eps)

        if rot :
            ipf = '{:s}{:d} {:7.2f} {:7.2f} {:7.2f} {:7.2f} {:7.2f} {:7.2f} {:7.2f} {:8.2f}'.format(
                   os.path.basename(name),i+1,vmicro,cm,nm,am,vrot,mh,logg,teff)
        else :
            ipf = '{:s}{:d} {:7.2f} {:7.2f} {:7.2f} {:7.2f} {:7.2f} {:7.2f} {:8.2f}'.format(
                   os.path.basename(name),i+1,vmicro,cm,nm,am,mh,logg,teff)
        fipf.write(ipf+'\n')

    f.close()
    fipf.close()

    # plots of sample
    allteff=np.array(allteff)
    alllogg=np.array(alllogg)
    allmh=np.array(allmh)
    allvmic=np.array(allvmic)
    allvrot=np.array(allvrot)
    allam=np.array(allam)
    allcm=np.array(allcm)
    allnm=np.array(allnm)
    t=[x[0] for x in grid]
    g=[x[1] for x in grid]
    m=[x[2] for x in grid]
    fig,ax=plots.multi(2,3,hspace=0.001,wspace=0.4)
    plots.plotc(ax[0,0],allteff+np.random.uniform(-30.,30.,size=len(allteff)),alllogg+np.random.uniform(-0.1,0.1,size=len(alllogg)),allmh,
                xr=[8000,2500],yr=[6.,-1],zr=mhlim,zt='[M/H]',colorbar=True,yt='log g')
    plots.plotc(ax[0,1],allteff+np.random.uniform(-30.,30.,size=len(allteff)),alllogg+np.random.uniform(-0.1,0.1,size=len(alllogg)),allam,
                xr=[8000,2500],yr=[6.,-1],zr=amlim,zt='[alpha/M]',colorbar=True,yt='log g')
    plots.plotc(ax[1,0],allteff+np.random.uniform(-30.,30.,size=len(allteff)),alllogg+np.random.uniform(-0.1,0.1,size=len(alllogg)),allvmic,
                xr=[8000,2500],yr=[6.,-1],zr=vmicrolim,zt='vmicro',colorbar=True,yt='log g')
    plots.plotc(ax[1,1],allteff+np.random.uniform(-30.,30.,size=len(allteff)),alllogg+np.random.uniform(-0.1,0.1,size=len(alllogg)),allvrot,
                xr=[8000,2500],yr=[6.,-1],zr=[0,30],zt='vrot',colorbar=True,yt='log g')
    plots.plotc(ax[2,0],allteff+np.random.uniform(-30.,30.,size=len(allteff)),alllogg+np.random.uniform(-0.1,0.1,size=len(alllogg)),allcm,
                xr=[8000,2500],yr=[6.,-1],zr=cmlim,zt='[C/M]',colorbar=True,yt='log g',xt='Teff')
    plots.plotc(ax[2,1],allteff+np.random.uniform(-30.,30.,size=len(allteff)),alllogg+np.random.uniform(-0.1,0.1,size=len(alllogg)),allnm,
                xr=[8000,2500],yr=[6.,-1],zr=nmlim,zt='[N/M]',colorbar=True,xt='Teff',yt='log g')
    fig.savefig(name+'.png')
    plt.close()
    fig,ax=plots.multi(2,2,hspace=0.4,wspace=0.4)
    plots.plotc(ax[0,0],allmh+np.random.uniform(-0.1,0.1,size=len(allmh)),allam+np.random.uniform(-0.1,0.1,size=len(allam)),allteff,
                xr=[-2.5,1],yr=[-0.75,1.],zr=[2500,8000],zt='Teff',colorbar=True,xt='[M/H]',yt='[alpha/M]')
    plots.plotc(ax[0,1],alllogg+np.random.uniform(-0.1,0.1,size=len(alllogg)),allcm+np.random.uniform(-0.1,0.1,size=len(allam)),allteff,
                xr=[6,-1],yr=[-1.5,2.],zr=[2500,8000],zt='Teff',colorbar=True,xt='log g',yt='[C/M]')
    plots.plotc(ax[1,1],alllogg+np.random.uniform(-0.1,0.1,size=len(alllogg)),allnm+np.random.uniform(-0.1,0.1,size=len(allam)),allteff,
                xr=[6,-1],yr=[-1.5,2.],zr=[2500,8000],zt='Teff',colorbar=True,xt='log g',yt='[N/M]')
    fig.tight_layout()
    fig.savefig(name+'_2.png')
    plt.close()

def dclip(d,lim=[-0.5,0.5]) :
    d[np.where(d < lim[0])]=lim[0]
    d[np.where(d > lim[1])]=lim[1]
    return d

def comp(file,true=None,truespec=None,hard=False,plot=False,minchi2=0.,testid=None,rot=False) :
    """ Compare input parameters with output results
    """
    if true is None: true=file+'.ipf'
    if rot :
        names=['id','vmicro','cm','nm','am','vrot','mh','logg','teff']
        names_spm=['id','vmicro','cm','nm','am','vrot','mh','logg','teff',
                   'evmicro','ecm','enm','eam','evrot','emh','elogg','eteff','a','b','chi2']
    else :
        names=['id','vmicro','cm','nm','am','mh','logg','teff']
        names_spm=['id','vmicro','cm','nm','am','mh','logg','teff',
                   'evmicro','ecm','enm','eam','emh','elogg','eteff','a','b','chi2']
    true=ascii.read(true,names=names)
    ##spec=np.loadtxt('test.dat')

    obs=ascii.read(file+'.spm')
    for i in range(len(names_spm) ):
        obs.rename_column('col{:d}'.format(i+1),names_spm[i])
    i1,i2=match.match(true['id'],obs['id'])

    # write out file with differences
    f=open(file+'.out','w')
    for i in range(len(i1)) :
        f.write(('{:<15s}'+'{:7.2f}'*14+'\n').format(true[i1[i]]['id'],
                 true[i1[i]]['teff'],true[i1[i]]['logg'],true[i1[i]]['mh'],true[i1[i]]['am'],true[i1[i]]['cm'],
                 true[i1[i]]['nm'],true[i1[i]]['vmicro'],
                 obs[i2[i]]['teff']-true[i1[i]]['teff'],obs[i2[i]]['logg']-true[i1[i]]['logg'],obs[i2[i]]['mh']-true[i1[i]]['mh'],
                 obs[i2[i]]['am']-true[i1[i]]['am'],obs[i2[i]]['cm']-true[i1[i]]['cm'],
                 obs[i2[i]]['nm']-true[i1[i]]['nm'],obs[i2[i]]['vmicro']-true[i1[i]]['vmicro']))
    f.close()

    # histogram of differences, in linear and log histograms
    fig,ax=plots.multi(8,2,wspace=0.001,hspace=0.001,figsize=(12,4),xtickrot=60)
    for iy,log in enumerate([False,True]) :
        ax[iy,0].hist(dclip(obs[i2]['teff']-true[i1]['teff'],lim=[-200,200]),bins=np.arange(-200,201,10),histtype='step',log=log)  
        ax[iy,0].set_xlabel('$\Delta$Teff')
        ax[iy,1].hist(dclip(obs[i2]['logg']-true[i1]['logg']),bins=np.arange(-0.5,0.51,0.01),histtype='step',log=log)  
        ax[iy,1].set_xlabel('$\Delta$logg')
        ax[iy,2].hist(dclip(obs[i2]['mh']-true[i1]['mh']),bins=np.arange(-0.5,0.51,0.01),histtype='step',log=log)  
        ax[iy,2].set_xlabel('$\Delta$[M/H]')
        ax[iy,3].hist(dclip(obs[i2]['am']-true[i1]['am']),bins=np.arange(-0.5,0.51,0.01),histtype='step',log=log)  
        ax[iy,3].set_xlabel('$\Delta$[alpha/M]')
        ax[iy,4].hist(dclip(obs[i2]['cm']-true[i1]['cm']),bins=np.arange(-0.5,0.51,0.01),histtype='step',log=log)  
        ax[iy,4].set_xlabel('$\Delta$[C/M]')
        ax[iy,5].hist(dclip(obs[i2]['nm']-true[i1]['nm']),bins=np.arange(-0.5,0.51,0.01),histtype='step',log=log)  
        ax[iy,5].set_xlabel('$\Delta$[N/M]')
        ax[iy,6].hist(dclip(obs[i2]['vmicro']-true[i1]['vmicro']),bins=np.arange(-0.5,0.51,0.01),histtype='step',log=log)  
        ax[iy,6].set_xlabel('$\Delta$vmicro')
        ax[iy,7].hist(dclip(10.**obs['chi2'],lim=[0,50]),bins=np.arange(0,51,0.1),log=log)  
        ax[iy,7].set_xlabel('chi2')
    fig.suptitle(file)
    if hard :
        fig.savefig(file+'.png')
        plt.close()

    # plots of differences vs Teff, color-coded by various quantities
    fig,ax=plots.multi(6,7,hspace=0.001,wspace=0.001,figsize=(16,8),xtickrot=60)
    x = true['teff'][i1] + np.random.uniform(-45.,45.,size=len(i1))
    for ix in range(6) :
      yt=''
      if ix == 0 : 
        z=true['logg'][i1]
        tit='color: logg'
        zr=[0,5]
      elif ix == 1 : 
        z=true['mh'][i1]
        tit='color: [M/H]'
        zr=[-1.5,0.5]
      elif ix == 2 : 
        z=true['mh'][i1]+true['am'][i1]
        tit='color: [alpha/H]'
        zr=[-1.5,1.5]
      elif ix == 3 : 
        z=true['mh'][i1]+true['cm'][i1]
        tit='color: [C/H]'
        zr=[-1.5,1.5]
      elif ix == 4 : 
        z=true['mh'][i1]+true['nm'][i1]
        tit='color: [N/H]'
        zr=[-1.5,1.5]
      elif ix == 5 : 
        z=obs['chi2'][i2]
        tit='color: chi2'
        zr=[0,10.]
      ax[0,ix].set_title(tit)
      if ix == 0 :
        plots.plotc(ax[0,ix],x,obs['teff'][i2]-true['teff'][i1],z,xt='Teff',yt=r'$\Delta$Teff',yr=[-1000,1000],zr=zr)
        plots.plotc(ax[1,ix],x,obs['logg'][i2]-true['logg'][i1],z,xt='Teff',yt=r'$\Delta$logg',yr=[-2.0,2.0],zr=zr)
        plots.plotc(ax[2,ix],x,obs['mh'][i2]-true['mh'][i1],z,xt='Teff',yt=r'$\Delta$[M/H]',yr=[-0.5,0.5],zr=zr)
        plots.plotc(ax[3,ix],x,obs['am'][i2]-true['am'][i1],z,xt='Teff',yt=r'$\Delta$[a/M]',yr=[-0.5,0.5],zr=zr)
        plots.plotc(ax[4,ix],x,obs['cm'][i2]-true['cm'][i1],z,xt='Teff',yt=r'$\Delta$[C/M]',yr=[-1.5,1.5],zr=zr)
        plots.plotc(ax[5,ix],x,obs['nm'][i2]-true['nm'][i1],z,xt='Teff',yt=r'$\Delta$[N/M]',yr=[-1.5,1.5],zr=zr)
        plots.plotc(ax[6,ix],x,10.**obs['vmicro'][i2]-10.**true['vmicro'][i1],z,xt='Teff',yt=r'$\Delta$vmicro',yr=[-1.0,1.0],zr=zr)
      else :
        plots.plotc(ax[0,ix],x,obs['teff'][i2]-true['teff'][i1],z,xt='Teff',yr=[-1000,1000],zr=zr)
        plots.plotc(ax[1,ix],x,obs['logg'][i2]-true['logg'][i1],z,xt='Teff',yr=[-2,2],zr=zr)
        plots.plotc(ax[2,ix],x,obs['mh'][i2]-true['mh'][i1],z,xt='Teff',yr=[-0.5,0.5],zr=zr)
        plots.plotc(ax[3,ix],x,obs['am'][i2]-true['am'][i1],z,xt='Teff',yr=[-0.5,0.5],zr=zr)
        plots.plotc(ax[4,ix],x,obs['cm'][i2]-true['cm'][i1],z,xt='Teff',yr=[-1.5,1.5],zr=zr)
        plots.plotc(ax[5,ix],x,obs['nm'][i2]-true['nm'][i1],z,xt='Teff',yr=[-1.5,1.5],zr=zr)
        plots.plotc(ax[6,ix],x,10.**obs['vmicro'][i2]-10.**true['vmicro'][i1],z,xt='Teff',yr=[-1.0,1.0],zr=zr)
    fig.suptitle(file)
    plt.show()
    if hard :
        fig.savefig(file+'_2.png')
        plt.close()
    plt.show()

    if plot :
        pdb.set_trace()
        obsspec=np.loadtxt(file+'.mdl')
        if truespec is None : truespec=file+'.frd'
        truespec=np.loadtxt(truespec)
        if testid is None : testid = range(1,len(i1)+1) 
        for tid in testid :
          i=np.where(np.core.defchararray.find(obs['id'][i2],'test'+str(tid))>=0)[0][0]
          if 10.**obs['chi2'][i2[i]] > minchi2 :
            plt.clf()
            plt.plot(truespec[i1[i],:],color='b')
            plt.plot(obsspec[i2[i],:],color='r')
            plt.plot(obsspec[i2[i],:]/truespec[i1[i],:]+0.1,color='g')
            plt.draw()
            print(true['id'][i1[i]])
            print('{:8.1f}{:7.2f}{:7.2f}{:7.2f}{:7.2f}{:7.2f}{:7.2f}'.format(
                  true['teff'][i1[i]],true['logg'][i1[i]],true['mh'][i1[i]],true['am'][i1[i]],true['cm'][i1[i]],true['nm'][i1[i]],true['vmicro'][i1[i]]))
            print('{:8.1f}{:7.2f}{:7.2f}{:7.2f}{:7.2f}{:7.2f}{:7.2f}{:10.2f}'.format(
                  obs['teff'][i2[i]],obs['logg'][i2[i]],obs['mh'][i2[i]],obs['am'][i2[i]],obs['cm'][i2[i]],obs['nm'][i2[i]],obs['vmicro'][i2[i]],obs['chi2'][i2[i]]))
            pdb.set_trace()

def clip(x,lim,eps=None) :
    """ Utility routine to clip values within limits, and move slightly off edges if requested
    """
    # set negative zero to zero
    if np.isclose(x,0.) : x=0.
    # clip to limits
    tmp=np.max([lim[0],np.min([lim[1],x])])
    # move off limit if requested
    if eps is not None :
        if np.isclose(tmp,lim[0]) : tmp+=eps
        if np.isclose(tmp,lim[1]) : tmp-=eps
    return tmp

#def __main__() :
#    for gridclass in ['GKg','Mg','Fd','GKd','Md'] :
#        thread.start_new_thread(sample,('test',gridclass))

