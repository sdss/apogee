from __future__ import print_function

import os
import pdb
import numpy as np
import multiprocessing as mp
import subprocess
from astropy.io import fits
from apogee.utils import yanny
from apogee.speclib import atmos
from apogee.aspcap import aspcap
from apogee.utils import spectra
from tools import plots
from tools import html
import matplotlib.pyplot as plt

def fill(planfile='tgGK_180625.par',dir='marcs/giantisotopes/tgGK_180625',
         cmrange=None, nmrange=None, vtrange=None,grid='GK',
         apstar=False,threads=30,fakehole=False,out='rbf_',r0=1.0) :
    """ routine to fill grid holes using RBF interpolation routine from Szabolcs
    """

    # Read planfile and set output file name
    if not os.path.isfile(planfile):
        print('{:s} does not exist'.format(planfile))
        return
    p=yanny.yanny(planfile,np=True)
    if p.get('r0') : r0 = float(p['r0'])

    # input directory 
    if dir is None :
        indir = os.environ['APOGEE_SPECLIB']+'/synth/'+p['specdir']+'/' if p.get('specdir') else './'
    else :
        indir=os.environ['APOGEE_SPECLIB']+'/synth/turbospec/'+dir+'/'
    print('indir: ', indir)


    if cmrange is None : cmrange=spectra.vector(p['cm0'],p['dcm'],p['ncm'])
    if nmrange is None : nmrange=spectra.vector(p['nm0'],p['dnm'],p['nnm'])
    if vtrange is None : 
        try:
            vtrange=10.**spectra.vector(p['vt0'],p['dvt'],p['nvt'])
        except :
            vtrange = [float(p['vmicro'])]

    # get configuration for grid
    # sizes of subgrids for RBF interpolation in [alpha/M], [M/H], logg, and Teff

    teffsize = subgrid(int(p['nteff']))
    amsize = subgrid(int(p['nam']))
    mhsize = subgrid(int(p['nmh']))
    loggsize = subgrid(int(p['nlogg']))

    if grid is None :
        if p['specdir'].find('GK_') >= 0 : grid = 'GK'
        elif p['specdir'].find('M_') >= 0 : grid = 'M'
        elif p['specdir'].find('F_') >= 0 : grid = 'F'

    # read holes file
    holefile='MARCS_'+grid+'_holefile.fits'
    holes=fits.open(os.environ['APOGEE_SPECLIB']+'/atmos/marcs/MARCS_v3_2016/'+holefile)[0]

    # total number of frequencies, and pixels to use
    if apstar :
        prefix=''
        nfreq=aspcap.nw_chip.sum()
        pix_apstar=aspcap.gridPix()
        pix_aspcap=aspcap.gridPix(apStar=False)
    else :
        prefix=''
        file=(prefix+'a{:s}c{:s}n{:s}v{:s}.fits').format(
               atmos.cval(0.),atmos.cval(cmrange[0]),atmos.cval(nmrange[0]),atmos.cval(vtrange[0]))
        grid=fits.open(indir+file)[0]
        nfreq=grid.data.shape[-1]

    # loop over [C/M], [N/M], and vmicro, and set up the subgrids for interpolation
    nrbf=0 
    for icm,cm in enumerate(cmrange) :
      hcm=int(round((cm - holes.header['CRVAL5'] ) / holes.header['CDELT5']))   
      if hcm< 0 :
          #print('off carbon grid edge!',hcm)
          hcm=0
      for inm,nm in enumerate(nmrange) :
       for ivt,vt in enumerate(vtrange) :
         # load into grids of [alpha,mh,logg,teff,wave]
         data=np.zeros([int(p['nam']),int(p['nmh']),int(p['nlogg']),int(p['nteff']),nfreq],dtype=np.float32)
         for iam,am in enumerate(spectra.vector(p['am0'],p['dam'],p['nam'])) :
           file=(prefix+'a{:s}c{:s}n{:s}v{:s}.fits').format(
                 atmos.cval(am),atmos.cval(cm),atmos.cval(nm),atmos.cval(vt))
           grid=fits.open(indir+file)[0]
           if apstar :
               # select out ASPCAP grid wavelengths
               for pasp,pap in zip(pix_aspcap,pix_apstar) :
                   data[iam,:,:,:,pasp[0]:pasp[1]]=grid.data[:,:,:,pap[0]:pap[1]]
           else :
               data[iam,:,:,:,:] = grid.data

         # within each of these grids, do rbf in subsets in parallel, so set up input for subgrids
         pars=[]
         npars=0
         sam=0 
         for nam in amsize :
           smh=0
           for nmh in mhsize :
             slogg=0
             for nlogg in loggsize :
               steff=0
               for nteff in teffsize :
                 # get the subgrid of holes

                 # try to extend subgrid by one to avoid edges
                 am1,am2,ham1,ham2 = extend(sam,sam+nam,spectra.vector(p['am0'],p['dam'],p['nam']),spectra.fits2vector(holes.header,4))
                 mh1,mh2,hmh1,hmh2 = extend(smh,smh+nmh,spectra.fits2vector(grid.header,4),spectra.fits2vector(holes.header,3))
                 logg1,logg2,hlogg1,hlogg2 = extend(slogg,slogg+nlogg,spectra.fits2vector(grid.header,3),spectra.fits2vector(holes.header,2))
                 teff1,teff2,hteff1,hteff2 = extend(steff,steff+nteff,spectra.fits2vector(grid.header,2),spectra.fits2vector(holes.header,1))

                 # get range in holes file
                 ham=int(round((float(p['am0'])+sam*float(p['dam']) - holes.header['CRVAL4'] ) / holes.header['CDELT4']))
                 hmh=int(round((grid.header['CRVAL4']+smh*grid.header['CDELT4'] - holes.header['CRVAL3'] ) / holes.header['CDELT3']))
                 hlogg=int(round((grid.header['CRVAL3']+slogg*grid.header['CDELT3'] - holes.header['CRVAL2'] ) / holes.header['CDELT2']))
                 hteff=int(round((grid.header['CRVAL2']+steff*grid.header['CDELT2'] - holes.header['CRVAL1'] ) / holes.header['CDELT1']))
                 nholes=len(np.where(holes.data[hcm,ham:ham+nam,hmh:hmh+nmh,hlogg:hlogg+nlogg,hteff:hteff+nteff]>0)[0])
                 # add "fake" holes if desired
                 if fakehole : 
                     holes.data[hcm,ham+nam/2,hmh+nmh/2,hlogg+nlogg/2,hteff] += 100.
                     holes.data[hcm,ham+nam/2,hmh+nmh/2,hlogg,hteff+nteff/2] += 100.
                     holes.data[hcm,ham+nam/2,hmh,hlogg+nlogg/2,hteff+nteff/2] += 100.
                     holes.data[hcm,ham,hmh+nmh/2,hlogg+nlogg/2,hteff+nteff/2] += 100.

                 # set up input for RBF:  (name, data, holes)
                 name = out+'c{:s}n{:s}v{:s}_{:02d}'.format(atmos.cval(cm),atmos.cval(nm),atmos.cval(vt),npars)
                 print(name,am1,am2,mh1,mh2,logg1,logg2,teff2,teff2,hcm,ham1,ham2,hmh1,hmh2,hlogg1,hlogg2,hteff1,hteff2)
                 pars.append((name,r0,data[am1:am2,mh1:mh2,logg1:logg2,teff1:teff2,:],
                              np.squeeze(holes.data[hcm,ham1:ham2,hmh1:hmh2,hlogg1:hlogg2,hteff1:hteff2])))
                 npars+=1
                 steff+=nteff
               slogg+=nlogg
             smh+=nmh
           sam+=nam

         # do the RBF in parallel for the subgrids
         pool = mp.Pool(threads)
         specs = pool.map_async(dorbf, pars).get()
         pool.close()
         pool.join()

         # fill in the holes in the full 4D grid with the interpolated spectra that are returned
         #   but don't use the "extended" regions
         ii=0
         sam=0 
         for nam in amsize :
           smh=0
           for nmh in mhsize :
             slogg=0
             for nlogg in loggsize :
               steff=0
               for nteff in teffsize :
                 am1,am2,ham1,ham2 = extend(sam,sam+nam,spectra.vector(p['am0'],p['dam'],p['nam']),spectra.fits2vector(holes.header,4))
                 mh1,mh2,hmh1,hmh2 = extend(smh,smh+nmh,spectra.fits2vector(grid.header,4),spectra.fits2vector(holes.header,3))
                 logg1,logg2,hlogg1,hlogg2 = extend(slogg,slogg+nlogg,spectra.fits2vector(grid.header,3),spectra.fits2vector(holes.header,2))
                 teff1,teff2,hteff1,hteff2 = extend(steff,steff+nteff,spectra.fits2vector(grid.header,2),spectra.fits2vector(holes.header,1))
                 # get range in holes file
                 ham=int(round((float(p['am0'])+sam*float(p['dam']) - holes.header['CRVAL4'] ) / holes.header['CDELT4']))
                 hmh=int(round((grid.header['CRVAL4']+smh*grid.header['CDELT4'] - holes.header['CRVAL3'] ) / holes.header['CDELT3']))
                 hlogg=int(round((grid.header['CRVAL3']+slogg*grid.header['CDELT3'] - holes.header['CRVAL2'] ) / holes.header['CDELT2']))
                 hteff=int(round((grid.header['CRVAL2']+steff*grid.header['CDELT2'] - holes.header['CRVAL1'] ) / holes.header['CDELT1']))
                 nholes=len(np.where(holes.data[hcm,ham:ham+nam,hmh:hmh+nmh,hlogg:hlogg+nlogg,hteff:hteff+nteff]>0)[0])
                 name = out+'c{:s}n{:s}v{:s}_{:02d}'.format(atmos.cval(cm),atmos.cval(nm),atmos.cval(vt),ii)
                 print('loading: ',name,am1,am2,mh1,mh2,logg1,logg2,teff2,teff2,hcm,ham1,ham2,hmh1,hmh2,hlogg1,hlogg2,hteff1,hteff2)
                 # replace data and holes with filled data
                 data[sam:sam+nam,smh:smh+nmh,slogg:slogg+nlogg,steff:steff+nteff,:] = \
                      specs[ii][0][sam-am1:sam-am1+nam,smh-mh1:smh-mh1+nmh,slogg-logg1:slogg-logg1+nlogg,steff-teff1:steff-teff1+nteff,:]
                 try :
                   holes.data[hcm,ham:ham+nam,hmh:hmh+nmh,hlogg:hlogg+nlogg,hteff:hteff+nteff] = \
                      specs[ii][1][ham-ham1:ham-ham1+nam,hmh-hmh1:hmh-hmh1+nmh,hlogg-hlogg1:hlogg-hlogg1+nlogg,hteff-hteff1:hteff-hteff1+nteff]
                 except :
                   pdb.set_trace()
                 ii+=1
                 steff+=nteff
               slogg+=nlogg
             smh+=nmh
           sam+=nam
         # write out the 4D grid across appropreiate output files (separate for each [alpha/m])
         for iam,am in enumerate(spectra.vector(p['am0'],p['dam'],p['nam'])) :
           file=('a{:s}c{:s}n{:s}v{:s}.fits').format(
                 atmos.cval(am),atmos.cval(cm),atmos.cval(nm),atmos.cval(vt))
           grid=fits.open(indir+prefix+file)
           if apstar :
               # fill in apStar wavelengths
               for pasp,pap in zip(pix_aspcap,pix_apstar) :
                   grid[0].data[:,:,:,pap[0]:pap[1]]=data[iam,:,:,:,pasp[0]:pasp[1]]
           else :
                   grid[0].data = data[iam,:,:,:,:]
           grid[0].header.add_comment('holes filled using RBF interpolation')
           grid[0].header.add_comment('APOGEE_VER:'+os.environ['APOGEE_VER'])
           grid[0].header.append(('R0',r0,'value of r0 used for RBF interpolation'))
           ham=int( round( (am-holes.header['CRVAL4']) / holes.header['CDELT4']) )
           # append the modified holes file for this subgrid
           hout=fits.ImageHDU(np.squeeze(holes.data[hcm,ham,:,:,:]))
           for idim in range(1,4) :
               spectra.add_dim(hout.header,holes.header['CRVAL'+str(idim)],holes.header['CDELT'+str(idim)],1,holes.header['CTYPE'+str(idim)],idim) 
           new = fits.HDUList()
           new.append(grid[0])
           new.append(hout)
           new.writeto(indir+out+file,overwrite=True)
           #grid.append(hout)
           #grid.writeto(indir+out+file,overwrite=True)
         #holes.writeto(out+holefile,overwrite=True)

def subgrid(n) :
    ''' return subgrid sizes for input dimension '''
    if n <=4 :
        return [n]
    elif n<=8 :
        return [4,n-4]
    elif n==9 :
        return [3,3,3]
    elif n==10 :
        return [3,4,3]
    elif n==11 :
        return [4,4,3]
    elif n==12 :
        return [4,4,4]
    elif n==13 :
        return [3,3,4,3]
    elif n==14 :
        return [3,4,4,3]
    elif n==15 :
        return [3,3,3,3,3]
    else :  
        print('unknown subgrid config')
        pdb.set_trace()

def dorbf(pars) :
    """ Routine that actually does one RBF interpolation
    """

    # input
    name=pars[0]
    r0=pars[1]
    data=pars[2]
    holes=pars[3]
    ndim=data.shape

    # do the rbf for this subgrid 
    print('doing rbf',name,ndim,holes.shape)

    # if we have no good models at the lowest logg, trim it off
    logg1 = 0
    #gd=np.where(holes[:,:,0,:] < 1.e-4)[0]
    #if len(gd) == 0 :
    #    print('Skipping lowest logg....')
    #    logg1 = 1
    #else :
    #    logg1 = 0

    nhole=0
    nsynhole=0
    ngood=0 
    fgood=open(name+'_good.dat','w')
    fhole=open(name+'_hole.dat','w')
    for iam in np.arange(ndim[0]) :
      for imh in np.arange(ndim[1]) :
        for ilogg in np.arange(logg1,ndim[2]) :
          for iteff in np.arange(ndim[3]) :
            spec=data[iam,imh,ilogg,iteff,:]
            dist=holes[iam,imh,ilogg,iteff]
            if spec.sum() == 0 :
              fhole.write(('{:8.2f} '*4+'\n').format(iam/(ndim[0]-1.),imh/(ndim[1]-1.),ilogg/(ndim[2]-1.),iteff/(ndim[3]-1.)))
              nsynhole+=1
            elif dist>0 : 
              fhole.write(('{:8.2f} '*4+'\n').format(iam/(ndim[0]-1.),imh/(ndim[1]-1.),ilogg/(ndim[2]-1.),iteff/(ndim[3]-1.)))
              nhole+=1
            else : 
              #normalize spectrum
              spec/=np.nanmean(spec)
              fgood.write(('{:8.2f} '*4).format(iam/(ndim[0]-1.),imh/(ndim[1]-1.),ilogg/(ndim[2]-1.),iteff/(ndim[3]-1.)))
              for i in range(ndim[-1]) : fgood.write('{:12.6f}'.format(spec[i]))
              fgood.write('\n')
              ngood+=1
    print('ngood, nsynhole, nhole: ', ngood,nsynhole,nhole)
    fgood.close()
    fhole.close()

    # if no holes, return original data
    if nhole+nsynhole == 0 : 
        os.remove(name+'_good.dat')
        os.remove(name+'_hole.dat')
        return data, holes

    # run rbf
    nd=len(ndim)-1
    cmd=['rbf',name+'_good.dat',name+'_hole.dat',str(r0),str(ngood),str(ngood),str(nd),'1000','1.e-5',str(ndim[-1]+nd),str(nhole+nsynhole),'3']
    for c in cmd : print('{:s} '.format(c),end='')
    print('\n')
    subprocess.call(cmd)

    # read and return rbf output
    filled=np.loadtxt(name+'_hole.dat.filled')
    ii=0
    for iam in np.arange(ndim[0]) :
      for imh in np.arange(ndim[1]) :
        for ilogg in np.arange(logg1,ndim[2]) :
          for iteff in np.arange(ndim[3]) :
            spec=data[iam,imh,ilogg,iteff,:]
            dist=holes[iam,imh,ilogg,iteff]
            if spec.sum() == 0 or dist > 0 :
              try :
                  data[iam,imh,ilogg,iteff,:]=np.atleast_2d(filled)[ii,:]
                  if dist > 0 : holes[iam,imh,ilogg,iteff]=-1.*dist
                  else : holes[iam,imh,ilogg,iteff]=-100.
              except :
                  print("failed to fill holes")
                  data[iam,imh,ilogg,iteff,:]=0.
              ii+=1
    # clean up files
    os.remove(name+'_good.dat.log')
    os.remove(name+'_good.dat')
    os.remove(name+'_hole.dat')
    os.remove(name+'_hole.dat.filled')
    return data, holes

def comp(planfile='tgGK_180625.par',dir='marcs/giantisotopes/tgGK_180625',grid='GK',fakehole=False,
         cmrange=None, nmrange=None, vtrange=None, apstar=False, hard=None,out='rbf_') :

    # Read planfile and set output file name
    if not os.path.isfile(planfile):
        print('{:s} does not exist'.format(planfile))
        return
    p=yanny.yanny(planfile,np=True)

    # input directory 
    if dir is None :
        indir = os.environ['APOGEE_SPECLIB']+'/synth/'+p['specdir']+'/' if p.get('specdir') else './'
    else :
        indir=os.environ['APOGEE_SPECLIB']+'/synth/turbospec/'+dir+'/'

    if grid is None :
        if p['name'].find('GK_') : grid = 'GK'
        elif p['name'].find('M_') : grid = 'M'
        elif p['name'].find('F_') : grid = 'F'

    if cmrange is None : cmrange=spectra.vector(p['cm0'],p['dcm'],p['ncm'])
    if nmrange is None : nmrange=spectra.vector(p['nm0'],p['dnm'],p['nnm'])
    if vtrange is None : 
        try:
            vtrange=10.**spectra.vector(p['vt0'],p['dvt'],p['nvt'])
        except :
            vtrange = [float(p['vmicro'])]

    if grid is None :
        if p['specdir'].find('GK_') >= 0 : grid = 'GK'
        elif p['specdir'].find('M_') >= 0 : grid = 'M'
        elif p['specdir'].find('F_') >= 0 : grid = 'F'

    holefile='MARCS_'+grid+'_holefile.fits'
    if fakehole : holes=fits.open(out+holefile)[0]
    else : holes=fits.open(os.environ['APOGEE_SPECLIB']+'/atmos/marcs/MARCS_v3_2016/'+holefile)[0]

    if apstar :
        prefix=''
    else :
        prefix=''

    ii=0
    nx=5
    ny=6
    figs=[]
    for icm,cm in enumerate(cmrange) :
      for inm,nm in enumerate(nmrange) :
       for ivt,vt in enumerate(vtrange) :
         for iam,am in enumerate(spectra.vector(p['am0'],p['dam'],p['nam'])) :
           print(am,cm,nm,vt)
           file=('a{:s}c{:s}n{:s}v{:s}.fits').format(
                 atmos.cval(am),atmos.cval(cm),atmos.cval(nm),atmos.cval(vt))
           raw=fits.open(indir+prefix+file)[0].data
           filled=fits.open(out+file)[0].data
           #holes=fits.open(out+file)[1].data

           for imh,mh in enumerate(spectra.vector(p['mh0'],p['dmh'],p['nmh'])) :
             for ilogg,logg in enumerate(spectra.vector(p['logg0'],p['dlogg'],p['nlogg'])) :
               for iteff,teff in enumerate(spectra.vector(p['teff0'],p['dteff'],p['nteff'])) :

                 hcm=int(round((cm - holes.header['CRVAL5'] ) / holes.header['CDELT5'])   )
                 if hcm< 0 :
                   #print('off carbon grid edge!',hcm)
                   hcm=0
                 ham=int(round((am - holes.header['CRVAL4'] ) / holes.header['CDELT4']))
                 hmh=int(round((mh - holes.header['CRVAL3'] ) / holes.header['CDELT3']))
                 hlogg=int(round((logg - holes.header['CRVAL2'] ) / holes.header['CDELT2']))
                 hteff=int(round((teff - holes.header['CRVAL1'] ) / holes.header['CDELT1']))

                 if ((fakehole and abs(holes.data[hcm,ham,hmh,hlogg,hteff]) > 99.) or
                     (not fakehole and abs(holes.data[hcm,ham,hmh,hlogg,hteff]) > 0) ): 
                     print(cm,nm,vt,am,mh,logg,teff,holes.data[hcm,ham,hmh,hlogg,hteff])
                     lab='{:6.2f}{:6.2f}{:6.2f}{:6.0f}{:7.2f}'.format(am,mh,logg,teff,holes.data[hcm,ham,hmh,hlogg,hteff])
                     lab1='{:6.2f}{:6.2f}{:6.2f}'.format(cm,nm,vt)
                     if hard :
                         if ii % (nx*ny) ==0 : fig,ax=plots.multi(nx,ny,hspace=0.001,wspace=0.001,figsize=(14,8))
                         iy=ii % ny
                         ix=(ii%(nx*ny))/ny
                         print(ii,iy,ix)
                         #ax[iy,ix].plot(raw[imh,ilogg,iteff,:]/np.nanmean(raw[imh,ilogg,iteff,:]))
                         ax[iy,ix].plot(filled[imh,ilogg,iteff,:],color='g')
                         ax[iy,ix].plot(raw[imh,ilogg,iteff,:]/np.nanmean(raw[imh,ilogg,iteff,:])/filled[imh,ilogg,iteff,:]+0.2)
                         ax[iy,ix].set_ylim([0.7,1.35])
                         if not np.isclose(holes.data[hcm,ham,hmh,hlogg,hteff],0.) and not np.isclose(holes.data[hcm,ham,hmh,hlogg,hteff],-100.) : color = 'r'
                         else : color='k'
                         ax[iy,ix].text(0.01,0.97,lab,transform=ax[iy,ix].transAxes,va='top',fontsize='x-small',color=color)
                         ax[iy,ix].text(0.01,0.9,lab1,transform=ax[iy,ix].transAxes,va='top',fontsize='x-small',color=color)
                         ii+=1
                         if ii % (nx*ny) == 0: 
                             fig.savefig(out+'{:02d}'.format(ii/(nx*ny))+'.png')
                             figs.append([out+'{:02d}'.format(ii/(nx*ny))+'.png'])
                             plt.close()
                     else :
                         plt.clf()
                         plt.plot(raw[imh,ilogg,iteff,:]/np.nanmean(raw[imh,ilogg,iteff,:])-0.2,color='b')
                         plt.plot(filled[imh,ilogg,iteff,:],color='g')
                         plt.plot(raw[imh,ilogg,iteff,:]/np.nanmean(raw[imh,ilogg,iteff,:])/filled[imh,ilogg,iteff,:]+0.2,color='r')
                         plt.ylim([0.7,1.3])
                         plt.draw() 
                         pdb.set_trace()

    html.htmltab(figs,file=out+'.html')

def mkhtml(n=24,r0s=['1.25','1.00','0.75','0.50','0.25']) :

    files=[]
    for i in range(1,n+1) :
      f=[]
      xtit=[]
      for r0 in r0s :
        f.append('rbf'+r0+'_{:02d}.png'.format(i))
        xtit.append(r0)
      files.append(f)

    html.htmltab(files,file='rbf.html',xtitle=xtit)


def extend(start,end,vector,holevector) :
    """ Extend the subgrids one point if possible, to avoid edges
    """
    s=np.max([0,start-1])
    e=np.min([len(vector),end+1])
    hs=np.max([0,np.where(np.isclose(holevector,vector[s]))[0][0]])
    he=np.min([len(holevector),np.where(np.isclose(holevector,vector[e-1]))[0][0]+1])
    if e-s != he-hs :
        print('ERROR: inconsistent size of data and hole grids',s,e,hs,he)
        pdb.set_trace()

    return s,e,hs,he

def mergeholes(planfile='tgGK_180625.par',dir='marcs/giantisotopes/tgGK_180625',grid='GK',fakehole=False,
         cmrange=None, nmrange=None, vtrange=None, apstar=False, hard=None,out='rbf_') :

    # Read planfile and set output file name
    if not os.path.isfile(planfile):
        print('{:s} does not exist'.format(planfile))
        return
    p=yanny.yanny(planfile,np=True)

    if dir is None :
        indir = os.environ['APOGEE_SPECLIB']+'/synth/'+p['specdir']+'/' if p.get('specdir') else './'
    else :
        indir=os.environ['APOGEE_SPECLIB']+'/synth/turbospec/'+dir+'/'

    if cmrange is None : cmrange=spectra.vector(p['cm0'],p['dcm'],p['ncm'])
    if nmrange is None : nmrange=spectra.vector(p['nm0'],p['dnm'],p['nnm'])
    if vtrange is None : 
        try:
            vtrange=10.**spectra.vector(p['vt0'],p['dvt'],p['nvt'])
        except :
            vtrange = [float(p['vmicro'])]

    if grid is None :
        if p['specdir'].find('GK_') >= 0 : grid = 'GK'
        elif p['specdir'].find('M_') >= 0 : grid = 'M'
        elif p['specdir'].find('F_') >= 0 : grid = 'F'

    holefile='MARCS_'+grid+'_holefile.fits'
    holes=fits.open(os.environ['APOGEE_SPECLIB']+'/atmos/marcs/MARCS_v3_2016/'+holefile)

    for icm,cm in enumerate(cmrange) :
      for inm,nm in enumerate(nmrange) :
       for ivt,vt in enumerate(vtrange) :
         for iam,am in enumerate(spectra.vector(p['am0'],p['dam'],p['nam'])) :
           print(am,cm,nm,vt)
           file=('a{:s}c{:s}n{:s}v{:s}.fits').format(
                 atmos.cval(am),atmos.cval(cm),atmos.cval(nm),atmos.cval(vt))
           filled=fits.open(out+file)[1].data
           hcm=int(round((cm - holes[0].header['CRVAL5'] ) / holes[0].header['CDELT5'])   )
           if hcm< 0 :
               #print('off carbon grid edge!',hcm)
               hcm=0
           ham=int(round((am - holes[0].header['CRVAL4'] ) / holes[0].header['CDELT4']))
           holes[0].data[hcm,ham,:,:,:] = filled
    holes.writeto(out+holefile,overwrite=True)
