# encoding: utf-8
#
# @Author: Jon Holtzman
# @Date: March 2018
# @Filename: synth.py
# @License: BSD 3-Clause
# @Copyright: Jon Holtzman

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import errno
import glob
import pdb
import shutil
import subprocess
import sys
import tempfile
import multiprocessing as mp
import numpy as np
import time
import struct
from apogee.aspcap import ferre
from apogee.speclib import atmos
from apogee.utils import atomic
from apogee.utils import spectra
#from sdss.utilities import yanny
from sdss import yanny
from astropy.io import ascii
from astropy.io import fits
from sklearn.decomposition import PCA
from sklearn.decomposition import IncrementalPCA
import matplotlib.pyplot as plt
from tools import plots

colors=['r','g','b','c','m','y']

def showtime(string) :
    print(string+' {:8.2f}'.format(time.time()))
    sys.stdout.flush()

def vector(header,axis) :
    caxis='{:1d}'.format(axis)
    return header['CRVAL'+caxis]+header['CDELT'+caxis]*np.arange(header['NAXIS'+caxis])

def kurucz2turbo(infile,outfile,trim=0) :
    """ Convert Kurucz model atmosphere for use by Turbospectrum 
        Allow for trimming of layers
    """
    try :
        fp=open(infile,'r')
    except :
        fp=open(infile+'.filled','r')
    fout=open(outfile,'w')
    layer=-1
    nlayers=-1
    for iline,line in enumerate(fp) :
        if iline == 0 : 
            # get parameters
            teff=line.split()[1]
            logg=float(line.split()[3])
        if line.split()[0] == 'READ' :
            # write header line once we know how many layers
            nlayers=int(line.split()[2])
            layer=0
            fout.write('KURUCZ {:d} 5000. {:8.2f} 0 0.\n'.format(nlayers-trim,logg))
        elif layer >= 0 and layer < nlayers:
            # write the layers
            layer+=1
            if layer > trim: fout.write(line)
    fp.close()
    fout.close()

def marcs2turbo(infile,outfile,trim=0) :
    """ Prepare MARCS input model for Turbospectrum, allowing for trimming of layers
    """
    try:
        fp=open(infile,'r')
    except :
        try :
            fp=open(infile+'.filled','r')
        except :
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filename)
            return

    fout=open(outfile,'w')
    layer = 0
    for iline,line in enumerate(fp) :
        if line.find('Number of depth') >= 0 :
            nlayers = int(line.split()[0])
            layer = 1
            fout.write(str(nlayers-trim)+' Number of depth points\n')
        else :
            try :
                i = int(line.split()[0])
                if i > trim : fout.write(line)
            except :
                fout.write(line)

def mkturbospec(teff,logg,mh,am,cm,nm,wrange=[15100.,17000],dw=0.05,vmicro=2.0,solarisotopes=False,elemgrid='',welem=None,
    els=None,atmod=None,kurucz=True,atmosroot=None,atmosdir=None,nskip=0,endskip=0,
    linelist='20150714',h2o=None,linelistdir=None,
    save=False,run=True,split=200,fluxcol=2) :

    """ Runs Turbospectrum for specified input parameters
    """

    # directory setup
    if atmosroot is None : atmosroot=os.environ['APOGEE_SPECLIB']+'/atmos/'
    if linelistdir is None : linelistdir=os.environ['APOGEE_SPECLIB']+'/linelists/' 
 
    # atmosphere and filename set up 
    geo = 'p'
    if kurucz : 
        atmoscode= 'k'
        model = 'Kurucz'
        if atmosdir is None : atmosdir = '/kurucz/'
    else : 
        atmoscode = 'm'
        model = 'MARCS'
        if atmosdir is None : atmosdir = '/marcs/MARCS_v3_2016/'
        if logg < 3 : geo = 's'
    atmosdir=atmosroot+'/'+atmosdir+'/'

    if solarisotopes: prefix=atmoscode+'d' 
    else : prefix=atmoscode+'g'

    # output directory and filename
    if save :
        workdir=(prefix+'m{:s}a{:s}c{:s}n{:s}v{:s}'+elemgrid).format(
                 atmos.cval(mh),atmos.cval(am),atmos.cval(cm),atmos.cval(nm),atmos.cval(vmicro))
        workdir=os.environ['APOGEE_LOCALDIR']+'/'+workdir
        try: os.mkdir(workdir)
        except: pass
    else :
        workdir=tempfile.mkdtemp(dir=os.environ['APOGEE_LOCALDIR'])

    root=workdir+'/'+(prefix+'t{:04d}g{:s}m{:s}a{:s}c{:s}n{:s}v{:s}'+elemgrid).format(teff, atmos.cval(logg), 
                      atmos.cval(mh), atmos.cval(am), atmos.cval(cm), atmos.cval(nm),atmos.cval(vmicro))
    if save :
        stdout = None
    else :
        stdout = open(os.devnull, 'w')

    # atmosphere: make local copy in Turbospectrum input format, allowing for trimmed layers
    # note that [N/M] is solar for atmospheres (since it doesn't have a big effect)
    if atmod is None :
        atmod=atmosdir+atmos.filename(teff,logg,mh,cm,am,model=model)
    if kurucz :
        if nskip == 0 : trim=0
        if nskip == 1 : trim=7
        if nskip == 2 : trim=15
        if nskip > 2 : return 0.
        kurucz2turbo(atmod,workdir+'/'+os.path.basename(atmod),trim=trim )
    else :
        try :
            marcs2turbo(atmod,workdir+'/'+os.path.basename(atmod),trim=nskip )
        except:
            return 0.

    # Turbospectrum setup
    try: os.symlink(os.environ['APOGEE_DIR']+'/src/turbospec/DATA',workdir+'/DATA')
    except: pass

    # individual element grid?
    nels = 3
    if elemgrid != '' :
        elemnum=atomic.periodic(elemgrid)[0]
        eabun=atomic.solar(elemgrid)[0]+np.arange(-0.75,1.20,0.25)
    else :
        eabun=np.array([0.])
    if els is not None :
        nels += len(els)

    # welem only computes in windows, but it is easier/faster to compute the whole range with a linelist that only 
    #   has lines in windows!
    wrange=np.array(wrange)
    if welem is not None :
        print('welem not yet implemented!')
        pdb.set_trace()
    else :
        welem=wrange

    sz=welem.shape
    if len(welem) > 1 : 
        # not tested!
        nrange = welem[0]
    else :
        drange=1

    # loop over individual element abundances
    for ielem,abun in enumerate(eabun) :
        file = root+'{:02d}'.format(ielem)
        # only compute opacities for a single nominal abundance
        if ielem == 0 :
            fout=open(root+'_babsma.csh','w')
            fout.write("#!/bin/csh -f\n")
            fout.write("{:s}/bin/babsma_lu << EOF\n".format(os.environ['APOGEE_DIR']))
            fout.write("'LAMBDA_MIN:'   '{:12.3f}'\n".format(welem.min()-dw))
            fout.write("'LAMBDA_MAX:'   '{:12.3f}'\n".format(welem.max()+dw))
            fout.write("'LAMBDA_STEP:'  '{:8.3f}'\n".format(dw))
            fout.write("'MODELINPUT:'  '{:s}'\n".format(os.path.basename(atmod)))
            if kurucz : fout.write("'MARCS-FILE:'  '.false.'\n")
            fout.write("'MODELOPAC:'  '{:s}opac'\n".format(os.path.basename(root)))
            fout.write("'METALLICITY:'  '{:8.3f}'\n".format(mh))
            fout.write("'ALPHA/Fe:'  '{:8.3f}'\n".format(am))
            fout.write("'HELIUM:'  '{:8.3f}'\n".format(0.00))
            fout.write("'R-PROCESS:'  '{:8.3f}'\n".format(0.00))
            fout.write("'S-PROCESS:'  '{:8.3f}'\n".format(0.00))
            fout.write("'INDIVIDUAL ABUNDANCES:'  '{:2d}'\n".format(nels))
            fout.write("    6  {:8.3f}\n".format(8.39+mh+cm))
            fout.write("    7  {:8.3f}\n".format(7.78+mh+nm))
            fout.write("    8  {:8.3f}\n".format(8.66+mh+am))
            if els is not None :
                for el in els :
                  num=atomic.periodic(el[0])[0]
                  ab=atomic.solar(el[0])[0]+el[1]
                  fout.write("{:5d}  {:8.3f}\n".format(num,ab+mh))
            if not solarisotopes :
              fout.write("'ISOTOPES:'  '2'\n")
              # adopt ratio of 12C/13C=15
              fout.write("   6.012 0.9375\n")
              fout.write("   6.013 0.0625\n")
            fout.write("'XIFIX:'  'T'\n")
            fout.write("{:8.3f}\n".format(vmicro))
            fout.write("EOF\n")
            fout.close()
            if run :
                os.chmod(root+'_babsma.csh', 0o777)
                cwd = os.getcwd()
                os.chdir(workdir)
                subprocess.call(['time','./'+os.path.basename(root)+'_babsma.csh'],stdout=stdout)
            if elemgrid != '' : nels+=1

        # create bsyn control file
        bsynfile = file+'bsyn{:02d}.inp'.format(ielem)
        fout = open(bsynfile,'w')
        fout.write("'LAMBDA_STEP:'  '{:8.3f}'\n".format(dw))
        fout.write("'LAMBDA_MIN:'   '{:12.3f}'\n".format(welem.min()))
        fout.write("'LAMBDA_MAX:'   '{:12.3f}'\n".format(welem.max()))
        fout.write("'INTENSITY/FLUX:'  'Flux'\n")
        fout.write("'COS(THETA):'  '1.00'\n")
        fout.write("'ABFIND:'  '.false'\n")
        fout.write("'MODELINPUT:'  '{:s}'\n".format(os.path.basename(atmod)))
        if kurucz : fout.write("'MARCS-FILE:'  '.false.'\n")
        fout.write("'MODELOPAC:'  '{:s}'\n".format(os.path.basename(root)+'opac'))
        fout.write("'RESULTFILE:'  '{:s}'\n".format(os.path.basename(file)))
        fout.write("'METALLICITY:'  '{:8.3f}\n".format(mh))
        fout.write("'ALPHA/Fe:'  '{:8.3f}'\n".format(am))
        fout.write("'HELIUM:'  '{:8.3f}'\n".format(0.00))
        fout.write("'R-PROCESS:'  '{:8.3f}'\n".format(0.00))
        fout.write("'S-PROCESS:'  '{:8.3f}'\n".format(0.00))
        fout.write("'INDIVIDUAL ABUNDANCES:'  '{:2d}'\n".format(nels))
        fout.write("    6  {:8.3f}\n".format(8.39+mh+cm))
        fout.write("    7  {:8.3f}\n".format(7.78+mh+nm))
        fout.write("    8  {:8.3f}\n".format(8.66+mh+am))
        if elemgrid != '' : fout.write("{:6d}  {:8.3f}\n".format(elemnum,abun+mh))
        if els is not None :
            for el in els :
              num=atomic.periodic(el[0])[0]
              ab=atomic.solar(el[0])[0]+el[1]
              fout.write("{:5d}  {:8.3f}\n".format(num,ab+mh))
        if  not solarisotopes :
          fout.write("'ISOTOPES:'  '2'\n")
          # adopt ratio of 12C/13C=15
          fout.write("   6.012 0.9375\n")
          fout.write("   6.013 0.0625\n")
        # do we need to add the H2O linelist?
        if h2o is None : 
          if teff < 4000 :
            if mh+am < -1.5 or teff > 3250 :
              h2o=1    
            else  :
              h2o=2
          else :
              h2o=0
        nlists=3
        # if no HI lines, don't use that list: it takes a while to read
        n_HI = len(open(linelistdir+'/turbospec.'+linelist+'.Hlinedata').readlines())
        if n_HI < 3 : nlists-=1
        # if we are using H2O, add that list
        if h2o > 0 : nlists+=1
        fout.write("'NFILES:'  '{:4d}'\n".format(nlists))
        if n_HI >= 3 : fout.write(linelistdir+'/turbospec.'+linelist+'.Hlinedata\n')
        fout.write(linelistdir+'/turbospec.'+linelist+'.atoms\n')
        fout.write(linelistdir+'/turbospec.'+linelist+'.molec\n')
        if h2o == 1 :
            fout.write(linelistdir+'/turbospec.h2o-BC8.5V'+'.molec\n')
        elif h2o == 2 :
            fout.write(linelistdir+'/turbospec.h2o-BC9.5V'+'.molec\n')
        if geo == 's' :
            fout.write("'SPHERICAL:'  'T'\n")
        else :
            fout.write("'SPHERICAL:'  'F'\n")
        fout.write("30\n")
        fout.write("300.00\n")
        fout.write("15\n")
        fout.write("1.3\n")
        fout.close()

        # control file, with special handling in case bsyn goes into infinite loop ...
        fout = open(root+"_bsyn.csh",'w')
        fout.write("#!/bin/csh -f\n")
        fout.write("{:s}/bin/bsyn_lu < {:s} &\n".format(os.environ['APOGEE_DIR'],os.path.basename(bsynfile)))
        fout.write('set bsynjob = $!\n')
        fout.write("set ok = 0\n")
        fout.write("set runtime = `ps -p $bsynjob -o cputime | tail -1 | awk -F: '{print ($1*3600)+($2*60)+$3}'`\n")
        tmax=120*int(.05/min([0.05,dw]))
        if h2o == 1 : tmax*=2
        if h2o == 2 : tmax*=5
        fout.write('while ( $runtime < {:d} && $ok == 0 )\n'.format(tmax))
        fout.write('  usleep 200000\n')
        fout.write("  set runtime = `ps -p $bsynjob -o cputime | tail -1 | awk -F: '{print ($1*3600)+($2*60)+$3}'`\n")
        fout.write('  if ( `ps -p $bsynjob -o comm=` == "" ) then\n')
        fout.write('    echo process done, exiting!\n')
        fout.write('    set ok = 1\n')
        fout.write('  endif\n')
        fout.write('end\n')
        fout.write('if ( $ok == 0 ) then\n')
        fout.write('  echo expired, killing job\n')
        fout.write('  kill $bsynjob\n')
        fout.write('endif\n')
        fout.close()
        if run :
            os.chmod(root+'_bsyn.csh', 0o777)
            os.chdir(workdir)
            subprocess.call(['time','./'+os.path.basename(root)+'_bsyn.csh'],stdout=stdout)
            os.chdir(cwd)
            try:
                out=np.loadtxt(file)
                if ielem == 0 : 
                    spec=out[:,fluxcol]
                else :
                    spec=np.vstack([spec,out[:,fluxcol]])
            except :
                print('failed...',file)
                return 0.

    if run :
        if not save : shutil.rmtree(workdir)
        return spec

def prange(start,delta,n) :
    return float(start)+np.arange(int(n))*float(delta)

def get_vmicro(vmicrofit,vmicro) :
    if vmicrofit == 0 :
        return float(vmicro)
    else :
        print('need to implement vmicrofit: ', vmicrofit)
        pdb.set_trace()
    return  float(vmicro)


def pca(planfile,dir='kurucz/giantisotopes/tgGK_150714_lsfcombo5',pcas=None,whiten=False,plot=False,writeraw=False,test=False, fz=False, incremental=False) :
    """ Read in grid of spectra
    """

    showtime('start:')
    # input directory 
    if dir is None : 
        dir='kurucz/giantisotopes/tgGK_150714_lsfcombo5'
    indir=os.environ['APOGEE_SPECLIB']+'/synth/turbospec/'+dir+'/'
    print('indir: ', indir)

    # Read planfile and set output file name
    if not os.path.isfile(indir+'/plan/'+planfile): 
        print('{:s} does not exist'.format(indir+'/plan/'+planfile))
        return
    p=yanny.yanny(indir+'/plan/'+planfile,np=True)
    outfile=os.path.basename(p['name'])
    if test :
        p['nvt'] = '1'
        p['nam'] = '1'
        p['ncm'] = '1'
        p['nnm'] = '1'

#   uncompress if needed
    if fz :
      for iam,am in enumerate(prange(p['am0'],p['dam'],p['nam'])) :
        for icm,cm in enumerate(prange(p['cm0'],p['dcm'],p['ncm'])) :
          for inm,nm in enumerate(prange(p['nm0'],p['dnm'],p['nnm'])) :
           for ivm,vm in enumerate(prange(p['vt0'],p['dvt'],p['nvt'])) :
            file=('a{:s}c{:s}n{:s}v{:s}.fits').format(
                   atmos.cval(am),atmos.cval(cm),atmos.cval(nm),atmos.cval(10**vm))
            subprocess.call(['funpack',indir+file+'.fz'])

    # Read reference spectrum to plot and to determine number of pixel wavelengths
    refhead=fits.open(indir+'ap00cp00np00vp20.fits')[1].header
    wave=10.**vector(refhead,1)
    ref=fits.open(indir+'ap00cp00np00vp20.fits')[1].data[10,5,0,:]
    for ichip in range(2,4) :
        refhead=fits.open(indir+'ap00cp00np00vp20.fits')[ichip].header
        wave=np.append(wave,10.**vector(refhead,1))
        ref=np.append(ref,fits.open(indir+'ap00cp00np00vp20.fits')[ichip].data[10,5,0,:])
    nwave=ref.shape[0]

    # loop over requested combinations of npieces and npca
    if pcas is None : pcas = (int(p['npart']),int(p['npca']))
    npiece,npca = pcas

    showtime('start config:')
    # determine number of pixels per piece
    nspec=int(np.ceil(nwave/npiece))
    print(npiece,npca,nspec)

    # initialize PCA object, output figure and file
    if incremental :
        pca = IncrementalPCA(n_components=npca,whiten=whiten,batch_size=1e5)
    else :
        pca = PCA(n_components=npca,whiten=whiten)
    if plot : fig,ax=plots.multi(1,3,hspace=0.001,wspace=0.001,figsize=(12,4))
    fout=open(outfile+'_{:03d}_{:03d}.txt'.format(npiece,npca),'w')

    # initialize eigenvector array
    eigen = np.zeros([npca,nwave])
    mean = np.zeros([nwave])

    # loop over pieces
    npixels = []
    for ipiece in range(npiece) :
      w1=ipiece*nspec
      w2=(ipiece+1)*nspec if ipiece < npiece -1 else nwave
      npix=w2-w1
      npixels.append(npix)
      print(ipiece,w1,w2)

      # load data for this piece
      pcadata=np.zeros([int(p['nam'])*int(p['ncm'])*int(p['nnm'])*int(p['nvt'])*int(p['nmh'])*int(p['nlogg'])*int(p['nteff']),npix],dtype=np.float32)
      t0=time.time()
      showtime('start piece:')
      nmod=0
      for iam,am in enumerate(prange(p['am0'],p['dam'],p['nam'])) :
        for icm,cm in enumerate(prange(p['cm0'],p['dcm'],p['ncm'])) :
          for inm,nm in enumerate(prange(p['nm0'],p['dnm'],p['nnm'])) :
           for ivm,vm in enumerate(prange(p['vt0'],p['dvt'],p['nvt'])) :
            file=('a{:s}c{:s}n{:s}v{:s}.fits').format(
                   atmos.cval(am),atmos.cval(cm),atmos.cval(nm),atmos.cval(10**vm))
            s1=fits.open(indir+file)[1].data
            s2=fits.open(indir+file)[2].data
            s3=fits.open(indir+file)[3].data
            for imh,mh in enumerate(prange(p['mh0'],p['dmh'],p['nmh'])) :
              for ilogg,logg in enumerate(prange(p['logg0'],p['dlogg'],p['nlogg'])) :
                for iteff,teff in enumerate(prange(p['teff0'],p['dteff'],p['nteff'])) :
                  s=np.append(s1[imh,ilogg,iteff,:],s2[imh,ilogg,iteff,:])
                  s=np.append(s,s3[imh,ilogg,iteff,:])
                  pcadata[nmod,:] = s[w1:w2]
                  nmod+=1
      del s1, s2, s3, s

      # do the PCA decomposition 
      print(pcadata.shape)
      t1=time.time()
      showtime('start pca:')
      model=pca.fit_transform(pcadata)
      print(pca.explained_variance_ratio_)
      eigen[:,w1:w2] = pca.components_
      mean[w1:w2] = pca.mean_
      t2=time.time()
      showtime('start inverse:')
      # do the PCA reconstruction
      fit=pca.inverse_transform(model)
      t3=time.time()
      print(t1-t0,t2-t1,t3-t2)
      fout.write('{:8.2f}{:8.2f}{:8.2f}\n'.format(t1-t0,t2-t1,t3-t2))
      rat=pcadata/fit
      # plot results
      if plot :
          showtime('start plot:')
          #for j in range(0,nmod,1000) :
          #    plots.plotl(ax[0],wave[w1:w2],rat[j,:],xr=[wave[0],wave[-1]],yr=[0.7,1.3])
          ax[0].text(wave[w1]+0.1*(wave[w2-1]-wave[w1]),1.3,'{:.2f}'.format(rat.min()),va='top',ha='left')
          ax[0].text(wave[w2-1]+0.1*(wave[w2-1]-wave[w1]),1.3,'{:.2f}'.format(rat.max()),va='top',ha='right')
          plots.plotl(ax[1],wave[w1:w2],ref[w1:w2],xr=[wave[0],wave[-1]],color=colors[ipiece%6],yr=[0.7,1.3])
          hist,bins=np.histogram(rat.flatten(),bins=np.arange(0.5,1.51,0.01),density=True)
          plots.plotl(ax[2],np.arange(0.5+0.005,1.5,0.01),hist,color=colors[ipiece%6],semilogy=True)
          #plt.draw()
          #plt.show()

      # write out uncompressed and compressed in binary to local file
      showtime('start write:')
      if writeraw: fraw=open(os.environ['APOGEE_LOCALDIR']+'/'+outfile+'_{:03d}_{:03d}_{:03d}.raw'.format(npiece,npca,ipiece),'wb')
      fpca=open(os.environ['APOGEE_LOCALDIR']+'/'+outfile+'_{:03d}_{:03d}_{:03d}.pca'.format(npiece,npca,ipiece),'wb')
      for i in range(nmod) :
          if writeraw: fraw.write(struct.pack('f'*npix,*pcadata[i,:]))
          fpca.write(struct.pack('f'*npca,*model[i,:]))
      if writeraw: fraw.close()
      fpca.close()
      showtime('done piece:')

    # save plot and close CPU time file
    if plot :
        fig.savefig(outfile+'_{:03d}_{:03d}.png'.format(npiece,npca))
        plt.close()
    fout.close()
    del pca

#   uncompress if needed
    if fz :
      for iam,am in enumerate(prange(p['am0'],p['dam'],p['nam'])) :
        for icm,cm in enumerate(prange(p['cm0'],p['dcm'],p['ncm'])) :
          for inm,nm in enumerate(prange(p['nm0'],p['dnm'],p['nnm'])) :
           for ivm,vm in enumerate(prange(p['vt0'],p['dvt'],p['nvt'])) :
            file=('a{:s}c{:s}n{:s}v{:s}.fits').format(
                   atmos.cval(am),atmos.cval(cm),atmos.cval(nm),atmos.cval(10**vm))
            os.remove(indir+file)

    # header file for PCA
    ferre.wrhead(p,'p_aps'+outfile+'_{:03d}_{:03d}.hdr'.format(npiece,npca),npca=npixels,npix=npiece*npca)
    # output eigenvectors
    fp = open('p_aps'+outfile+'_{:03d}_{:03d}.hdr'.format(npiece,npca),'a')
    out=np.append(mean.reshape((1,nwave)),mean.reshape((1,nwave))*0.,axis=0)
    out=np.append(out,eigen,axis=0)
    np.savetxt(fp,out)
    
    # bundle the files into a single file
    allpca=open('p_aps'+outfile+'_{:03d}_{:03d}.unf'.format(npiece,npca),'wb')
    fpca=[]
    if writeraw: 
        ferre.wrhead(p,'f_aps'+outfile+'.hdr',npix=nwave)
        allraw=open('f_aps'+outfile+'.unf','wb')
        fraw=[]
    for ipiece in range(npiece) :
        if writeraw: fraw.append(open(os.environ['APOGEE_LOCALDIR']+'/'+outfile+'_{:03d}_{:03d}_{:03d}.raw'.format(npiece,npca,ipiece),'rb'))
        fpca.append(open(os.environ['APOGEE_LOCALDIR']+'/'+outfile+'_{:03d}_{:03d}_{:03d}.pca'.format(npiece,npca,ipiece),'rb'))
    for i in range(nmod) :
        for ipiece in range(npiece) :
            w1=ipiece*nspec
            w2=(ipiece+1)*nspec if ipiece < npiece -1 else nwave
            npix=w2-w1
            pca=fpca[ipiece].read(npca*4)
            allpca.write(pca)
            if writeraw: 
                raw=fraw[ipiece].read(npix*4)
                allraw.write(raw)
    allpca.close()
    if writeraw: allraw.close()
    for ipiece in range(npiece) :
        fpca[ipiece].close()
        os.remove(os.environ['APOGEE_LOCALDIR']+'/'+outfile+'_{:03d}_{:03d}_{:03d}.pca'.format(npiece,npca,ipiece))
        if writeraw: 
            fraw[ipiece].close()
            os.remove(os.environ['APOGEE_LOCALDIR']+'/'+outfile+'_{:03d}_{:03d}_{:03d}.raw'.format(npiece,npca,ipiece))

def mkgrid(planfile,clobber=False,resmooth=False,renorm=False,save=False,run=True,split=None,highres=9) :
    """ Create a grid of synthetic spectra 
    """

    # Read planfile
    if not os.path.isfile(planfile): 
        print('{:s} does not exist'.format(planfile))
        return
    p=yanny.yanny(planfile,np=True)

    # header information
    wrange=[float(x) for x in p['wrange'].split()]
    dw=float(p['dw'])
    vacuum = int(p['vacuum']) if p.get('vacuum') else 0
    kurucz = True if p['atmos'] == 'kurucz' else False
    marcsdir = p['marcsdir'] if p.get('marcsdir') else None
    solarisotopes = p['solarisotopes'] if p.get('solarisotopes') else 0
    elem = p['elem'] if p.get('elem') else ''
    maskdir = p['maskdir'] if p.get('maskdir') else None
    vmicrofit = p['vmicrofit'] if p.get('vmicrofit') else 0
    vmicro = p['vmicro'] if p.get('vmicro') else 0
    vmacrofit = p['vmacrofit'] if p.get('vmacrofit') else 0
    vmacro = p['vmacro'] if p.get('vmacro') else 0
    specdir = os.environ['APOGEE_SPECLIB']+'/synth/'+p['specdir'] if p.get('specdir') else './'
    linelistdir=os.environ['APOGEE_SPECLIB']+'/linelists/' 
    linelist = p['linelist'] if p.get('linelist') else None

    # wavelength array
    nspec=int((wrange[1]-wrange[0])/dw)+1
    rawwave=wrange[0]+np.arange(nspec)*dw

    # if element minigrid, create mini linelist and get wavelengths to store
    if elem == '' :
        nelem=1
        gd=range(nspec)
    else :
        nelem=8
        wvac = mini_linelist(elem,linelist,maskdir)
        nwind=wvac.shape[0]
        gd=[]
        for iwind in range(nwind) :
            gd.extend(np.where( (rawwave >= wvac[iwind,0]) & (rawwave <= wvac[iwind,1]) )[0])
        nspec=len(gd)

    # make the grid(s)
    for am in prange(p['am0'],p['dam'],p['nam']) :
      for cm in prange(p['cm0'],p['dcm'],p['ncm']) :
        for nm in prange(p['nm0'],p['dnm'],p['nnm']) :
          specdata=np.zeros([nelem,int(p['nmh']),int(p['nlogg']),int(p['nteff']),nspec],dtype=np.float32)
          for imh,mh in enumerate(prange(p['mh0'],p['dmh'],p['nmh'])) :
            for ilogg,logg in enumerate(prange(p['logg0'],p['dlogg'],p['nlogg'])) :
              for iteff,teff in enumerate(prange(p['teff0'],p['dteff'],p['nteff'])) :

                print(teff, logg, mh)
                sys.stdout.flush()
                nskip=0 
                dskip = 1 if kurucz else 2
                vout = get_vmicro(vmicrofit,vmicro)
                while nskip >= 0 and nskip < 10 :
                  spec=mkturbospec(int(teff),logg,mh,am,cm,nm,
                    wrange=wrange,dw=dw,atmosdir=marcsdir,
                    elemgrid=elem,linelistdir=linelistdir+'/'+elem+'/',linelist=linelist,vmicro=vout,
                    solarisotopes=solarisotopes,
                    nskip=nskip,kurucz=kurucz,run=run,save=save,split=split) 
                  nskip = nskip+dskip if isinstance(spec,float) else -1
                if elem == '' :
                    specdata[0,imh,ilogg,iteff,:]=spec
                else :
                    specdata[:,imh,ilogg,iteff,:]=spec[:,gd]

          # FITS header and output
          hdu=fits.PrimaryHDU(np.squeeze(specdata))
          idim=1
          if elem == '' :
              add_dim(hdu.header,rawwave[0],rawwave[1]-rawwave[0],1,'WAVELENGTH',idim)
          else :
              hdu.header.append('CDELT1',rawwave[1]-rawwave[0])
              for iwind in range(nwind) :
                  hdu.header.append(('WIND0_{:d}'.format(iwind),wvac[iwind,0]))
                  hdu.header.append(('WIND1_{:d}'.format(iwind),wvac[iwind,1]))
          if int(p['nteff']) > 1 :
              idim+=1
              add_dim(hdu.header,float(p['teff0']),float(p['dteff']),1,'TEFF',idim)
          if int(p['nlogg']) > 1 :
              idim+=1
              add_dim(hdu.header,float(p['logg0']),float(p['dlogg']),1,'LOGG',idim)
          if int(p['nmh']) > 1 :
              idim+=1
              add_dim(hdu.header,float(p['mh0']),float(p['dmh']),1,'M_H',idim)
          hdu.header['LOGW'] = 0
          if p.get('width') : hdu.header['width'] = p['width']
          if p.get('linelist') : hdu.header['linelist'] = p['linelist']
          if p['synthcode'] == 'asset'  : hdu.header.add_comment('ASSET generated synthetic spectra')
          if p['synthcode'] == 'turbospec' : hdu.header.add_comment('Turbospec generated synthetic spectra')
          if p['synthcode'] == 'moog ' : hdu.header.add_comment('MOOG generated synthetic spectra')
          try : os.mkdir(specdir)
          except: pass
          hdu.writeto(specdir+'/'+p['name']+'.fits',overwrite=True)

def add_dim(header,crval,cdelt,crpix,ctype,idim) :
    """ Add a set of CRVAL/CDELT,CRPIX,CTYPE cards to header
    """

    header.append(('CRVAL{:d}'.format(idim),crval))
    header.append(('CDELT{:d}'.format(idim),cdelt))
    header.append(('CRPIX{:d}'.format(idim),crpix))
    header.append(('CTYPE{:d}'.format(idim),ctype))


def mkgriddirs(configfile) :

    """ Script to create output directories and plan and batch queue files for all grids listed in master grid configuration file
    """

    # Read grid configuration file
    if not os.path.isfile(configfile): 
        print('{:s} does not exist'.format(configfile))
        return
    p=yanny.yanny(configfile,np=True)

    # loop over each grid
    for i in range(len(p['GRID']['specdir'])) :

        # construct name and create output directory
        name = p['GRID']['specdir'][i]+'_'+p['GRID']['smooth'][i]

        if abs(p['GRID']['solarisotopes'][i]) == 1 :
            iso = 'solarisotopes'
        else :
            iso = 'giantisotopes'
        if p['GRID']['solarisotopes'][i] < 0 :
            iso = 'tests/'+iso
        dir = os.getenv('APOGEE_SPECLIB')+'/synth/'+p['synthcode'].strip("'")+'/'+p['GRID']['atmos'][i]+'/'+iso+'/'+name+'/plan/'
        print(dir)
        try: os.makedirs(dir)
        except: pass

        # remove any old plan files
        os.chdir(dir)
        for filePath in glob.glob("*.par"):
            if os.path.isfile(filePath): os.remove(filePath)

        # write the master planfile and one for each minigrid
        elems=['']
        elems.extend(p['GRID']['elem'][i])
        f = open(dir+name+'.par','w')
        f.write('{:20s}{:20s}\n'.format('name', name))
        for key in p.keys() :
            if ( key != 'GRID' ) & (key != 'symbols') :
                f.write('{:30s}{:30s}\n'.format(key, p[key]))
        for key in p['GRID'].dtype.names :
            out='{:30s}'.format(str(p['GRID'][key][i])).strip("'").strip('[]')
            out='{:30s}{:30s}\n'.format(key,out.replace(']',''))
            f.write(out.strip('[]'))
        f.close()

        # make all of the individual planfiles from the master planfile
        subprocess.call(['idl','-e',"speclib_allplan,'"+name+".par'"])

        # make pbs scripts
        os.chdir('..')
        specdir = p['synthcode'].strip("'")+'/'+p['GRID']['atmos'][i]+'/'+iso+'/'+name
        os.environ['NO_NODES'] = 'yes'
        subprocess.call(['mkslurm','mkgrid','"plan/'+name+'_a[mp]*vp20.par"','"plan/'+name+'_a[mp]*vp??.par"'])
        subprocess.call(['mkslurm','mkgridlsf','"plan/'+name+'_a[mp]*vp??.par"'])
        subprocess.call(['mkslurm','bundle','"plan/'+name+'_??.par"'])

def mkspec(pars) :
    """ Makes a single spectrum given input pars
    """
    teff=pars[0].astype('int')
    logg=pars[1]
    mh=pars[2]
    vmicro=pars[3]
    vrot=pars[4]
    cm=round(pars[5]/0.25)*0.25
    nm=round(pars[6]/0.5)*0.5
    am=round(pars[7]/0.25)*0.25
    elems=[]
    els = ['C','N','O','Na','Mg','Al','Si','P','S','K','Ca','Ti','V','Cr','Mn','Co','Ni','Cu','Ge','Rb','Ce','Nd']
    for j,el in enumerate(els) :
        elems.append([el,pars[5+j]])
    print(teff,logg,mh,vmicro,am,cm,nm)
    spec=mkturbospec(teff,logg,mh,am,cm,nm,vmicro=vmicro,els=elems,kurucz=False)
    return pars,spec
    

def mksynth(file,wrange=[15100,17000],threads=8) :
    """ Make a series of spectra from parameters in an input file
    """
    pars=np.loadtxt(file)
    out=[]
    outpar=[]

    pool = mp.Pool(threads)
    specs = pool.map_async(mkspec, pars).get()
    pool.close()
    pool.join()

    for spec in specs :
        if isinstance(spec[1],np.ndarray) :
            out.append(spec[1])
            outpar.append(spec[0])

    # write the spectra out
    hdu=fits.HDUList()
    hdu.append(fits.ImageHDU(out))
    hdu.append(fits.ImageHDU(outpar))
    hdu.writeto('synth.fits',overwrite=True)
   
def filter_lines(infile,outfile,wind,nskip=0) :
    """ Read from infile, output comments and lines falling in windows of [w1,w2] to outfile
    """
    fin=open(infile,'r')
    fout=open(outfile,'w')
    nout = 0
    nwind=wind.shape[0]
    for iline,line in enumerate(fin) :
        if iline >= nskip : 
            w=line.split()[0]
            if w[0] == '#' :
                    fout.write(line)
            else :
              for i in range(nwind) :
                if (float(w) >= wind[i,0]) and (float(w) <=wind[i,1]) :
                    fout.write(line)
                    nout+=1
    fout.close()
    return nout
 
def mini_linelist(elem,linelist,maskdir) :
    """ Produce an abbreviated line list for minigrid construction given mask file and lineist file IN AIR
        Return array of vacuum wavelength ranges
    """

    wind=np.loadtxt(os.environ['APOGEE_DIR']+'/data/windows/'+maskdir+'/'+elem+'.wave')
    wair=spectra.vactoair(wind)
   
    outdir = os.environ['APOGEE_SPECLIB']+'/linelists/'+elem+'/'
    try: os.mkdir(outdir)
    except: pass
    nout=filter_lines(os.environ['APOGEE_SPECLIB']+'/linelists/linelist.'+linelist,outdir+linelist,wair/10.)
    subprocess.call(['turboscript',outdir+linelist])

    lists=['turbospec.20170418.Hlinedata','turbospec.h2o-BC8.5V.molec','turbospec.h2o-BC9.5V.molec']
    code=['01.000000','010108.000000000','010108.00000000']
    comment=['HI culled','Barber culled','Barber culled']
    for i,list in enumerate(lists) :
        nout=filter_lines(os.environ['APOGEE_SPECLIB']+'/linelists/'+list,outdir+list+'.tmp',wair,nskip=2)
        fin=open(outdir+list+'.tmp','r')
        fout=open(outdir+list,'w')
        fout.write("'"+code[i]+" '  1 "+'{:d}\n'.format(nout))
        fout.write("'"+comment[i]+"' \n")
        for line in fin :
            fout.write(line)
        fout.close()
        fin.close()
        os.remove(outdir+list+'.tmp')
    return wind
