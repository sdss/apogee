# encoding: utf-8
#
# @Author: Jon Holtzman
# @Date: March 2018
# @Filename: pca.py
# @License: BSD 3-Clause
# @Copyright: Jon Holtzman

# routines for PCA compression of synthetic spectral grids, and
# output into FERRE style library files. Also test routine
# for calculating and comparing results from PCA library

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
from apogee.aspcap import aspcap
from apogee.aspcap import ferre
from apogee.speclib import atmos
from apogee.speclib import lsf
from apogee.speclib import sample
from apogee.utils import atomic
from apogee.utils import spectra
from apogee.plan import mkslurm
#from sdss.utilities import yanny
from sdss import yanny
from astropy.io import ascii
from astropy.io import fits
from sklearn.decomposition import PCA
from sklearn.decomposition import IncrementalPCA
import matplotlib.pyplot as plt
from tools import plots
from tools import match
from tools import html

colors=['r','g','b','c','m','y']

def pca(planfile,dir='kurucz/giantisotopes/tgGK_150714_lsfcombo5',pcas=None,whiten=False,writeraw=False,test=False, 
        incremental=False, threads=4, rawsynth=False, prefix='',piece=None) :
    """ Read in grid of spectra and do PCA compression

    Args :
        planfile (str) : input file name with grid parameteris
        dir (str) :  input directory relative to #$APOGEE_SPECLIB/synth/turbospec (default='kurucz/giantisotopes/tgGK_150714_lsfcombo5')
        pcas (2-tuple) : input number of PCA components and number of pieces (default=None --> get from input planfile)
        pieces (list) : process only specified pieces, useful if want to split across processors, use pieces=[] to bundle existing output (default=None --> all pieces)
        whiten (bool) : pre-whiten data? (default=False)
        writeraw (bool ) : output uncompressed grid in FERRE format ? (default=False)
        test (bool) : used abreviated grid for testing (speed) (default=False)
        incremental (bool) : use incremental PCA routline (default=False) (WARNING: if not incremental, default PCA is not perfectly repeatable!)
        threads (int) : number of threads to use for parallel calculation (default=4), large grids may require threads=0
        rawsynth (bool) : work on raw highres synthesis output (default=False), UNTESTED??

    Output: compressed grid in FERRE format
    """

    # Read planfile and set output file name
    if not os.path.isfile(planfile):
        print('{:s} does not exist'.format(planfile))
        return
    p=yanny.yanny(planfile,np=True)

    if dir is None :
        if int(p['solarisotopes']) == 1 : isodir = 'solarisotopes'
        else : isodir = 'giantisotopes'
        dir = p['atmos'] + '/' + isodir + '/' + p['name']+'/' if p.get('name') else './'

    showtime('start:')
    # input directory 
    indir=os.environ['APOGEE_SPECLIB']+'/synth/turbospec/'+dir+'/'
    print('indir: ', indir)

    outfile=os.path.basename(p['name'])
    if test :
        p['nvt'] = '1'
        p['nam'] = '1'
        p['ncm'] = '1'
        p['nnm'] = '1'
        outfile = 'test'+outfile

    nmod = int(p['nvt'])*int(p['ncm'])*int(p['nnm'])*int(p['nam'])*int(p['nrot'])*int(p['nmh'])*int(p['nlogg'])*int(p['nteff'])

    # loop over requested combinations of npieces and npca
    if pcas is None : pcas = (int(p['npart']),int(p['npca']))
    npiece,npca = pcas

    indata={}
    indata['npiece'] = npiece
    indata['npca'] = npca
    indata['whiten'] = whiten
    indata['indir'] = indir
    indata['outfile'] = outfile
    indata['writeraw'] = writeraw
    indata['rawsynth'] = rawsynth
    indata['prefix'] = prefix
    indata['incremental'] = incremental
    for key in ['am0','dam','nam','cm0','dcm','ncm','nm0','dnm','nnm','vt0','dvt','nvt','mh0','dmh','nmh','logg0','dlogg','nlogg','teff0','dteff','nteff','rot0','drot','nrot'] :
        indata[key] = p[key]

    # determine number of pixels per piece
    if rawsynth :
        wave=np.arange(15100.,17000.01,0.05)
        nwave=len(wave)
    else :
        nwave=aspcap.nw_chip.sum()
    nspec=int(np.ceil(nwave/npiece))
    print(npiece,npca,nspec,nwave)

    # loop over pieces
    pars = []
    if piece == None : 
        piece=range(npiece)
        outdir = os.environ['APOGEE_LOCALDIR']
    else :
        outdir = './'
    indata['outdir'] = outdir
    for ipiece in piece :
        w1=ipiece*nspec
        w2=(ipiece+1)*nspec if ipiece < npiece -1 else nwave
        npix=w2-w1
        print(ipiece,w1,w2)
        pars.append((ipiece,indata,[w1,w2]))

    if threads == 0 :
        outputs=[]
        for par in pars :
            wrange=par[2]
            outputs.append(dopca(par))
    else :
        # do the PCA in parallel for the subgrids
        pool = mp.Pool(threads)
        outputs = pool.map_async(dopca, pars).get()
        pool.close()
        pool.join()

    # header file informationfor output files
    wchip=[ [aspcap.nw_chip[0],aspcap.logw0_chip[0],aspcap.dlogw], 
            [aspcap.nw_chip[1],aspcap.logw0_chip[1],aspcap.dlogw], 
            [aspcap.nw_chip[2],aspcap.logw0_chip[2],aspcap.dlogw] ]
    cont=[(1,1,0.,0.),(1,1,0.,0.),(1,1,0.,0)]

    # if we've just done some of the pieces (e.g. on different nodes), stop here
    if len(piece) != npiece and len(piece) != 0 : return

    # put together output eigenvectors and means (components are written to disk in pieces)
    if npca > 0 :
        eigen = np.zeros([npca,nwave])
        mean = np.zeros([nwave])
        npixels = []
        for ipiece in range(npiece) :
            w1=ipiece*nspec
            w2=(ipiece+1)*nspec if ipiece < npiece -1 else nwave
            npix=w2-w1
            npixels.append(npix)
            if npca > 0 : 
                feigen=open(outdir+'/'+outfile+'_{:03d}_{:03d}_{:03d}.eigen'.format(npiece,npca,ipiece),'rb')
                for i in range(npca) : eigen[i,w1:w2] = struct.unpack('f'*npix,feigen.read(npix*4))
                mean[w1:w2] = struct.unpack('f'*npix,feigen.read(npix*4))
                feigen.close()
        #for ipiece,output in enumerate(outputs) :
        #    wrange=pars[ipiece][2]
        #    eigen[:,wrange[0]:wrange[1]] = output[0]
        #    mean[wrange[0]:wrange[1]] = output[1]

        ferre.wrhead(p,'p_aps'+outfile+'_{:03d}_{:03d}.hdr'.format(npiece,npca),npca=npixels,npix=npiece*npca,wchip=wchip,cont=cont)
        # output eigenvectors
        fp = open('p_aps'+outfile+'_{:03d}_{:03d}.hdr'.format(npiece,npca),'a')
        out=np.append(mean.reshape((1,nwave)),mean.reshape((1,nwave))*0.,axis=0)
        out=np.append(out,eigen,axis=0)
        np.savetxt(fp,out)
        fp.close()
    
        # bundle the component files into a single file
        allpca=open('p_aps'+outfile+'_{:03d}_{:03d}.unf'.format(npiece,npca),'wb')

    if npca > 0 : fpca=[]
    if writeraw: 
        ferre.wrhead(p,'f_aps'+outfile+'.hdr',npix=nwave,wchip=wchip,cont=cont)
        allraw=open('f_aps'+outfile+'.unf','wb')
        fraw=[]
    for ipiece in range(npiece) :
        if writeraw: fraw.append(open(outdir+'/'+outfile+'_{:03d}_{:03d}_{:03d}.raw'.format(npiece,npca,ipiece),'rb'))
        if npca > 0 : fpca.append(open(outdir+'/'+outfile+'_{:03d}_{:03d}_{:03d}.pca'.format(npiece,npca,ipiece),'rb'))
    for i in range(nmod) :
        for ipiece in range(npiece) :
            w1=ipiece*nspec
            w2=(ipiece+1)*nspec if ipiece < npiece -1 else nwave
            npix=w2-w1
            if npca > 0 : pca=fpca[ipiece].read(npca*4)
            if npca > 0 : allpca.write(pca)
            if writeraw: 
                raw=fraw[ipiece].read(npix*4)
                allraw.write(raw)
    if npca > 0 : allpca.close()
    if writeraw: allraw.close()
    for ipiece in range(npiece) :
        if npca > 0 : 
            fpca[ipiece].close()
            os.remove(outdir+'/'+outfile+'_{:03d}_{:03d}_{:03d}.pca'.format(npiece,npca,ipiece))
        if writeraw: 
            fraw[ipiece].close()
            os.remove(outdir+'/'+outfile+'_{:03d}_{:03d}_{:03d}.raw'.format(npiece,npca,ipiece))

def dopca(pars) :
    """ do a single PCA decomposition for a limited number of pixels

    Args:
        pars : input tuple giving required input
    """
    ipiece = pars[0]
    showtime('start piece: '+str(ipiece))
    p=pars[1]
    wrange=pars[2]
    w1=wrange[0]
    w2=wrange[1]
    npix=w2-w1

    npiece=p['npiece']
    npca=p['npca']
    whiten=p['whiten']
    indir=p['indir']
    outfile=p['outfile']
    outdir=p['outdir']
    writeraw=p['writeraw']
    rawsynth=p['rawsynth']
    prefix=p['prefix']

    if p['incremental'] :
        print('using incremental PCA')
        pca = IncrementalPCA(n_components=npca,whiten=whiten,batch_size=1000)
    else :
        pca = PCA(n_components=npca,whiten=whiten)

    # load data for this piece
    pcadata=np.zeros([int(p['nam'])*int(p['ncm'])*int(p['nnm'])*int(p['nvt'])*
                      int(p['nrot'])*int(p['nmh'])*int(p['nlogg'])*int(p['nteff']),npix],dtype=np.float32)
    nmod=0
    pix_apstar=aspcap.gridPix()
    pix_aspcap=aspcap.gridPix(apStar=False)
    if rawsynth:
        wave=np.arange(15100.,17000.01,0.05)
        nwave=len(wave)
    else :
        nwave=aspcap.nw_chip.sum()

    for ivm,vm in enumerate(spectra.vector(p['vt0'],p['dvt'],p['nvt'])) :
      for icm,cm in enumerate(spectra.vector(p['cm0'],p['dcm'],p['ncm'])) :
        for inm,nm in enumerate(spectra.vector(p['nm0'],p['dnm'],p['nnm'])) :
          for iam,am in enumerate(spectra.vector(p['am0'],p['dam'],p['nam'])) :
            file=('a{:s}c{:s}n{:s}v{:s}.fits').format(
                   atmos.cval(am),atmos.cval(cm),atmos.cval(nm),atmos.cval(10**vm))
            # read file and pack into ASPCAP grid size
            sap=fits.open(indir+prefix+file)[0].data
            sap=sap.reshape((int(p['nrot']),int(p['nmh']),int(p['nlogg']),int(p['nteff']),sap.shape[-1]))
            s=np.zeros([int(p['nrot']),int(p['nmh']),int(p['nlogg']),int(p['nteff']),nwave])
            for pasp,pap in zip(pix_aspcap,pix_apstar) :
                s[:,:,:,:,pasp[0]:pasp[1]]=sap[:,:,:,:,pap[0]:pap[1]]
            #s1=fits.open(indir+file)[1].data
            #s2=fits.open(indir+file)[2].data
            #s3=fits.open(indir+file)[3].data
            for irot,rot in enumerate(spectra.vector(p['rot0'],p['drot'],p['nrot'])) :
             for imh,mh in enumerate(spectra.vector(p['mh0'],p['dmh'],p['nmh'])) :
              for ilogg,logg in enumerate(spectra.vector(p['logg0'],p['dlogg'],p['nlogg'])) :
                for iteff,teff in enumerate(spectra.vector(p['teff0'],p['dteff'],p['nteff'])) :
                  #s=np.append(s1[imh,ilogg,iteff,:],s2[imh,ilogg,iteff,:])
                  #s=np.append(s,s3[imh,ilogg,iteff,:])
                  #s=sall[imh,ilogg,iteff,:]
                  if s[irot,imh,ilogg,iteff,w1:w2].sum() == 0. or not np.isfinite(np.nanmean(s[irot,imh,ilogg,iteff,:]))  : 
                     print('!!ZERO or NaN MODEL',am,cm,nm,vm,rot,mh,logg,teff,w1,w2)
                     s[irot,imh,ilogg,iteff,w1:w2] = 1.
            #      s/=np.nanmean(s)
            #      if len(np.where(np.isfinite(s) is False)[0]) > 0 : print(imh,ilogg,iteff,np.where(np.isfinite(s) is False) )
                  pcadata[nmod,:] = s[irot,imh,ilogg,iteff,w1:w2]/np.nanmean(s[irot,imh,ilogg,iteff,:])
                  nmod+=1
      #del s1, s2, s3, s
    del s

    # do the PCA decomposition 
    if npca > 0 :
        print(pcadata.shape)
        showtime('start pca: '+str(ipiece))
        model=pca.fit_transform(pcadata)
        print(pca.explained_variance_ratio_)
        eigen = pca.components_
        mean = pca.mean_
    else :
        eigen = 0.
        mean = 0.

    # write out uncompressed and compressed in binary to local file
    showtime('start write: '+str(ipiece))
    if writeraw: fraw=open(outdir+'/'+outfile+'_{:03d}_{:03d}_{:03d}.raw'.format(npiece,npca,ipiece),'wb')
    if npca > 0 : 
        fpca=open(outdir+'/'+outfile+'_{:03d}_{:03d}_{:03d}.pca'.format(npiece,npca,ipiece),'wb')
        feigen=open(outdir+'/'+outfile+'_{:03d}_{:03d}_{:03d}.eigen'.format(npiece,npca,ipiece),'wb')
        for i in range(npca) : feigen.write(struct.pack('f'*npix,*eigen[i,:]))
        feigen.write(struct.pack('f'*npix,*mean))
        feigen.close()
    for i in range(nmod) :
        if writeraw: fraw.write(struct.pack('f'*npix,*pcadata[i,:]))
        if npca > 0 : fpca.write(struct.pack('f'*npca,*model[i,:]))
    if writeraw: fraw.close()
    if npca > 0 : fpca.close()
    showtime('done piece: '+str(ipiece))
    del pca 
    del pcadata
    return eigen,mean

def showtime(string) :
    """ Utiltiy routine to print a string and clock time
    """
    print(string+' {:8.2f}'.format(time.time()))
    sys.stdout.flush()

def liblink(lib,outdir) :
    """ Utility routine to create links in output directory to FERRE library files
       
    Args:
        lib (str) : root name of library
        outdir (str) : name of output directory to create links in
    Returns :
        prefix (str) : prefix giving relative directory of library files 
    """
    try : os.mkdir(outdir)
    except : pass
    if outdir[-1] != '/' : outdir=outdir+'/'
    nsub = len(outdir.split('/'))-1
    prefix=''
    for i in range(nsub) :
        prefix=prefix+'../'
    for suffix in ['.hdr','.unf'] :
        try: os.remove(outdir+'/'+lib+suffix)
        except: pass
        os.symlink(prefix+lib+suffix,outdir+'/'+lib+suffix)
    return prefix

def test(planfile,grid='GKg',npiece=12,npca=75,runraw=True,runpca=True,fit=True, fast=False, niso=None, sns=[1000], rot=False) :
    """ Routine to set up and run a series of tests for a PCA library
        Includes comparison of raw and PCA spectra for a sample suite,
        and several FERRE runs to recover parameters for the sample suite

        Requires existing test.ipf file with sample suite parameters

    Args :
        planfile (str) : name of input plan file that includes filename
        npiece (int) : number of PCA pieces (for directory/file names)
        npca (int) : number of PCA pieces (for directory/file names)
        runraw (bool) : create the sample and input FERRE spectra (default=True)
        runpca (bool) : create the PCA-derived spectra (default=True)
        fit (bool) : do the raw-PCA FERRE runs and comparison (default=True)
        fit (bool) : submit the FERRE test runs to queue (default=True)
        fast (bool) : use the sdss-fast queue (default=False)

    """
    # Read planfile and set output file name
    if not os.path.isfile(planfile):
        print('{:s} does not exist'.format(planfile))
        return
    p=yanny.yanny(planfile,np=True)
    outfile=os.path.basename(p['name'])

    # create test directory
    try : os.mkdir('test')
    except : pass

    # create test sample for this grid
    if runraw: sample.sample('test/test',gridclass=grid,niso=niso)

    # raw and PCA library file root names
    flib='f_aps'+outfile
    plib='p_aps'+outfile+('_{:03d}_{:03d}').format(npiece,npca)

    # make raw directory and setup to produce test spectra
    prefix=liblink(flib,'test/raw/')
    try: os.remove('test/raw/test.ipf')
    except: pass
    os.symlink(prefix+'test/test_'+grid+'.ipf','test/raw/test.ipf')
    l=ferre.rdlibhead('f_aps'+outfile+'.hdr')[0]
    ferre.writenml('test/raw/input.nml','test',l,nov=0,ncpus=1,f_access=1)
    mkslurm.write('ferre.x',outdir='test/raw/',runplans=False,cwd=os.getcwd()+'/test/raw',fast=fast)
    if runraw : 
        print('running ferre in raw to create spectra')
        subprocess.call(['test/raw/ferre.x'],shell=False)
        # create uncertainty spectra
        true=np.loadtxt('test/raw/test.mdl')
        ferre.writespec('test/raw/test.err',true*0.+0.001)
    else :
        true=np.loadtxt('test/raw/test.mdl')
    # add noise
    for sn in sns :
        if sn < 1000 :
            name = 'test/raw/testsn{:d}'.format(sn)
            if not os.path.isfile(name+'.mdl') :
                dim=true.shape
                truen=true.flatten()+np.random.normal(0.,1./sn,true.flatten().shape)
                truen=np.reshape(truen,(dim[0],dim[1]))
                ferre.writespec(name+'.mdl',truen)
                ferre.writespec(name+'.err',truen*0.+1./sn)
 
    # list of all the different tests and their input.nml files; latter are set below
    # first fdir is just for creating PCA
    root='test/pca_{:d}_{:d}'.format(npiece,npca)
    fdirs = [root+'/',root+'/algor3/',root+'/algor1/',root+'/algor1_12/',root+'/algor3_12/',
             root+'/algor3_001/',root+'/algor3_01/',root+'/algor3_1',
             root+'/algor3_12_01/',root+'/algor3_12_1',
             root+'/algor3_indi1234567/',root+'/algor3_indi3142567/',
             root+'/algor3_4d' ] 

    # general setup for all FERRE directories using PCA library
    l=ferre.rdlibhead('p_aps'+outfile+('_{:03d}_{:03d}.hdr').format(npiece,npca))[0]
    print("set up...")
    for fdir in fdirs :
        # create library links and return relative directory
        print(fdir)
        prefix=liblink(plib,fdir)
        # create input ipf
        cmds=[]
        for sn in sns :
            if sn < 1000 : name = fdir+'/testsn{:d}'.format(sn)
            else : name = fdir+'/test'
            # create input ipf
            try: os.remove(name+'.ipf')
            except: pass
            os.symlink(prefix+'test/test_'+grid+'.ipf',name+'.ipf')
            # create obs file with true spectra
            try: os.remove(name+'.obs')
            except: pass
            os.symlink(prefix+'test/raw/'+os.path.basename(name)+'.mdl',name+'.obs')
            # create uncertainty file
            try: os.remove(name+'.err')
            except: pass
            os.symlink(prefix+'test/raw/'+os.path.basename(name)+'.err',name+'.err')
            cmds.extend(['cp '+os.path.basename(name)+'.nml input.nml','ferre.x'])
        # create the command file to run FERRE
        mkslurm.write(cmds,name='ferre.x',outdir=fdir+'/',runplans=False,cwd=os.getcwd()+'/'+fdir,time='48:00:00',fast=fast)

    # test-specific input.nml files
    ferre.writenml(root+'/input.nml','test',l,nov=0,ncpus=1)
    for sn in sns :
        if sn < 1000 : name = 'testsn{:d}'.format(sn)
        else : name = 'test'
        ferre.writenml(root+'/algor3/'+name+'.nml',name,l,algor=3,renorm=4,obscont=1,ncpus=32,init=1)
        ferre.writenml(root+'/algor1_12/'+name+'.nml',name,l,algor=1,renorm=4,obscont=1,ncpus=32,indini=[1,1,1,1,2,2,3],init=1)
        ferre.writenml(root+'/algor1/'+name+'.nml',name,l,algor=1,renorm=4,obscont=1,ncpus=32,init=1)
        ferre.writenml(root+'/algor3_12/'+name+'.nml',name,l,algor=3,renorm=4,obscont=1,ncpus=32,indini=[1,1,1,1,2,2,3],init=1)
        ferre.writenml(root+'/algor3_001/'+name+'.nml',name,l,algor=3,renorm=4,obscont=1,ncpus=32,stopcr=0.001,init=1)
        ferre.writenml(root+'/algor3_01/'+name+'.nml',name,l,algor=3,renorm=4,obscont=1,ncpus=32,stopcr=0.01,init=1)
        ferre.writenml(root+'/algor3_12_01/'+name+'.nml',name,l,algor=3,renorm=4,obscont=1,ncpus=32,stopcr=0.01,indini=[1,1,1,1,2,2,3],init=1)
        ferre.writenml(root+'/algor3_1/'+name+'.nml',name,l,algor=3,renorm=4,obscont=1,ncpus=32,stopcr=0.1,init=1)
        ferre.writenml(root+'/algor3_12_1/'+name+'.nml',name,l,algor=3,renorm=4,obscont=1,ncpus=32,stopcr=0.1,indini=[1,1,1,1,2,2,3],init=1)
        ferre.writenml(root+'/algor3_indi1234567/'+name+'.nml',name,l,algor=3,renorm=4,obscont=1,ncpus=32,indi=[1,2,3,4,5,6,7],init=1)
        ferre.writenml(root+'/algor3_indi3142567/'+name+'.nml',name,l,algor=3,renorm=4,obscont=1,ncpus=32,indi=[3,1,4,2,5,6,7],init=1)
        ferre.writenml(root+'/algor3_4d/'+name+'.nml',name,l,algor=3,renorm=4,obscont=1,ncpus=32,indv=[4,5,6,7],init=1)

    # produce PCA version of test spectra
    if runpca : 
        print('running ferre in 12_75 to create spectra')
        subprocess.call([root+'/ferre.x'],shell=False)
    pca=np.loadtxt(root+'/test.mdl')
    # histogram of ratio of pca to true
    print("making pca/raw comparison histogram ...")
    fig,ax=plots.multi(1,2)
    hist,bins=np.histogram((pca/true).flatten(),bins=np.linspace(0.9,1.1,4001))
    plots.plotl(ax[0],np.linspace(0.8005,1.2,4000),hist/hist.sum(),semilogy=True,xt='pca/true')
    ax[1].hist(np.abs((pca-true).flatten()),bins=np.logspace(-7,3,50),histtype='step',normed=True,cumulative=True,color='k')
    ax[1].set_xlim(0.,0.01)
    ax[1].set_xlabel('|pca-true|')
    ax[1].set_ylabel('Cumulative fraction')
    fig.tight_layout()
    fig.savefig(root+'/'+outfile+'_pca.png')
    plt.close()

    # use the PCA libraries to fit raw spectra
    print("running/making plots...")
    tab=[]
    yt=[]
    # first fdir is just for creating PCA
    for fdir in fdirs[1:] :
        print(fdir)
        if fit : subprocess.call(['sbatch','ferre.x'],shell=False,cwd=fdir)
        xtab=[]
        xt=[]
        tmp=''
        for sn in sns :
            if sn < 1000 : name = fdir+'/testsn{:d}'.format(sn)
            else : name = fdir+'/test'
            try : 
                sample.comp(name,hard=True,true='test/raw/test.ipf',rot=rot)
                dt=subprocess.check_output("grep ellapsed "+fdir+"/ferre.x.out | tail -1 | awk '{print $3}'",shell=True)
            except : 
                print('failed comp: ',name)
                dt='-1.'
            xtab.extend(['../'+name+'.png','../'+name+'_2.png'])
            xt.extend(['S/N = '+str(sn),'S/N = '+str(sn)])
            tmp=tmp+dt+','
        tab.append(xtab)
        yt.append('<a href=../'+fdir+'>'+fdir+'</a><br>FERRE time: '+tmp+' s')

    header = '<h3>'+outfile+'</h3>\n'
    header = header + '<br> Test sample: <br><img src=test_'+grid+'.png width=40%> <img src=test_'+grid+'_2.png width=40%>'
    header = header + '<br> Histogram of pixel ratios between raw and PCA for test sample:<br>'
    header = header + '<img src=../'+root+'/'+outfile+'_pca.png width=40%>\n'

    html.htmltab(tab,file='test/test.html',header=header,ytitle=yt,xtitle=xt)
