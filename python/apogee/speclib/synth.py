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
from apogee.aspcap import aspcap
from apogee.aspcap import ferre
from apogee.speclib import atmos
from apogee.speclib import lsf
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
from tools import match

colors=['r','g','b','c','m','y']

def showtime(string) :
    """ Utiltiy routine to print a string and clock time
    """
    print(string+' {:8.2f}'.format(time.time()))
    sys.stdout.flush()

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

def marcs2turbo(infile,outfile,trim=0,fill=True) :
    """ Prepare MARCS input model for Turbospectrum, allowing for trimming of layers
    """
    try:
        fp=open(infile,'r')
    except :
        if not fill : return -1
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
    return 0

def mkturbospec(teff,logg,mh,am,cm,nm,wrange=[15100.,17000],dw=0.05,vmicro=2.0,solarisotopes=False,elemgrid='',welem=None,
    els=None,atmod=None,kurucz=True,atmosroot=None,atmosdir=None,nskip=0,endskip=0,fill=True,
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
        if logg <= 3.001 : geo = 's'
    atmosdir=atmosroot+'/'+atmosdir+'/'

    if abs(solarisotopes)==1: prefix=atmoscode+'d' 
    else : prefix=atmoscode+'g'

    # output directory and filename
    if save :
        workdir=(prefix+'m{:s}a{:s}c{:s}n{:s}v{:s}'+elemgrid).format(
                 atmos.cval(mh),atmos.cval(am),atmos.cval(cm),atmos.cval(nm),atmos.cval(vmicro))
        workdir=os.environ['APOGEE_LOCALDIR']+'/'+workdir
        try: os.mkdir(workdir)
        except: pass
    else :
        workdir=tempfile.mkdtemp(dir=os.environ['APOGEE_LOCALDIR'],suffix=os.environ['HOSTNAME'])

    root=workdir+'/'+(prefix+'t{:04d}g{:s}m{:s}a{:s}c{:s}n{:s}v{:s}'+elemgrid).format(int(teff), atmos.cval(logg), 
                      atmos.cval(mh), atmos.cval(am), atmos.cval(cm), atmos.cval(nm),atmos.cval(vmicro))
    print(root)
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
        if nskip > 2 : return 0.,0.
        kurucz2turbo(atmod,workdir+'/'+os.path.basename(atmod),trim=trim )
    else :
        try :        
            ret = marcs2turbo(atmod,workdir+'/'+os.path.basename(atmod),trim=nskip,fill=fill )
            if not fill and ret<0 :
                print('HOLE NOT SYNTHESIZED: ', atmod)
                return 0.,0.
        except:
            print('PROBLEM: ',atmod)
            fail('mkturbospec problem: '+atmod)
            return 0.,0.

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
        if 'C' in np.array(els)[:,0] : nels-=1
        if 'N' in np.array(els)[:,0] : nels-=1
        if 'O' in np.array(els)[:,0] : nels-=1

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
            if els is None or not 'C' in np.array(els)[:,0] :
                fout.write("    6  {:8.3f}\n".format(8.39+mh+cm))
            if els is None or not 'N' in np.array(els)[:,0] :
                fout.write("    7  {:8.3f}\n".format(7.78+mh+nm))
            if els is None or not 'O' in np.array(els)[:,0] :
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
        fout.write("'ABFIND:'  '.false.'\n")
        fout.write("'MODELINPUT:'  '{:s}'\n".format(os.path.basename(atmod)))
        if kurucz : fout.write("'MARCS-FILE:'  '.false.'\n")
        fout.write("'MODELOPAC:'  '{:s}'\n".format(os.path.basename(root)+'opac'))
        fout.write("'RESULTFILE:'  '{:s}'\n".format(os.path.basename(file)))
        fout.write("'METALLICITY:'  '{:8.3f}'\n".format(mh))
        fout.write("'ALPHA/Fe:'  '{:8.3f}'\n".format(am))
        fout.write("'HELIUM:'  '{:8.3f}'\n".format(0.00))
        fout.write("'R-PROCESS:'  '{:8.3f}'\n".format(0.00))
        fout.write("'S-PROCESS:'  '{:8.3f}'\n".format(0.00))
        fout.write("'INDIVIDUAL ABUNDANCES:'  '{:2d}'\n".format(nels))
        if els is None or not 'C' in np.array(els)[:,0] :
            fout.write("    6  {:8.3f}\n".format(8.39+mh+cm))
        if els is None or not 'N' in np.array(els)[:,0] :
            fout.write("    7  {:8.3f}\n".format(7.78+mh+nm))
        if els is None or not 'O' in np.array(els)[:,0] :
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
                    spec=out[:,2]
                    specnorm=out[:,1]
                else :
                    spec=np.vstack([spec,out[:,2]])
                    specnorm=np.vstack([specnorm,out[:,1]])
            except :
                print('failed...',file)
                fail('mkturbospec array problem: {:8d} {:8.2f} {:8.2f} {:8.2f}  {:8.2f} {:8.2f} {:8.2f}'.format(
                               teff,logg,mh,am,cm,nm,vmicro))
                return 0.,0.

    if run :
        if not save : shutil.rmtree(workdir)
        return spec, specnorm

def fail(out) :
    """ Routine to log to FAILURE file
    """
    ferr = open('FAILURE','a+')
    ferr.write(out+'\n')
    ferr.close()

def prange(start,delta,n) :
    """ Routine to return vector of values given start, delta, n
    """
    return float(start)+np.arange(int(n))*float(delta)

def get_vmicro(vmicrofit,vmicro) :
    """ Placeholder routine to return vmicro given a fit type
    """ 
    if vmicrofit == 0 :
        return float(vmicro)
    else :
        print('need to implement vmicrofit: ', vmicrofit)
        pdb.set_trace()
    return  float(vmicro)

def mkgrid(planfile,clobber=False,save=False,run=True,split=None,highres=9) :
    """ Create a grid of synthetic spectra using Turbospectrum given input parameter file
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
    solarisotopes = int(p['solarisotopes']) if p.get('solarisotopes') else 0
    solarisotopes = True if abs(solarisotopes) == 1 else False
    elem = p['elem'] if p.get('elem') else ''
    maskdir = p['maskdir'] if p.get('maskdir') else None
    vmicrofit = int(p['vmicrofit']) if p.get('vmicrofit') else 0
    vmicro = p['vmicro'] if p.get('vmicro') else 0
    vmacrofit = int(p['vmacrofit']) if p.get('vmacrofit') else 0
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
          specnormdata=np.zeros([nelem,int(p['nmh']),int(p['nlogg']),int(p['nteff']),nspec],dtype=np.int16)
          # allow for restart after certain number of [M/H] have been completed
          if clobber :
              nmh = 0
          else :
              try :
                  # does output file exist?
                  old=fits.open(specdir+'/'+p['name']+'.fits')[0]
                  specdata=old.data
                  if len(old.shape) < 5 : specdata=np.expand_dims(specdata,axis=0)
                  # is it a partially completed file with nmh card, or a completed file?
                  try:
                      nmh =old.header['nmh']
                  except :
                      # file is completed
                      nmh = p['nmh']
                      print('file already done!')
                      return
              except :
                  nmh = 0
          for imh,mh in enumerate(prange(p['mh0'],p['dmh'],p['nmh'])) :
           if imh >= nmh :
            for ilogg,logg in enumerate(prange(p['logg0'],p['dlogg'],p['nlogg'])) :
              for iteff,teff in enumerate(prange(p['teff0'],p['dteff'],p['nteff'])) :

                print(teff, logg, mh)
                sys.stdout.flush()
                nskip=0 
                dskip = 1 if kurucz else 2
                vout = get_vmicro(vmicrofit,vmicro)
                while nskip >= 0 and nskip < 10 :
                  spec,specnorm=mkturbospec(int(teff),logg,mh,am,cm,nm,
                    wrange=wrange,dw=dw,atmosdir=marcsdir,
                    elemgrid=elem,linelistdir=linelistdir+'/'+elem+'/',linelist=linelist,vmicro=vout,
                    solarisotopes=solarisotopes,
                    nskip=nskip,kurucz=kurucz,run=run,save=save,split=split) 
                  nskip = nskip+dskip if isinstance(spec,float) else -1
                if nskip > 0 : 
                    print('FAILED Turbospec',nskip)
                    fail('failed Turbospec convergence: {:8d} {:8.2f} {:8.2f} {:8.2f}  {:8.2f} {:8.2f} {:8.2f} {:d}'.format(
                               int(teff),logg,mh,am,cm,nm,vout,nskip))
                try:
                    if elem == '' :
                        specdata[0,imh,ilogg,iteff,:]=spec
                        specnormdata[0,imh,ilogg,iteff,:]=np.round((specnorm-0.5)*65534.).astype(int)
                    else :
                        specdata[:,imh,ilogg,iteff,:]=spec[:,gd]
                        specnormdata[:,imh,ilogg,iteff,:]=np.round((specnorm[:,gd]-0.5)*65534).astype(int)
                except :
                    print(specdata.shape)
                    specdata[:,imh,ilogg,iteff,:]=0.
                    specnormdata[:,imh,ilogg,iteff,:]=-32767
                    fail('error loading specdata: {:8d} {:8.2f} {:8.2f} {:8.2f}  {:8.2f} {:8.2f} {:8.2f} {:d}'.format(
                               int(teff),logg,mh,am,cm,nm,vout,len(spec)))

            # FITS header and output after each metallicity subgrid
            hdu=fits.PrimaryHDU(np.squeeze(specdata))
            idim=1
            if elem == '' :
                spectra.add_dim(hdu.header,rawwave[0],rawwave[1]-rawwave[0],1,'WAVELENGTH',idim)
            else :
                hdu.header.append(('CDELT1',rawwave[1]-rawwave[0]))
                for iwind in range(nwind) :
                    hdu.header.append(('WIND0_{:d}'.format(iwind),wvac[iwind,0]))
                    hdu.header.append(('WIND1_{:d}'.format(iwind),wvac[iwind,1]))
            if int(p['nteff']) > 1 :
                idim+=1
                spectra.add_dim(hdu.header,float(p['teff0']),float(p['dteff']),1,'TEFF',idim)
            if int(p['nlogg']) > 1 :
                idim+=1
                spectra.add_dim(hdu.header,float(p['logg0']),float(p['dlogg']),1,'LOGG',idim)
            if int(p['nmh']) > 1 :
                idim+=1
                spectra.add_dim(hdu.header,float(p['mh0']),float(p['dmh']),1,'M_H',idim)
            hdu.header['LOGW'] = 0
            if p.get('width') : hdu.header['width'] = p['width']
            if p.get('linelist') : hdu.header['linelist'] = p['linelist']
            if p['synthcode'] == 'asset'  : hdu.header.add_comment('ASSET generated synthetic spectra')
            if p['synthcode'] == 'turbospec' : hdu.header.add_comment('Turbospec generated synthetic spectra')
            if p['synthcode'] == 'moog ' : hdu.header.add_comment('MOOG generated synthetic spectra')
            hdu.header.add_comment('APOGEE_VER:'+os.environ['APOGEE_VER'])
            # following card for partial completion output
            if imh+1 < int(p['nmh']) : hdu.header['nmh'] = imh+1
            try : os.mkdir(specdir)
            except: pass
            hdunorm=fits.ImageHDU(np.squeeze(specnormdata))
            hdunorm.header.extend(hdu.header.copy(strip=True))
            hdunorm.header['BZERO'] = 0.5
            hdunorm.header['BSCALE'] = 1./65534.
            hdulist=fits.HDUList()
            hdulist.append(hdu)
            hdulist.append(hdunorm)
            hdulist.writeto(specdir+'/'+p['name']+'.fits',overwrite=True)

def mkgridlsf(planfile,clobber=False,highres=9,fiber=None,ls=None,comp=False,apred='r8',threads=32,rbf=False) :
    """ Create a grid of synthetic spectra using Turbospectrum given input parameter file
    """

    # Read planfile
    if not os.path.isfile(planfile): 
        print('{:s} does not exist'.format(planfile))
        return
    p=yanny.yanny(planfile,np=True)
    # header information
    specdir = os.environ['APOGEE_SPECLIB']+'/synth/'+p['specdir'] if p.get('specdir') else './'
    if fiber is None : fiber=np.array(p.get('lsffiber').split()).astype(int).tolist()
    if isinstance(fiber,int): fiber= [fiber]
    lsfid=int(p.get('lsfid'))
    waveid=int(p.get('waveid'))

    if ls is None :
        lsfile = 'lsf_{:08d}_{:08d}.fits'.format(lsfid,waveid)
        while os.path.isfile(lsfile+'.lock') : 
            print('waiting for lock: ',lsfile+'.lock')
            time.sleep(10)

        if os.path.isfile(lsfile) :
            x=fits.open(lsfile)[1].data
            ls=fits.open(lsfile)[2].data
        else :
            fp = open(lsfile+'.lock','w')
            fp.close()
            x,ls = lsf.get(lsfid,waveid,fiber,highres=highres,apred=apred)
            hdu=fits.HDUList()
            hdu.append(fits.PrimaryHDU())
            hdu[0].header['APRED'] = apred
            hdu[0].header['LSFID'] = lsfid
            hdu[0].header['WAVEID'] = waveid
            hdu[0].header['HIGHRES'] = highres
            for i,f in enumerate(fiber) :
                hdu[0].header['FIBER{:d}'.format(i)] = f
            hdu.append(fits.ImageHDU(x))
            hdu.append(fits.ImageHDU(ls))
            hdu.writeto(lsfile,overwrite=True)
            os.remove(lsfile+'.lock') 

    if rbf : prefix = 'rbf_'
    else : prefix = ''
    specdata = fits.open(specdir+'/'+prefix+p['name']+'.fits')[0]
    npix = specdata.data.shape[-1]
    nspec=1
    for i in range(len(specdata.data.shape)-1) :
        nspec*=specdata.data.shape[i]
    print('nspec: ', nspec)
    specdata.data=np.reshape(specdata.data,(nspec,npix))

    # synthesis is in air, we want vacuum
    ws=spectra.fits2vector(specdata.header,1)
    ws=spectra.airtovac(ws)

    # output wavelength grid
    wa=aspcap.apStarWave()
    nout=wa.shape[0]

    # create vmacro array
    vmacro=[]
    dlam=np.log10(wa[1])-np.log10(wa[0])
    for k,mh in enumerate(prange(p['mh0'],p['dmh'],p['nmh'])) :
      for j,logg in enumerate(prange(p['logg0'],p['dlogg'],p['nlogg'])) :
        for i,teff in enumerate(prange(p['teff0'],p['dteff'],p['nteff'])) :
            vm = 10.**(0.470794-0.254*mh)
            vm = vm if vm<15 else 15.
            vmacro.append(vm)
    vmacro=np.array(vmacro)

    # LSF and rotation convolution all spectra at the same time
    vrot=10.**.176
    smoothdata=lsf.convolve(ws,specdata.data,lsf=ls,xlsf=x,vrot=vrot,vmacro=vmacro)

    nmh=int(p['nmh'])
    nlogg=int(p['nlogg'])
    nteff=int(p['nteff'])
    specdata.data=np.reshape(specdata.data,(nmh,nlogg,nteff,npix))
    smoothdata=np.reshape(smoothdata,(nmh,nlogg,nteff,nout)).astype(np.float32)

    hdu=fits.PrimaryHDU(np.squeeze(smoothdata))
    hdu.header.extend(specdata.header.copy(strip=True))
    hdu.header['CRVAL1'] = aspcap.logw0
    hdu.header['CDELT1'] = aspcap.dlogw
    hdu.header['CTYPE1'] = 'LOG(WAVELENGTH)'
    hdu.writeto(p['name']+'.fits',overwrite=True)

    return smoothdata

def complsf(name) :
    new = fits.open(name+'.fits')[0]
    a = fits.open('old_'+name+'.fits')[1].data
    b = fits.open('old_'+name+'.fits')[2].data
    c = fits.open('old_'+name+'.fits')[3].data

    for k in range(new.header['NAXIS4']) :
      for j in range(new.header['NAXIS3']) :
        for i in range(new.header['NAXIS2']) :
          old = np.append(a[k,j,i,:],b[k,j,i,:])
          old = np.append(old,c[k,j,i,:])

          plt.clf()
          plt.plot(aspcap.aspcap2apStar(old))
          plt.plot(new.data[k,j,i,:]/600000.)
          plt.plot(aspcap.aspcap2apStar(old)/new.data[k,j,i,:]*600000.)
          plt.draw()
          pdb.set_trace()

def mkspec(pars) :
    """ Makes a single spectrum given input pars
    """
    teff=pars[0].astype('int')
    logg=pars[1]
    mh=pars[2]
    am=round(pars[3]/0.25)*0.25
    cm=round(pars[4]/0.25)*0.25
    nm=round(pars[5]/0.5)*0.5
    vmicro=pars[6]
    vrot=pars[7]
    elems=[]
    els = ['O','Na','Mg','Al','Si','P','S','K','Ca','Ti','V','Cr','Mn','Co','Ni','Cu','Ge','Rb','Ce','Nd']
    for j,el in enumerate(els) :
        elems.append([el,pars[8+j]])
    print(teff,logg,mh,vmicro,am,cm,nm)
    spec,specnorm=mkturbospec(teff,logg,mh,am,cm,nm,vmicro=vmicro,els=elems,kurucz=False,fill=False,linelist='20180721',wrange=[15100.,17000.],save=True)
    return pars,spec
    

def mksynth(file,threads=8,highres=9,waveid=2420038,lsfid=5440020,apred='r10',fiber='combo',plot=False,lines=None) :
    """ Make a series of spectra from parameters in an input file, with parallel processing for turbospec
    """
    pars=np.loadtxt(file)
    if lines is not None :
        pars = pars[lines[0]:lines[1]]
        suffix = '_{:d}'.format(lines[0])
    else :
        suffix = ''

    pool = mp.Pool(threads)
    specs = pool.map_async(mkspec, pars).get()
    pool.close()
    pool.join()

    # convolved and bundle output spectra into output fits file
    wa=aspcap.apStarWave()
    x,ls=lsf.get(lsfid,waveid,fiber,highres=highres,apred=apred)

    out=[]
    conv=[]
    outpar=[]
    if plot : plt.clf()
    ws=np.linspace(15100.,17000., len(specs[0][1]))
    # synthesis is in air, we want vacuum
    ws=spectra.airtovac(ws)
    for spec in specs :
        if isinstance(spec[1],np.ndarray) :
          if spec[1].sum() > 0.001 :
            mh=spec[0][2]
            vmacro = 10.**(0.470794-0.254*mh)
            vmacro = vmacro if vmacro<15 else 15.
            vrot=spec[0][7]
            if vrot < 0.5 : vrot=None
            # convolve one at a time because we have different vrot for each
            z=lsf.convolve(ws,spec[1],lsf=ls,xlsf=x,vmacro=vmacro,vrot=vrot)
            out.append(spec[1])
            conv.append(np.squeeze(z))
            outpar.append(spec[0])
            if plot : plt.plot(wa,np.squeeze(z))

    # write the spectra out
    hdu=fits.HDUList()
    hdu.append(fits.ImageHDU(outpar))
    hdu.append(fits.ImageHDU(out))
    hdu.append(fits.ImageHDU(conv))
    hdu.writeto(file+suffix+'.fits',overwrite=True)
   
def filter_lines(infile,outfile,wind,nskip=0) :
    """ Read from input linelist file, output comments and lines falling in windows of [w1,w2] to outfile
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

