# encoding: utf-8
#
# @Author: Jon Holtzman
# @Date: March 2018
# @Filename: synth.py
# @License: BSD 3-Clause
# @Copyright: Jon Holtzman

# routines related to creation of synthetic spectra

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
from synple import synple

colors=['r','g','b','c','m','y']

def showtime(string) :
    """ Utility routine to print a string and clock time, with flush to stdout

    Args:
        string(str) : string to print
    """
    print(string+' {:8.2f}'.format(time.time()))
    sys.stdout.flush()

def kurucz2turbo(infile,outfile,trim=0) :
    """ Convert Kurucz model atmosphere for use by Turbospectrum 
        Allows for trimming of layers

    Args:
        infile (str)  : name of input Kurucz atmosphere file
        outfile (str) : name of output file
        trim  (int)   : number of layers to remove (default=0)
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

    Args :
        infile (str)  : name of input MARCS file
        outfile (str) : name of output MARCS fie
        trim (int)    : number of layers to trim (default=0)

    Returns :
        0 if successful, -1 if not
    """
    try:
        fp=open(infile,'r')
    except :
        if not fill : return -1
        try :
            fp=open(infile+'.filled','r')
        except :
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filename)
            return -1

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

def get_atmod_file(teff,logg,mh,am,cm,nm,atmos_type='marcs',atmosroot=None,nskip=0,fill=True,workdir='./',atmosdir=None) :
    """ Copy atmosphere file to working directory in Turbospec format,  trimming as requested
    """
    # directory setup
    if atmosroot is None : atmosroot=os.environ['APOGEE_SPECLIB']+'/atmos/'
 
    # atmosphere and filename set up 
    geo = 'p'
    if atmos_type == 'kurucz' : 
        atmoscode= 'k'
        model = 'Kurucz'
        if atmosdir is None : atmosdir = '/kurucz/'
    elif atmos_type == 'marcs' : 
        atmoscode = 'm'
        model = 'MARCS'
        if atmosdir is None : atmosdir = '/marcs/MARCS_v3_2016/'
        if logg <= 3.001 : geo = 's'
    else :
        print('unknown atmos_type: ', atmos_type)
        pdb.set_trace()
    atmosdir=atmosroot+'/'+atmosdir+'/'
    atmod=atmosdir+atmos.filename(teff,logg,mh,cm,am,model=model)

    # atmosphere: make local copy in Turbospectrum input format, allowing for trimmed layers
    # note that [N/M] is solar for atmospheres (since it doesn't have a big effect)
    outmod = workdir+'/'+os.path.basename(atmod)
    if atmos_type == 'kurucz' :
        if nskip == 0 : trim=0
        if nskip == 1 : trim=7
        if nskip == 2 : trim=15
        if nskip > 2 : return 0.,0.
        kurucz2turbo(atmod,outmod,trim=trim )
    else :
        try :        
            ret = marcs2turbo(atmod,outmod,trim=nskip,fill=fill )
            if not fill and ret<0 : return ret
        except:
            return -2
    return outmod

def get_workdir(teff,logg,mh,am,cm,nm,atmos_type='marcs',solarisotopes=False,save=False,elemgrid=None,vmicro=1.) :
    """ Create work directory name for synthesis
    """
    if atmos_type == 'kurucz' : 
        atmoscode= 'k'
    elif atmos_type == 'marcs' : 
        atmoscode = 'm'
    else :
        print('unknown atmosphere type: ', atmos_type)
        pdb.set_trace()

    if abs(solarisotopes)==1: prefix=atmoscode+'d' 
    else : prefix=atmoscode+'g'

    if save :
        workdir=(prefix+'m{:s}a{:s}c{:s}n{:s}v{:s}'+elemgrid).format(
                 atmos.cval(mh),atmos.cval(am),atmos.cval(cm),atmos.cval(nm),atmos.cval(vmicro))
        workdir=os.environ['APOGEE_LOCALDIR']+'/'+workdir
        try: os.mkdir(workdir)
        except: pass
    else :
        workdir=tempfile.mkdtemp(dir=os.environ['APOGEE_LOCALDIR'],suffix=os.environ['HOSTNAME'])

    return workdir

def mk_synthesis(code,teff,logg,mh,am,cm,nm,wrange=[15100.,17000],dw=0.05,vmicro=2.0,solarisotopes=False,elemgrid='',welem=None,
    els=None,atmod=None,atmos_type='marcs',atmosroot=None,atmosdir=None,nskip=0,fill=True,
    linelist='20180901',h2o=None,linelistdir=None,atoms=True,molec=True, msuffix='', save=False,run=True) :
    """ Do synthesis

    Args:
        teff (int) : Effective temperature
        logg (float) : log(surface gravity)
        mh (float ) : [M/H]
        am (float ) : [alpha/M]
        cm (float ) : [C/M]
        nm (float ) : [N/M]
        wrange (list) : wavelength range (default [15100,17000.])
        dw (float)  : wavelength spacing (default 0.05 A)
        vmicro (float ) : microturbulent velocity (default 2.0 km/s)
        solarisotopes (bool) : use solar isotope ratios, else "giant" isotope ratios ( default False )
        elemgrid (str) :  name of element for minigrid calculation
        welem (str) :  NOT IMPLEMENTED set of wavelength ranges
        els (list of pairs) : list of [element name, abundance] pairs
        atmod (str) : name of atmosphere model (default=None, model is determined from input parameters)
        atmosroot (str) : root atmosphere directory (default=None --> $APOGEE_SPECLIB/atmos)
        atmosdir (str) :  atmosphere directory under atmosroot (default=None --> kurucz or marcs/MARCS_v3_2016, depending on kurucz boolean)
        nskip (int)  : number of layers to strip (for kurucz, nskip*7), default=0
        fill (bool)  : use filled atmosphere hole if model requests it ( default = True), else skip synthesis
        linelist (str) : name of linelist (default='20150714')
        linelistdir (str) : directory for linelist (default=None --> $APOGEE_SPECLIB/linelists)
        h2o (int or None) : H2O linelist to use: 0=none, 1=8.5V, 2=9.5V, None-->chooses by Teff and [alpha/H]
        save (bool ) : save temporary directory and files for synthesis (default=False)
        run (bool) :  actually do the synthesis (default=True)

    Returns:
        spec (np.array) : synthetic spectrum
        specnorm (np.array) : synthetic spectrum, normalized
    """
 
    # output directory and filename
    workdir = get_workdir(teff,logg,mh,am,cm,nm,atmos_type='marcs',solarisotopes=solarisotopes,save=save,elemgrid=elemgrid,vmicro=vmicro)
    root=(atmos_type+'_t{:04d}g{:s}m{:s}a{:s}c{:s}n{:s}v{:s}'+elemgrid).format(int(teff), atmos.cval(logg), 
                      atmos.cval(mh), atmos.cval(am), atmos.cval(cm), atmos.cval(nm),atmos.cval(vmicro))
    print(root)

    # atmosphere: make local copy in Turbospectrum input format, allowing for trimmed layers
    # note that [N/M] is solar for atmospheres (since it doesn't have a big effect)
    atmod = get_atmod_file(teff,logg,mh,am,cm,nm,atmos_type=atmos_type,nskip=nskip, atmosroot=atmosroot,atmosdir=atmosdir,workdir=workdir)
    if type(atmod) is int :
        if atmod == -1 :        
            print('HOLE NOT SYNTHESIZED: ', atmod)
            return 0.,0.
        else :
            print('PROBLEM: ',atmod)
            fail('mkturbospec problem: '+atmod)
            return 0.,0.

    # linelists
    if linelistdir is None : linelistdir=os.environ['APOGEE_SPECLIB']+'/linelists/' 
    if code == 'turbospec' :
        linelists = []
        n_HI = len(open(linelistdir+'/turbospec.'+linelist+'.Hlinedata').readlines())
        if n_HI > 2 : linelists.append(linelistdir+'/turbospec.'+linelist+'.Hlinedata')
        if atoms : linelists.append(linelistdir+'/turbospec.'+linelist+'.atoms')
        if molec : linelists.append(linelistdir+'/turbospec.'+linelist+msuffix+'.molec')
    elif code == 'synspec' :
        linelists = [linelistdir+'/synspec/synspec.'+linelist+'.atoms']
        if molec :
            #if solarisotopes : linelists.append(linelistdir+'/synspec/synspec.'+linelist+'sun_nofeh_noh2o_noc2.molec')
            #else : linelists.append(linelistdir+'/synspec/synspec.'+linelist+'giant_nofeh_noh2o_noc2.molec')
            if solarisotopes : linelists.append(linelistdir+'/synspec/synspec.'+linelist+'sun'+msuffix+'.molec')
            else : linelists.append(linelistdir+'/synspec/synspec.'+linelist+'giant'+msuffix+'.molec')
        h2o=0  
        #linelists = ['apogeeDR16.20180901.19','apogeeDR16_arc.20']
    else :
        print('unknown code!')
        pdb.set_trace()
    print('linelists: ',linelists)

    tfactor=1
    if h2o is None : 
      if teff < 4000 :
        if mh+am < -1.5 or teff > 3250 :
            linelists.append(linelistdir+'/turbospec.h2o-BC8.5V'+'.molec')
            tfactor=2
        else  :
            linelists.append(linelistdir+'/turbospec.h2o-BC9.5V'+'.molec')
            tfactor=5
    elif h2o == 1 :
        linelists.append(linelistdir+'/turbospec.h2o-BC8.5V'+'.molec')
        tfactor=2
    elif h2o == 2 :
        linelists.append(linelistdir+'/turbospec.h2o-BC9.5V'+'.molec')
        tfactor=5

    # default abundances
    abundances = atomic.solar()
    abundances[2:] += mh
    abundances[6-1] += cm
    abundances[7-1] += nm
    for i in [8,10,12,14,16,18,20,22] : abundances[i-1] += am

    # abundance overrides from els, given as [X/M]
    if els is not None :
        for el in els :
            atomic_num = atomic.periodic(el[0])
            abundances[atomic_num-1] = atomic.solar(el[0]) + mh + el[1] 

    if elemgrid != '' :
        elemnum=atomic.periodic(elemgrid)[0]
        elem0=abundances[elemnum-1]
        eabun=np.arange(-0.75,1.20,0.25)
    else :
        eabun=np.array([0.])

    cwd = os.getcwd()
    os.chdir(workdir)

    # if turbospec do the opacity calculations first, then synthesis for all eabun
    if code == 'turbospec' :
        out = do_turbospec(root,atmod,linelists,mh,am,abundances,wrange,dw,save=save,run=run,
                           solarisotopes=solarisotopes,bsyn=False,atmos_type=atmos_type,vmicro=vmicro)
    elif code == 'synspec' :
        abundances = 10.**(np.array(abundances)-abundances[0])

    # if this is an elemgrid, loop over abundances
    for ielem,abun in enumerate(eabun) :
        file = root+'{:02d}'.format(ielem)
        if elemgrid != '' :
            abundances[elemnum-1] = elem0+abun
     
        if code == 'turbospec' :
            if atmos_type == 'marcs' and logg <= 3.001 :  spherical= True
            else : spherical = False
            if ielem == 0 : 
                wave,flux,fluxnorm = do_turbospec(file,atmod,linelists,mh,am,abundances,wrange,dw,
                                                  save=save,run=run,solarisotopes=solarisotopes,
                                                  babsma=root+'opac',atmos_type=atmos_type,spherical=spherical,tfactor=tfactor)
            else :
                out = do_turbospec(file,atmod,linelists,mh,am,abundances,wrange,dw,
                                   save=save,run=run,solarisotopes=solarisotopes,
                                   babsma=root+'opac',atmos_type=atmos_type,spherical=spherical,tfactor=tfactor)
        elif code == 'synspec' :
            wave,flux,cont = synple.syn(atmod,wrange,linelist=linelists,dw=dw,vmicro=vmicro,save=save,clean=not save)
            #wave,flux,cont = synple.syn(atmod,wrange,linelist=linelists,dw=dw,abu=abundances,vmicro=vmicro,save=save)
            fluxnorm = flux/cont
        else :
            print('unknown synth code!')
            pdb.set_trace()

        # load into final arrays
        if ielem == 0 :
            spec = flux
            specnorm = fluxnorm
        else :
            spec=np.vstack([spec,flux])
            specnorm=np.vstack([specnorm,fluxnorm])

    os.chdir(cwd)
    if run :
        if not save : shutil.rmtree(workdir)
        return spec, specnorm
    pdb.set_trace()


def do_turbospec(file,atmod,linelists,mh,am,abundances,wrange,dw,save=False,run=True,solarisotopes=False,babsma=None,bsyn=True,atmos_type='marcs',spherical=True,vmicro=1.0,tfactor=1.) :
    """ Runs Turbospectrum for specified input parameters

    Args:

    Returns:
    """

    # Turbospectrum setup
    try: os.symlink(os.environ['APOGEE_DIR']+'/src/turbospec/DATA','./DATA')
    except: pass

    if save : stdout = None
    else : stdout = open(os.devnull, 'w')

    # individual element grid?
    nels = len(abundances)

    welem=np.array(wrange)
    # only compute opacities for a single nominal abundance
    if babsma is None :
        fout=open(file+'_babsma.csh','w')
        fout.write("#!/bin/csh -f\n")
        fout.write("{:s}/bin/babsma_lu << EOF\n".format(os.environ['APOGEE_DIR']))
        fout.write("'LAMBDA_MIN:'   '{:12.3f}'\n".format(welem.min()-dw))
        fout.write("'LAMBDA_MAX:'   '{:12.3f}'\n".format(welem.max()+dw))
        fout.write("'LAMBDA_STEP:'  '{:8.3f}'\n".format(dw))
        fout.write("'MODELINPUT:'  '{:s}'\n".format(os.path.basename(atmod)))
        if atmos_type != 'marcs' : fout.write("'MARCS-FILE:'  '.false.'\n")
        fout.write("'MODELOPAC:'  '{:s}opac'\n".format(os.path.basename(file)))
        fout.write("'METALLICITY:'  '{:8.3f}'\n".format(mh))
        fout.write("'ALPHA/Fe:'  '{:8.3f}'\n".format(am))
        fout.write("'HELIUM:'  '{:8.3f}'\n".format(0.00))
        fout.write("'R-PROCESS:'  '{:8.3f}'\n".format(0.00))
        fout.write("'S-PROCESS:'  '{:8.3f}'\n".format(0.00))
        fout.write("'INDIVIDUAL ABUNDANCES:'  '{:2d}'\n".format(nels))
        for iel,abun in enumerate(abundances) :
            fout.write("{:5d}  {:8.3f}\n".format(iel+1,abun))
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
            os.chmod(file+'_babsma.csh', 0o777)
            subprocess.call(['time','./'+os.path.basename(file)+'_babsma.csh'],stdout=stdout)
        babsma = os.path.basename(file)+'opac'

    if not bsyn : return

    # create bsyn control file
    bsynfile = file
    fout = open(bsynfile+'.inp','w')
    fout.write("'LAMBDA_STEP:'  '{:8.3f}'\n".format(dw))
    fout.write("'LAMBDA_MIN:'   '{:12.3f}'\n".format(welem.min()))
    fout.write("'LAMBDA_MAX:'   '{:12.3f}'\n".format(welem.max()))
    fout.write("'INTENSITY/FLUX:'  'Flux'\n")
    fout.write("'COS(THETA):'  '1.00'\n")
    fout.write("'ABFIND:'  '.false.'\n")
    fout.write("'MODELINPUT:'  '{:s}'\n".format(os.path.basename(atmod)))
    if atmos_type != 'marcs' : fout.write("'MARCS-FILE:'  '.false.'\n")
    fout.write("'MODELOPAC:'  '{:s}'\n".format(babsma))
    fout.write("'RESULTFILE:'  '{:s}'\n".format(os.path.basename(file)))
    fout.write("'METALLICITY:'  '{:8.3f}'\n".format(mh))
    fout.write("'ALPHA/Fe:'  '{:8.3f}'\n".format(am))
    fout.write("'HELIUM:'  '{:8.3f}'\n".format(0.00))
    fout.write("'R-PROCESS:'  '{:8.3f}'\n".format(0.00))
    fout.write("'S-PROCESS:'  '{:8.3f}'\n".format(0.00))
    fout.write("'INDIVIDUAL ABUNDANCES:'  '{:2d}'\n".format(len(abundances)))
    for iel,abun in enumerate(abundances) :
        fout.write("{:5d}  {:8.3f}\n".format(iel+1,abun))
    if  not solarisotopes :
        fout.write("'ISOTOPES:'  '2'\n")
        # adopt ratio of 12C/13C=15
        fout.write("   6.012 0.9375\n")
        fout.write("   6.013 0.0625\n")
    fout.write("'NFILES:'  '{:4d}'\n".format(len(linelists)))
    for linelist in linelists: 
        fout.write(linelist+"\n")
    if spherical: fout.write("'SPHERICAL:'  'T'\n")
    else : fout.write("'SPHERICAL:'  'F'\n")
    fout.write("30\n")
    fout.write("300.00\n")
    fout.write("15\n")
    fout.write("1.3\n")
    fout.close()

    # control file, with special handling in case bsyn goes into infinite loop ...
    fout = open(file+"_bsyn.csh",'w')
    fout.write("#!/bin/csh -f\n")
    fout.write("{:s}/bin/bsyn_lu < {:s} &\n".format(os.environ['APOGEE_DIR'],os.path.basename(bsynfile)+'.inp'))
    fout.write('set bsynjob = $!\n')
    fout.write("set ok = 0\n")
    fout.write("set runtime = `ps -p $bsynjob -o cputime | tail -1 | awk -F: '{print ($1*3600)+($2*60)+$3}'`\n")
    tmax=120*int(.05/min([0.05,dw]))*tfactor
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
        os.chmod(file+'_bsyn.csh', 0o777)
        subprocess.call(['time','./'+os.path.basename(file)+'_bsyn.csh'],stdout=stdout)
        try:
            out=np.loadtxt(file)
            wave=out[:,1]
            specnorm=out[:,1]
            spec=out[:,2]
        except :
            print('failed...',file,atmod,mh,am)
            return 0.,0.,0.
        return wave,spec,specnorm

def mkturbospec(teff,logg,mh,am,cm,nm,wrange=[15100.,17000],dw=0.05,vmicro=2.0,solarisotopes=False,elemgrid='',welem=None,
    els=None,atmod=None,kurucz=True,atmosroot=None,atmosdir=None,nskip=0,fill=True,
    linelist='20150714',h2o=None,linelistdir=None,atoms=True,molec=True,
    save=False,run=True) :
    """ Runs Turbospectrum for specified input parameters

    Args:
        teff (int) : Effective temperature
        logg (float) : log(surface gravity)
        mh (float ) : [M/H]
        am (float ) : [alpha/M]
        cm (float ) : [C/M]
        nm (float ) : [N/M]
        wrange (list) : wavelength range (default [15100,17000.])
        dw (float)  : wavelength spacing (default 0.05 A)
        vmicro (float ) : microturbulent velocity (default 2.0 km/s)
        solarisotopes (bool) : use solar isotope ratios, else "giant" isotope ratios ( default False )
        elemgrid (str) :  name of element for minigrid calculation
        welem (str) :  NOT IMPLEMENTED set of wavelength ranges
        els (list of pairs) : list of [element name, abundance] pairs
        atmod (str) : name of atmosphere model (default=None, model is determined from input parameters)
        atmosroot (str) : root atmosphere directory (default=None --> $APOGEE_SPECLIB/atmos)
        atmosdir (str) :  atmosphere directory under atmosroot (default=None --> kurucz or marcs/MARCS_v3_2016, depending on kurucz boolean)
        nskip (int)  : number of layers to strip (for kurucz, nskip*7), default=0
        fill (bool)  : use filled atmosphere hole if model requests it ( default = True), else skip synthesis
        linelist (str) : name of linelist (default='20150714')
        linelistdir (str) : directory for linelist (default=None --> $APOGEE_SPECLIB/linelists)
        h2o (int or None) : H2O linelist to use: 0=none, 1=8.5V, 2=9.5V, None-->chooses by Teff and [alpha/H]
        save (bool ) : save temporary directory and files for synthesis (default=False)
        run (bool) :  actually do the synthesis (default=True)

    Returns:
        spec (np.array) : synthetic spectrum
        specnorm (np.array) : synthetic spectrum, normalized
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
        nlists=0
        # if no HI lines, don't use that list: it takes a while to read/process
        n_HI = len(open(linelistdir+'/turbospec.'+linelist+'.Hlinedata').readlines())
        if n_HI > 2 : nlists+=1
        if atoms : nlists+=1
        if molec : nlists+=1
        # if we are using H2O, add that list
        if h2o > 0 : nlists+=1
        fout.write("'NFILES:'  '{:4d}'\n".format(nlists))
        if n_HI >= 3 : fout.write(linelistdir+'/turbospec.'+linelist+'.Hlinedata\n')
        if atoms: fout.write(linelistdir+'/turbospec.'+linelist+'.atoms\n')
        if molec: fout.write(linelistdir+'/turbospec.'+linelist+'.molec\n')
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

def fail(out,file='FAILURE') :
    """ Routine to log to FAILURE file

    Args:
        out (str) : string to write to file
        file (str) : name of failure file (default='FAILURE')
    """
    ferr = open(file,'a+')
    ferr.write(out+'\n')
    ferr.close()

def prange(start,delta,n) :
    """ Routine to return vector of values given start, delta, n as strings

    Args :
        start (str) : starting vlue
        delta (str) : delta
        n (str) : number of values

    Returns :
        np.array : array of n values, from start incremented by delta
    """
    return float(start)+np.arange(int(n))*float(delta)

def get_vmicro(vmicrofit,vmicro,teff=None,logg=None,mh=None) :
    """ Return vmicro from input functional type and coefficients
   
    Args :
        vmicrofit (int) : input vmicro code to set fitting function
        vmicro (float)  : input vmicro fitting function coefficients

    Returns:
        vmicro (float) : microturbulent velocity

    """ 
    if vmicrofit == 0 :
        vm = vmicro[0]
    elif vmicrofit == 1 :
        #cubic in logg
        vm = 10.**(vmicro[0]+vmicro[1]*logg+vmicro[2]*logg**2+vmicro[3]*logg**3)
    elif vmicrofit == 2 :
        #linear with terms in  teff, logg, mh
        vm = 10.**(vmicro[0]+vmicro[1]*teff+vmicro[2]*logg+vmicro[3]*mh)
    else :
        print('need to implement vmicrofit: ', vmicrofit)
        pdb.set_trace()
    return vm

def mkgrid(planfile,code=None,clobber=False,save=False,run=True) :
    """ Create a grid of synthetic spectra using Turbospectrum given specifications in  input parameter file
        Outputs results in FITS file

    Args :
        planfile (str) : name of input Yanny planfile with grid specifications
        clobber ( bool ) : skip trying to use previously saved syntheses (by [M/H]) (default = False)
        save (bool) :  save all temporary files for syntheses (default = False)
        run (bool)  : actually run syntheses (default = True)
    """

    # Read planfile
    if not os.path.isfile(planfile): 
        print('{:s} does not exist'.format(planfile))
        return
    p=yanny.yanny(planfile,np=True)

    # grid specifications from input planfile
    wrange=[float(x) for x in p['wrange'].split()]
    dw=float(p['dw'])
    vacuum = int(p['vacuum']) if p.get('vacuum') else 0
    kurucz = True if p['atmos'] == 'kurucz' else False
    marcsdir = p['marcsdir'] if p.get('marcsdir') else None
    solarisotopes = int(p['solarisotopes']) if p.get('solarisotopes') else 0
    solarisotopes = True if abs(solarisotopes) == 1 else False
    enhanced_o = p['enhanced_o'] if p.get('enhanced_o') else 0
    elem = p['elem'] if p.get('elem') else ''
    maskdir = p['maskdir'] if p.get('maskdir') else None
    vmicrofit = int(p['vmicrofit']) if p.get('vmicrofit') else 0
    vmicro = np.array(p['vmicro'].split()).astype(float) if p.get('vmicro') else 0.
    vmacrofit = int(p['vmacrofit']) if p.get('vmacrofit') else 0
    vmacro = p['vmacro'] if p.get('vmacro') else 0
    specdir = os.environ['APOGEE_SPECLIB']+'/synth/'+p['specdir'] if p.get('specdir') else './'
    linelistdir=os.environ['APOGEE_SPECLIB']+'/linelists/' 
    linelist = p['linelist'] if p.get('linelist') else None
    oa0 = p['oa0'] if p.get('oa0') else 0.
    doa = p['doa'] if p.get('doa') else 0.
    noa = p['noa'] if p.get('noa') else 1

    # wavelength array
    nspec=int((wrange[1]-wrange[0])/dw)+1
    rawwave=wrange[0]+np.arange(nspec)*dw

    # if element minigrid, create mini linelist and get wavelengths to store
    if elem == '' :
        nelem=1
        gdspec=range(nspec)
        nwind = 1
        pixels=[[0,nspec]]
    else :
        # number of minigrid abundances
        nelem=8
        wvac,wair = mini_linelist(elem,linelist,maskdir=maskdir)
        nwind=wair.shape[0]
        gdspec=[]
        pixels=[]
        for iwind in range(nwind) :
            pix=np.where( (rawwave >= wair[iwind,0]) & (rawwave <= wair[iwind,1]) )[0]
            gdspec.extend(pix)
            pixels.append([pix[0],pix[-1]+1])
        nspec=len(gdspec)
        # for elem with all waves, use next line and comment out previous 8
        #gdspec=range(nspec)
        #nwind = 1
        #pixels=[[0,nspec]]

    # make the grid(s)
    for oa in prange(oa0,doa,noa) :
     for am in prange(p['am0'],p['dam'],p['nam']) :
      if enhanced_o : oa = [('O',2*am)]
      else : oa = [('O',oa+am)]
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
                  old=fits.open(specdir+'/'+p['name']+elem+'.fits')[0]
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
                vout = get_vmicro(vmicrofit,vmicro,teff=teff,logg=logg,mh=mh)
                print(teff,logg,mh,am,cm,nm,oa,vout)
                while nskip >= 0 and nskip < 10 :
                  if code is None :
                      spec,specnorm=mkturbospec(int(teff),logg,mh,am,cm,nm,els=oa,
                        wrange=wrange,dw=dw,atmosdir=marcsdir,
                        elemgrid=elem,linelistdir=linelistdir+'/'+elem+'/',linelist=linelist,vmicro=vout,
                        solarisotopes=solarisotopes,
                        nskip=nskip,kurucz=kurucz,run=run,save=save) 
                  else :
                      spec,specnorm=mk_synthesis(code,int(teff),logg,mh,am,cm,nm,els=oa,
                        wrange=wrange,dw=dw,atmosdir=marcsdir,
                        elemgrid=elem,linelistdir=linelistdir+'/'+elem+'/',linelist=linelist,vmicro=vout,
                        solarisotopes=solarisotopes,
                        nskip=nskip,atmos_type=p['atmos'],run=run,save=save,h2o=0) 
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
                        specdata[:,imh,ilogg,iteff,:]=spec[:,gdspec]
                        specnormdata[:,imh,ilogg,iteff,:]=np.round((specnorm[:,gdspec]-0.5)*65534).astype(int)
                except :
                    print(specdata.shape)
                    specdata[:,imh,ilogg,iteff,:]=0.
                    specnormdata[:,imh,ilogg,iteff,:]=-32767
                    try: lspec=len(spec)
                    except: lspec=1
                    fail('error loading specdata: {:8d} {:8.2f} {:8.2f} {:8.2f}  {:8.2f} {:8.2f} {:8.2f} {:d}'.format(
                               int(teff),logg,mh,am,cm,nm,vout,lspec))

            # FITS header and output after each metallicity subgrid
            # for minigrids, each section is output in a separate HDU
            p1=0
            hdulist=fits.HDUList()
            for iwind in range(nwind) :
                p2 = p1 + pixels[iwind][1]-pixels[iwind][0]
                print(iwind,p1,p2)
                if iwind == 0 :
                    hdu=fits.PrimaryHDU(np.squeeze(specdata[:,:,:,:,p1:p2]))
                else :
                    hdu=fits.ImageHDU(np.squeeze(specdata[:,:,:,:,p1:p2]))
                p1 += pixels[iwind][1]-pixels[iwind][0]
                idim=1
                # for elem with all waves, use next line and comment out following 7
                #spectra.add_dim(hdu.header,rawwave[0],rawwave[1]-rawwave[0],1,'WAVELENGTH',idim)
                if elem == '' :
                    spectra.add_dim(hdu.header,rawwave[0],rawwave[1]-rawwave[0],1,'WAVELENGTH',idim)
                else :
                    hdu.header.append(('ELEM',elem))
                    gd = np.where( (rawwave >= wair[iwind,0]) & (rawwave <= wair[iwind,1]) )[0]
                    spectra.add_dim(hdu.header,rawwave[gd[0]],rawwave[1]-rawwave[0],1,'WAVELENGTH',idim)
                    hdu.header.append(('CDELT1',rawwave[1]-rawwave[0]))
                if int(p['nteff']) > 1 :
                    idim+=1
                    spectra.add_dim(hdu.header,float(p['teff0']),float(p['dteff']),1,'TEFF',idim)
                if int(p['nlogg']) > 1 :
                    idim+=1
                    spectra.add_dim(hdu.header,float(p['logg0']),float(p['dlogg']),1,'LOGG',idim)
                if int(p['nmh']) > 1 :
                    idim+=1
                    spectra.add_dim(hdu.header,float(p['mh0']),float(p['dmh']),1,'M_H',idim)
                if elem != '' :
                    idim+=1
                    spectra.add_dim(hdu.header,-0.75,0.25,1,elem,idim)
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
                hdulist.append(hdu)
            hdunorm=fits.ImageHDU(np.squeeze(specnormdata))
            hdunorm.header.extend(hdu.header.copy(strip=True))
            hdunorm.header['BZERO'] = 0.5
            hdunorm.header['BSCALE'] = 1./65534.
            hdulist.append(hdunorm)
            hdulist.writeto(specdir+'/'+p['name']+elem+'.fits',overwrite=True)

def mkgridlink(planfile,suffix=None) :
    """  DEVELOPMENT : create coarse grid by merging syntheses from multiple grids
    """

    # Read planfile
    if not os.path.isfile(planfile): 
        print('{:s} does not exist'.format(planfile))
        return
    p=yanny.yanny(planfile,np=True)

    linelist=p['linelist'][2:]
    if suffix is None :
        suffix='_'+p['smooth']

    for ivm,vm in enumerate(spectra.vector(p['vt0'],p['dvt'],p['nvt'])) :
      for icm,cm in enumerate(spectra.vector(p['cm0'],p['dcm'],p['ncm'])) :
        for inm,nm in enumerate(spectra.vector(p['nm0'],p['dnm'],p['nnm'])) :
          for iam,am in enumerate(spectra.vector(p['am0'],p['dam'],p['nam'])) :
            file=('a{:s}c{:s}n{:s}v{:s}.fits').format(
                   atmos.cval(am),atmos.cval(cm),atmos.cval(nm),atmos.cval(10**vm))
            print(file)

            GKg= fits.open('../../giantisotopes/tgGK_'+linelist+suffix+'/rbf_'+file)[0]
            Mg= fits.open('../../giantisotopes/tgM_'+linelist+suffix+'/rbf_'+file)[0]
            Fg= fits.open('../../giantisotopes/tgF_'+linelist+suffix+'/'+file)[0]
            Fd= fits.open('../../solarisotopes/tdF_'+linelist+suffix+'/rbf_'+file)[0]
            GKd= fits.open('../../solarisotopes/tdGK_'+linelist+suffix+'/rbf_'+file)[0]
            Md= fits.open('../../solarisotopes/tdM_'+linelist+suffix+'/rbf_'+file)[0]
            nwave = GKg.shape[-1]
            s=np.zeros([int(p['nmh']),int(p['nlogg']),int(p['nteff']),nwave],dtype=np.float32)
            grids = [Fg, GKg, Mg, Fd, GKd, Md]
            for imh,mh in enumerate(spectra.vector(p['mh0'],p['dmh'],p['nmh'])) :
              for ilogg,logg in enumerate(spectra.vector(p['logg0'],p['dlogg'],p['nlogg'])) :
                for iteff,teff in enumerate(spectra.vector(p['teff0'],p['dteff'],p['nteff'])) :
                    igrid=0
                    i,j,k=getindex(Fg.header,(2,3,4),(teff,logg,mh))
                    if i<0 or j<0 or k<0 :
                        igrid=1
                        i,j,k=getindex(GKg.header,(2,3,4),(teff,logg,mh))
                    if i<0 or j<0 or k<0 :
                        igrid=2
                        i,j,k=getindex(Mg.header,(2,3,4),(teff,logg,mh))
                    if i<0 or j<0 or k<0 :
                        igrid=3
                        i,j,k=getindex(Fd.header,(2,3,4),(teff,logg,mh))
                    if i<0 or j<0 or k<0 :
                        igrid=4
                        i,j,k=getindex(GKd.header,(2,3,4),(teff,logg,mh))
                    if i<0 or j<0 or k<0 :
                        igrid=5
                        i,j,k=getindex(Md.header,(2,3,4),(teff,logg,mh))
                    if i<0 or j<0 or k<0 :
                        print("can't find model to fill!",mh,logg,teff)
                    else :
                        #print(mh,logg,teff,igrid,i,j,k)
                        try:
                            #print(grids[igrid].data.shape)
                            s[imh,ilogg,iteff,:]=grids[igrid].data[k,j,i,:]
                        except:
                            pdb.set_trace()
            hdulist=fits.HDUList()
            hdu=fits.ImageHDU(s)
            idim=1
            hdu.header['CRVAL1'] = GKg.header['CRVAL1']
            hdu.header['CDELT1'] = GKg.header['CDELT1']
            hdu.header['CTYPE1'] = GKg.header['CTYPE1']
            if int(p['nteff']) > 1 :
                idim+=1
                spectra.add_dim(hdu.header,float(p['teff0']),float(p['dteff']),1,'TEFF',idim)
            if int(p['nlogg']) > 1 :
                idim+=1
                spectra.add_dim(hdu.header,float(p['logg0']),float(p['dlogg']),1,'LOGG',idim)
            if int(p['nmh']) > 1 :
                idim+=1
                spectra.add_dim(hdu.header,float(p['mh0']),float(p['dmh']),1,'M_H',idim)
            hdulist.append(hdu)
            hdulist.writeto(file,overwrite=True)

def getindex(header,axes,vals) :
    """ DEVELOPMENT : get index of requested model from input grid header

    """
    out=[]
    for ax,val in zip(axes,vals) :
        i = (val-header['CRVAL'+str(ax)])/header['CDELT'+str(ax)]
        if i+1 > header['NAXIS'+str(ax)] : i=-1.
        out.append(int(round(i)))
    return out

def mkgridlsf(planfile,highres=9,fiber=None,ls=None,apred=None,prefix=None,telescope=None) :
    """ Create a grid of LSF-convolved spectra given specifications in input parameter file and existing raw syntheses

    Args :
        planfile (str) : name of Yanny parameter file with specifications
        highres (int ) : number of subpixels to use for convolution (default=9)
        fiber (int or None ) : fiber number for LSF, if None use parameter file (default=None)
        ls () : lsf convolution kernel (?), if precomputed (default=None)
        apred (str) : name of reduction version to get LSF and WAVE from (default='r8')
        prefix (str) : string to prepend for input syntheses file, e.g. 'rbf_' to use rbf-filled grids (default='')
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
    if apred is None :apred = p['apred'] if p.get('apred') else 'r10'
    if telescope is None : telescope = p['telescope'] if p.get('telescope') else 'apo25m'
    lsfid=int(p.get('lsfid'))
    waveid=int(p.get('waveid'))
    vmacrofit = int(p['vmacrofit']) if p.get('vmacrofit') else 0
    vmacro_arr=np.array(p['vmacro'].split()).astype(float)
    kernel=p['kernel'] if p.get('kernel') else 'rot'

    if ls is None :
        x, ls = getlsf(lsfid,waveid,apred=apred,telescope=telescope,fiber=fiber,highres=highres)
    else :
        x = ls[0]
        ls = ls[1]

    # read in raw spectra, from RBF interpolated if we have it
    if prefix is None :  
        if p.get('r0') and float(p['r0']) >= -0.001 : prefix = 'rbf_'
        else : prefix=''

    # get the synthesis
    # for regular grids, we just want to process the first extension,
    # for minigrids, we want to process all except the last (which is the normalized synthesis)
    speclist = fits.open(specdir+'/'+prefix+p['name']+'.fits')
    nexten = len(speclist) - 1

    hdulist=fits.HDUList()
    for exten in range(nexten) :

        specdata = speclist[exten]
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
    #wa=aspcap.apStarWave()
    #nout=wa.shape[0]

        # is this a minigrid?
        try :
            if p['elem'] == '' : nelem = 1
            else : nelem=8
        except : nelem=1

        # create vmacro array
        vmacro=[]
        #dlam=np.log10(wa[1])-np.log10(wa[0])
        for l in range(nelem) :
          for k,mh in enumerate(prange(p['mh0'],p['dmh'],p['nmh'])) :
            for j,logg in enumerate(prange(p['logg0'],p['dlogg'],p['nlogg'])) :
              for i,teff in enumerate(prange(p['teff0'],p['dteff'],p['nteff'])) :
                if vmacrofit == 1 :
                    vm = 10.**(0.470794-0.254*mh)
                    vm = 10.**(vmacro_arr[0]+vmacro_arr[1]*teff+vmacro_arr[2]*logg+vmacro_arr[3]*mh)
                    vm = vm if vm<15 else 15.
                elif vmacrofit == 0 :
                    vm = 0.
                vmacro.append(vm)
        vmacro=np.array(vmacro)

        # LSF and rotation convolution all spectra at the same time
        nmh=int(p['nmh'])
        nlogg=int(p['nlogg'])
        nteff=int(p['nteff'])
        nrot = int(p['nrot'])

        if nrot == 1 :
            vrot=10.**.176
            smoothdata,waveout=lsf.convolve(ws,specdata.data,lsf=ls,xlsf=x,vrot=vrot,vmacro=vmacro)
            nout=smoothdata.shape[-1]
            smoothdata=np.reshape(smoothdata,(nelem,nmh,nlogg,nteff,nout)).astype(np.float32)
        else :
            smoothdata=np.zeros([nrot,nelem,nmh,nlogg,nteff,nout],dtype=np.float32)
            for irot,vrot in enumerate(prange(p['rot0'],p['drot'],p['nrot'])) :
                if kernel == 'rot' :
                    smooth,waveout=lsf.convolve(ws,specdata.data,lsf=ls,xlsf=x,vrot=10.**vrot,vmacro=vmacro)
                elif kernel == 'gauss' :
                    smooth,waveout=lsf.convolve(ws,specdata.data,lsf=ls,xlsf=x,vmacro=10.**vrot)
                else :
                    print('Unknown kernel!')
                    pdb.set_trace()
                smoothdata[irot,:,:,:,:,:]=np.reshape(smooth,(nelem,nmh,nlogg,nteff,nout)).astype(np.float32)

        specdata.data=np.reshape(specdata.data,(nelem,nmh,nlogg,nteff,npix))
        if exten == 0 : hdu=fits.PrimaryHDU(np.squeeze(smoothdata))
        else : hdu=fits.ImageHDU(np.squeeze(smoothdata))
        hdu.header.extend(specdata.header.copy(strip=True))
        #hdu.header['CRVAL1'] = aspcap.logw0
        hdu.header['CRVAL1'] = np.log10(waveout[0])
        hdu.header['CDELT1'] = aspcap.dlogw
        hdu.header['CTYPE1'] = 'LOG(WAVELENGTH)'
        if nrot > 1 :
            hdu.header.insert('CTYPE4',('CRVAL5',float(p['rot0']),''),after=True)
            hdu.header.insert('CRVAL5',('CDELT5',float(p['drot']),''),after=True)
            hdu.header.insert('CDELT5',('CRPIX5',1,''),after=True)
            hdu.header.insert('CRPIX5',('CTYPE5','LOG(VSINI)',''),after=True)
        hdu.header['INFILE'] = p['specdir']+'/'+prefix+p['name']+'.fits'
        hdu.header['APRED'] = apred
        hdu.header['LSFID'] = lsfid
        hdu.header['WAVEID'] = waveid
        hdu.header['HIGHRES'] = highres
        hdu.header.add_comment('LSF convolved spectra')
        hdu.header.add_comment('APOGEE_VER:'+os.environ['APOGEE_VER'])
        hdulist.append(hdu)

    hdulist.writeto(p['name']+'.fits',overwrite=True)

    return smoothdata

def complsf(name) :
    """ TEST : routine to compare calculated LSFs between Python and IDL
    """
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

def mkspec(input) :
    """ Makes a single spectrum given input pars
        Used by mksynth for multi-processor calculations

    Args :
        pars (list ) : list of spectrum parameters: [Teff, logg, [M/H], [alpha/M], [C/M], [N/M], vmicro, vrot, 20 abundances]

    Returns :
        pars : input parameters
        spec : synthetic spectrum
    """
    pars=input[0]
    indata=input[1]

    teff=pars[0].astype('int')
    logg=pars[1]
    mh=pars[2]
    am=round(pars[3]/0.25)*0.25
    cm=round(pars[4]/0.25)*0.25
    #nm=round(pars[5]/0.5)*0.5
    nm=pars[5]
    vmicro=pars[6]
    vrot=pars[7]
    elems=[]
    els = ['O','Na','Mg','Al','Si','P','S','K','Ca','Ti','V','Cr','Mn','Co','Fe','Ni','Cu','Ge','Rb','Ce','Nd']
    for j,el in enumerate(els) :
        elems.append([el,pars[8+j]])
    print(teff,logg,mh,vmicro,am,cm,nm)
    spec,specnorm=mkturbospec(teff,logg,mh,am,cm,nm,vmicro=vmicro,els=elems,kurucz=indata['kurucz'],fill=False,
                              linelist=indata['linelist'],linelistdir=indata['linelistdir'],
                              wrange=indata['wrange'],dw=indata['dw'],h2o=indata['h2o'],atoms=indata['atoms'],save=False)
    return pars,specnorm
    

def mksynth(file,threads=8,highres=9,waveid=2420038,lsfid=5440020,apred='r10',telescope='apo25m',
            fiber='combo',linelist='20180901',linelistdir=None,kurucz=False,h2o=None,atoms=True,plot=False,lines=None,ls=None) :
    """ Make a series of spectra from parameters in an input file, with parallel processing for turbospec
        Outputs to FITS file {file}.fits

    Args:
        file (str) : name of input file with parameters of spectra to calculate
        threads (int) : number of parallel processes
        highres (int) : number of subpixels for LSF convolution (default=9)
        waveid (int) : ID for wavelength calibration file (default=2420038)
        lsfid (int) : ID for LSF calibration file (default=5440020)
        apred (str) : reduction version for waveid,lsfid (default='r10')
        fiber (int or str) : fiber for LSF (default='combo')
        plot (bool) : plot each spectrum (default=False)
        lines (list) : specify limited range of input lines to calculate (default=None, i.e. all lines)
    """
    names=ascii.read(file).colnames
    pars=np.loadtxt(file)
    if lines is not None :
        pars = pars[lines[0]:min([len(pars),lines[1]])]
        suffix = '_{:d}'.format(lines[0])
    else :
        suffix = ''

    indata={}
    indata['linelist']=linelist
    indata['linelistdir']=linelistdir
    indata['kurucz']=kurucz
    indata['h2o']=h2o
    indata['atoms']=atoms
    indata['wrange']=[15100.,17000.]
    indata['dw']=0.05
    nspec=38001
    inputs=[]
    for par in pars :
        inputs.append((par,indata))

    if threads == 0 :
        specs=[]
        for input in inputs :
            specs.append(mkspec(input))
    else :
        pool = mp.Pool(threads)
        specs = pool.map_async(mkspec, inputs).get()
        pool.close()
        pool.join()

    # convolved and bundle output spectra into output fits file
    wa=aspcap.apStarWave()
    prefix='lsf_'
    if ls is None :
        x, ls = getlsf(lsfid,waveid,prefix=prefix,apred=apred,telescope=telescope,fiber=fiber,highres=highres)
    else :
        x=ls[0]
        ls=ls[1]

    out=[]
    conv=[]
    outpar=[]
    if plot : plt.clf()
    ws=np.linspace(15100.,17000., nspec)
    # synthesis is in air, we want vacuum
    ws=spectra.airtovac(ws)
    print('convolving...')
    for spec in specs :
        if isinstance(spec[1],np.ndarray) :
          if spec[1].sum() > 0.001 :
            mh=spec[0][2]
            vmacro = 10.**(0.470794-0.254*mh)
            vmacro = vmacro if vmacro<15 else 15.
            vrot=spec[0][7]
            if vrot < 0.5 : vrot=None
            # convolve one at a time because we have different vrot for each
            z,waveout=lsf.convolve(ws,spec[1],lsf=ls,xlsf=x,vmacro=vmacro,vrot=vrot)
            out.append(spec[1])
            conv.append(np.squeeze(z))
            outpar.append(spec[0])
            if plot : plt.plot(wa,np.squeeze(z))

    # write the spectra out
    hdu=fits.HDUList()
    h=fits.ImageHDU(outpar)
    h.header['NPAR' ] = len(names)
    for ipar,par in enumerate(names) : h.header['PAR{:d}'.format(ipar)] = par
    hdu.append(h)
    hdu.append(fits.ImageHDU(out))
    h=fits.ImageHDU(conv)
    h.header['CRVAL1'] = np.log10(wa[0])
    h.header['CDELT1'] = np.log10(wa[1])-np.log10(wa[0])
    h.header['CTYPE1'] = 'LOG(WAVELENGTH)'
    hdu.append(h)
    hdu.writeto(file+suffix+'.fits',overwrite=True)
    return file+suffix+'.fits' 

def getlsf(lsfid,waveid,apred='r10',telescope='apo25m',highres=9,prefix='lsf_',fiber='combo',clobber=False,fill=False) :
    """ Create LSF FITS file or read if already created
    """
    lsfile = prefix+'{:08d}_{:08d}.fits'.format(lsfid,waveid)
    print(lsfile)
    while os.path.isfile(lsfile+'.lock') :
        # if another process is creating LSF wait until done
        print('waiting for lock: ',lsfile+'.lock')
        time.sleep(10)

    if os.path.isfile(lsfile) and not clobber :
        # if file exists, read it
        x=fits.open(lsfile)[1].data
        ls=fits.open(lsfile)[2].data
    else :
        fp = open(lsfile+'.lock','w')
        fp.close()
        # lsf.get does the real work
        x,ls = lsf.get(lsfid,waveid,fiber,highres=highres,apred=apred,telescope=telescope)
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

    if fill :
        # for all non-finite pixels, fill in LSF from nearest good pixel
        gd = np.where(np.isfinite(ls[:,0]))[0]
        mask = np.zeros(ls.shape[0],dtype=bool)
        mask[gd] = True
        bd = np.where(mask == False)[0]
        for i in bd:
            j = np.argmin(np.abs(i-gd))
            ls[i,:] = ls[gd[j],:]
    return x, ls

   
def elemsens(files=None,outfile='elemsens',highres=9,waveid=13140000,lsfid=14600018,apred='r12',telescope='apo25m',fiber='combo',
             calc=False,plot=True,ls=None,htmlfile='elemsens.html',filt=None,filtdir=None,linelist='20180901') :
    """ Create spectra at a range of paramters with individual abundances varied independently
    """
    if files==None : files=sample.elemsens()
    hdu=fits.HDUList()
    hdumask=fits.HDUList()
    if ls is None :
        xls, ls = getlsf(lsfid,waveid,apred=apred,telescope=telescope,fiber=fiber,highres=highres)
    else :
        xls = ls[0]
        ls = ls[1]
    grid=[]
    ytit=[]
    for i,name in enumerate(files[3:]) :
        elem=name.replace('.dat','')
        print(name)
        if name == 'C.dat' or name == 'N.dat' : continue
        if calc : 
            #out=mksynth(name,threads=16,ls=(x,ls))
            mini_linelist(elem,linelist,only=True)
            out0=mksynth('ref.dat',threads=16,ls=(xls,ls),linelistdir=os.environ['APOGEE_SPECLIB']+'/linelists/'+elem+'_only',h2o=0)
            os.rename('ref.dat.fits',elem+'_ref.fits')
            out1=mksynth(name,threads=16,ls=(xls,ls),linelistdir=os.environ['APOGEE_SPECLIB']+'/linelists/'+elem+'_only',h2o=0)
            os.rename(name+'.fits',elem+'.fits')
        #ehdu=fits.open(out)[2]
        ehdu=fits.open(elem+'.fits')[2]
        ref=fits.open(elem+'_ref.fits')[2].data
        ehdu.header['ELEM'] = elem
        ehdu.header['CRVAL1'] = aspcap.logw0
        ehdu.header['CDELT1'] = aspcap.dlogw
        ehdu.header['CTYPE1'] = 'LOG(WAVELENGTH)'
        if plot :
            if i==-1 :
                ref=ehdu.data
            else :
                if filtdir is not None: 
                    filt=aspcap.aspcap2apStar(ascii.read(filtdir+ehdu.header['ELEM']+'.mask',format='fixed_width_no_header')['col1'])
                fig,ax=plots.multi(3,3,hspace=0.001,wspace=0.001,xtickrot=60,figsize=(12,8))
                nspec = ehdu.data.shape[0]
                ispec=0
                gd=[]
                for ix in range(3) :
                    te=3500+ix*1000
                    for iy in range(3) :
                        logg=1.+iy*2.
                        x=10.**spectra.fits2vector(ehdu.header,1)
                        y=ehdu.data[ispec,:]/ref[ispec,:]
                        plots.plotl(ax[iy,ix],x,y,yr=[0.9,1.1],xr=[15100,16950])
                        ax[iy,ix].text(0.05,0.9,'Teff:{:6.0f} logg:{:6.1f}'.format(te,logg),transform=ax[iy,ix].transAxes)
                        if filt is not None: plots.plotl(ax[iy,ix],x,filt*0.1+1.005)
                        j=np.where(y < 0.99)[0]
                        gd.extend(j)
                        ispec+=1
                mask=np.zeros(ehdu.data.shape[-1])
                mask[list(set(gd))] = 1.
                figname = ehdu.header['ELEM'].strip()+'.png'
                fig.savefig(figname)
                plt.close()
                grid.append([figname])
                ytit.append(ehdu.header['ELEM'])
                mhdu=fits.ImageHDU(mask)
                mhdu.header['ELEM'] = name.replace('.dat','')
                mhdu.header['CRVAL1'] = aspcap.logw0
                mhdu.header['CDELT1'] = aspcap.dlogw
                mhdu.header['CTYPE1'] = 'LOG(WAVELENGTH)'
                hdumask.append(mhdu)
        ehdu.data -= ref
        hdu.append(ehdu)
    # output single file
    if outfile is not None: 
        hdu.writeto(outfile+'.fits',overwrite=True)
        hdumask.writeto(outfile+'_mask.fits',overwrite=True)
    html.htmltab(grid,ytitle=ytit,file=htmlfile)

def mkmask(file='elemsens')  :

    mask=fits.open(file+'.fits')
    els=[]
    for i in range(len(mask)) : els.append(mask[i].header['ELEM'])
    els = np.array(els)
    alphas=np.array(['O','Mg','Si','S','Ca','Ti'])
    metals=np.array(['Na','Al','P','K','V','Cr','Mn','Co','Fe','Ni','Cu','Ge','Rb','Ce','Nd'])
    fig,ax=plots.multi(1,2,hspace=0.001,sharex=True)
    x=10.**spectra.fits2vector(mask[0].header,1)
    for i,el in enumerate(els) :
        print(el)
        gd=[]
        ax[0].cla()
        for k in range(9) : 
            plots.plotl(ax[0],x,mask[i].data[k,:])
            gd.extend(np.where(mask[i].data[k,:] < -0.01)[0])
        if el in alphas :
            bd=[]
            for al in alphas :
                ax[1].cla()
                if el != al :
                    print(el,al)
                    j=np.where(els == al)[0][0]
                    for k in range(9) : 
                        plots.plotl(ax[1],x,mask[j].data[k,:])
                        #plots.plotl(ax[1],x,mask[j].data[k,:]/mask[i].data[k,:] )
                        bd.extend(np.where((mask[j].data[k,:]/mask[i].data[k,:] < 0.2) & (mask[j].data[k,:]<-0.01) )[0])
                    plt.show()
                    #gd=np.where(abs(mask[j].data) > 0.)[0]
                    #mask[i].data[gd] = -1.*mask[i].data[gd]
        elif el in metals :
            for al in metals :
                ax[1].cla()
                if el != al :
                    print(el,al)
                    j=np.where(els == al)[0][0]
                    for k in range(9) : 
                        plots.plotl(ax[1],x,mask[j].data[k,:])
                        #plots.plotl(ax[1],x,mask[j].data[k,:]/mask[i].data[k,:] )
                        bd.extend(np.where((mask[j].data[k,:]/mask[i].data[k,:] < 0.2) & (mask[j].data[k,:]<-0.01) )[0])
                    #gd=np.where(abs(mask[j].data) > 0.)[0]
                    #print(el,al,j,len(gd))
                    #mask[i].data[gd] = -1.*abs(mask[i].data[gd])
                    plt.show()

        new=np.zeros(mask[i].data.shape[-1])
        new[gd] = 1.
        new[bd] = -1*new[bd]
        ax[1].cla()
        plots.plotl(ax[1],x,new)
          
        plt.draw()
        plt.show()
        pdb.set_trace()

def filter_lines(infile,outfile,wind,nskip=0) :
    """ Read from input linelist file, output comments and lines falling in windows of [w1,w2] to outfile

    Args :
        infile (str) : name of input file
        outfile (str) : name of output file
        wind (list of [w1,w2] pairs) : list of window ranges
        nskip (int) : number of header lines to skip (default=0)

    Returns :
        nout (int) : number of windows
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
 
def mini_linelist(elem,linelist,maskdir=None,only=False,clobber=False) :
    """ Produce abbreviated Turbospec linelists, e.g. for minigrid construction, given mask file and linelist file IN AIR
        With only, produce linelist with only lines from input element 
        Return arrays of wavelength ranges wind,wair
    """

    # get window ranges in vacuum and convert to air
    if maskdir is not None :
        wind=np.loadtxt(os.environ['APOGEE_DIR']+'/data/windows/'+maskdir+'/'+elem+'.wave')
        nwind=wind.shape[0]
        wair=spectra.vactoair(wind)
    else : 
        nwind = 1
        wair=np.zeros([2,2])
        wair[0,0] = -1.
        wair[0,1] = 1.e10
        wind = wair

    # setup output directory
    if only :
        outdir = os.environ['APOGEE_SPECLIB']+'/linelists/'+elem+'_only/'
    else :
        outdir = os.environ['APOGEE_SPECLIB']+'/linelists/'+elem+'/'
    try: os.mkdir(outdir)
    except: pass

    # if files are being create by other process, wait until done
    while os.path.isfile(outdir+elem+'.lock') : 
        print('waiting for lock: ',outdir+elem+'.lock')
        time.sleep(10)

    # if files are already created, return, otherwise open .lock file and create
    if not clobber and os.path.isfile(outdir+elem+'.done') : return wind,wair
    fp = open(outdir+elem+'.lock','w')
    fp.close()

    # convert Turbospectrum files to filtered Turbospectrum files
    # Turbospectrum files are in air wavelengths
    lists=['turbospec.'+linelist+'.atoms','turbospec.'+linelist+'.molec',
           'turbospec.'+linelist+'.Hlinedata','turbospec.h2o-BC8.5V.molec','turbospec.h2o-BC9.5V.molec']
    for i,linelist in enumerate(lists) :
        filepath=os.environ['APOGEE_SPECLIB']+'/linelists/'+linelist
        fout=open(outdir+linelist,'w')
        with open(filepath) as fp:  
            out = ''
            nelem = 0
            n = 0
            line = fp.readline()
            while line :
                if line[0] == "'" :
                    # we have a new element
                    if nelem > 0 :
                        # if it's not the first element, write out the previous one!
                        if n > 0 :
                            if 'molec' in linelist :
                                tmp = int(float(head.split("'")[1]))
                                elemcode = [tmp//100,tmp%100]
                            else :
                                elemcode = [int(float(head.split("'")[1]))]
                            if not only or atomic.periodic(elem) in elemcode :
                                j=head.split("'")[2].split()[0]
                                # for the header line, include the new number of lines
                                fout.write("'"+head.split("'")[1]+"'   "+j+'{:10d}\n'.format(n))
                                # write the accumlated data output
                                fout.write(out)
                            n=0
                    head = line
                    # start the line data output with the comment line
                    out = fp.readline()
                    nelem += 1
                else :
                    # accumulate the linelist for this element if it's within the desired range
                    w = line.split()[0]
                    for i in range(nwind) :
                      if (float(w) >= wair[i,0]) and (float(w) <=wair[i,1]) : 
                          out=out+line
                          n+=1
                line = fp.readline()          
            if n > 0 :
                # last element
                if 'molec' in linelist :
                    tmp = int(float(head.split("'")[1]))
                    elemcode = [tmp//100,tmp%100]
                else :
                    elemcode = [int(float(head.split("'")[1]))]
                if not only or atomic.periodic(elem) in elemcode :
                    j=head.split("'")[2].split()[0]
                    fout.write("'"+head.split("'")[1]+"'   "+j+'{:10d}\n'.format(n))
                    fout.write(out)
        fout.close()

    # write .done file and remove .lock
    fp = open(outdir+elem+'.done','w')
    fp.close()
    os.remove(outdir+elem+'.lock') 
    return wind,wair

def plotcross(a,val=[0,0,0],hard=None,sum=True) :
    """ plot cross sections of input 3D grid of synthetic spectra

        Args:
            a : input HDU with synthetic spectral grid
           val : 3 ([M/H], logg, Teff) grid indices to use when varying other dimensions (default=[0,0,0])
    """
    dim=a.data.shape
    # if we have a rotation dimension, take the lowest rotation
    if len(dim) == 5 : 
        data=np.squeeze(a.data[0,:,:,:,:])
        dim=data.shape
    else : data = a.data
    n=dim[3]
    x=10.**spectra.fits2vector(a.header,1)
    mh=spectra.fits2vector(a.header,4)
    logg=spectra.fits2vector(a.header,3)
    teff=spectra.fits2vector(a.header,2)
    colors=['r','g','b','c','m','y','black']*3
    for i in range(dim[0]) :
        if i == 0 : fig,ax=aspcap.plot(x,data[i,val[1],val[2],:],color=colors[i],sum=sum)
        else : aspcap.plot(x,data[i,val[1],val[2],:],ax=ax,color=colors[i],sum=sum)
    fig.suptitle('[M/H] varied from {:6.2} to {:6.2f} at logg {:6.1f}, Teff {:6.0f}'.format(mh[0],mh[-1],logg[val[1]],teff[val[2]]))
    if hard is not None : 
        fig.savefig(hard+'_mh.pdf')
        plt.close()
    for i in range(dim[1]) :
        if i == 0 : fig,ax=aspcap.plot(x,data[val[0],i,val[2],:],color=colors[i],sum=sum)
        else : aspcap.plot(x,data[val[0],i,val[2],:],ax=ax,color=colors[i],sum=sum)
    fig.suptitle('log g varied from {:6.1} to {:6.2f} at [M/H] {:6.2f}, Teff {:6.0f}'.format(logg[0],logg[-1],mh[val[0]],teff[val[2]]))
    if hard is not None : 
        fig.savefig(hard+'_logg.pdf')
        plt.close()
    for i in range(dim[2]) :
        if i == 0 : fig,ax=aspcap.plot(x,data[val[0],val[1],i,:],color=colors[i],sum=sum)
        else : aspcap.plot(x,data[val[0],val[1],i,:],ax=ax,color=colors[i],sum=sum)
    fig.suptitle('Teff varied from {:6.0f} to {:6.0f} at logg {:6.1f}, [M/H] {:6.2f}'.format(teff[0],teff[-1],logg[val[1]],mh[val[0]]))
    if hard is not None : 
        fig.savefig(hard+'_teff.pdf')
        plt.close()

    dim=a.data.shape
    if len(dim) == 5 :
        for i in range(dim[0]) :
            if i == 0 : fig,ax=aspcap.plot(x,data[i,val[0],val[1],val[2],:],color=colors[i],sum=sum)
            else : aspcap.plot(x,data[i,val[0],val[1],val[2],:],ax=ax,color=colors[i],sum=sum)
        #fig.suptitle('vsini varied from {:6.0f} to {:6.0f} at logg {:6.1f}, [M/H] {:6.2f}'.format(teff[0],teff[-1],logg[val[1]],mh[val[0]]))
        if hard is not None : 
            fig.savefig(hard+'_vsini.pdf')
            plt.close()

