# encoding: utf-8
#
# @Author: Jon Holtzman
# @Date: March 2018
# @Filename: synth.py
# @License: BSD 3-Clause
# @Copyright: Jon Holtzman

import os
import glob
import pdb
import shutil
import subprocess
import sys
import numpy as np
from apogee.speclib import atmos
from apogee.utils import atomic
from sdss.utilities import yanny
from astropy.io import ascii
from astropy.io import fits

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
        fp=open(infile+'.filled','r')
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

def mkturbospec(teff,logg,mh,am,cm,nm,wrange=[15100.,17000],dw=0.05,vmicro=2.0,solarisotopes=False,elem='',welem=None,
    atmod=None,kurucz=True,atmosroot=None,atmosdir=None,nskip=0,endskip=0,
    linelist='20150714',h2o=None,linedir=None,
    save=False,run=True,split=200) :

    """ Runs Turbospectrum for specified input parameters
    """

    # directory setup
    if atmosroot is None : atmosroot=os.environ['APOGEE_SPECLIB']+'/atmos/'
    if linedir is None : linelistdir=os.environ['APOGEE_SPECLIB']+'/linelists/' 
 
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
    workdir=(prefix+'m{:s}a{:s}c{:s}n{:s}v{:s}'+elem).format(
             atmos.cval(mh),atmos.cval(am),atmos.cval(cm),atmos.cval(nm),atmos.cval(vmicro))
    try: os.mkdir(workdir)
    except: pass
    root=workdir+'/'+(prefix+'t{:04d}g{:s}m{:s}a{:s}c{:s}n{:s}v{:s}'+elem).format(teff, atmos.cval(logg), 
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
        marcs2turbo(atmod,workdir+'/'+os.path.basename(atmod),trim=nskip )

    # Turbospectrum setup
    try: os.symlink(os.environ['APOGEE_DIR']+'/src/turbospec/DATA',workdir+'/DATA')
    except: pass

    # individual element grid?
    if elem != '' :
        elemnum=atomic.periodic(elem)
        eabun=np.arange(atomic.solar(elem)-0.75,0.25,10)
        nelem=1
        linelistdir=linelistdir+elem+'/'
    else :
        eabun=np.array([0.])
        nelem=0

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
        nrange=1

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
            fout.write("'INDIVIDUAL ABUNDANCES:'  '3'\n")
            fout.write("    6  {:8.3f}\n".format(8.39+mh+cm))
            fout.write("    7  {:8.3f}\n".format(7.78+mh+nm))
            fout.write("    8  {:8.3f}\n".format(8.66+mh+am))
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
                os.chmod(root+'_babsma.csh', 0777)
                os.chdir(workdir)
                subprocess.call(['time','./'+os.path.basename(root)+'_babsma.csh'],stdout=stdout)
                os.chdir('..')

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
        fout.write("'INDIVIDUAL ABUNDANCES:'  '3'\n")
        fout.write("    6  {:8.3f}\n".format(8.39+mh+cm))
        fout.write("    7  {:8.3f}\n".format(7.78+mh+nm))
        fout.write("    8  {:8.3f}\n".format(8.66+mh+am))
        if elem != '' : fout.write("{:6d}  {:8.3f}\n".format(elemnum,abun+mh))
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
        n_HI = len(open(linelistdir+'turbospec.'+linelist+'.Hlinedata').readlines())
        if n_HI < 3 : nlists-=1
        # if we are using H2O, add that list
        if h2o > 0 : nlists+=1
        fout.write("'NFILES:'  '{:4d}'\n".format(nlists))
        if n_HI >= 3 : fout.write(linelistdir+'turbospec.'+linelist+'.Hlinedata\n')
        fout.write(linelistdir+'turbospec.'+linelist+'.atoms\n')
        fout.write(linelistdir+'turbospec.'+linelist+'.molec\n')
        if h2o == 1 :
            fout.write(linelistdir+'turbospec.h2o-BC8.5V'+'.molec\n')
        elif h2o == 2 :
            fout.write(linelistdir+'turbospec.h2o-BC9.5V'+'.molec\n')
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
        fout.write("set runtime = `ps -q $bsynjob -o cputime | tail -1 | awk -F: '{print ($1*3600)+($2*60)+$3}'`\n")
        tmax=120*int(.05/min([0.05,dw]))
        if h2o == 1 : tmax*=2
        if h2o == 2 : tmax*=5
        fout.write('while ( $runtime < {:d} && $ok == 0 )\n'.format(tmax))
        fout.write('  usleep 200000\n')
        fout.write("  set runtime = `ps -q $bsynjob -o cputime | tail -1 | awk -F: '{print ($1*3600)+($2*60)+$3}'`\n")
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
            os.chmod(root+'_bsyn.csh', 0777)
            os.chdir(workdir)
            subprocess.call(['time','./'+os.path.basename(root)+'_bsyn.csh'],stdout=stdout)
            os.chdir('..')
            try:
                out=ascii.read(file)
                if ielem == 0 : 
                    spec=out['col3']
                else :
                    spec=np.vstack(spec,out['col3'])
            except :
                return 0.

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
    vmicrofit = p['vmicrofit'] if p.get('vmicrofit') else 0
    vmicro = p['vmicro'] if p.get('vmicro') else 0
    vmacrofit = p['vmacrofit'] if p.get('vmacrofit') else 0
    vmacro = p['vmacro'] if p.get('vmacro') else 0
    specdir = os.environ['APOGEE_SPECLIB']+'/synth/'+p['specdir'] if p.get('specdir') else './'
    linelistdir = None
    linelist = p['linelist'] if p.get('linelist') else None

    if elem == '' :
        nelem=1

    # wavelength array
    nspec=int((wrange[1]-wrange[0])/dw)+1
    rawwave=wrange[0]+np.arange(nspec)*dw

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
                    elem=elem,linedir=linelistdir,linelist=linelist,vmicro=vout,
                    solarisotopes=solarisotopes,
                    nskip=nskip,kurucz=kurucz,run=run,save=save,split=split) 
	          nskip = nskip+dskip if isinstance(spec,float) else -1
                specdata[:,imh,ilogg,iteff,:]=spec

          # FITS header and output
          hdu=fits.PrimaryHDU(np.squeeze(specdata))
          idim=1
          add_dim(hdu.header,rawwave[0],rawwave[1]-rawwave[0],1,'WAVELENGTH',idim)
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
        subprocess.call(['mkslurm','mkgrid','"plan/'+name+'_a[mp]*vp??.par"'])
        subprocess.call(['mkslurm','bundle','"plan/'+name+'_??.par"'])

def mini_linelist(elem,linelist,maskdir) :
    """ Produce an abbreviated line list for minigrid construction given mask file
    """

    wind=ascii.read(os.environ['APOGEE_DIR']+'/lib/'+maskdir+'/'+elem+'./wave',names=['w1','w2'])
"""    
pro speclib_wline,elem,linelist,maskdir,wair=wair,wvac=wvac
;
; takes a linelist and element, and produces modified linelists only with lines in the
; windows for that element in a subdirectory of the linelist directory
; handles master linelist + 2 H2O linelists + HI linelist
;
;elem='Ce'
;linelist='linelist.20170418'
readcol,getenv('SPECLIB_DIR')+'/lib/'+maskdir+'/'+elem+'.wave',w1,w2,format='(d,d)'
vactoair,w1,w1air
vactoair,w2,w2air
wair=[]
wvac=[]
awk='$1<0'
tawk='NR>2&&$1<0'
for i=0,n_elements(w1)-1 do begin
  wair=[[wair],[w1air[i],w2air[i]]]
  wvac=[[wvac],[w1[i],w2[i]]]
  awk=awk+'||($1>='+string(w1air[i]/10.,format='(f9.4)')+'&&$1<='+string(w2air[i]/10.,format='(f9.4)')+')'
  tawk=tawk+'||($1>='+string(w1air[i],format='(f10.4)')+'&&$1<='+string(w2air[i],format='(f10.4)')+')'
endfor
linelistdir=getenv('APOGEE_SPECLIB')+'/linelists/'
file_mkdir,elem
openw,lun,elem+'.csh',/get_lun
printf,lun,'#!/bin/csh -f'
printf,lun,'awk '+"'"+awk+"' "+linelistdir+linelist+' >'+elem+'/'+linelist
printf,lun,'turboscript '+elem+'/'+linelist

lists=['turbospec.20170418.Hlinedata','turbospec.h2o-BC8.5V.molec','turbospec.h2o-BC9.5V.molec']
code=['01.000000','010108.000000000','010108.00000000']
comment=['HI culled','Barber culled','Barber culled']
for ilist=0,n_elements(lists)-1 do begin
  printf,lun,'awk '+"'"+tawk+"' "+linelistdir+lists[ilist]+' >'+elem+'/'+lists[ilist]+'.tmp'
  printf,lun,"set n=`wc -l "+elem+"/"+lists[ilist]+".tmp | awk '{print $1}'`"
  printf,lun,"echo \'"+code[ilist]+"                 \'    1   $n>"+elem+'/'+lists[ilist]
  printf,lun,"echo \'"+comment[ilist]+"\'>>"+elem+'/'+lists[ilist]
  printf,lun,'cat '+elem+'/'+lists[ilist]+'.tmp >>'+elem+'/'+lists[ilist]
endfor
free_lun,lun
spawn,'csh '+elem+'.csh'
end
"""
