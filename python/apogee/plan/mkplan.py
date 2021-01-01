import copy
import numpy as np
import os
import glob
import pdb
import subprocess
import yaml
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from sdss import yanny
from apogee.speclib import atmos
from apogee.utils import spectra
from apogee.plan import mkslurm

def mkgriddirs(configfile,nosynth=False,synthonly=False,writeraw=False,queryport=1052,digits=3,py2=False) :
    """ Script to create output directories and plan and batch queue files for all grids listed in master grid configuration file
    """

    # Read grid configuration file
    if not os.path.isfile(configfile+'.yml'):
        print('{:s} does not exist'.format(configfile+'.yml'))
        return
    p=yaml.safe_load(open(configfile+'.yml','r'))

    # loop over each grid
    for i in range(len(p['GRID'])) :

      # do both "raw" directory and final directory: former may be repeated!
      specdir=p['GRID'][i]['specdir']
      smooth=p['GRID'][i]['smooth']
      synthcode=p['GRID'][i]['synthcode']
      atmos=p['GRID'][i]['atmos']
      if synthonly : names = [ specdir ]
      elif nosynth : names = [ specdir+'_'+smooth ]
      else : names = [ specdir+'_'+smooth, specdir ]
      for igrid,name in enumerate(names) :
        # construct name and create output directory

        if abs(p['GRID'][i]['solarisotopes']) == 1 :
            iso = 'solarisotopes'
        elif abs(p['GRID'][i]['solarisotopes']) == 2 :
            iso = 'giantisotopes'
        elif p['GRID'][i]['solarisotopes'] < 0 :
            iso = 'tests/'+iso
        dir = os.getenv('APOGEE_SPECLIB')+'/synth/'+synthcode.strip("'")+'/'+atmos+'/'+iso+'/'+name+'/plan/'
        print(dir)
        try: os.makedirs(dir)
        except: pass

        # remove any old plan files
        os.chdir(dir)
        for filePath in glob.glob("*.yml"):
            if os.path.isfile(filePath): os.remove(filePath)

        # move GRID keys up one level
        #out = copy.deepcopy(p)
        out = p['GRID'][i]
        out['name'] = name
        #for key in out['GRID'][i].keys() : out[key] = out['GRID'][i][key]
        #out.pop('GRID')

        fp = open(dir+name+'.yml','w')
        fp.write(yaml.dump(out,sort_keys=False))
        fp.close()
        for elem in p['GRID'][i]['elem'] : speclib_split(dir+name,el=elem,digits=digits,py2=py2)

        # make pbs scripts
        os.chdir('..')
        specdir = synthcode.strip("'")+'/'+atmos+'/'+iso+'/'+name

        if name == p['GRID'][i]['specdir'] :
            speclib_split(dir+name,amsplit=False,digits=digits,py2=py2)
            mkslurm.write('mkgrid plan/'+name+'_a[mp]*vp20.yml plan/'+name+'_a[mp]*vp48.yml plan/'+name+'_a[mp]*vp??.yml',queryhost=os.uname()[1],queryport=queryport,maxrun=32)
            mkslurm.write('mkrbf plan/'+name+'_c[mp]*vp??.yml',queryhost=os.uname()[1],queryport=queryport,maxrun=1,time='72:00:00')
            mkslurm.write('mkrbf --nofill plan/'+name+'.yml',name='mkrbfholes',runplans=False,time='72:00:00')
        else :
            if writeraw : raw = '--writeraw'
            else : raw = ''
            mkslurm.write('mkgridlsf plan/'+name+'_a[mp]*vp??.yml',queryhost=os.uname()[1],queryport=queryport,maxrun=12,time='24:00:00')
            #mkslurm.write('bundle plan/'+name+'_??.yml',queryhost=os.uname()[1],queryport=queryport,maxrun=32)
            mkslurm.write('"pca --pcas 12 75 --incremental --threads 0" '+raw+' plan/'+name+'.yml',maxrun=1,time='72:00:00',queryhost=os.uname()[1],queryport=queryport)
            mkslurm.write('mkgridlsf plan/'+name+'_a[mp]*vp??.yml',queryhost=os.uname()[1],queryport=queryport,maxrun=12,time='72:00:00',
                          postcmd='pca --pcas 12 75 --incremental --threads 0 '+raw+' plan/'+name+'.yml',name='mkgridlsf_pca')

def speclib_split(planfile,amsplit=True,cmsplit=True,nmsplit=True,oasplit=True,vtsplit=True,el='',digits=3,py2=False) :
    """ Make a bunch of individual plan files from master, splitting [alpha/M],[C/M],[N/M],vt
    """
    # read master plan file
    print('splitting: ', planfile)
    p=yaml.safe_load(open(planfile+'.yml','r'))

    # some cards removed in split par files
    p.pop('npart',None)
    p.pop('npca',None)
    p.pop('vmsuffix',None)
    if el is not '' : p['elem'] = el 
    else: p.pop('elem',None)

    # make specdir the full path relative to $APOGEE_SPECLIB/synth
    if int(p['solarisotopes']) == 1 : isodir='solarisotopes'
    elif int(p['solarisotopes']) == 2 : isodir='giantisotopes'

    #for key in ['synthcode','atmos','specdir','linelist','config'] : 
    #    p[key] = p[key].decode().strip("'")
 
    p['specdir'] = p['synthcode']+'/'+p['atmos']+'/'+isodir+'/'+p['specdir']

    # get ranges in [alpha/M], [C/M], [N/M], and vt
    if amsplit :
        amrange=spectra.vector(p['am0'],p['dam'],p['nam'])
        p['nam'] = 1
    else :
        amrange = [0]
    if cmsplit :
        cmrange=spectra.vector(p['cm0'],p['dcm'],p['ncm'])
        p['ncm'] = 1
    else :
        cmrange = [0]
    if nmsplit :
        nmrange=spectra.vector(p['nm0'],p['dnm'],p['nnm'])
        p['nnm'] = 1
    else :
        nmrange = [0]
    if oasplit :
        try :
            oarange=spectra.vector(p['oa0'],p['doa'],p['noa'])
        except :
            oasplit=False
            oarange=[0.]
        p['noa'] = 1
    else :
        oarange = [0]
    if int(p['vmicrofit']) == 0 :
        vtrange=spectra.vector(p['vt0'],p['dvt'],p['nvt'])
        p.pop('vt0') 
        p.pop('dvt') 
        p.pop('nvt')
    else :
        vtrange = [0]
        vtsplit = False

    # loop through all and make individual plan files
    dw = float(p['dw'])
    for am in amrange :
        if amsplit : p['am0'] = float(am)
        for cm in cmrange :
            if cmsplit : p['cm0'] = float(cm)
            for nm in nmrange :
                if nmsplit : p['nm0'] = float(nm)
                for oa in oarange :
                    if oasplit : p['oa0'] = float(oa)
                    for vt in vtrange :
                        # vmicro handled differently
                        if int(p['vmicrofit']) == 0 : 
                            p['vmicro'] = [float(10.**vt)]
                            # special handling for dw
                            if np.isclose(dw,-1.) :
                                if p['vmicro'] < 3.99 : p['dw'] = 0.05
                                else : p['dw'] = 0.10   
    
                        suffix=''
                        if amsplit : suffix+='a'+atmos.cval(am,digits=digits,py2=py2)
                        if cmsplit : suffix+='c'+atmos.cval(cm,digits=digits,py2=py2)
                        if nmsplit : suffix+='n'+atmos.cval(nm,digits=digits,py2=py2)
                        if oasplit : suffix+='o'+atmos.cval(oa,digits=digits,py2=py2)
                        if vtsplit : suffix+='v'+atmos.cval(10.**vt)
                        p['name'] = suffix

                        with open(planfile+'_'+suffix+el+'.yml', 'w') as fp:
                            fp.write(yaml.dump(p,sort_keys=False,Dumper=Dumper))

def aspcap(field,apred='r13',telescope='apo25m',aspcap_vers='l33',aspcap_config='l33cnmask',ncpus=16, minmjdlast=None) :

    plan={}
    plan['apogee_ver'] = os.environ['APOGEE_VER']
    plan['apvisit'] = 0
    plan['apred_vers'] = apred
    plan['telescope'] = telescope
    if telescope == 'lco25m' : plan['instrument'] = 'apogee-s'
    else : plan['instrument'] = 'apogee-n'
    plan['apstar_vers'] = 'stars'
    plan['aspcap_vers'] = aspcap_vers
    plan['aspcap_config'] = aspcap_config
    plan['ncpus'] =  ncpus
    plan['queue'] =  0
    plan['qname'] =  'apogee'
    plan['qgroup'] =  'apogee'
    plan['caldir'] = 'cal'
    if minmjdlast is not None : plan['minmjdlast'] = 58814
    plan['field'] = field

    outdir=os.environ['APOGEE_ASPCAP']+'/'+apred+'/'+aspcap_vers+'/plan'
    os.makedirs(outdir,exist_ok=True)

    with open(outdir+'/'+field+'_'+telescope+'.yml', 'w') as fp:
        fp.write(yaml.dump(plan,sort_keys=False,Dumper=Dumper))

