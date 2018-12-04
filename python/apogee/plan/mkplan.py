import numpy as np
import os
import glob
import pdb
import subprocess
from sdss import yanny
from apogee.speclib import atmos
from apogee.utils import spectra
from apogee.plan import mkslurm


def mkgriddirs(configfile) :
    """ Script to create output directories and plan and batch queue files for all grids listed in master grid configuration file
        Calls IDL routine to make the individual subplan files
    """

    # Read grid configuration file
    if not os.path.isfile(configfile):
        print('{:s} does not exist'.format(configfile))
        return
    p=yanny.yanny(configfile,np=True)

    # loop over each grid
    for i in range(len(p['GRID']['specdir'])) :

      # do both "raw" directory and final directory: former may be repeated!
      for igrid,name in enumerate([ p['GRID']['specdir'][i]+'_'+p['GRID']['smooth'][i], p['GRID']['specdir'][i] ]) :
        # construct name and create output directory
        #name = p['GRID']['specdir'][i]+'_'+p['GRID']['smooth'][i]

        if abs(p['GRID']['solarisotopes'][i]) == 1 :
            iso = 'solarisotopes'
        elif abs(p['GRID']['solarisotopes'][i]) == 2 :
            iso = 'giantisotopes'
        elif p['GRID']['solarisotopes'][i] < 0 :
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
                f.write('{:30s}{:30s}\n'.format(key, p[key].strip("'")))
        for key in p['GRID'].dtype.names :
            if igrid == 0 or key != 'smooth' :
              out='{:30s}'.format(str(p['GRID'][key][i])).strip("'").strip('[]')
              out='{:30s}{:30s}\n'.format(key,out.replace(']',''))
              f.write(out.strip('[]'))
        f.close()

        # make all of the individual planfiles from the master planfile
        speclib_split(dir+name)
        #subprocess.call(['idl','-e',"speclib_allplan,'"+name+".par'"])

        # make pbs scripts
        os.chdir('..')
        specdir = p['synthcode'].strip("'")+'/'+p['GRID']['atmos'][i]+'/'+iso+'/'+name
        #os.environ['NO_NODES'] = 'yes'
        #subprocess.call(['mkslurm.csh','mkgrid','"plan/'+name+'_a[mp]*vp20.par"','"plan/'+name+'_a[mp]*vp48.par"','"plan/'+name+'_a[mp]*vp??.par"'],shell=False)
        #subprocess.call(['mkslurm.csh','mkrbf','"plan/'+name+'_c[mp]*vp??.par"'],shell=False)
        #subprocess.call(['mkslurm.csh','mkgridlsf','"plan/'+name+'_a[mp]*vp??.par"'],shell=False)
        #subprocess.call(['mkslurm.csh','bundle','"plan/'+name+'_??.par"'],shell=False)

        if name == p['GRID']['specdir'][i] :
            speclib_split(dir+name,amsplit=False)
            mkslurm.write('mkgrid plan/'+name+'_a[mp]*vp20.par plan/'+name+'_a[mp]vp48.par plan/'+name+'_a[mp]*vp??.par',queryhost=os.uname()[1],queryport=1052,maxrun=32)
            mkslurm.write('mkrbf plan/'+name+'_c[mp]*vp??.par',queryhost=os.uname()[1],queryport=1052,maxrun=1,time='72:00:00')
            mkslurm.write('mkrbf --nofill plan/'+name+'.par',name='mkrbfholes',runplans=False,time='72:00:00')
        else :
            mkslurm.write('mkgridlsf plan/'+name+'_a[mp]*vp??.par',queryhost=os.uname()[1],queryport=1052,maxrun=12,time='24:00:00')
            #mkslurm.write('bundle plan/'+name+'_??.par',queryhost=os.uname()[1],queryport=1052,maxrun=32)
            mkslurm.write('pca --pcas 12 75 --incremental --threads 0 --writeraw plan/'+name+'.par',runplans=False,time='72:00:00')

def speclib_split(planfile,amsplit=True,cmsplit=True,nmsplit=True,vtsplit=True,el=None) :
    """ Make a bunch of individual plan files from master, splitting [alpha/M],[C/M],[N/M],vt
    """
    # read master plan file
    p=yanny.yanny(planfile+'.par')

    # some cards removed in split par files
    p.pop('npart',None)
    p.pop('npca',None)
    p.pop('vmsuffix',None)
    if el is not None : p['el'] = el 
    else: p.pop('elem',None)

    # make specdir the full path relative to $APOGEE_SPECLIB/synth
    if int(p['solarisotopes']) == 1 : isodir='solarisotopes'
    elif int(p['solarisotopes']) == 2 : isodir='giantisotopes'

    for key in ['synthcode','atmos','specdir','linelist','config'] : 
        p[key] = p[key].strip("'")
 
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
    if vtsplit :
        vtrange=spectra.vector(p['vt0'],p['dvt'],p['nvt'])
        p.pop('vt0') 
        p.pop('dvt') 
        p.pop('nvt')
    else :
        vtrange = [0]

    # loop through all and make individual plan files
    dw = float(p['dw'])
    for am in amrange :
        if amsplit : p['am0'] = am
        for cm in cmrange :
            if cmsplit : p['cm0'] = cm
            for nm in nmrange :
                if nmsplit : p['nm0'] = nm
                for vt in vtrange :
                    # vmicro handled differently
                    if vtsplit : 
                        p['vmicrofit'] = 0
                        p['vmicro'] = 10.**vt
                        # special handling for dw
                        if np.isclose(dw,-1.) :
                            if p['vmicro'] < 3.99 : p['dw'] = 0.05
                            else : p['dw'] = 0.10   
    
                    suffix=''
                    if amsplit : suffix+='a'+atmos.cval(am)
                    if cmsplit : suffix+='c'+atmos.cval(cm)
                    if nmsplit : suffix+='n'+atmos.cval(nm)
                    if vtsplit : suffix+='v'+atmos.cval(10.**vt)
                    p['name'] = suffix
                    p.write(planfile+'_'+suffix+'.par')

