#!/usr/bin/env python

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import argparse
import glob
import numpy as np
import os
import sys
import subprocess
import pdb
from apogee.plan import mkslurm
import yaml

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='Creates apStar plan files')

    parser.add_argument("config",type=str, help='configfile')
    parser.add_argument("--fields",type=str,nargs='+',help='list of fields',default=[])
    args=parser.parse_args()

    cfg = yaml.safe_load(open(args.config,'r'))

    outdir = os.environ['APOGEE_REDUX']+'/'+cfg['apred_vers']+'/'+cfg['apstar_vers']+'/plan'
    os.makedirs(outdir,exist_ok=True)

    for i,field in enumerate(args.fields) :
        comp=field.split('/')
        name=comp[-1]
        telescope=comp[-2]
        if telescope == 'lco25m' : prefix = 'as'
        else : prefix = 'ap'
        fp = open(outdir+'/{:s}Star-{:s}.yml'.format(prefix,name),'w')
        fp.write('---\n')
        fp.write('apogee_ver : {:s}\n'.format(cfg['apogee_ver']))
        fp.write('apred_vers : {:s}\n'.format(cfg['apred_vers']))
        fp.write('apstar_vers : {:s}\n'.format(cfg['apstar_vers']))
        fp.write('mjdstart : {:d}\n'.format(cfg['mjdstart']))
        fp.write('mjdend : {:d}\n'.format(cfg['mjdend']))
        fp.write('telescope : {:s}\n'.format(telescope))
        #fp.write('survey : {:s}'.format(survey))
        fp.write('field : {:s}\n'.format(name))
        fp.close()


    cmd = 'rv --clobber --threads=32'
    mkslurm.write('"'+cmd+'" plan/a?Star*.yml' ,maxrun=1,idlthreads=16,queryport=1051,queryhost=os.uname()[1],pythreads=1)
    #sort=np.argsort(nstars)[::-1]
    #fp=open(topdir+'/slurm/fields.sort','w')
    #for i in range(len(sort)) : 
    #    tel=args.fields[sort[i]].split('/')[0]
    #    field=args.fields[sort[i]].split('/')[1]
    #    fp.write('{:s}/plan/aspcapStar-{:s}.par {:d}\n'.format(tel,field,nstars[sort[i]]))
    #fp.close()
    #print('Modify slurm/aspcap to use fields.sort if desired...')

