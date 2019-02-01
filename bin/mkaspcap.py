#!/usr/bin/env python

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import argparse
import os
import sys
import subprocess
import pdb
from apogee.plan import mkslurm

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='Creates ASPCAP directories and plan files')

    parser.add_argument("apred",type=str, help='apred version')
    parser.add_argument("aspcap",type=str,help='aspcap version')
    parser.add_argument("config",type=str,help='aspcap configuration')
    parser.add_argument("--apstar",type=str, help='apred version',default='stars')
    parser.add_argument("--queue",type=int,default=0)
    parser.add_argument("--ncpus",type=int,default=16)
    parser.add_argument("--noplot",type=int,default=0)
    parser.add_argument("--noelem",type=int,default=0)
    parser.add_argument("--nstars",type=int,default=0)
    parser.add_argument("--commiss",type=int,default=0)
    parser.add_argument("--nored",type=int,default=0)
    parser.add_argument("--visits",type=int,default=0)
    parser.add_argument("--caldir",type=str,default='0')
    parser.add_argument("--npar",type=int,default=0)
    parser.add_argument("--renorm",type=int,default=0)
    parser.add_argument("--maxwind",type=int,default=0)
    parser.add_argument("--fields",type=str,nargs='+',help='list of fields')
    args=parser.parse_args()

    for field in args.fields :
        print('field: ', field)
        cmd=["idl","-e","aspcap_mkplan,'"+field+"'"+
             ",apred_vers='{:s}'".format(args.apred)+
             ",apstar_vers='{:s}'".format(args.apstar)+
             ",aspcap_vers='{:s}'".format(args.aspcap)+
             ",aspcap_config='{:s}'".format(args.config)+
             ",ncpus={:d}".format(args.ncpus)+
             ",queue={:d}".format(args.queue)+
             ",noplot={:d}".format(args.noplot)+
             ",caldir='{:s}'".format(args.caldir)+
             ",noelem={:d}".format(args.noelem)+
             ",nstars={:d}".format(args.nstars)+
             ",commiss={:d}".format(args.commiss)+
             ",nored={:d}".format(args.nored)+
             ",visits={:d}".format(args.visits)+
             ",npar={:d}".format(args.npar)+
             ",renorm={:d}".format(args.renorm)+
             ",maxwind={:d}".format(args.maxwind)]
        print(cmd)
        subprocess.call(cmd,shell=False)
    for inst in ['apogee-n','apogee-s'] :
        outdir=os.environ['APOGEE_ASPCAP']+'/'+args.apred+'/'+args.aspcap+'/config/'+inst+'/'
        cmd=["idl","-e","aspcap_mklib,'"+args.config+"'"+
             ",outdir='{:s}'".format(outdir)+
             ",apred='{:s}'".format(args.apred)+
             ",renorm={:d}".format(args.renorm)+
             ",maxwind={:d}".format(args.maxwind)+
             ",inst='"+inst+"'"]
        print(cmd)
        subprocess.call(cmd,shell=False)

        if args.maxwind > 0 :
            cmd=["idl","-e","readcol,'"+outdir+"/elem.list"+"'"+
                 ",format='(a)',els,stringskip='#',/silent","&","print,els"]
            print(cmd)
            subprocess.call(cmd,shell=False)

        f=open(outdir+'/done','w')
        f.close()
    topdir=os.environ['APOGEE_ASPCAP']+'/'+args.apred+'/'+args.aspcap
    os.chdir(topdir)
    cmd='aspcap '
    if args.noelem != 0 : cmd+=' --noelem'
    mkslurm.write('"'+cmd+'" apo*/plan/aspcapStar*.par lco*/plan/aspcapStar*.par',maxrun=2,idlthreads=16,queryport=1051,queryhost=os.uname()[1])

#r=aspcap_root+apred_vers+'/'+aspcap_vers+'/config/'+dirs.instrument+'/'
#file_mkdir,configdir
#while file_test(configdir+'lockfile') do apwait,configdir+'lockfile',10
#if ~file_test(configdir+'/done') then begin
#  openw,lock,/get_lun,configdir+'/lockfile'
#  free_lun,lock
#  aspcap_mklib,aspcap_config,outdir=configdir,maskdir=getenv('APOGEE_DIR')+'/data/windows/filters_26042016/'
#  ; create individual windows if requested
#  if keyword_set(maxwind) then  begin
#    readcol,configdir+'/elem.list',format='(a)',els,stringskip='#',/silent
#    for i=0,n_elements(els)-1 do n=filtsplit(els[i],maxwind=maxwind,outdir=configdir)
#  endif
#  openw,done,/get_lun,configdir+'/done'
#  free_lun,done
#  file_delete,configdir+'/lockfile'
#endif


#cd $APOGEE_ASPCAP/$apred_vers/$aspcap_vers
#mkslurm "aspcap apo*/plan/aspcapStar*.par lco*/plan/aspcapStar*.par" --maxrun=2 --idlthreads=16 --queryport=1051
#
#echo "\n\nREMEMBER to modify the aspcapStar*.par files if non-default options,"
##echo "  e.g., noplot, noelem, fits, are desired"
