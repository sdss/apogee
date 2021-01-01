#!/usr/bin/env python
import sys
import os
import argparse

def write(cmd,outdir='slurm/',cwd=None,queryhost=None,queryport=None,maxrun=None,idlthreads=1,runplans=True,
          time='240:00:00',name=None,fast=False,tacc=False, notchpeak=False, np=False, postcmd=None, flag=None, pythreads=None) :
    """ Routine to write a slurm file with specified input command and parameteres
    """

    try :
        os.mkdir(outdir)
    except:
        pass

    if name is None : name=cmd.split()[0].strip('"')
    file = outdir+name
    f=open(file,'w')

    f.write('#!/bin/csh\n')
    if tacc :
        f.write('#SBATCH --partition=normal\n')
        f.write('#SBATCH --ntasks-per-node=24\n')
    elif notchpeak :
        f.write('#SBATCH --account=sdss\n')
        f.write('#SBATCH --partition=notchpeak\n')
        f.write('#SBATCH --ntasks=16\n')
        f.write('#SBATCH -C rom\n')
        f.write('#SBATCH --mem=48G\n')
        f.write('#SBATCH --time=48:00:00\n')
        f.write('#SBATCH --nodes=1\n')
    elif np :
        f.write('#SBATCH --account=sdss-np\n')
        f.write('#SBATCH --partition=sdss-np\n')
        f.write('#SBATCH --ntasks=64\n')
        f.write('#SBATCH --time='+time+'\n')
        f.write('#SBATCH --nodes=1\n')
    else :
        if fast : f.write('#SBATCH --account=sdss-kp-fast\n')
        else : f.write('#SBATCH --account=sdss-kp\n')
        f.write('#SBATCH --partition=sdss-kp\n')
        f.write('#SBATCH --ntasks=16\n')
        f.write('#SBATCH --time='+time+'\n')
        f.write('#SBATCH --nodes=1\n')
    f.write('#SBATCH -o '+os.path.basename(file)+'.out\n' )
    f.write('#SBATCH -e '+os.path.basename(file)+'.out\n' )
    if runplans :
        f.write('setenv QUERYHOST '+queryhost+'\n' )
        f.write('setenv QUERYPORT '+str(queryport)+'\n')
        f.write('setenv APOGEE_MAXRUN '+str(maxrun)+'\n' )
        if flag is not None :f.write('setenv APOGEE_FLAG '+flag+'\n' )
    f.write('setenv IDL_CPU_TPOOL_NTHREADS '+str(idlthreads)+'\n' )
    if pythreads is not None :
        f.write('setenv OMP_NUM_THREADS '+str(pythreads)+'\n')
        f.write('setenv OPENBLAS_NUM_THREADS '+str(pythreads)+'\n')
        f.write('setenv MKL_NUM_THREADS '+str(pythreads)+'\n')
        f.write('setenv VECLIB_MAXIMUM_THREADS '+str(pythreads)+'\n')
        f.write('setenv NUMEXPR_NUM_THREADS '+str(pythreads)+'\n')

    if cwd is None : cwd=os.getcwd()
    f.write('cd '+cwd+'\n' )
    if runplans : 
        f.write('runplans ' + cmd+'\n' )
    elif isinstance(cmd,list) :
        for c in cmd : f.write(c+'\n' )
    else : 
        f.write(cmd+'\n' )
    f.write('wait\n' )
    if postcmd is not None : f.write(postcmd+'\n')
    f.write('echo DONE\n' )
    f.close()
    os.chmod(file,0o770)

def main(args) :
    """ Main routine to create SLURM file from command line
    """
    parser = argparse.ArgumentParser(
        prog=os.path.basename(args[0]),
        description='Makes a SLURM batch file')
    parser.add_argument('cmd', type=str, help='cmd')
    parser.add_argument('--outdir', type=str, help='directory',default='slurm/')
    parser.add_argument('--name', type=str, help='output file name',default=None)
    parser.add_argument('--norunplans', help='directory',action="store_false")
    parser.add_argument('--queryport', type=int, help='port to use for queue manager',default=1050)
    parser.add_argument('--queryhost', type=str, help='host to use for queue manager',default=os.uname()[1])
    parser.add_argument('--maxrun', type=int, help='maximum jobs to run at a time',default=1)
    parser.add_argument('--flag', type=str, help='value for APOGEE_FLAG',default='1111111')
    parser.add_argument('--time', type=str, help='maximum wall clock time',default='240:00:00')
    parser.add_argument('--idlthreads', type=int, help='maximum IDL threads',default=1)
    parser.add_argument('--postcmd', type=str, help='post cmd command',default=None)
    args=parser.parse_args(args)

    write(args.cmd,outdir=args.outdir,name=args.name,queryhost=args.queryhost,queryport=args.queryport,maxrun=args.maxrun,idlthreads=args.idlthreads,runplans=args.norunplans,time=args.time,postcmd=args.postcmd,flag=args.flag)

