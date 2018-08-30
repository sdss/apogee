#!/usr/bin/env python
import sys
import os
import argparse


def write(cmd,outdir='slurm/',cwd=None,queryhost=None,queryport=None,maxrun=None,idlthreads=1,runplans=True,time='240:00:00',name=None) :

    try :
        os.mkdir(outdir)
    except:
        pass

    if name is None : name=cmd.split()[0].strip('"')
    file = outdir+name
    f=open(file,'w')

    f.write('#!/bin/csh\n')
    f.write('#SBATCH --account=sdss-kp\n')
    f.write('#SBATCH --partition=sdss-kp\n')
    f.write('#SBATCH --time='+time+'\n')
    f.write('#SBATCH --ntasks=16\n')
    f.write('#SBATCH --nodes=1\n')
    f.write('#SBATCH -o '+os.path.basename(file)+'.out\n' )
    f.write('#SBATCH -e '+os.path.basename(file)+'.out\n' )
    if runplans :
        f.write('setenv QUERYHOST '+queryhost+'\n' )
        f.write('setenv QUERYPORT '+str(queryport)+'\n')
        f.write('setenv APOGEE_MAXRUN '+str(maxrun)+'\n' )
    f.write('setenv IDL_CPU_TPOOL_NTHREADS '+str(idlthreads)+'\n' )

    if cwd is None : cwd=os.getcwd()
    f.write('cd '+cwd+'\n' )
    if runplans : f.write('runplans ' + cmd+'\n' )
    else : f.write(cmd+'\n' )
    f.write('wait\n' )
    f.write('echo DONE\n' )
    f.close()
    os.chmod(file,0o770)

def main(args) :

    parser = argparse.ArgumentParser(
        prog=os.path.basename(args[0]),
        description='Makes a SLURM batch file')
    parser.add_argument('cmd', type=str, help='cmd')
    parser.add_argument('--outdir', type=str, help='directory',default='slurm/')
    parser.add_argument('--runplans', type=str, help='directory',default='runplans ')
    parser.add_argument('--queryport', type=int, help='port to use for queue manager',default=1050)
    parser.add_argument('--queryhost', type=str, help='host to use for queue manager',default=os.uname()[1])
    parser.add_argument('--maxrun', type=int, help='maximum jobs to run at a time',default=1)
    parser.add_argument('--time', type=str, help='maximum wall clock time',default='240:00:00')
    parser.add_argument('--idlthreads', type=int, help='maximum IDL threads',default=1)
    args=parser.parse_args(args)

    write(args.cmd,outdir=args.outdir,queryhost=args.queryhost,queryport=args.queryport,maxrun=args.maxrun,idlthreads=args.idlthreads,runplans=args.runplans,time=args.time)

