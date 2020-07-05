#!/usr/bin/env python
#
# JOB_DAEMON.PY - Simple batch job manager.
#

from __future__ import print_function

__authors__ = 'David Nidever <dnidever@noao.edu>'
__version__ = '20181010'  # yyyymmdd

import os
import sys
import numpy as np
import warnings
import socket
import time
#import shutil
import subprocess
#import logging
import tempfile
from . import utils as dln

# Ignore these warnings, it's a bug
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
        

def mkstatstr(n=None):
    """ This returns the stat structure schema or an instance of the stat structure."""
    dtype = np.dtype([('jobid',np.str,20),('name',np.str,100),('user',np.str,100),('timeuse',np.str,100),('status',np.str,10),('queue',np.str,20)])
    if n is None:
        return dtype
    else:
        statstr = np.zeros(n,dtype=dtype)    
        return statstr


def mkjobstr(n=None):
    """ This returns the job structure schema or an instance of the job structure."""
    dtype = np.dtype([('host',np.str,20),('jobid',(np.str,100)),('input',np.str,1000),('dir',np.str,500),('name',np.str,100),('scriptname',np.str,200),
                      ('logfile',np.str),('submitted',np.bool),('done',np.bool),('begtime',np.float64),('endtime',np.float64),('duration',float)])
    if n is None:
        return dtype
    else:
        jobstr = np.zeros(n,dtype=dtype)
        jobstr['submitted'] = False
        jobstr['done'] = False
        return jobstr


def mkrunbatch():
    curdir = os.getcwd()
    batchfile = os.path.join(curdir,'runbatch')
    if os.path.exists(batchfile) is False:
        lines = []
        lines.append("if test $# -eq 0\n")
        lines.append("then\n")
        lines.append("  echo 'Syntax - runbatch program'\n")
        lines.append("else\n")
        lines.append("  echo 'Log file: '$1'.log'\n")
        lines.append("  ( nohup  $1 > $1.log 2>&1 ) &\n")
        lines.append("  echo 'JOBID='$!\n")
        lines.append("fi\n")
        dln.writelines(batchfile,lines,overwrite=True,raw=True)
        os.chmod(batchfile,0o755)
    return batchfile


def mkidlbatch():
    curdir = os.getcwd()
    batchfile = os.path.join(curdir,'idlbatch')
    if os.path.exists(batchfile) is False:
        try:
            out = subprocess.check_output(['which','idl'],stderr=subprocess.STDOUT,shell=False)
        except subprocess.CalledProcessError:
            raise Exception("IDL program not available")
        idlprog = out.decode().strip()
        if os.path.exists(idlprog) is False:
            raise Exception("IDL program "+idlprog+" not found")
        lines = []
        lines.append("if test $# -eq 0\n")
        lines.append("then\n")
        lines.append("  echo 'Syntax - idlbatch idl.batch'\n")
        lines.append("else\n")
        lines.append("  echo 'Log file: '$1'.log'\n")
        lines.append("  ( nohup "+idlprog+" < $1 > $1.log 2>&1 ) &\n")
        lines.append("  echo 'JOBID='$!\n")
        lines.append("fi\n")
        dln.writelines(batchfile,lines,overwrite=True,raw=True)
        os.chmod(batchfile,0o755)
    return batchfile


def check_diskspace(indir=None,updatestatus=False):
    """ This checks if there's enough free disk space. """
    if indir is None: raise ValueError("Must give directory")
    statvfs = os.statvfs(indir)
    available = statvfs.f_frsize * statvfs.f_bavail / 1e6  # in MB
    if updatestatus: print('Disk Space: '+str(available[0]),' MB available')
    # Not enough disk space available
    if available<100.:
        raise Exception('NOT enough disk space available')


def check_killfile(jobs=None,hyperthread=True):
    """ This checks for a kill file and if found kills all the active jobs. """
    if jobs is None: raise ValueError("No jobs structure input")
    killfile = 'killjobs'
    if os.path.exists(killfile) is True:
        sub, = np.where((jobs['submitted']==1) & (jobs['done']==0))
        nsub = len(sub)
        print('Kill file found.  Killing all '+str(nsub)+' job(s)')
        for i in range(nsub):
            # Killing the job
            print('Killing '+jobs['name'][sub[i]]+'  JobID='+jobs['jobid'][sub[i]])
            if hyperthread is False:
                out = subprocess.check_output(['qdel',jobs['jobid'][sub[i]]],stdout=sf,stderr=subprocess.STDOUT,shell=False)
            else:
                out = subprocess.check_output(['kill','-9',jobs['jobid'][sub[i]]],stderr=subprocess.STDOUT,shell=False)
        # Remove the kill file
        print('Deleting kill file "'+killfile+'"')
        os.remove(killfile)
        return True
    return False


def makescript(inp=None,indir=None,name=None,prefix=None,hyperthread=True,idle=False):
    """This makes job scripts for the job_daemon program.

    Parameters
    ----------
    inp : string list or array
          The command to execute.  Can be an array.  Must input absolute path names.
          idlbatch will be used if it's an IDL command.  If the command is a
          series of unix commands separated by commas then these
          will be put on separate lines.
    indir : string list or array
          The directory to put the job script in.
    name : string list or array, optional
          The name to call the job script (without the '.sh' or '.batch' ending).
          If this is not provided then an autogenerated name will
          be made and returned.
    prefix : string, optional
          The prefix to use for the job script name. "pr" by default.
    hyperthread : bool, optional
          Not on a job server but one with multiple hyperthreaded
                   processors.
    idle : bool, optional
          This is an IDL command, otherwise a SHELL command.

    Returns
    -------
    scriptname : string array
               The absolute names of the scripts.
    Job scripts output to the directories and with the names specified.


    Example
    -------

    .. code-block:: python

        scriptname = makescript(inp,dir,name,hyperthread=True)

    """
  
    # Not enough inputs
    if (inp is None):
        raise ValueError('No input given')

    ninp = dln.size(inp)
    ndir = dln.size(indir)
    nname = dln.size(name)

    # Not enough directories input
    if (ndir>0) & (ndir!=ninp):
        raise ValueError('INPUT and DIRECTORIES are of different size')

    # Current directory
    curdir = os.getcwd()

    # No directories input
    if ndir==0: indir = np.repeat(curdir,ninp)
    if ndir==1: indir = np.repeat(indir,ninp)    # multiple commands in same dir

    # Construct names
    if (nname==0):
        name = np.zeros(ninp,dtype=(np.str,200))
        if prefix is not None:
            pre = dln.first_el(prefix)
        else:
            pre = 'job'
        for i in range(ninp):
            tid,tfile = tempfile.mkstemp(prefix=pre,dir=indir[i])
            os.close(tid)   # mkstemp opens the file, close it
            name[i] = os.path.basename(tfile)

    # Make scriptnames
    scriptname = np.array(dln.strjoin(dln.pathjoin(indir,name),'.sh'),ndmin=1)
    # Script loop
    for i,input1 in enumerate(np.array(inp,ndmin=1)):
        base = str(name[i])
        sname = str(indir[i])+'/'+base+'.sh'

        #------
        # PBS
        #------
        if hyperthread is False:
            # IDL command
            if idle is True:
                # Making an IDL batch file
                bname = str(indir[i])+'/'+base+'.batch'
                dln.writelines(bname,input1,overwrite=True,raw=True)
                # The execution command
                cmd = 'idl < '+base+'.batch'
            # SHELL command
            else:
                # The execution command
                cmd = input1
                # If there are commas in the line then break it up into multiple lines
                if cmd.find(';') != -1:
                    cmd = cmd.replace(';','\n;')
                    cmd = cmd.split('j')

            # Make the command
            #----------------------
            lines = []
            lines.append('#!/bin/sh\n')
            lines.append('#PBS -l nodes=1:ppn=1\n')
            lines.append('#PBS -l walltime=96:00:00\n')
            lines.append('#PBS -o '+base+'.report.out\n')
            lines.append('#PBS -e '+base+'.error.out\n')
            #lines.append('#PBS -M dln5q@virginia.edu\n')
            lines.append('#PBS -V\n')
            lines.append('\n')
            lines.append('echo Running on host `hostname`\n')
            lines.append('echo Time is `date`\n')
            lines.append('echo "Nodes used for this job:"\n')
            lines.append('echo "------------------------"\n')
            lines.append('cat $PBS_NODEFILE\n')
            lines.append('echo "------------------------"\n')
            lines.append('\n')
            lines.append('cd '+indir[i]+'\n')
            for j in range(len(cmd)): lines.append(cmd[j]+'\n')
            lines.append('\n')
            lines.append('# print end time\n')
            lines.append('echo\n')
            lines.append('echo "Job Ended at `date`"\n')
            lines.append('echo\n')
            # Writing the file
            dln.writelines(scriptname[i],lines,overwrite=True,raw=True)
            # Print info
            print('PBS script written to: '+str(scriptname[i]))

        #----------------
        # Hyperthreaded
        #----------------
        else:
            # Just make batch file
            # treat shell and idl the same
            # The execution command
            cmd = input1
            # If there are commas in the line then break it up into multiple lines
            if cmd.find(';') != -1:
                cmd = cmd.replace(';','\n;')
                cmd = cmd.split(';')
            # IDL files should end in .batch
            if idle is True: scriptname[i] = str(indir[i])+'/'+base+'.batch'
            # Writing the file
            dln.writelines(scriptname[i],cmd,overwrite=True,raw=True)
            # Make SHELL scripts executable
            if idle is False: os.chmod(scriptname[i],0o755)
            # Print info
            print('HYPERTHREAD script written to: '+str(scriptname[i]))

    # Erase the temporary files that mkstemp makes
    dln.remove(dln.pathjoin(indir,name),allow=True)

    if dln.size(scriptname)==1: scriptname=scriptname[0]
    return scriptname


def submitjob(scriptname=None,indir=None,hyperthread=True,idle=False):
    """ This submits new jobs.

    Parameters
    ----------
    scriptname : string
           Name of the script to run.
    indir : string, optional
           Directory to change to before script is executed.
    hyperthread : bool, optional
          Not on a job server but one with multiple hyperthreaded
             processors.  By default, this is True.
    idle : bool, optional
         This is an IDL program.  By default his is False.

    Returns
    -------
    jobid : string
         The job ID for this job.
    logfile : string
         Name of the output log file for this job/script.
    The script `scriptname` will be run.

    Example
    -------

    .. code-block:: python

        jobid,logfile = submitjob(scriptname,indir,hyperthread=True)

    """
    if scriptname is None: raise ValueError("scriptname must be input")

    curdir = os.getcwd()
    # Submitting the job
    if hyperthread is False:
        try:
            out = subprocess.check_output('qsub '+scriptname,stderr=subprocess.STDOUT,shell=True)
        except subprocess.CalledProcessError:
            raise Exception("Problem submitting PBS job")
        jobid = dln.first_el(out)
        logfile = scriptname+'.log'
    else:
        if idle is True:
            batchprog = mkidlbatch()
        else:
            batchprog = mkrunbatch()
        if indir is not None: os.chdir(indir)
        try:
            out = subprocess.check_output(batchprog+' '+scriptname,stderr=subprocess.STDOUT,shell=True)
        except:
            raise Exception("Problem submitting shell job")
        if indir is not None: os.chdir(curdir)
        # Get the JOBID
        out = out.decode().split('\n')
        jobid_ind = dln.grep(out,'^JOBID=',index=True)
        njobid_ind = len(jobid_ind)
        jobid = out[jobid_ind[0]].split('=')[1]
        jobid = jobid.strip()
        # Get the logfile
        logfile_ind = dln.grep(out,'^Log file: ',index=True)
        nlogfile_ind = len(logfile_ind)
        logfile = out[logfile_ind[0]].split(':')[1]
        logfile = logfile.strip()
    # Printing info
    print('Submitted '+scriptname+'  JobID='+jobid)

    return jobid, logfile


def getstat(jobid=None,hyperthread=True):
    """Get the status of jobs

    This checks the status of jobs for the job_daemon program.
    If no jobs are found in the queue then an empty
    statstr structure is returned.

    Parameters
    ----------
    jobid : string or int
          Specific JOBID to check.
    hyperthread : bool, optional
          Not on a PBS machine but one with multipe hyperthreaded
          processors running simultaneously.  Default is True.

    Results
    -------
    statstr : numpy structured array
          Status structure.

    Example
    -------

    .. code-block:: python

        stat = getstat(jobid)

    """
    if jobid is None: raise ValueError("Must input jobid")
    njobid = dln.size(jobid)

    # PBS
    #--------
    if hyperthread is False:
        if njobid>0: 
            out = subprocess.check_output(['qstat',jobid[0]],stderr=subprocess.STDOUT,shell=False)
        else:
            out = subprocess.check_output(['qstat'],stderr=subprocess.STDOUT,shell=False)
        nout = np.sum(out.strip() != '')
        gd = dln.grep(out,'^'+jobid,index=True)
        ngd = len(gd)
        if ngd>0:
            statlines = out[gd[0]]
            nstat = len(statlines)
        else:
            statlines = None
        # Some jobs in queue
        if statlines is not None:
            arr = statlines
            arr = arr.split(' ')
            statstr = mkstatstr(nstat)
            statstr['jobid'] = arr[0]
            statstr['name'] = arr[1]
            statstr['user'] = arr[2]
            statstr['timeuse'] = arr[3]
            statstr['status'] = arr[4]
            statstr['queue'] = arr[5]
        # No jobs in queue
        else:
            statstr = mkstatstr(1)

    # Hyperthreaded.  Need a jobid
    #-----------------------------
    else:
        # No JOBID input
        if njobid==0:
            print('Need JOBID with /hyperthread')
            return mkstatstr(1)

        try:
            out = subprocess.check_output(['ps','-o','pid,user,etime,command','-p',str(dln.first_el(jobid))],
                                          stderr=subprocess.STDOUT,shell=False)
        except:
            statstr = mkstatstr(1)
            statstr['jobid'] = jobid
            statstr['queue'] = 'hyperthread'
            return statstr
        # can put in the column that you want
        # ps -o etime -p jobid
        out = dln.strsplit(out.decode(),'\n')
        out = dln.strip(out)
        gd = dln.grep(out,'^'+str(jobid),index=True)
        ngd = len(gd)
        if ngd>0:
            statlines = np.array(out,ndmin=1)[gd[0]]
        else:
            statlines = None
        # Some jobs in queue
        if statlines is not None:
            arr = statlines.split()
            statstr = mkstatstr(1)
            statstr['jobid'] = arr[0]
            statstr['user'] = arr[1]
            # CAN'T get the name.
            statstr['timeuse'] = arr[2]
            statstr['status'] = 'R'
            statstr['queue'] = 'hyperthread'
        # No jobs in queue
        else:
            statstr = mkstatstr(1)
            statstr['jobid'] = jobid
            statstr['queue'] = 'hyperthread'
    return statstr


def checkjobs(jobs=None,hyperthread=True):
    """ Check on the status of the running jobs"""
    if jobs is None: raise ValueError("jobs must be input")
    sub, = np.where((jobs['submitted']==True) & (jobs['done']==False))
    nsub = len(sub)
    nfinished = 0
    for i in range(nsub):
        # Checking status
        jobid = jobs[sub[i]]['jobid']
        statstr = getstat(jobid,hyperthread=hyperthread)
        # Job done
        if statstr['status']=='':
            if nfinished==0: print('')
            print(time.ctime()+'  Input '+str(sub[i]+1)+' '+jobs['name'][sub[i]]+' JobID='+jobs['jobid'][sub[i]]+' FINISHED')
            jobs['done'][sub[i]] = True
            jobs['endtime'][sub[i]] = time.time()/3600/24   # in days
            jobs['duration'][sub[i]] = (jobs['endtime'][sub[i]] - jobs['begtime'][sub[i]])*3600*24   # in sec
            nfinished += 1
            # Check for errors as well!! and put in jobs structure
    return jobs


def status_update(jobs=None):
    """ Print out the status update."""
    if jobs is None: raise ValueError("jobs must be input")
    njobs = dln.size(jobs)                                                   # Number of total jobs
    n_inqueue = np.sum((jobs['submitted']==True) & (jobs['done']==False))    # Number of jobs still in queue  
    n_nosubmit = np.sum(jobs['submitted']==False)                            # Number of jobs left to do 
    n_finished = np.sum(jobs['done']==True)                                  # Number of jobs finished 
    # Print the status
    print('')
    print(time.ctime())
    print(('Jobs Summary: %d total, %d finished, %d running, %d left') % (njobs,n_finished,n_inqueue,n_nosubmit))


def job_daemon(inp=None,dirs=None,inpname=None,nmulti=4,prefix="job",hyperthread=True,idle=False,
               waittime=0.2,statustime=60):
    """This program is a simple batch job manager

    NOTE:  If you want to "kill" all of the jobs create a file "killjobs"
    in the same directory that JOB_DAEMON is being run in and all of the
    PBS jobs will be killed.

    Parameters
    ----------
    inp : string array/list
          A string array with the commands (i.e. jobs) to be run.
    dirs : string array/list
          The directories in which the commands are to be run.
    inpname : string array/list, optional
           The name to be used for the script (without the path or
           the .sh/.batch ending.  This is normally autogenerated.
    nmulti : int, optional
           How many nodes to run these jobs on.  Default is 4.
    prefix : string
          The prefix for the script names.  Default is "job".
    hyperthread : bool, optional
           Not on a PBS server but one that has multiple processors
           hyperthreaded.  Run multiple jobs at the same time on
           the same server.  Default is True.
    idle : bool, optional
         This is an IDL command, otherwise a SHELL command.  Default is False.
    statustime : float or int, optional
           The time between status updates.  However, the status
           will always be updated if something has actually changed.
           Default is 60.
    waittime : float or int, optional
           Time to wait between checking the running jobs.  Default is 0.2 sec.

    Results
    -------
    jobs : numpy structured array
          The jobs structure with information on the JOBID, script name, etc.
    Jobs are run in a batch mode.

    Example
    -------

    .. code-block:: python

        jobs = job_daemon(input,dirs,hyperthread=True,nmulti=5)

    """

    # How many input lines
    if inp is None:
        raise ValueError("Nothing input")
    ninp = dln.size(inp)

    # Current directory
    curdir = os.getcwd()

    # Checking DIRS array
    if dirs is None:
        dirs = np.repeat(curdir,ninp)
    else:
        ndirs = dln.size(dirs)
        if ndirs!=ninp:
            raise ValueError('DIRS array must be same size as INPUT')            
        if ndirs==1: dirs = np.repeat(dirs,ninp)

    # Check INPNAME array
    if inpname is not None:
        ninpname = len(inpname)
        if ninpname!=ninp:
            raise ValueError('INPNAME array must be same size as INPUT')

    # Defaults
    if waittime<0.1: waittime=0.1
    if statustime<1: statustime=1

    # Host name
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    # Which IDL are we using?
    if idle is True:
        try:
            out = subprocess.check_output(['which','idl'],stderr=subprocess.STDOUT,shell=False)
        except subprocess.CalledProcessError:
            raise Exception("IDL program not available")
        idlprog = out.strip()
        if os.path.exists(idlprog) is False:
            raise Exception("IDL program "+idlprog+" not found")

    print('---------------------------------')
    print(' RUNNING JOB_DAEMON for '+str(ninp)+' JOB(S)')
    print('---------------------------------')
    print('Host='+host)
    print('Nmulti='+str(nmulti))
    
    t0 = time.time()
    timesincelaststatus = time.time()  # initializing the update time

    #--------
    # DAEMON
    #--------
    # -Keep submitting jobs until nmulti is reached
    # -Check every minute or so to see how many jobs are still
    #  running.  If it falls below nmulti and more jobs are left then
    #  submit more jobs
    # -Don't return until all jobs are done.

    # Initialize the "jobs" structure
    # id will be the ID from Pleione
    jobs = mkjobstr(ninp)
    jobs['host'] = host
    jobs['input'] = inp
    njobs = ninp

    # Loop until all jobs are done
    # On each loop check the pleione queue and figure out what to do
    count = np.longlong(0)
    endflag = 0
    while (endflag==0):
        # Status update
        dtstatus_sec = time.time()-timesincelaststatus
        if (dtstatus_sec>statustime):
            updatestatus = True
            timesincelaststatus = time.time()
        else:
            updatestatus = False
  
        # Check disk space
        check_diskspace(dln.first_el(dirs))
        # Check for kill file
        if check_killfile(jobs) is True: return jobs

        # Check status of running jobs
        #-----------------------------
        jobs = checkjobs(jobs,hyperthread=hyperthread)

        # Status update
        #---------------
        if updatestatus is True: status_update(jobs)

        # Submit new jobs
        #----------------
        n_inqueue = np.sum((jobs['submitted']==True) & (jobs['done']==False))  # Number of jobs still in queue  
        n_nosubmit = np.sum(jobs['submitted']==False)                            # Number of jobs left to do 
        nnew = dln.limit(nmulti-n_inqueue,0,n_nosubmit)
        if (nnew>0):
            # Get the indices of new jobs to be submitted
            nosubmit, = np.where(jobs['submitted']==False)
            newind = nosubmit[0:nnew]
            # Update immediately if there are new jobs to submit
            print('')
            print(time.ctime())
            print('Updating Queue: '+str(n_inqueue)+' JOB(S) running, out of '+str(nmulti)+' Maximum. Submitting '+str(nnew)+' more job(s)')

            # Loop through the new submits
            for i in range(nnew):
                print('')
                cmd = jobs[newind[i]]['input']
                if idle is True: cmd = 'IDL>'+str(cmd)
                print('Input '+str(newind[i]+1)+'  Command: >>'+str(cmd)+'<<')
                # Make script
                name = None
                if inpname is not None: name = inpname[newind[i]]  # use input name      
                scriptname = makescript(jobs[newind[i]]['input'],indir=dirs[newind[i]],name=name,
                                        prefix=prefix,hyperthread=hyperthread,idle=idle)
                name = os.path.basename(os.path.splitext(scriptname)[0])
                # Submitting the job
                jobid, logfile = submitjob(scriptname,dirs[newind[i]],hyperthread=hyperthread,idle=idle)
                # Updating the jobs structure
                jobs['submitted'][newind[i]] = True
                jobs['jobid'][newind[i]] = jobid
                jobs['name'][newind[i]] = name
                jobs['dir'][newind[i]] = dirs[newind[i]]
                jobs['scriptname'][newind[i]] = scriptname
                jobs['logfile'][newind[i]] = logfile
                jobs['begtime'][newind[i]] = time.time()/3600/24   # in days

        # Are we done?
        #-------------
        ndone = np.sum(jobs['done']==True)
        if (ndone==njobs): endflag=1
        # Wait a bit
        #--------------
        if (endflag==0): time.sleep(waittime)
        # Increment the counter
        count += 1

    print('DONE')
    print('dt = '+str(time.time()-t0)+' sec')
    return jobs
