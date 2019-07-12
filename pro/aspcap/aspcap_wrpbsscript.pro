pro aspcap_wrpbsscript,name,path_src,ncpus,jobsid,libsize,email=email,qname=qname,group=group,queue=queue,workdir=workdir

; This procedure creates a PBS job script file called name.pbs. This
; file is used for running FERRE in a queue system. The name of the
; queue and group can be passed to the routine using the keywords as
; well as the e-mail address. The path to FERRE needs to be specified
; as well as the number of nodes to be used. The output of the
; procedure is the file and a job id (the name of the file with the -
; removed). 
;
;
; INPUT:   
; 
;    name          -string       the jobs script filename
;    path_src      -string       the directory where ferre.out is
;    ncpus           -int          number of nodes to be used
; OUTPUT:
;    jobsid        -string       the job id as specified by the user
;     
; KEYWORDS:
;
;    email         -string       the e-mail address where PBS send the
;                                message once the jobs left the queue
;                                system
;    qname         -string       the name of the PBS queue (apogee by default)
;    group         -string       the name of the PBS group (apogee by default)
;
;
;
; By Ana Elia Garcia Perez - July 2010
;
;
;


  if n_elements(qname) eq 0 then qname='apogee'
  if n_elements(group) eq 0 then group='apogee'
 
  openw,numl,name+'.pbs',/get_lun
  printf,numl,'#!/bin/bash'
  printf,numl,'#PBS -l select=1:ncpus='+strcompress(ncpus,/remove_all)+':mem='+libsize+'GB'
  printf,numl,'#PBS -l walltime=99:00:00'
  printf,numl,'#PBS -W group_list='+group
  printf,numl,'#PBS -q '+qname
  printf,numl,'#PBS -W umask=002'
  date=(systime(/julian,/utc)-2455562.500000)*1d6      ; that date is january 2011 0 ut
  date=string(date,format='(I11)')
;  limit=strpos(date,'.') 
;  if limit gt -1 then begin date=strmid(date,0,limit)
;  jobsid='aspc'+strcompress(date,/remove_all)
;trim(string(date,format='(I10)'),2)
  jobsid=name
  posit=strpos(jobsid,'-')  
  while posit ge 0 do begin 
    jobsid=strmid(jobsid,0,posit)+strmid(jobsid,posit+1,strlen(jobsid)-(posit+1))
    posit=strpos(jobsid,'-')
  endwhile
  if strlen(jobsid) gt 15 then begin 
     xerr='Warning: the name of the jobs has more than 15 allowed number of characters'
     print,xerr
  endif
  printf,numl,'#PBS -N '+jobsid
  if n_elements(email) gt 0 then begin
     printf,numl,'#PBS -m e'                  ; send an e-mail when finished 
     printf,numl,'#PBS -M '+email
  endif
  printf,numl,'#PBS -V'                   ; do I need to change directory?
  printf,numl,'#PBS -z'
  printf,numl,'#PBS -j oe'
  printf,numl,''
  if queue eq 1 then begin
     printf,numl,'echo pbs_o_workdir: $PBS_O_WORKDIR'
     printf,numl,'cd $PBS_O_WORKDIR'
     if keyword_set(workdir) then printf,numl,'workdir: ',workdir
     if keyword_set(workdir) then printf,numl,'cd ',workdir
  endif
  printf,numl,path_src+'/ferre.x 1> '+name+'.out 2> '+name+'.err.out' 
  printf,numl,''
  ; clean up temporary files
  ;printf,numl,'rm input.nml lib'
  printf,numl,'rm input.nml'
  close,numl 
  free_lun,numl
  file_chmod,name+'.pbs','775'o

   
end
