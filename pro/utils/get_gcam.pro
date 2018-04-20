function get_gcam,cmjd

dirs=getdir()

; get guider summary information
file=dirs.expdir+'/'+cmjd+'/gcam-'+cmjd+'.fits'

if file_test(file+'.lock') then apwait,file+'.lock',10
openw,lock,file+'.lock',/get_lun
free_lun,lock

; if guider summary file doesn't exist, create it
if ~file_test(file) then spawn,['gcam_process','-m',cmjd,'-i',dirs.instrument,'-o',file],/noshell
gcam=mrdfits(file,1)

file_delete,file+'.lock',/allow
return,gcam
end
