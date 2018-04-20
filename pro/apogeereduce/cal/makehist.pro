pro makehist,mjd,clobber=clobber,dark=dark

cmjd=string(format='(i5.5)',mjd)

dirs=getdir(apodir,caldir,specdir,apovers,libdir,datadir=datadir,apred_vers=apred_vers)

outdir=dirs.expdir+cmjd+'/'
file_mkdir,outdir
file=apogee_filename('Hist',mjd=cmjd,chip='c')

; Is file already being created?
lockfile=file+'.lock'
while file_test(lockfile) do apwait,lockfile,10

; does file already exist?
if file_test(file) and not keyword_set(clobber) then begin
  print,' Hist file: ', file, ' already made'
  return
endif

; create lock file
openw,lock,/get_lun,lockfile
free_lun,lock

; call python script to make apHist file
if keyword_set(dark) then  $
spawn,'makehist '+cmjd+' --apred ' +apred_vers + ' --darkid '+ string(dark) else $
spawn,'makehist '+cmjd+' --apred ' +apred_vers

; remove lock file
file_delete,lockfile,/allow_non


end
