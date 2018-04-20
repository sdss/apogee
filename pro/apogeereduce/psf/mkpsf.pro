;======================================================================
pro mkpsf,psfid,darkid=darkid,flatid=flatid,sparseid=sparseid,fiberid=fiberid,littrowid=littrowid,average=average,clobber=clobber

  dirs=getdir(apodir,caldir,spectrodir,vers)
  caldir=dirs.caldir

  file=dirs.prefix+string(format='("PSF-c-",i8.8)',psfid[0])
  ;if another process is alreadying make this file, wait!
  while file_test(caldir+'psf/'+file+'.lock') do apwait,file,10
  ; does product already exist?
  if file_test(caldir+'/psf/'+file+'.fits') and not keyword_set(clobber) then begin
    print,' PSF file: ', file+'.fits', ' already made'
    return
  endif
  if not keyword_set(fiberid) then fiberid=0
  if not keyword_set(sparseid) then sparseid=0

  print,'making PSF: ', psfid[0]
  ; open .lock file
  openw,lock,/get_lun,caldir+'psf/'+file+'.lock'
  free_lun,lock

  cmjd=getcmjd(psfid)
  ;d=approcess(psfid,darkid=darkid,flatid=flatid,/nocr,nfs=1,/clobber,/doap3dproc)
print,'mkpsf approcess...'
  d=approcess(psfid,darkid=darkid,flatid=flatid,littrowid=littrowid,/nocr,nfs=1,/doap3dproc)
  psffile=dirs.expdir+cmjd+'/'+string(format='(i8.8)',psfid)
  apmkpsf,psffile,caldir+'psf/',sparseid=sparseid,fiberid=fiberid,average=average,clobber=clobber

  file_delete,caldir+'psf/'+file+'.lock'
end

