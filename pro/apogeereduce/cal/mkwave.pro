;======================================================================
pro mkwave,waveid,darkid=darkid,flatid=flatid,psfid=psfid,clobber=clobber

  dirs=getdir(apodir,caldir,spectrodir,vers)
  caldir=dirs.caldir
  file=dirs.prefix+string(format='("Wave-",i8.8)',waveid[0])
  ;if another process is alreadying make this file, wait!
  while file_test(caldir+'wave/'+file+'.lock') do apwait,file,10
  ; does product already exist?
  if file_test(caldir+'/wave/'+file+'.dat') and not keyword_set(clobber) then begin
    print,' Wavecal file: ', file+'.dat', ' already made'
    return
  endif

  print,'making wave: ', waveid
  ; open .lock file
  openw,lock,/get_lun,caldir+'wave/'+file+'.lock'
  free_lun,lock

  cmjd=getcmjd(psfid)
  mkpsf,psfid,darkid=darkid,flatid=flatid
;  ntrace=mkepsf(psfid,dark=darkid,flat=flatid)
  ;w=approcess(waveid,dark=darkid,flat=flatid,psf=psfid,flux=0,/clobber,/doproc)
  w=approcess(waveid,dark=darkid,flat=flatid,psf=psfid,flux=0,/doproc)
  psffile = caldir+'psf/'+string(format='(i8.8)',psfid)
  wavefile = dirs.expdir+cmjd+'/'+string(format='(i8.8)',waveid)
  apwavecal,wavefile,psfid=psffile,/shortname,/save;,/pl

  file_delete,caldir+'wave/'+file+'.lock'

end
