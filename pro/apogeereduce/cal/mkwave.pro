;======================================================================
pro mkwave,waveid,name=name,darkid=darkid,flatid=flatid,psfid=psfid,clobber=clobber,nowait=nowait,nofit=nofit

  if ~keyword_set(name) then name=string(waveid[0])

  dirs=getdir(apodir,caldir,spectrodir,vers)
  caldir=dirs.caldir
  file=dirs.prefix+string(format='("Wave-",i8.8)',name)
  ;if another process is alreadying make this file, wait!
  while file_test(caldir+'wave/'+file+'.lock') do begin
    if keyword_set(nowait) then return
    apwait,file,10
  endwhile
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
  w=approcess(waveid,dark=darkid,flat=flatid,psf=psfid,flux=0,/doproc)

  ; new Python version! 
  if keyword_set(nofit) then nofit='--nofit' else nofit=''
  cmd=['apmultiwavecal','--name',name,nofit,'--plot','--hard','--inst',dirs.instrument,'--verbose']
  for i=0,n_elements(waveid)-1 do cmd=[cmd,string(waveid[i])]
  spawn,cmd,/noshell

  ;psffile = caldir+'psf/'+string(format='(i8.8)',psfid)
  ;wavefile = dirs.expdir+cmjd+'/'+string(format='(i8.8)',waveid)
  ;apmultiwavecal,wavefile,psfid=psffile,name=name,/save,clobber=clobber ;,/pl

  file_delete,caldir+'wave/'+file+'.lock'

end
