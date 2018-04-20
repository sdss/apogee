;======================================================================

pro mkbpm,bpmid,darkid=darkid,flatid=flatid,badrow=badrow,clobber=clobber

 dirs=getdir()
 caldir=dirs.caldir

 file=apogee_filename('BPM',num=bpmid,chip='c')
 ;if another process is alreadying make this file, wait!
 while file_test(file+'.lock') do apwait,file,10
 ; does product already exist?
 if file_test(file) and not keyword_set(clobber) then begin
   print,' BPM file: ', file, ' already made'
   return
 endif

 print,'making BPM: ', bpmid
 ; open .lock file
 openw,lock,/get_lun,file+'.lock'
 free_lun,lock

 chips=['a','b','c']

 for ichip=0,2 do begin
  chip=chips[ichip]

  mask=intarr(2048,2048)

  ; bad pixels from dark frame
  file=apogee_filename("Dark",chip=chip,num=darkid)
  darkmask=mrdfits(file,3)
  bad=where(darkmask gt 0)
  mask[bad]=mask[bad] or maskval('BADDARK')

  ; bad pixels from flat frame
  if flatid gt 0 then begin
    file=apogee_filename("Flat",chip=chip,num=flatid)
    flatmask=mrdfits(file,3)
    bad=where(flatmask gt 0)
    mask[bad]=mask[bad] or maskval('BADFLAT')
  endif else flatmask=darkmask*0

  ; flag them both as bad pixel (for historical compatibility?)
  bad=where((darkmask or flatmask) gt 0)
  mask[bad]=mask[bad] or maskval('BADPIX')
  if keyword_set(badrow) then begin
    for i=0,n_elements(badrow)-1 do begin
      if badrow.chip eq ichip then mask[*,badrow.row]=mask[*,badrow.row] or maskval('BADPIX')
    endfor
  endif
  file=apogee_filename('BPM',chip=chip,num=bpmid)
  mwrfits,mask,file,/create
 endfor

 file_delete,file+'.lock'

end
